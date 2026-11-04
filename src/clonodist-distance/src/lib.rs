//! Deterministic row tiles. Bounded output blocks; no dense repertoire allocation.
use clonodist_core::{add, distance, germline, mul, PreparedSequence, Receptor, Result, CDR3};
use rayon::prelude::*;
use std::collections::BTreeMap;
#[derive(Clone, Debug)]
pub struct Config {
    pub threads: usize,
    pub tile: usize,
    pub memory: usize,
    pub optimized: bool,
    pub radius: Option<u32>,
}
/// Precompute only used V-loop combinations; cache is capped by the run budget.
pub struct Engine<'a> {
    pub q: &'a [Receptor],
    pub t: &'a [Receptor],
    pub config: Config,
    cache: Vec<Vec<u32>>,
    q_ids: Vec<Vec<usize>>,
    t_ids: Vec<Vec<usize>>,
    widths: Vec<usize>,
    q_prepared: Vec<Vec<PreparedSequence>>,
    t_prepared: Vec<Vec<PreparedSequence>>,
    pool: rayon::ThreadPool,
    radius_order: Vec<usize>,
    radius_groups: Vec<(usize, usize)>,
}
impl<'a> Engine<'a> {
    pub fn new(q: &'a [Receptor], t: &'a [Receptor], config: Config) -> Result<Self> {
        if q.is_empty() || t.is_empty() {
            return Err("empty query/target".into());
        }
        if config.threads == 0 || config.threads > 256 || config.tile == 0 {
            return Err("threads must be 1..256; tile must be positive".into());
        }
        let rowbytes = t.len().checked_mul(8).ok_or("row size overflow")?;
        let blockbytes = rowbytes
            .checked_mul(config.tile)
            .ok_or("tile size overflow")?;
        let input = q.iter().chain(t).try_fold(0usize, |sum, r| {
            let n = serde_size(r)?;
            sum.checked_add(n).ok_or("input size overflow".to_owned())
        })?;
        let row_scratch = t.len().checked_mul(256).ok_or("row scratch overflow")?;
        let planned = input
            .checked_add(row_scratch)
            .ok_or("row scratch overflow")?
            .checked_add(blockbytes)
            .and_then(|x| x.checked_add(config.threads.checked_mul(2 * 1024 * 1024)?))
            .ok_or("memory plan overflow")?;
        if planned > config.memory / 2 {
            return Err(format!(
                "memory plan {planned} exceeds conservative half-budget {}",
                config.memory / 2
            ));
        }
        if q[0].chains.is_empty() || q[0].chains.len() > 2 {
            return Err("chain configuration mismatch".into());
        }
        if q.iter()
            .chain(t)
            .any(|r| r.chains.len() != q[0].chains.len())
        {
            return Err("chain configuration mismatch".into());
        }
        let mut identities = BTreeMap::new();
        for r in q.iter().chain(t) {
            for (c, ch) in r.chains.iter().enumerate() {
                if ch.cdr3.len() > 127
                    || ch.cdr3.iter().any(|b| *b > 23)
                    || ch
                        .loops
                        .iter()
                        .any(|s| s.len() > 127 || s.iter().any(|b| *b > 23))
                {
                    return Err("invalid encoded sequence".into());
                }
                if let Some(previous) = identities.insert((c, ch.v.clone()), ch.loops.clone()) {
                    if previous != ch.loops {
                        return Err("inconsistent germline loops for V key".into());
                    }
                }
            }
        }
        let mut cache = vec![];
        let mut q_ids = vec![];
        let mut t_ids: Vec<Vec<usize>> = vec![];
        let mut widths = vec![];
        if config.optimized {
            for c in 0..q[0].chains.len() {
                let a: BTreeMap<_, _> = q.iter().map(|r| (&r.chains[c].v, &r.chains[c])).collect();
                let b: BTreeMap<_, _> = t.iter().map(|r| (&r.chains[c].v, &r.chains[c])).collect();
                let entries = a.len().checked_mul(b.len()).ok_or("cache size overflow")?;
                if entries.checked_mul(512).ok_or("cache size overflow")? > config.memory / 8 {
                    return Err("germline cache exceeds budget; use --scalar".into());
                }
                let av: Vec<_> = a.values().copied().collect();
                let bv: Vec<_> = b.values().copied().collect();
                let ai: BTreeMap<_, _> = a
                    .keys()
                    .enumerate()
                    .map(|(i, v)| ((*v).clone(), i))
                    .collect();
                let bi: BTreeMap<_, _> = b
                    .keys()
                    .enumerate()
                    .map(|(i, v)| ((*v).clone(), i))
                    .collect();
                q_ids.push(q.iter().map(|r| ai[&r.chains[c].v]).collect());
                t_ids.push(t.iter().map(|r| bi[&r.chains[c].v]).collect());
                widths.push(b.len());
                let mut values = Vec::with_capacity(entries);
                for x in &av {
                    for y in &bv {
                        let z = germline(x, y)?;
                        values.push(add(add(z[0], z[1])?, z[2])?);
                    }
                }
                cache.push(values);
            }
        }
        let prepare = |rs: &[Receptor]| -> Result<Vec<Vec<PreparedSequence>>> {
            rs.iter()
                .map(|r| {
                    r.chains
                        .iter()
                        .map(|c| PreparedSequence::new(&c.cdr3))
                        .collect()
                })
                .collect()
        };
        // Linear grouping storage <=32 bytes/target plus small constant scratch.
        // This is <= one eighth of the already reserved 256-byte/target scratch;
        // the unchanged half-budget admission leaves that additional reserve.
        let mut radius_order = Vec::new();
        let mut radius_groups = Vec::new();
        if config.optimized && config.radius.is_some() {
            radius_order = Vec::with_capacity(t.len());
            radius_order.extend(0..t.len());
            let compare = |a: usize, b: usize| {
                #[allow(clippy::needless_range_loop)]
                // Preserve the original indexing and refusal behavior.
                for c in 0..q[0].chains.len() {
                    let order = t_ids[c][a]
                        .cmp(&t_ids[c][b])
                        .then_with(|| t[a].chains[c].cdr3.len().cmp(&t[b].chains[c].cdr3.len()));
                    if !order.is_eq() {
                        return order;
                    }
                }
                std::cmp::Ordering::Equal
            };
            radius_order.sort_unstable_by(|a, b| compare(*a, *b));
            let mut start = 0;
            for end in 1..=t.len() {
                if end == t.len() || !compare(radius_order[start], radius_order[end]).is_eq() {
                    radius_groups.push((start, end));
                    start = end;
                    if radius_groups.len() >= t.len() / 2 {
                        radius_groups = Vec::new();
                        radius_order = Vec::new();
                        break;
                    }
                }
            }
        }
        let q_prepared = prepare(q)?;
        let t_prepared = prepare(t)?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .stack_size(2 * 1024 * 1024)
            .build()
            .map_err(|e| e.to_string())?;
        Ok(Self {
            q,
            t,
            config,
            cache,
            q_ids,
            t_ids,
            widths,
            q_prepared,
            t_prepared,
            pool,
            radius_order,
            radius_groups,
        })
    }
    fn score(&self, i: usize, j: usize) -> Result<u32> {
        if !self.config.optimized {
            return distance(&self.q[i], &self.t[j]);
        }
        // All encodings are validated before any nonnegative lower-bound rejection.
        if let Some(radius) = self.config.radius {
            if self.q[i].chains.len() == 2 {
                let mut lower = 0u64;
                for (c, (a, b)) in self.q[i].chains.iter().zip(&self.t[j].chains).enumerate() {
                    lower += u64::from(
                        self.cache[c][self.q_ids[c][i] * self.widths[c] + self.t_ids[c][j]],
                    );
                    lower += a.cdr3.len().abs_diff(b.cdr3.len()) as u64 * 12;
                }
                if lower > u64::from(radius) {
                    return Ok(u32::MAX);
                }
            }
        }
        let mut total = 0;
        for (c, (a, b)) in self.q[i].chains.iter().zip(&self.t[j].chains).enumerate() {
            let v = self.cache[c][self.q_ids[c][i] * self.widths[c] + self.t_ids[c][j]];
            total = add(total, v)?;
            let length = mul(a.cdr3.len().abs_diff(b.cdr3.len()) as u32, 12)?;
            if self
                .config
                .radius
                .is_some_and(|r| u64::from(total) + u64::from(length) > u64::from(r))
            {
                return Ok(u32::MAX);
            }
            total = add(
                total,
                mul(
                    self.q_prepared[i][c].score(&self.t_prepared[j][c], CDR3)?,
                    3,
                )?,
            )?;
        }
        Ok(total)
    }
    fn group_pruned(&self, i: usize, group: &(usize, usize)) -> bool {
        let Some(radius) = self.config.radius else {
            return false;
        };
        let j = self.radius_order[group.0];
        let mut lower = 0u64;
        for (c, (a, b)) in self.q[i].chains.iter().zip(&self.t[j].chains).enumerate() {
            lower += u64::from(self.cache[c][self.q_ids[c][i] * self.widths[c] + self.t_ids[c][j]]);
            lower += a.cdr3.len().abs_diff(b.cdr3.len()) as u64 * 12;
        }
        lower > u64::from(radius)
    }
    fn row(&self, i: usize) -> Result<Vec<u32>> {
        if self.config.optimized && self.config.radius.is_some() && !self.radius_groups.is_empty() {
            let skipped: usize = self
                .radius_groups
                .iter()
                .filter(|g| self.group_pruned(i, g))
                .map(|g| g.1 - g.0)
                .sum();
            if skipped >= self.t.len() / 2 {
                let mut row = vec![u32::MAX; self.t.len()];
                for group in &self.radius_groups {
                    if !self.group_pruned(i, group) {
                        for &j in &self.radius_order[group.0..group.1] {
                            row[j] = self.score(i, j)?;
                        }
                    }
                }
                return Ok(row);
            }
        }
        (0..self.t.len()).map(|j| self.score(i, j)).collect()
    }
    pub fn run<F>(&self, consume: F) -> Result<()>
    where
        F: FnMut(usize, &[Vec<u32>]) -> Result<()>,
    {
        self.run_profile(consume).map(|_| ())
    }
    pub fn run_profile<F>(&self, mut consume: F) -> Result<(f64, f64)>
    where
        F: FnMut(usize, &[Vec<u32>]) -> Result<()>,
    {
        let (mut compute, mut consume_seconds) = (0., 0.);
        for start in (0..self.q.len()).step_by(self.config.tile) {
            let end = start.saturating_add(self.config.tile).min(self.q.len());
            let clock = std::time::Instant::now();
            let rows: Result<Vec<Vec<u32>>> = self
                .pool
                .install(|| (start..end).into_par_iter().map(|i| self.row(i)).collect());
            let rows = rows?;
            compute += clock.elapsed().as_secs_f64();
            let clock = std::time::Instant::now();
            consume(start, &rows)?;
            consume_seconds += clock.elapsed().as_secs_f64();
        }
        Ok((compute, consume_seconds))
    }
}
fn serde_size(r: &Receptor) -> Result<usize> {
    let mut n = std::mem::size_of::<Receptor>();
    for part in std::iter::once(r.id.len())
        .chain(r.metadata.iter().flat_map(|(k, v)| [k.len(), v.len(), 96]))
        .chain(r.chains.iter().flat_map(|c| {
            [
                c.v.len(),
                c.cdr3.len(),
                c.loops[0].len(),
                c.loops[1].len(),
                c.loops[2].len(),
                256,
            ]
        }))
    {
        n = n.checked_add(part).ok_or("input size overflow")?;
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use clonodist_core::{encode, Chain};
    fn record() -> Receptor {
        Receptor {
            id: "x".into(),
            count: 1,
            metadata: BTreeMap::new(),
            chains: vec![Chain {
                v: "V*01".into(),
                loops: [vec![0; 10], vec![0; 10], vec![0; 6]],
                cdr3: encode("AAAACAAAAA", true).unwrap(),
            }],
        }
    }
    fn config(optimized: bool) -> Config {
        Config {
            threads: 1,
            tile: 2,
            memory: 128 * 1024 * 1024,
            optimized,
            radius: Some(0),
        }
    }
    #[test]
    fn empty_and_overflow_plans() {
        assert!(Engine::new(&[], &[], config(true)).is_err());
        let q = vec![record()];
        let mut c = config(true);
        c.tile = usize::MAX;
        assert!(Engine::new(&q, &q, c).is_err());
    }
    #[test]
    fn malformed_refused_before_pruning() {
        let a = vec![record()];
        let mut bad = record();
        bad.chains[0].cdr3[0] = 24;
        let b = vec![bad];
        for opt in [true, false] {
            assert!(Engine::new(&a, &b, config(opt)).is_err());
        }
    }
    #[test]
    fn inconsistent_cache_refused() {
        let a = vec![record()];
        let mut bad = record();
        bad.chains[0].loops[0][0] = 1;
        let b = vec![bad];
        for opt in [true, false] {
            assert!(Engine::new(&a, &b, config(opt)).is_err());
        }
    }
    #[test]
    fn native_zero_and_pruned_values() {
        let a = vec![record()];
        let e = Engine::new(&a, &a, config(true)).unwrap();
        e.run(|_, rows| {
            assert_eq!(rows, &[vec![0]]);
            Ok(())
        })
        .unwrap();
    }
    #[test]
    fn grouped_radius_matches_original_rows_and_scalar_edges() {
        for chain_count in [1, 2] {
            let mut targets = Vec::new();
            for v in 0..2 {
                for length in [5, 8, 40, 127] {
                    for variant in 0..8 {
                        let mut r = record();
                        r.id = format!("{v}-{length}-{variant}");
                        r.chains = (0..chain_count)
                            .map(|c| {
                                let mut cdr3 = vec![(variant % 5) as u8; length];
                                if variant == 7 {
                                    cdr3[length / 2] = 23;
                                }
                                Chain {
                                    v: format!("V-{v}-{c}"),
                                    loops: [vec![v as u8; 10], vec![v as u8; 10], vec![v as u8; 6]],
                                    cdr3,
                                }
                            })
                            .collect();
                        targets.push(r);
                    }
                }
            }
            let queries: Vec<_> = [0usize, 8, 16, 24, 32, 40, 48, 56]
                .iter()
                .map(|i| targets[*i].clone())
                .collect();
            for radius in [0u32, 1, 35, 36, 37, 100, 1000] {
                for (threads, tile) in [(1, 1), (2, 7), (4, 32), (8, 64)] {
                    let mut cfg = config(true);
                    cfg.radius = Some(radius);
                    cfg.threads = threads;
                    cfg.tile = tile;
                    let engine = Engine::new(&queries, &targets, cfg.clone()).unwrap();
                    assert!(!engine.radius_groups.is_empty());
                    let mut original = Engine::new(&queries, &targets, cfg).unwrap();
                    original.radius_groups.clear();
                    original.radius_order.clear();
                    let collect = |e: &Engine<'_>| {
                        let mut rows = Vec::new();
                        e.run(|_, block| {
                            rows.extend_from_slice(block);
                            Ok(())
                        })
                        .unwrap();
                        rows
                    };
                    let rows = collect(&engine);
                    assert_eq!(rows, collect(&original));
                    for (i, q) in queries.iter().enumerate() {
                        for (j, t) in targets.iter().enumerate() {
                            let score = clonodist_core::distance(q, t).unwrap();
                            assert_eq!(rows[i][j] <= radius, score <= radius);
                            if score <= radius {
                                assert_eq!(rows[i][j], score);
                            }
                        }
                    }
                }
            }
            let mut bad = targets.clone();
            bad.last_mut().unwrap().chains[0].cdr3[0] = 24;
            assert!(Engine::new(&queries, &bad, config(true)).is_err());
        }
    }
}
