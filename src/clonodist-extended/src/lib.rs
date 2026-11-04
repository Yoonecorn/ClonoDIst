//! Additive Extended contracts. All residue scoring uses the frozen Core kernel.
use clonodist_core::{add, encode, mul, region, Chain, Params, Receptor, Result, CDR3, GERMLINE};
use clonodist_distance::{Config, Engine};
use clonodist_germline::Reference;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};
pub mod meta;
pub mod mouse_trd;

pub fn checksum<T: Serialize>(v: &T) -> Result<String> {
    let value = serde_json::to_value(v).map_err(|e| e.to_string())?;
    let bytes = serde_json::to_vec(&value).map_err(|e| e.to_string())?;
    Ok(format!("{:x}", Sha256::digest(bytes)))
}
fn strict_params<'de, D: serde::Deserializer<'de>>(d: D) -> std::result::Result<Params, D::Error> {
    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    struct P {
        dist_weight: u32,
        gap_penalty: u32,
        ntrim: usize,
        ctrim: usize,
        fixed_gappos: bool,
    }
    let p = P::deserialize(d)?;
    Ok(Params {
        dist_weight: p.dist_weight,
        gap_penalty: p.gap_penalty,
        ntrim: p.ntrim,
        ctrim: p.ctrim,
        fixed_gappos: p.fixed_gappos,
    })
}
pub fn matrix_checksum() -> String {
    let mut bytes = Vec::new();
    for i in 0..24 {
        for j in 0..24 {
            bytes.extend_from_slice(&clonodist_core::cost(i, j).to_le_bytes());
        }
    }
    format!("{:x}", Sha256::digest(bytes))
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Region {
    pub name: String,
    pub weight: u32,
    #[serde(deserialize_with = "strict_params")]
    pub parameters: Params,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Profile {
    pub schema: String,
    pub id: String,
    pub status: String,
    pub substitution: String,
    pub substitution_sha256: String,
    pub reference_sha256: String,
    pub encoding: String,
    pub regions: Vec<Region>,
    pub chain_weights: Vec<u32>,
    pub chain_aggregation: String,
    pub missing_chain: String,
    pub normalization: String,
    pub evidence: String,
}
impl Profile {
    pub fn pinned(n: usize) -> Self {
        Self {
            schema: "clonodist.extended.profile.v1".into(),
            id: "tcrdist3-pinned-0e6efa3".into(),
            status: "reference-compatible".into(),
            substitution: "core-tcr-nb-24x24-v1".into(),
            substitution_sha256: matrix_checksum(),
            reference_sha256: clonodist_germline::DB_SHA256.into(),
            encoding: "canonical-CDR3; aligned-germline; Core-23-symbol-unknown-zero".into(),
            regions: ["cdr1", "cdr2", "pmhc", "cdr3"]
                .iter()
                .enumerate()
                .map(|(i, n)| Region {
                    name: n.to_string(),
                    weight: if i == 3 { 3 } else { 1 },
                    parameters: if i == 3 { CDR3 } else { GERMLINE },
                })
                .collect(),
            chain_weights: vec![1; n],
            chain_aggregation: "weighted_sum".into(),
            missing_chain: "reject".into(),
            normalization: "none-integer".into(),
            evidence:
                "pinned repertoire.py:565-682; local component/workflow differential required"
                    .into(),
        }
    }
    pub fn cdr3(n: usize) -> Self {
        let mut p = Self::pinned(n);
        p.id = "cdr3-focused-v1".into();
        p.status = "experimental".into();
        p.evidence =
            "Ablation hypothesis: remove germline contributions; no predictive-validity claim"
                .into();
        for r in &mut p.regions[..3] {
            r.weight = 0;
        }
        p
    }
    pub fn validate(&self, n: usize, experimental: bool) -> Result<()> {
        if self.substitution_sha256 != matrix_checksum()
            || self.reference_sha256 != clonodist_germline::DB_SHA256
            || self.encoding != "canonical-CDR3; aligned-germline; Core-23-symbol-unknown-zero"
        {
            return Err("matrix/reference/encoding mismatch".into());
        }
        if self.schema != "clonodist.extended.profile.v1"
            || self.id.is_empty()
            || self.id.len() > 128
            || self.evidence.is_empty()
            || self.regions.len() != 4
            || !(1..=2).contains(&n)
            || self.chain_weights.len() != n
            || self.substitution != "core-tcr-nb-24x24-v1"
            || self.chain_aggregation != "weighted_sum"
            || self.missing_chain != "reject"
            || self.normalization != "none-integer"
        {
            return Err("unsupported profile contract".into());
        }
        for (r, name) in self.regions.iter().zip(["cdr1", "cdr2", "pmhc", "cdr3"]) {
            let p = r.parameters;
            if r.name != name
                || r.weight > 16
                || p.dist_weight > 16
                || p.gap_penalty > 16
                || p.ntrim > 127
                || p.ctrim > 127
            {
                return Err("invalid region or parameter outside admitted bounds".into());
            }
        }
        if self.regions.iter().all(|r| {
            r.weight == 0
                || (r.parameters.gap_penalty == 0
                    && (r.parameters.dist_weight == 0
                        || r.parameters.ntrim + r.parameters.ctrim >= 127))
        }) || self.chain_weights.iter().any(|w| *w == 0 || *w > 16)
        {
            return Err("zero metric or invalid chain weight".into());
        }
        if self.status == "reference-compatible" {
            if checksum(self)? != checksum(&Self::pinned(n))? {
                return Err("reference profile is immutable".into());
            }
        } else if self.status != "experimental" || !experimental {
            return Err("custom profiles require status experimental and --experimental".into());
        }
        Ok(())
    }
    pub fn is_pinned(&self) -> bool {
        self.status == "reference-compatible"
    }
    pub fn score(&self, a: &Receptor, b: &Receptor) -> Result<u32> {
        if a.chains.len() != self.chain_weights.len() || b.chains.len() != a.chains.len() {
            return Err("chain configuration mismatch".into());
        }
        let mut total = 0;
        for ((a, b), w) in a.chains.iter().zip(&b.chains).zip(&self.chain_weights) {
            let mut sum = 0;
            for (i, r) in self.regions.iter().enumerate() {
                if r.weight == 0 {
                    continue;
                }
                let (x, y) = if i == 3 {
                    (&a.cdr3, &b.cdr3)
                } else {
                    (&a.loops[i], &b.loops[i])
                };
                sum = add(sum, mul(region(x, y, r.parameters)?, r.weight)?)?;
            }
            total = add(total, mul(sum, *w)?)?;
        }
        Ok(total)
    }
}
pub fn chain_config(chains: &str, organism: &str) -> Result<Vec<String>> {
    let allowed = matches!(
        chains,
        "alpha" | "beta" | "alpha,beta" | "gamma" | "delta" | "gamma,delta"
    );
    if !allowed
        || !matches!(organism, "human" | "mouse")
        || (organism == "mouse" && chains.contains("delta"))
    {
        return Err("unsupported organism/chain qualification scope".into());
    }
    Ok(chains.split(',').map(str::to_owned).collect())
}
pub fn read(
    path: &Path,
    db: &Reference,
    org: &str,
    chains: &[String],
    memory: usize,
) -> Result<Vec<Receptor>> {
    read_inner(path, db, org, chains, memory, None)
}
pub fn read_with_mouse_trd(
    path: &Path,
    db: &Reference,
    org: &str,
    chains: &[String],
    memory: usize,
    supplement: &mouse_trd::MouseTrdReference,
) -> Result<Vec<Receptor>> {
    if org != "mouse" || !(chains == ["delta"] || chains == ["gamma", "delta"]) {
        return Err("mouse TRD reference requires mouse delta or gamma,delta".into());
    }
    read_inner(path, db, org, chains, memory, Some(supplement))
}
fn read_inner(
    path: &Path,
    db: &Reference,
    org: &str,
    chains: &[String],
    memory: usize,
    supplement: Option<&mouse_trd::MouseTrdReference>,
) -> Result<Vec<Receptor>> {
    use std::io::Read;
    let meta = std::fs::metadata(path).map_err(|e| e.to_string())?;
    if !meta.is_file() || meta.len() > memory as u64 / 64 {
        return Err("input exceeds conservative admission".into());
    }
    let f = std::fs::File::open(path).map_err(|e| e.to_string())?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(f.take(meta.len() + 1));
    let headers = rdr.headers().map_err(|e| e.to_string())?.clone();
    let mut seen = BTreeSet::new();
    if headers
        .iter()
        .any(|h| h.is_empty() || !seen.insert(h.to_owned()))
    {
        return Err("empty/duplicate columns".into());
    }
    let mut records = vec![];
    let mut ids = BTreeSet::new();
    let mut planned = 0usize;
    let mut count_sum = 0u64;
    for row in rdr.records() {
        let row = row.map_err(|e| e.to_string())?;
        planned = planned
            .checked_add(
                row.as_slice()
                    .len()
                    .checked_mul(64)
                    .and_then(|v| v.checked_add(16384))
                    .ok_or("input plan overflow")?,
            )
            .ok_or("input plan overflow")?;
        if planned > memory / 4 {
            return Err("normalized input exceeds admission".into());
        }
        let get = |key: &str| -> Result<&str> {
            let i = headers
                .iter()
                .position(|h| h == key)
                .ok_or_else(|| format!("missing column {key}"))?;
            if row[i].is_empty() {
                Err(format!("empty {key}"))
            } else {
                Ok(&row[i])
            }
        };
        let id = get("id")?.to_owned();
        if id.chars().any(|c| c.is_control() || c == '"') || !ids.insert(id.clone()) {
            return Err("invalid/duplicate id".into());
        }
        let count = if headers.iter().any(|h| h == "count") {
            get("count")?.parse::<u64>().map_err(|_| "invalid count")?
        } else {
            1
        };
        if count == 0 {
            return Err("zero count".into());
        }
        count_sum = count_sum.checked_add(count).ok_or("count overflow")?;
        let mut resolved = vec![];
        for ch in chains {
            let (letter, label) = match ch.as_str() {
                "alpha" => ("a", "A"),
                "beta" => ("b", "B"),
                "gamma" => ("g", "G"),
                "delta" => ("d", "D"),
                _ => return Err("unsupported chain".into()),
            };
            let v = get(&format!("v_{letter}_gene"))?;
            let j = get(&format!("j_{letter}_gene"))?;
            if org == "mouse" && ch == "delta" {
                if let Some(reference) = supplement {
                    resolved.push(reference.resolve(v, j, get(&format!("cdr3_{letter}_aa"))?)?);
                    continue;
                }
            }
            let vg = db.gene(org, v)?;
            let jg = db.gene(org, j)?;
            if !v.contains('*')
                || !j.contains('*')
                || vg.chain != label
                || vg.region != "V"
                || jg.chain != label
                || jg.region != "J"
                || vg.loops.len() < 3
            {
                return Err("exact V/J locus, region or allele mismatch; shared alpha/delta alleles unsupported".into());
            }
            resolved.push(Chain {
                v: v.into(),
                loops: [
                    encode(&vg.loops[0], false)?,
                    encode(&vg.loops[1], false)?,
                    encode(&vg.loops[2], false)?,
                ],
                cdr3: encode(get(&format!("cdr3_{letter}_aa"))?, true)?,
            });
        }
        let metadata: BTreeMap<String, String> = headers
            .iter()
            .zip(row.iter())
            .filter(|(h, _)| !["id", "count"].contains(h))
            .map(|(h, v)| (h.into(), v.into()))
            .collect();
        if metadata.values().any(|v| v.is_empty()) {
            return Err("empty metadata".into());
        }
        records.push(Receptor {
            id,
            count,
            metadata,
            chains: resolved,
        });
    }
    if rdr.position().byte() != meta.len() || records.is_empty() {
        return Err("empty or changed input".into());
    }
    Ok(records)
}
/// Profile orchestration cache; residue/gap scoring remains Core PreparedSequence.
struct PreparedProfile {
    cache: Vec<Vec<u32>>,
    q_ids: Vec<Vec<usize>>,
    t_ids: Vec<Vec<usize>>,
    widths: Vec<usize>,
    q: Vec<Vec<clonodist_core::PreparedSequence>>,
    t: Vec<Vec<clonodist_core::PreparedSequence>>,
}
impl PreparedProfile {
    fn new(q: &[Receptor], t: &[Receptor], p: &Profile, memory: usize) -> Result<Self> {
        let mut out = Self {
            cache: vec![],
            q_ids: vec![],
            t_ids: vec![],
            widths: vec![],
            q: vec![],
            t: vec![],
        };
        let mut bytes = 0usize;
        for c in 0..p.chain_weights.len() {
            let av: BTreeMap<_, _> = q.iter().map(|r| (&r.chains[c].v, &r.chains[c])).collect();
            let bv: BTreeMap<_, _> = t.iter().map(|r| (&r.chains[c].v, &r.chains[c])).collect();
            let entries = av
                .len()
                .checked_mul(bv.len())
                .ok_or("profile cache overflow")?;
            bytes = bytes
                .checked_add(entries.checked_mul(128).ok_or("profile cache overflow")?)
                .ok_or("profile cache overflow")?;
            if bytes > memory / 16 {
                return Err("profile germline cache exceeds admission; use --scalar".into());
            }
            let ai: BTreeMap<_, _> = av
                .keys()
                .enumerate()
                .map(|(i, v)| ((*v).clone(), i))
                .collect();
            let bi: BTreeMap<_, _> = bv
                .keys()
                .enumerate()
                .map(|(i, v)| ((*v).clone(), i))
                .collect();
            out.q_ids
                .push(q.iter().map(|r| ai[&r.chains[c].v]).collect());
            out.t_ids
                .push(t.iter().map(|r| bi[&r.chains[c].v]).collect());
            out.widths.push(bv.len());
            let mut scores = Vec::with_capacity(entries);
            for a in av.values() {
                for b in bv.values() {
                    let mut sum = 0;
                    for (i, r) in p.regions[..3].iter().enumerate() {
                        if r.weight > 0 {
                            sum = add(
                                sum,
                                mul(region(&a.loops[i], &b.loops[i], r.parameters)?, r.weight)?,
                            )?;
                        }
                    }
                    scores.push(mul(sum, p.chain_weights[c])?);
                }
            }
            out.cache.push(scores);
        }
        let prepare = |rs: &[Receptor]| -> Result<Vec<Vec<clonodist_core::PreparedSequence>>> {
            rs.iter()
                .map(|r| {
                    r.chains
                        .iter()
                        .map(|c| clonodist_core::PreparedSequence::new(&c.cdr3))
                        .collect()
                })
                .collect()
        };
        out.q = prepare(q)?;
        out.t = prepare(t)?;
        Ok(out)
    }
    fn score(&self, i: usize, j: usize, p: &Profile) -> Result<u32> {
        let mut total = 0;
        for c in 0..p.chain_weights.len() {
            total = add(
                total,
                self.cache[c][self.q_ids[c][i] * self.widths[c] + self.t_ids[c][j]],
            )?;
            let r = &p.regions[3];
            if r.weight > 0 {
                total = add(
                    total,
                    mul(
                        mul(self.q[i][c].score(&self.t[j][c], r.parameters)?, r.weight)?,
                        p.chain_weights[c],
                    )?,
                )?;
            }
        }
        Ok(total)
    }
}
/// Default profiles dispatch through the unchanged qualified engine. Custom profiles
/// orchestrate its same checked regional kernel; no second alignment implementation.
pub fn rows<F>(
    q: &[Receptor],
    t: &[Receptor],
    p: &Profile,
    c: Config,
    mut consume: F,
) -> Result<(f64, f64)>
where
    F: FnMut(usize, &[Vec<u32>]) -> Result<()>,
{
    p.validate(q.first().ok_or("empty query")?.chains.len(), true)?;
    if p.is_pinned() {
        return Engine::new(q, t, c)?.run_profile(consume);
    }
    use rayon::prelude::*;
    // Invoke unchanged engine admission and validation before custom scheduling.
    let admission = Engine::new(
        q,
        t,
        Config {
            optimized: false,
            radius: None,
            ..c.clone()
        },
    )?;
    drop(admission);
    let prepared = if c.optimized {
        Some(PreparedProfile::new(q, t, p, c.memory)?)
    } else {
        None
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(c.threads)
        .stack_size(2 * 1024 * 1024)
        .build()
        .map_err(|e| e.to_string())?;
    let (mut compute, mut output) = (0., 0.);
    for start in (0..q.len()).step_by(c.tile) {
        let timer = std::time::Instant::now();
        let end = start.saturating_add(c.tile).min(q.len());
        let block: Result<Vec<Vec<u32>>> = pool.install(|| {
            (start..end)
                .into_par_iter()
                .map(|i| {
                    t.iter()
                        .enumerate()
                        .map(|(j, b)| {
                            if let Some(prepared) = &prepared {
                                prepared.score(i, j, p)
                            } else {
                                p.score(&q[i], b)
                            }
                        })
                        .collect()
                })
                .collect()
        });
        let block = block?;
        compute += timer.elapsed().as_secs_f64();
        let timer = std::time::Instant::now();
        consume(start, &block)?;
        output += timer.elapsed().as_secs_f64();
    }
    Ok((compute, output))
}
pub fn signature(r: &Receptor) -> Result<String> {
    checksum(&r.chains)
}
