//! Prespecified empirical neighborhoods; discovery and untouched qualification are separate.
use super::*;
// Decimal strings preserve exact f64 round trips independently of serde_json's
// optional float parser feature, without changing any Core dependency features.
mod decimal {
    use serde::{Deserialize, Deserializer, Serializer};
    pub fn serialize<S: Serializer>(v: &f64, s: S) -> std::result::Result<S::Ok, S::Error> {
        s.serialize_str(&v.to_string())
    }
    pub fn deserialize<'de, D: Deserializer<'de>>(d: D) -> std::result::Result<f64, D::Error> {
        #[derive(Deserialize)]
        #[serde(untagged)]
        enum N {
            Text(String),
            Number(f64),
        }
        match N::deserialize(d)? {
            N::Text(s) => s.parse().map_err(serde::de::Error::custom),
            N::Number(v) => Ok(v),
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Rules {
    pub radii: Vec<u32>,
    #[serde(with = "decimal")]
    pub background_cap: f64,
    pub minimum_hits: usize,
    #[serde(with = "decimal")]
    pub maximum_jaccard: f64,
    pub same_vj: bool,
    pub motifs: Vec<String>,
    pub resamples: usize,
    pub seed: u64,
}
impl Rules {
    pub fn validate(&self, n: usize) -> Result<()> {
        if self.radii.is_empty()
            || self.radii.len() > 1024
            || self.radii.windows(2).any(|w| w[0] >= w[1])
            || !self.background_cap.is_finite()
            || !(0.0..=1.0).contains(&self.background_cap)
            || !self.maximum_jaccard.is_finite()
            || !(0.0..=1.0).contains(&self.maximum_jaccard)
            || self.minimum_hits == 0
            || self.resamples < 2
            || self.resamples > 1000
            || (!self.motifs.is_empty() && self.motifs.len() != n)
        {
            return Err("invalid discovery rules".into());
        }
        if self.motifs.iter().any(|m| {
            m.is_empty()
                || m.len() > 127
                || m.bytes()
                    .any(|b| b != b'.' && !clonodist_core::ALPHABET[..20].contains(&b))
        }) {
            return Err("motifs must be full-length canonical AA/dot patterns per chain".into());
        }
        Ok(())
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Candidate {
    pub center: Receptor,
    pub radius: Option<u32>,
    pub discovery_hits: usize,
    pub background_hits: usize,
    #[serde(with = "decimal")]
    pub enrichment: f64,
    #[serde(with = "decimal")]
    pub stability: f64,
    pub radius_range: Option<[u32; 2]>,
    #[serde(with = "decimal")]
    pub membership_jaccard: f64,
    pub reason: String,
    pub heldout_hits: usize,
    pub heldout_background_hits: usize,
    #[serde(with = "decimal")]
    pub adjusted_p: f64,
    pub qualified: bool,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Model {
    pub schema: String,
    pub phase: String,
    pub organism: String,
    pub chains: Vec<String>,
    pub profile: Profile,
    pub profile_sha256: String,
    pub rules: Rules,
    pub candidates: Vec<Candidate>,
    pub discovery_subjects: Vec<String>,
    pub discovery_cohorts: Vec<String>,
    pub discovery_signatures: Vec<String>,
    pub dataset_hashes: BTreeMap<String, String>,
    pub qualification: serde_json::Value,
    pub seal: String,
}
impl Model {
    pub fn reseal(&mut self) -> Result<()> {
        self.seal.clear();
        self.seal = checksum(self)?;
        Ok(())
    }
    pub fn verify(&self) -> Result<()> {
        let mut m = self.clone();
        let want = m.seal.clone();
        m.seal.clear();
        if checksum(&m)? != want || checksum(&m.profile)? != m.profile_sha256 {
            return Err("model seal/profile mismatch".into());
        }
        m.profile.validate(m.chains.len(), true)?;
        m.rules.validate(m.chains.len())?;
        chain_config(&m.chains.join(","), &m.organism)?;
        if !["discovered", "qualified"].contains(&m.phase.as_str()) {
            return Err("invalid model phase".into());
        }
        let mut ids = BTreeSet::new();
        for c in &m.candidates {
            if !ids.insert(&c.center.id)
                || c.center.chains.len() != m.chains.len()
                || c.radius.is_some_and(|r| !m.rules.radii.contains(&r))
                || c.radius_range.is_some_and(|[a, b]| {
                    a > b || !m.rules.radii.contains(&a) || !m.rules.radii.contains(&b)
                })
                || ![
                    c.enrichment,
                    c.stability,
                    c.membership_jaccard,
                    c.adjusted_p,
                ]
                .iter()
                .all(|v| v.is_finite())
                || !(0.0..=1.0).contains(&c.stability)
                || !(0.0..=1.0).contains(&c.membership_jaccard)
                || !(0.0..=1.0).contains(&c.adjusted_p)
                || c.enrichment < 0.
                || (c.qualified
                    && (m.phase != "qualified"
                        || c.reason != "discovery_selected"
                        || c.radius.is_none()))
            {
                return Err("invalid candidate structure".into());
            }
            if c.reason != "discovery_selected"
                && c.reason != "no_radius_satisfies_background_and_support"
                && !c.reason.starts_with("redundant_with:")
            {
                return Err("invalid candidate reason".into());
            }
        }
        if m.schema != "clonodist.extended.meta.v1" {
            return Err("unsupported model schema".into());
        }
        Ok(())
    }
}
pub fn fields(rs: &[Receptor], key: &str) -> Result<BTreeSet<String>> {
    rs.iter()
        .map(|r| {
            r.metadata
                .get(key)
                .filter(|v| !matches!(v.as_str(), "NA" | "unknown" | "unreported" | "None"))
                .cloned()
                .ok_or_else(|| format!("meta workflows require known {key}"))
        })
        .collect()
}
fn validate_units(fg: &[Receptor], bg: &[Receptor]) -> Result<()> {
    let mut seen = BTreeSet::new();
    for r in fg.iter().chain(bg) {
        if !seen.insert(signature(r)?) {
            return Err("duplicate/overlapping inference receptor signatures".into());
        }
    }
    Ok(())
}
fn union_field(a: &[Receptor], b: &[Receptor], key: &str) -> Result<Vec<String>> {
    Ok(fields(a, key)?.union(&fields(b, key)?).cloned().collect())
}
/// V/J identity is restricted to selected validated receptor chains.
pub fn gene_identity(center: &Receptor, r: &Receptor) -> bool {
    center.chains.iter().all(|c| {
        let letter = match c.v.get(..4) {
            Some("TRAV") => "a",
            Some("TRBV") => "b",
            Some("TRGV") => "g",
            Some("TRDV") => "d",
            _ => return false,
        };
        [format!("v_{letter}_gene"), format!("j_{letter}_gene")]
            .iter()
            .all(|key| {
                center.metadata.contains_key(key) && center.metadata.get(key) == r.metadata.get(key)
            })
    })
}
pub fn constraints(center: &Receptor, r: &Receptor, rule: &Rules) -> (bool, bool) {
    let gene = !rule.same_vj || gene_identity(center, r);
    let motif = rule.motifs.is_empty()
        || rule.motifs.iter().zip(&r.chains).all(|(m, c)| {
            m.len() == c.cdr3.len()
                && m.bytes().zip(&c.cdr3).all(|(x, y)| {
                    x == b'.' || Some(&x) == clonodist_core::ALPHABET.get(*y as usize)
                })
        });
    (gene, motif)
}
pub fn matching(
    center: &Receptor,
    rs: &[Receptor],
    distances: &[u32],
    radius: u32,
    rules: &Rules,
) -> BTreeSet<usize> {
    rs.iter()
        .enumerate()
        .filter(|(i, r)| {
            let (g, m) = constraints(center, r, rules);
            distances[*i] <= radius && g && m
        })
        .map(|(i, _)| i)
        .collect()
}
fn enrichment(f: usize, n: usize, b: usize, m: usize) -> f64 {
    ((f + 1) as f64 / (n + 2) as f64) / ((b + 1) as f64 / (m + 2) as f64)
}
fn sample(subject: &str, seed: u64, rep: usize) -> bool {
    let mut h = Sha256::new();
    h.update(seed.to_le_bytes());
    h.update((rep as u64).to_le_bytes());
    h.update(subject.as_bytes());
    h.finalize()[0] < 205
}
#[allow(clippy::too_many_arguments)] // Preserve the existing public API.
pub fn discover(
    fg: &[Receptor],
    bg: &[Receptor],
    profile: Profile,
    rules: Rules,
    org: String,
    chains: Vec<String>,
    config: Config,
    hashes: BTreeMap<String, String>,
) -> Result<Model> {
    rules.validate(chains.len())?;
    validate_units(fg, bg)?;
    let subjects = union_field(fg, bg, "subject")?;
    let cohorts = union_field(fg, bg, "cohort")?;
    let signatures = fg
        .iter()
        .chain(bg)
        .map(signature)
        .collect::<Result<BTreeSet<_>>>()?;
    let mut candidates = vec![];
    let mut memberships = vec![];
    // One center row at a time bounds distance scratch independently of center count.
    for center in fg {
        let mut a = vec![];
        let mut b = vec![];
        rows(
            std::slice::from_ref(center),
            fg,
            &profile,
            Config {
                radius: None,
                ..config.clone()
            },
            |_, r| {
                a = r[0].clone();
                Ok(())
            },
        )?;
        rows(
            std::slice::from_ref(center),
            bg,
            &profile,
            Config {
                radius: None,
                ..config.clone()
            },
            |_, r| {
                b = r[0].clone();
                Ok(())
            },
        )?;
        let mut selected = None;
        for &radius in &rules.radii {
            let f = matching(center, fg, &a, radius, &rules);
            let g = matching(center, bg, &b, radius, &rules);
            if f.len() >= rules.minimum_hits
                && g.len() as f64 / bg.len() as f64 <= rules.background_cap
            {
                selected = Some((radius, f, g));
            }
        }
        let mut candidate = Candidate {
            center: center.clone(),
            radius: None,
            discovery_hits: 0,
            background_hits: 0,
            enrichment: 0.,
            stability: 0.,
            radius_range: None,
            membership_jaccard: 0.,
            reason: "no_radius_satisfies_background_and_support".into(),
            heldout_hits: 0,
            heldout_background_hits: 0,
            adjusted_p: 1.,
            qualified: false,
        };
        let mut membership = BTreeSet::new();
        if let Some((radius, f, g)) = selected {
            candidate.radius = Some(radius);
            candidate.discovery_hits = f.len();
            candidate.background_hits = g.len();
            candidate.enrichment = enrichment(f.len(), fg.len(), g.len(), bg.len());
            candidate.reason = "discovery_selected".into();
            // Refit radius on deterministic 80% donor-block subsamples. A missing center counts as unstable.
            let mut stable = 0;
            let mut sampled_radii = vec![];
            let mut jaccard = 0.;
            for rep in 0..rules.resamples {
                if !sample(&center.metadata["subject"], rules.seed, rep) {
                    continue;
                }
                let fi: BTreeSet<_> = fg
                    .iter()
                    .enumerate()
                    .filter(|(_, r)| sample(&r.metadata["subject"], rules.seed, rep))
                    .map(|(i, _)| i)
                    .collect();
                let bi: BTreeSet<_> = bg
                    .iter()
                    .enumerate()
                    .filter(|(_, r)| sample(&r.metadata["subject"], rules.seed, rep))
                    .map(|(i, _)| i)
                    .collect();
                if bi.is_empty() {
                    continue;
                }
                let mut chosen = None;
                for &r in &rules.radii {
                    let x = matching(center, fg, &a, r, &rules);
                    let y = matching(center, bg, &b, r, &rules);
                    if x.intersection(&fi).count() >= rules.minimum_hits
                        && y.intersection(&bi).count() as f64 / bi.len() as f64
                            <= rules.background_cap
                    {
                        chosen = Some(r);
                    }
                }
                if chosen == Some(radius) {
                    stable += 1;
                }
                if let Some(r) = chosen {
                    sampled_radii.push(r);
                    let sampled = matching(center, fg, &a, r, &rules);
                    let expected: BTreeSet<_> = f.intersection(&fi).copied().collect();
                    let observed: BTreeSet<_> = sampled.intersection(&fi).copied().collect();
                    let union = expected.union(&observed).count();
                    if union > 0 {
                        jaccard += expected.intersection(&observed).count() as f64 / union as f64;
                    }
                }
            }
            candidate.stability = stable as f64 / rules.resamples as f64;
            candidate.membership_jaccard = jaccard / rules.resamples as f64;
            if let (Some(a), Some(b)) = (sampled_radii.iter().min(), sampled_radii.iter().max()) {
                candidate.radius_range = Some([*a, *b]);
            }
            membership = f;
        }
        candidates.push(candidate);
        memberships.push(membership);
        // Stored memberships are bounded explicitly; reject before excessive growth.
        let used: usize = memberships.iter().map(|s| s.len().saturating_mul(64)).sum();
        if used > config.memory / 8 {
            return Err("discovery membership admission exceeded".into());
        }
    }
    let mut order: Vec<_> = (0..candidates.len()).collect();
    order.sort_by(|&a, &b| {
        candidates[b]
            .enrichment
            .total_cmp(&candidates[a].enrichment)
            .then_with(|| candidates[a].center.id.cmp(&candidates[b].center.id))
    });
    let mut retained: Vec<usize> = vec![];
    for i in order {
        if candidates[i].radius.is_none() {
            continue;
        }
        let redundant = retained.iter().find(|&&j| {
            let intersection = memberships[i].intersection(&memberships[j]).count();
            let union = memberships[i].union(&memberships[j]).count();
            intersection as f64 / union as f64 > rules.maximum_jaccard
        });
        if let Some(&j) = redundant {
            candidates[i].reason = format!("redundant_with:{}", candidates[j].center.id);
        } else {
            retained.push(i);
        }
    }
    let mut model = Model {
        schema: "clonodist.extended.meta.v1".into(),
        phase: "discovered".into(),
        organism: org,
        chains,
        profile_sha256: checksum(&profile)?,
        profile,
        rules,
        candidates,
        discovery_subjects: subjects,
        discovery_cohorts: cohorts,
        discovery_signatures: signatures.into_iter().collect(),
        dataset_hashes: hashes
            .into_iter()
            .map(|(k, v)| (format!("discovery:{k}"), v))
            .collect(),
        qualification: serde_json::Value::Null,
        seal: String::new(),
    };
    model.reseal()?;
    Ok(model)
}
/// One-sided Fisher exact tail, computed with log-factorials and log-sum-exp.
pub fn fisher(a: usize, n: usize, b: usize, m: usize) -> f64 {
    let total = n + m;
    let hits = a + b;
    let mut lf = vec![0.; total + 1];
    for i in 1..=total {
        lf[i] = lf[i - 1] + (i as f64).ln();
    }
    let choose = |x: usize, y: usize| lf[x] - lf[y] - lf[x - y];
    let mut logs = vec![];
    for x in a..=n.min(hits) {
        if hits - x <= m {
            logs.push(choose(n, x) + choose(m, hits - x) - choose(total, hits));
        }
    }
    if logs.is_empty() {
        return 1.;
    }
    let peak = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    (peak.exp() * logs.iter().map(|v| (*v - peak).exp()).sum::<f64>()).min(1.)
}
pub fn qualify(
    mut model: Model,
    fg: &[Receptor],
    bg: &[Receptor],
    config: Config,
    hashes: BTreeMap<String, String>,
) -> Result<Model> {
    model.verify()?;
    if model.phase != "discovered" {
        return Err("qualification requires discovered artifact".into());
    }
    validate_units(fg, bg)?;
    let subjects = union_field(fg, bg, "subject")?;
    let cohorts = union_field(fg, bg, "cohort")?;
    if subjects
        .iter()
        .any(|s| model.discovery_subjects.contains(s))
        || cohorts.iter().any(|s| model.discovery_cohorts.contains(s))
    {
        return Err("held-out subject/cohort leakage".into());
    }
    if fields(fg, "subject")?.len() < 2 || fields(bg, "subject")?.len() < 2 {
        return Err("held-out evaluation requires two foreground and background subjects".into());
    }
    let old: BTreeSet<_> = model.discovery_signatures.iter().collect();
    let f: Result<BTreeSet<_>> = fg.iter().map(signature).collect();
    let f = f?;
    let b: Result<BTreeSet<_>> = bg.iter().map(signature).collect();
    let b = b?;
    if f.iter().chain(&b).any(|s| old.contains(s)) || !f.is_disjoint(&b) {
        return Err("held-out exact receptor leakage".into());
    }
    let tests = model
        .candidates
        .iter()
        .filter(|c| c.reason == "discovery_selected")
        .count()
        .max(1);
    let mut heldout_reports = vec![];
    for c in &mut model.candidates {
        if c.reason != "discovery_selected" {
            continue;
        }
        let radius = c.radius.ok_or("selected candidate lacks radius")?;
        let mut a = vec![];
        let mut b = vec![];
        rows(
            std::slice::from_ref(&c.center),
            fg,
            &model.profile,
            Config {
                radius: None,
                ..config.clone()
            },
            |_, r| {
                a = r[0].clone();
                Ok(())
            },
        )?;
        rows(
            std::slice::from_ref(&c.center),
            bg,
            &model.profile,
            Config {
                radius: None,
                ..config.clone()
            },
            |_, r| {
                b = r[0].clone();
                Ok(())
            },
        )?;
        let matched = matching(&c.center, fg, &a, radius, &model.rules);
        let donors: BTreeSet<_> = matched
            .iter()
            .map(|i| &fg[*i].metadata["subject"])
            .collect();
        c.heldout_hits = matched.len();
        c.heldout_background_hits = matching(&c.center, bg, &b, radius, &model.rules).len();
        c.adjusted_p = (fisher(
            c.heldout_hits,
            fg.len(),
            c.heldout_background_hits,
            bg.len(),
        ) * tests as f64)
            .min(1.);
        c.qualified = c.heldout_hits >= model.rules.minimum_hits
            && c.heldout_background_hits as f64 / bg.len() as f64 <= model.rules.background_cap
            && (c.heldout_hits as f64 / fg.len() as f64
                > c.heldout_background_hits as f64 / bg.len() as f64)
            && donors.len() >= 2
            && c.stability >= 0.6;
        heldout_reports.push(serde_json::json!({"center":c.center.id,"supporting_subjects":donors.len(),"foreground_prevalence":(c.heldout_hits as f64/fg.len()as f64).to_string(),"background_prevalence":(c.heldout_background_hits as f64/bg.len()as f64).to_string(),"smoothed_prevalence_ratio":enrichment(c.heldout_hits,fg.len(),c.heldout_background_hits,bg.len()).to_string()}));
    }
    model.phase = "qualified".into();
    model.dataset_hashes.extend(
        hashes
            .into_iter()
            .map(|(k, v)| (format!("qualification:{k}"), v)),
    );
    model.qualification = serde_json::json!({"heldout_candidates":heldout_reports,"foreground_n":fg.len(),"background_n":bg.len(),"subjects":subjects,"cohorts":cohorts,"tests":tests,"statistic":"EXPLORATORY one-sided Fisher, Bonferroni across retained discovery centers; receptor independence unproved; not a qualification gate","stability_min":0.6,"exploratory_independent_receptor_binomial_zero_event_95_bound":(1.-0.05f64.powf(1./bg.len()as f64)).to_string(),"qualification_gate":"operational support/background cap, positive raw foreground enrichment, >=2 heldout supporting donors, stability>=0.6; no inferential p-value claim","interpretation":"empirical enrichment against supplied background; receptor independence unproved; no antigen specificity or clinical claim"});
    model.reseal()?;
    Ok(model)
}
