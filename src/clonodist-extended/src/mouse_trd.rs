//! Explicit current-IMGT mouse delta context. Never changes the base reference.
use clonodist_core::{encode, Chain, Result, ALPHABET};
use clonodist_germline::Reference;
use serde::Deserialize;
use sha2::{Digest, Sha256};
use std::{collections::BTreeMap, io::Read, path::Path};
pub const REFERENCE_SHA256: &str =
    "482280a9dfacda0d45acf5243b68e088e7649819a05bc24b90a21f776536f649";
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Source {
    url: String,
    sha256: String,
    bytes: usize,
    license: String,
    retrieved_utc: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Entry {
    id: String,
    kind: String,
    regions: [String; 3],
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Excluded {
    id: String,
    reason: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Contract {
    schema: String,
    organism: String,
    version: String,
    base_sha256: String,
    sources: Vec<Source>,
    coordinate_policy: String,
    intervals: Vec<[usize; 2]>,
    alias_policy: String,
    entries: Vec<Entry>,
    j_alleles: Vec<String>,
    excluded: Vec<Excluded>,
}
pub struct MouseTrdReference {
    entries: BTreeMap<String, [String; 3]>,
    j_alleles: Vec<String>,
}
impl MouseTrdReference {
    pub fn load(path: &Path, base: &Reference) -> Result<Self> {
        let size = std::fs::metadata(path).map_err(|e| e.to_string())?;
        if !size.is_file() || size.len() > 65536 {
            return Err("mouse TRD reference must be a bounded regular file".into());
        }
        let mut bytes = vec![];
        std::fs::File::open(path)
            .map_err(|e| e.to_string())?
            .take(65537)
            .read_to_end(&mut bytes)
            .map_err(|e| e.to_string())?;
        if format!("{:x}", Sha256::digest(&bytes)) != REFERENCE_SHA256 {
            return Err("mouse TRD reference SHA-256 mismatch".into());
        }
        let c: Contract = serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
        if c.schema != "clonodist.extended.mouse-trd-reference.v1"
            || c.organism != "mouse"
            || c.base_sha256 != base.sha256
            || c.version != "IMGT-VQUEST-content-snapshot-2026-10-29"
            || c.intervals != [[28, 39], [57, 66], [82, 88]]
            || c.entries.len() != 51
            || c.sources.len() != 3
            || c.excluded.len() != 2
            || c.coordinate_policy.is_empty()
            || c.alias_policy.is_empty()
        {
            return Err("mouse TRD reference contract mismatch".into());
        }
        for s in &c.sources {
            if s.url.is_empty()
                || s.sha256.len() != 64
                || s.bytes == 0
                || s.license.is_empty()
                || s.retrieved_utc.is_empty()
            {
                return Err("incomplete reference provenance".into());
            }
        }
        for x in &c.excluded {
            if x.id.is_empty() || x.reason.is_empty() {
                return Err("invalid exclusion ledger".into());
            }
        }
        let mut entries = BTreeMap::new();
        let shared = [
            "TRAV13-4/DV7",
            "TRAV14D-3/DV8",
            "TRAV15-1/DV6-1",
            "TRAV15-2/DV6-2",
            "TRAV15D-1/DV6D-1",
            "TRAV15D-2/DV6D-2",
            "TRAV16D/DV11",
            "TRAV21/DV12",
            "TRAV4-4/DV10",
            "TRAV6-7/DV9",
        ];
        for e in c.entries {
            let (gene, allele) = e.id.split_once('*').ok_or("missing reference allele")?;
            if allele.len() != 2
                || !allele.bytes().all(|b| b.is_ascii_digit())
                || !match e.kind.as_str() {
                    "shared" => shared.contains(&gene),
                    "dedicated" => {
                        ["TRDV1", "TRDV2-1", "TRDV2-2", "TRDV4", "TRDV5"].contains(&gene)
                    }
                    _ => false,
                }
            {
                return Err("invalid mouse delta membership".into());
            }
            for (region, len) in e.regions.iter().zip([12, 10, 7]) {
                if region.len() != len
                    || region
                        .bytes()
                        .any(|b| b != b'.' && !ALPHABET[..20].contains(&b))
                    || region.bytes().all(|b| b == b'.')
                {
                    return Err("incomplete/ambiguous mouse delta region".into());
                }
            }
            if entries.insert(e.id, e.regions).is_some() {
                return Err("duplicate reference allele".into());
            }
        }
        for j in &c.j_alleles {
            let g = base.gene("mouse", j)?;
            if g.chain != "D" || g.region != "J" {
                return Err("delta J reference mismatch".into());
            }
        }
        Ok(Self {
            entries,
            j_alleles: c.j_alleles,
        })
    }
    pub fn resolve(&self, v: &str, j: &str, cdr3: &str) -> Result<Chain> {
        let loops=self.entries.get(v).ok_or_else(||format!("unsupported exact current mouse TRD V allele {v}; aliases, ordinary TRA, ORF and missing-loop records are not admitted"))?;
        if !self.j_alleles.iter().any(|x| x == j) {
            return Err("unsupported exact mouse TRDJ allele".into());
        }
        Ok(Chain {
            v: v.into(),
            loops: [
                encode(&loops[0], false)?,
                encode(&loops[1], false)?,
                encode(&loops[2], false)?,
            ],
            cdr3: encode(cdr3, true)?,
        })
    }
}
