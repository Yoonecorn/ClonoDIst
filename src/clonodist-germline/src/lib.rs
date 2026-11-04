//! External, checksum-pinned reference. No germline dataset is bundled in the executable.
use clonodist_core::{encode, Chain, Result};
use sha2::{Digest, Sha256};
use std::{collections::BTreeMap, path::Path};
pub const DB_SHA256: &str = "1d7fcac37ad6d955b3a3f693526e36623186464aef2c528cf8e0322dc89706cd";
#[derive(Clone, Debug)]
pub struct Gene {
    pub chain: String,
    pub region: String,
    pub loops: Vec<String>,
}
pub struct Reference {
    pub genes: BTreeMap<(String, String), Gene>,
    pub sha256: String,
}
impl Reference {
    pub fn load(path: &Path) -> Result<Self> {
        let meta = std::fs::metadata(path).map_err(|e| e.to_string())?;
        if !meta.is_file() || meta.len() != 1234156 {
            return Err("reference must be regular file of pinned byte size".into());
        }
        let mut bytes = Vec::new();
        let f = std::fs::File::open(path).map_err(|e| e.to_string())?;
        std::io::Read::read_to_end(&mut std::io::Read::take(f, 1234157), &mut bytes)
            .map_err(|e| e.to_string())?;
        let sha256 = format!("{:x}", Sha256::digest(&bytes));
        if sha256 != DB_SHA256 {
            return Err("reference SHA-256 mismatch".into());
        }
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(bytes.as_slice());
        let headers = rdr.headers().map_err(|e| e.to_string())?.clone();
        let col = |s: &str| {
            headers
                .iter()
                .position(|h| h == s)
                .ok_or_else(|| format!("reference missing {s}"))
        };
        let (id, org, ch, reg, cdr, prot, cols) = (
            col("id")?,
            col("organism")?,
            col("chain")?,
            col("region")?,
            col("cdrs")?,
            col("aligned_protseq")?,
            col("cdr_columns")?,
        );
        let mut genes = BTreeMap::new();
        for r in rdr.records() {
            let r = r.map_err(|e| e.to_string())?;
            let loops: Vec<String> = r[cdr].split(';').map(str::to_owned).collect();
            if !r[cdr].is_empty() {
                for (s, interval) in loops.iter().zip(r[cols].split(';')) {
                    let (a, b) = interval
                        .split_once('-')
                        .ok_or("invalid reference interval")?;
                    let a: usize = a.parse().map_err(|_| "invalid interval")?;
                    let b: usize = b.parse().map_err(|_| "invalid interval")?;
                    if a == 0
                        || r[prot].get((a - 1).min(r[prot].len())..b.min(r[prot].len()))
                            != Some(s.as_str())
                    {
                        return Err(format!(
                            "reference loop extraction mismatch {}:{} {} {} {:?} {:?}",
                            &r[org],
                            &r[id],
                            a,
                            b,
                            s,
                            r[prot].get((a - 1).min(r[prot].len())..b.min(r[prot].len()))
                        ));
                    }
                }
            }
            genes.insert(
                (r[org].to_owned(), r[id].to_owned()),
                Gene {
                    chain: r[ch].to_owned(),
                    region: r[reg].to_owned(),
                    loops,
                },
            );
        }
        Ok(Self { genes, sha256 })
    }
    /// Educational reference containing only invented identifiers and invented loops.
    /// This opt-in path accepts one exact project-authored marker file, never a real reference.
    pub fn load_synthetic_demo(path: &Path) -> Result<Self> {
        const MARKER: &[u8] = b"ClonoDist synthetic demonstration reference v1\n";
        let meta = std::fs::metadata(path).map_err(|e| e.to_string())?;
        if !meta.is_file() || meta.len() != MARKER.len() as u64 {
            return Err("synthetic demonstration reference marker mismatch".into());
        }
        let data = std::fs::read(path).map_err(|e| e.to_string())?;
        if data != MARKER {
            return Err("synthetic demonstration reference marker mismatch".into());
        }
        let mut genes = BTreeMap::new();
        for (label, name) in [("A", "A"), ("B", "B"), ("G", "G"), ("D", "D")] {
            for region in ["V", "J"] {
                let gene = Gene {
                    chain: label.into(),
                    region: region.into(),
                    loops: if region == "V" {
                        vec!["AAAA".into(), "RRRR".into(), "NNNN".into()]
                    } else {
                        vec![]
                    },
                };
                genes.insert(("human".into(), format!("FAKE_{name}_{region}*01")), gene);
            }
        }
        Ok(Self {
            genes,
            sha256: format!("{:x}", Sha256::digest(MARKER)),
        })
    }
    pub fn gene(&self, organism: &str, v: &str) -> Result<&Gene> {
        self.genes
            .get(&(organism.to_owned(), v.to_owned()))
            .ok_or_else(|| format!("unknown exact gene/allele {organism}:{v}"))
    }
    pub fn resolve(
        &self,
        organism: &str,
        chain: &str,
        v: &str,
        j: Option<&str>,
        cdr3: &str,
        strict: bool,
    ) -> Result<Chain> {
        if !["human", "mouse"].contains(&organism) {
            return Err("organism must be human or mouse".into());
        }
        let expected = match chain {
            "alpha" => "A",
            "beta" => "B",
            _ => return Err("chain must be alpha or beta".into()),
        };
        let g = self.gene(organism, v)?;
        if g.loops.len() < 3 {
            return Err("V gene lacks three germline loops".into());
        }
        if strict {
            if g.chain != expected || g.region != "V" || !v.contains('*') {
                return Err("V chain/region/allele mismatch".into());
            }
            let j = j.ok_or("missing J allele")?;
            let jg = self.gene(organism, j)?;
            if jg.chain != expected || jg.region != "J" || !j.contains('*') {
                return Err("J chain/region/allele mismatch".into());
            }
        }
        Ok(Chain {
            v: v.into(),
            loops: [
                encode(&g.loops[0], false)?,
                encode(&g.loops[1], false)?,
                encode(&g.loops[2], false)?,
            ],
            cdr3: encode(cdr3, strict)?,
        })
    }
}
