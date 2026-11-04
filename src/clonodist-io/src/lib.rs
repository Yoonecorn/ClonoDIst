//! Strict normalization, conservation, and atomic directory transactions.
use clonodist_core::{Receptor, Result};
use clonodist_germline::Reference;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
    path::{Path, PathBuf},
};
#[derive(Debug, Serialize, Deserialize)]
pub struct Ledger {
    pub input_rows: usize,
    pub accepted_rows: usize,
    pub output_records: usize,
    pub input_count: u64,
    pub output_count: u64,
    pub mapping: Vec<(usize, String)>,
}
pub fn hash(path: &Path) -> Result<String> {
    use std::io::Read;
    let meta = std::fs::metadata(path).map_err(|e| e.to_string())?;
    if !meta.is_file() {
        return Err("hash input must be a regular file".into());
    }
    let expected = meta.len();
    let mut total = 0u64;
    let mut f = std::fs::File::open(path)
        .map_err(|e| e.to_string())?
        .take(expected.checked_add(1).ok_or("file size overflow")?);
    let mut h = Sha256::new();
    let mut b = [0; 65536];
    loop {
        let n = f.read(&mut b).map_err(|e| e.to_string())?;
        if n == 0 {
            break;
        }
        total = total.checked_add(n as u64).ok_or("file size overflow")?;
        h.update(&b[..n]);
    }
    if total != expected {
        return Err("file size changed during hashing".into());
    }
    Ok(format!("{:x}", h.finalize()))
}
pub fn read_tsv(
    path: &Path,
    reference: &Reference,
    organism: &str,
    chains: &[String],
    dedup: bool,
    memory: usize,
) -> Result<(Vec<Receptor>, Ledger)> {
    let meta = std::fs::metadata(path).map_err(|e| e.to_string())?;
    if !meta.is_file() {
        return Err("input must be a regular file".into());
    }
    let bytes = meta.len();
    if bytes > memory as u64 / 16 {
        return Err("input exceeds memory budget (conservative 16x text allowance)".into());
    }
    let file = std::fs::File::open(path).map_err(|e| e.to_string())?;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(std::io::Read::take(
            file,
            bytes.checked_add(1).ok_or("input byte size overflow")?,
        ));
    let headers = reader.headers().map_err(|e| e.to_string())?.clone();
    let mut seen = BTreeSet::new();
    for h in &headers {
        if h.is_empty() || !seen.insert(h) {
            return Err("empty/duplicate column name".into());
        }
    }
    let pos = |h: &str| {
        headers
            .iter()
            .position(|x| x == h)
            .ok_or_else(|| format!("missing column {h}"))
    };
    let mut planned = 0usize;
    let mut records: Vec<Receptor> = Vec::new();
    let mut groups: BTreeMap<Vec<(String, String)>, usize> = BTreeMap::new();
    let mut ids = BTreeSet::new();
    let mut ledger = Ledger {
        input_rows: 0,
        accepted_rows: 0,
        output_records: 0,
        input_count: 0,
        output_count: 0,
        mapping: vec![],
    };
    for row in reader.records() {
        let row = row.map_err(|e| e.to_string())?;
        ledger.input_rows += 1;
        let rownum = ledger.input_rows;
        let row_plan = row
            .as_slice()
            .len()
            .checked_mul(16)
            .and_then(|n| n.checked_add(headers.len().checked_mul(1024)?))
            .and_then(|n| n.checked_add(8192))
            .ok_or("input allocation overflow")?;
        planned = planned
            .checked_add(row_plan)
            .ok_or("input allocation overflow")?;
        if planned > memory {
            return Err("normalized input allocation plan exceeds budget".into());
        }
        let get = |h: &str| -> Result<&str> {
            let v = &row[pos(h)?];
            if v.is_empty() {
                Err(format!("row {rownum}: empty {h}"))
            } else {
                Ok(v)
            }
        };
        let count = match headers.iter().position(|x| x == "count") {
            Some(i) => row[i]
                .parse::<u64>()
                .map_err(|_| format!("row {rownum}: count must be positive integer"))?,
            None => 1,
        };
        if count == 0 {
            return Err(format!("row {rownum}: count must be positive"));
        }
        ledger.input_count = ledger
            .input_count
            .checked_add(count)
            .ok_or("count overflow")?;
        let id = match headers.iter().position(|x| x == "id") {
            Some(i) => {
                if row[i].is_empty() {
                    return Err("empty id".into());
                }
                row[i].to_owned()
            }
            None => format!("row-{rownum:09}"),
        };
        if id.chars().any(|c| c.is_control() || c == '"') {
            return Err("id contains TSV framing character".into());
        }
        if !ids.insert(id.clone()) {
            return Err(format!("duplicate id {id}"));
        }
        let mut resolved = vec![];
        for ch in chains {
            let l = match ch.as_str() {
                "alpha" => "a",
                "beta" => "b",
                _ => return Err("unsupported chain".into()),
            };
            resolved.push(
                reference
                    .resolve(
                        organism,
                        ch,
                        get(&format!("v_{l}_gene"))?,
                        Some(get(&format!("j_{l}_gene"))?),
                        get(&format!("cdr3_{l}_aa"))?,
                        true,
                    )
                    .map_err(|e| format!("row {rownum}: {e}"))?,
            );
        }
        let metadata: BTreeMap<String, String> = headers
            .iter()
            .zip(row.iter())
            .filter(|(h, _)| !["id", "count", "cell_id", "clone_id"].contains(h))
            .map(|(h, v)| (h.into(), v.into()))
            .collect();
        if metadata.values().any(|v| v.is_empty()) {
            return Err(format!("row {rownum}: empty grouping metadata"));
        }
        let key = metadata
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect::<Vec<_>>();
        if let Some(i) = if dedup {
            groups.get(&key).copied()
        } else {
            None
        } {
            records[i].count = records[i]
                .count
                .checked_add(count)
                .ok_or("count aggregation overflow")?;
            ledger.mapping.push((rownum, records[i].id.clone()));
        } else {
            if dedup {
                groups.insert(key, records.len());
            }
            ledger.mapping.push((rownum, id.clone()));
            records.push(Receptor {
                id,
                count,
                metadata,
                chains: resolved,
            });
        }
        ledger.accepted_rows += 1;
    }
    if reader.position().byte() != bytes {
        return Err("input size changed during parse".into());
    }
    if records.is_empty() {
        return Err("empty repertoire".into());
    }
    ledger.output_records = records.len();
    ledger.output_count = records
        .iter()
        .try_fold(0u64, |s, r| s.checked_add(r.count).ok_or("count overflow"))?;
    if ledger.input_count != ledger.output_count {
        return Err("count conservation failure".into());
    }
    Ok((records, ledger))
}
/// Commit all output artifacts together by renaming an adjacent new directory.
pub struct Transaction {
    pub temporary: PathBuf,
    destination: PathBuf,
    committed: bool,
}
impl Transaction {
    pub fn new(destination: &Path) -> Result<Self> {
        if destination.exists() {
            return Err("output directory already exists".into());
        }
        let parent = destination
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        if !parent.is_dir() {
            return Err("output parent does not exist".into());
        }
        let name = destination
            .file_name()
            .ok_or("invalid output directory")?
            .to_string_lossy();
        let temporary = parent.join(format!(".{name}.tmp-{}", std::process::id()));
        std::fs::create_dir(&temporary).map_err(|e| e.to_string())?;
        Ok(Self {
            temporary,
            destination: destination.into(),
            committed: false,
        })
    }
    pub fn json<T: Serialize>(&self, name: &str, value: &T) -> Result<()> {
        let p = self.temporary.join(name);
        let mut f = std::io::BufWriter::with_capacity(
            65536,
            std::fs::File::create(p).map_err(|e| e.to_string())?,
        );
        serde_json::to_writer_pretty(&mut f, value).map_err(|e| e.to_string())?;
        f.write_all(b"\n").map_err(|e| e.to_string())?;
        f.flush().map_err(|e| e.to_string())?;
        f.get_ref().sync_all().map_err(|e| e.to_string())
    }
    pub fn commit(mut self) -> Result<()> {
        if self.destination.exists() {
            return Err("output appeared during run".into());
        }
        std::fs::rename(&self.temporary, &self.destination).map_err(|e| e.to_string())?;
        self.committed = true;
        Ok(())
    }
}
impl Drop for Transaction {
    fn drop(&mut self) {
        if !self.committed {
            let _ = std::fs::remove_dir_all(&self.temporary);
        }
    }
}

pub mod compat;
