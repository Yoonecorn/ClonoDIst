//! Native CSR stores explicit zero distances; absence means no retained edge.
use clonodist_core::Result;
use serde::{Deserialize, Serialize};
#[derive(Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Csr {
    pub rows: usize,
    pub cols: usize,
    pub indptr: Vec<usize>,
    pub indices: Vec<usize>,
    pub data: Vec<u32>,
}
impl Csr {
    pub fn from_dense(rows: &[Vec<u32>], radius: u32, budget: usize) -> Result<Self> {
        let pointers = rows
            .len()
            .checked_add(1)
            .and_then(|n| n.checked_mul(std::mem::size_of::<usize>()))
            .ok_or("CSR pointer overflow")?;
        if pointers > budget {
            return Err("CSR pointer budget exceeded".into());
        }
        let cols = rows.first().map_or(0, Vec::len);
        let mut pointers_vec = Vec::new();
        pointers_vec
            .try_reserve_exact(rows.len().checked_add(1).ok_or("CSR pointer overflow")?)
            .map_err(|_| "CSR allocation failed")?;
        pointers_vec.push(0);
        let mut out = Self {
            rows: rows.len(),
            cols,
            indptr: pointers_vec,
            indices: vec![],
            data: vec![],
        };
        for row in rows {
            if row.len() != cols {
                return Err("ragged dense input".into());
            }
            for (j, d) in row.iter().enumerate() {
                if *d <= radius {
                    let bytes = out
                        .data
                        .len()
                        .checked_add(1)
                        .and_then(|n| n.checked_mul(32))
                        .ok_or("CSR size overflow")?;
                    if bytes > budget - pointers {
                        return Err("CSR budget exceeded".into());
                    }
                    out.indices.push(j);
                    out.data.push(*d);
                }
            }
            out.indptr.push(out.data.len());
        }
        Ok(out)
    }
    pub fn validate(&self) -> Result<()> {
        if self.indptr.len() != self.rows.checked_add(1).ok_or("CSR overflow")?
            || self.indptr.first() != Some(&0)
            || self.indptr.last() != Some(&self.data.len())
            || self.indices.len() != self.data.len()
        {
            return Err("invalid CSR shape".into());
        }
        for w in self.indptr.windows(2) {
            if w[0] > w[1] || w[1] > self.data.len() {
                return Err("invalid CSR pointer".into());
            }
            let ix = &self.indices[w[0]..w[1]];
            if ix.iter().any(|x| *x >= self.cols) || ix.windows(2).any(|w| w[0] >= w[1]) {
                return Err("invalid CSR columns".into());
            }
        }
        Ok(())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn explicit_zero() {
        assert!(Csr::from_dense(&vec![vec![100]; 10], 0, 0).is_err());
        let c = Csr::from_dense(&[vec![0, 4, 10], vec![4, 0, 0]], 4, 1024).unwrap();
        assert_eq!(c.data, vec![0, 4, 4, 0, 0]);
        c.validate().unwrap();
        assert!(Csr::from_dense(&[vec![0; 10]], 0, 1).is_err());
    }
}
