//! Independent checked implementation of the frozen integer scientific profile.
use serde::{Deserialize, Serialize};
pub type Result<T> = std::result::Result<T, String>;
pub const PROFILE: &str = "m0-tcrdist3-0e6efa3-aligned";
pub const ALPHABET: &[u8] = b"ARNDCQEGHILKMFPSTWYVBZX";
pub const MAX_LEN: usize = 127;
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Params {
    pub dist_weight: u32,
    pub gap_penalty: u32,
    pub ntrim: usize,
    pub ctrim: usize,
    pub fixed_gappos: bool,
}
pub const CDR3: Params = Params {
    dist_weight: 1,
    gap_penalty: 4,
    ntrim: 3,
    ctrim: 2,
    fixed_gappos: false,
};
pub const GERMLINE: Params = Params {
    dist_weight: 1,
    gap_penalty: 4,
    ntrim: 0,
    ctrim: 0,
    fixed_gappos: true,
};
pub fn add(a: u32, b: u32) -> Result<u32> {
    a.checked_add(b).ok_or_else(|| "distance overflow".into())
}
pub fn mul(a: u32, b: u32) -> Result<u32> {
    a.checked_mul(b).ok_or_else(|| "distance overflow".into())
}
pub fn encode(s: &str, strict: bool) -> Result<Vec<u8>> {
    if !s.is_ascii() || s.len() > MAX_LEN {
        return Err("sequence must be ASCII and at most 127 residues".into());
    }
    if strict && (s.is_empty() || s.bytes().any(|b| !ALPHABET[..20].contains(&b))) {
        return Err("CDR3 requires nonempty canonical uppercase amino acids".into());
    }
    Ok(s.bytes()
        .map(|b| ALPHABET.iter().position(|x| *x == b).unwrap_or(23) as u8)
        .collect())
}
const fn raw_cost(a: u8, b: u8) -> u32 {
    if a >= 20 || b >= 20 || a == b {
        return 0;
    }
    let (x, y) = (ALPHABET[a as usize], ALPHABET[b as usize]);
    let p = if x < y { [x, y] } else { [y, x] };
    match &p {
        b"AS" | b"DN" | b"EK" | b"FW" | b"HN" | b"IM" | b"KQ" | b"LV" | b"MV" | b"NS" | b"QR"
        | b"ST" => 3,
        b"DE" | b"EQ" | b"HY" | b"IL" | b"KR" | b"LM" | b"WY" => 2,
        b"FY" | b"IV" => 1,
        _ => 4,
    }
}
const fn cost_table() -> [[u32; 24]; 24] {
    let mut out = [[0; 24]; 24];
    let mut i = 0;
    while i < 24 {
        let mut j = 0;
        while j < 24 {
            out[i][j] = raw_cost(i as u8, j as u8);
            j += 1;
        }
        i += 1;
    }
    out
}
const COSTS: [[u32; 24]; 24] = cost_table();
#[inline]
pub fn cost(a: u8, b: u8) -> u32 {
    if a > 23 || b > 23 {
        0
    } else {
        COSTS[a as usize][b as usize]
    }
}
/// Immutable validated sequence; enables repeated scoring without repeated validation.
#[derive(Clone, Debug)]
pub struct PreparedSequence(Vec<u8>);
impl PreparedSequence {
    pub fn new(bytes: &[u8]) -> Result<Self> {
        if bytes.len() > MAX_LEN || bytes.iter().any(|b| *b > 23) {
            return Err("invalid encoded sequence".into());
        }
        Ok(Self(bytes.to_vec()))
    }
    pub fn score(&self, other: &Self, p: Params) -> Result<u32> {
        region_inner(&self.0, &other.0, p)
    }
}
pub fn region(a: &[u8], b: &[u8], p: Params) -> Result<u32> {
    if a.len() > MAX_LEN || b.len() > MAX_LEN || a.iter().chain(b).any(|x| *x > 23) {
        return Err("invalid encoded sequence".into());
    }
    region_inner(a, b, p)
}
fn region_inner(a: &[u8], b: &[u8], p: Params) -> Result<u32> {
    if p.ntrim > MAX_LEN || p.ctrim > MAX_LEN {
        return Err("trim outside supported domain".into());
    }
    let short = a.len().min(b.len());
    let mut best;
    if a.len() == b.len() {
        let mut sum = 0;
        for i in p.ntrim..a.len().saturating_sub(p.ctrim) {
            sum = add(sum, cost(a[i], b[i]))?;
        }
        best = sum;
    } else {
        let (lo, hi) = if p.fixed_gappos {
            let g = 6isize.min(3 + (short as isize - 5).div_euclid(2));
            (g, g)
        } else {
            let (mut lo, mut hi) = (5isize, short as isize - 5);
            while lo > hi {
                lo -= 1;
                hi += 1;
            }
            (lo, hi)
        };
        // The gap coordinates are within [0, short] for both admitted formulas.
        // Score the first placement, then replace its moving prefix/suffix boundary.
        let mut sum = 0u32;
        for i in p.ntrim as isize..lo {
            let i = usize::try_from(i).map_err(|_| "unsafe gap coordinate")?;
            sum = add(
                sum,
                cost(
                    *a.get(i).ok_or("unsafe gap coordinate")?,
                    *b.get(i).ok_or("unsafe gap coordinate")?,
                ),
            )?;
        }
        for offset in p.ctrim as isize..short as isize - lo {
            let i = a.len() as isize - 1 - offset;
            let j = b.len() as isize - 1 - offset;
            if i < 0 || j < 0 {
                return Err("unsafe gap coordinate".into());
            }
            sum = add(sum, cost(a[i as usize], b[j as usize]))?;
        }
        best = sum;
        for g in lo..hi {
            if best == 0 {
                break;
            }
            let prefix = g as usize;
            if prefix >= p.ntrim {
                sum = add(sum, cost(a[prefix], b[prefix]))?;
            }
            let offset = short - prefix - 1;
            if offset >= p.ctrim {
                sum -= cost(a[a.len() - 1 - offset], b[b.len() - 1 - offset]);
            }
            best = best.min(sum);
        }
    }
    add(
        mul(best, p.dist_weight)?,
        mul(a.len().abs_diff(b.len()) as u32, p.gap_penalty)?,
    )
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Chain {
    pub v: String,
    pub loops: [Vec<u8>; 3],
    pub cdr3: Vec<u8>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Receptor {
    pub id: String,
    pub count: u64,
    pub metadata: std::collections::BTreeMap<String, String>,
    pub chains: Vec<Chain>,
}
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Breakdown {
    pub cdr1: u32,
    pub cdr2: u32,
    pub pmhc: u32,
    pub cdr3: u32,
    pub total: u32,
}
pub fn germline(a: &Chain, b: &Chain) -> Result<[u32; 3]> {
    Ok([
        region(&a.loops[0], &b.loops[0], GERMLINE)?,
        region(&a.loops[1], &b.loops[1], GERMLINE)?,
        region(&a.loops[2], &b.loops[2], GERMLINE)?,
    ])
}
pub fn chain(a: &Chain, b: &Chain) -> Result<Breakdown> {
    let [cdr1, cdr2, pmhc] = germline(a, b)?;
    let cdr3 = region(&a.cdr3, &b.cdr3, CDR3)?;
    let total = add(add(add(cdr1, cdr2)?, pmhc)?, mul(cdr3, 3)?)?;
    Ok(Breakdown {
        cdr1,
        cdr2,
        pmhc,
        cdr3,
        total,
    })
}
pub fn distance(a: &Receptor, b: &Receptor) -> Result<u32> {
    if a.chains.len() != b.chains.len() || a.chains.is_empty() {
        return Err("chain configuration mismatch".into());
    }
    a.chains
        .iter()
        .zip(&b.chains)
        .try_fold(0, |sum, (x, y)| add(sum, chain(x, y)?.total))
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn checked_overflow() {
        assert!(add(u32::MAX, 1).is_err());
        assert!(mul(u32::MAX, 3).is_err());
        let a = encode("AAAAAAAA", false).unwrap();
        let b = encode("WWWWWWWW", false).unwrap();
        assert!(region(
            &a,
            &b,
            Params {
                dist_weight: u32::MAX,
                ..GERMLINE
            }
        )
        .is_err());
    }
    #[test]
    fn malformed() {
        for s in ["", "A?", "é", "cas"] {
            assert!(encode(s, true).is_err());
        }
        assert!(encode(&"A".repeat(128), false).is_err());
        assert!(region(&[24], &[0], GERMLINE).is_err());
    }
    #[test]
    fn symmetry_exhaustive() {
        let letters = ["A", "S", "W"];
        for n in 0..8 {
            for m in 0..8 {
                for x in letters {
                    for y in letters {
                        let a = encode(&x.repeat(n), false).unwrap();
                        let b = encode(&y.repeat(m), false).unwrap();
                        assert_eq!(region(&a, &b, CDR3).unwrap(), region(&b, &a, CDR3).unwrap());
                    }
                }
            }
        }
    }
}
