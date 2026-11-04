//! Descriptive neighborhoods. Similarity is not antigen specificity evidence.
use clonodist_core::{Receptor, Result};
use std::collections::BTreeSet;
pub struct Summary {
    pub neighbors: usize,
    pub weighted_count: u64,
    pub subjects: usize,
    pub samples: usize,
}
/// Exact partial selection for inputs whose target IDs have been validated unique.
pub fn ranked_unique(
    row: &[u32],
    targets: &[Receptor],
    radius: u32,
    k: Option<usize>,
    exclude: Option<&str>,
) -> Vec<(usize, u32)> {
    let mut out: Vec<_> = row
        .iter()
        .enumerate()
        .filter(|(j, d)| **d <= radius && exclude != Some(targets[*j].id.as_str()))
        .map(|(j, d)| (j, *d))
        .collect();
    let compare = |(i, d): &(usize, u32), (j, e): &(usize, u32)| {
        d.cmp(e).then_with(|| targets[*i].id.cmp(&targets[*j].id))
    };
    if let Some(k) = k {
        if k < out.len() {
            out.select_nth_unstable_by(k, compare);
            out.truncate(k);
        }
    }
    out.sort_unstable_by(compare);
    out
}
pub fn ranked(
    row: &[u32],
    targets: &[Receptor],
    radius: u32,
    k: Option<usize>,
    exclude: Option<&str>,
) -> Vec<(usize, u32)> {
    let mut out: Vec<_> = row
        .iter()
        .enumerate()
        .filter(|(j, d)| **d <= radius && exclude != Some(targets[*j].id.as_str()))
        .map(|(j, d)| (j, *d))
        .collect();
    out.sort_unstable_by(|(i, d), (j, e)| {
        d.cmp(e).then_with(|| targets[*i].id.cmp(&targets[*j].id))
    });
    if let Some(k) = k {
        out.truncate(k);
    }
    out
}
pub fn summarize(edges: &[(usize, u32)], targets: &[Receptor]) -> Result<Summary> {
    let mut count = 0u64;
    let mut subjects = BTreeSet::new();
    let mut samples = BTreeSet::new();
    for (j, _) in edges {
        let t = &targets[*j];
        count = count
            .checked_add(t.count)
            .ok_or("neighborhood count overflow")?;
        if let Some(s) = t.metadata.get("subject") {
            subjects.insert(s);
        }
        if let Some(s) = t.metadata.get("sample") {
            samples.insert(s);
        }
    }
    Ok(Summary {
        neighbors: edges.len(),
        weighted_count: count,
        subjects: subjects.len(),
        samples: samples.len(),
    })
}

#[cfg(test)]
mod selection_tests {
    use super::*;
    #[test]
    fn unique_selection_matches_original_sort_at_ties_and_boundaries() {
        for n in [0usize, 1, 2, 5, 32, 100] {
            let targets: Vec<_> = (0..n)
                .map(|i| Receptor {
                    id: format!("id-{:03}", n - i),
                    count: 1,
                    metadata: Default::default(),
                    chains: vec![],
                })
                .collect();
            let row: Vec<_> = (0..n)
                .map(|i| [0, 36, 37, 5, 5, 1, u32::MAX][i % 7])
                .collect();
            for radius in [0, 36, u32::MAX] {
                for k in [None, Some(0), Some(1), Some(5), Some(n), Some(n + 1)] {
                    for exclude in [None, targets.first().map(|r| r.id.as_str())] {
                        assert_eq!(
                            ranked_unique(&row, &targets, radius, k, exclude),
                            ranked(&row, &targets, radius, k, exclude)
                        );
                    }
                }
            }
        }
    }
}
