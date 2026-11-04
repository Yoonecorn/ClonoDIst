//! Project-authored, reference-free probes. These do not replace scientific qualification.
use clonodist_core::{distance, encode, Chain, Receptor};
use clonodist_distance::{Config, Engine};
use clonodist_extended::{chain_config, Profile};
use clonodist_neighbors::{ranked, summarize};
use std::collections::BTreeMap;

fn invented(id: &str, seq: &str, chains: usize) -> Receptor {
    let make = |i| Chain {
        v: format!("invented-v-{i}"),
        loops: [
            encode("AAAA", true).unwrap(),
            encode("RRRR", true).unwrap(),
            encode("NNNN", true).unwrap(),
        ],
        cdr3: encode(seq, true).unwrap(),
    };
    Receptor {
        id: id.into(),
        count: 1,
        metadata: BTreeMap::new(),
        chains: (0..chains).map(make).collect(),
    }
}

#[test]
fn synthetic_profiles_and_streaming_agree() {
    for chains in [1, 2] {
        let records = vec![
            invented("a", "AAAACAAAAA", chains),
            invented("b", "AAAADAAAAA", chains),
        ];
        let profile = Profile::pinned(chains);
        profile.validate(chains, false).unwrap();
        let expected = if chains == 1 { 12 } else { 24 };
        assert_eq!(distance(&records[0], &records[1]).unwrap(), expected);
        assert_eq!(profile.score(&records[0], &records[1]).unwrap(), expected);
        for optimized in [false, true] {
            let config = Config {
                threads: 1,
                tile: 1,
                memory: 64 * 1024 * 1024,
                optimized,
                radius: None,
            };
            let engine = Engine::new(&records, &records, config).unwrap();
            let mut rows = Vec::new();
            engine
                .run(|_, batch| {
                    rows.extend_from_slice(batch);
                    Ok(())
                })
                .unwrap();
            assert_eq!(rows, vec![vec![0, expected], vec![expected, 0]]);
        }
    }
}

#[test]
fn synthetic_edges_preserve_inclusive_radius_and_zero() {
    let targets = vec![
        invented("b", "AAAACAAAAA", 1),
        invented("a", "AAAADAAAAA", 1),
    ];
    let edges = ranked(&[0, 12], &targets, 12, None, None);
    assert_eq!(edges, vec![(0, 0), (1, 12)]);
    let summary = summarize(&edges, &targets).unwrap();
    assert_eq!((summary.neighbors, summary.weighted_count), (2, 2));
    assert_eq!(
        chain_config("gamma,delta", "human").unwrap(),
        ["gamma", "delta"]
    );
    assert!(chain_config("gamma,delta", "mouse").is_err());
    assert!(encode("invalid!", true).is_err());
}
