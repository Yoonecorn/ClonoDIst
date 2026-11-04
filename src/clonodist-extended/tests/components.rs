use clonodist_extended::*;
#[test]
fn strict_profile_and_identity() {
    let p = Profile::pinned(2);
    p.validate(2, false).unwrap();
    assert_eq!(
        checksum(&p).unwrap(),
        checksum(&serde_json::to_value(&p).unwrap()).unwrap()
    );
    let mut x = serde_json::to_value(&p).unwrap();
    x["regions"][0]["parameters"]["typo"] = 1.into();
    assert!(serde_json::from_value::<Profile>(x).is_err());
    let mut p = Profile::cdr3(1);
    assert!(p.validate(1, false).is_err());
    p.validate(1, true).unwrap();
    p.chain_weights[0] = 0;
    assert!(p.validate(1, true).is_err());
}
#[test]
fn fisher_known() {
    assert!((meta::fisher(3, 3, 0, 3) - 0.05).abs() < 1e-12);
    assert!((meta::fisher(0, 4, 4, 4) - 1.).abs() < 1e-12);
}
#[test]
fn scope() {
    assert!(chain_config("gamma,delta", "human").is_ok());
    assert!(chain_config("gamma", "mouse").is_ok());
    assert!(chain_config("delta", "mouse").is_err());
    assert!(chain_config("gamma", "rhesus").is_err());
}
#[test]
fn model_decimal_roundtrip() {
    use clonodist_extended::meta::*;
    use std::collections::BTreeMap;
    for i in 1..1000 {
        let rules = Rules {
            radii: vec![0, 12],
            background_cap: 0.01,
            minimum_hits: 2,
            maximum_jaccard: 0.8,
            same_vj: false,
            motifs: vec![],
            resamples: 64,
            seed: 2,
        };
        let mut model = Model {
            schema: "clonodist.extended.meta.v1".into(),
            phase: "discovered".into(),
            organism: "human".into(),
            chains: vec!["beta".into()],
            profile: Profile::pinned(1),
            profile_sha256: checksum(&Profile::pinned(1)).unwrap(),
            rules,
            candidates: vec![],
            discovery_subjects: vec![],
            discovery_cohorts: vec![],
            discovery_signatures: vec![],
            dataset_hashes: BTreeMap::new(),
            qualification: serde_json::Value::Null,
            seal: String::new(),
        };
        let c = Candidate {
            center: clonodist_core::Receptor {
                id: "x".into(),
                count: 1,
                metadata: BTreeMap::new(),
                chains: vec![clonodist_core::Chain {
                    v: "v".into(),
                    loops: [vec![], vec![], vec![]],
                    cdr3: vec![0],
                }],
            },
            radius: Some(12),
            discovery_hits: 3,
            background_hits: 0,
            enrichment: i as f64 / 37.,
            stability: 0.734375,
            radius_range: Some([0, 12]),
            membership_jaccard: i as f64 / 1000.,
            reason: "discovery_selected".into(),
            heldout_hits: 0,
            heldout_background_hits: 0,
            adjusted_p: 1. / i as f64,
            qualified: false,
        };
        model.candidates.push(c);
        model.reseal().unwrap();
        let text = serde_json::to_string_pretty(&model).unwrap();
        let restored: Model = serde_json::from_str(&text).unwrap();
        restored.verify().unwrap();
    }
}
