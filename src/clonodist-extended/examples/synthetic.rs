//! Educational API demonstration with invented sequences; no germline reference is loaded.
use clonodist_core::{distance, encode, Chain, Receptor};
use clonodist_distance::Config;
use clonodist_extended::{
    chain_config,
    meta::{discover, matching, qualify, Rules},
    Profile,
};
use clonodist_neighbors::{ranked, summarize};
use std::collections::BTreeMap;

fn chain(name: &str, cdr3: &str) -> Chain {
    Chain {
        v: name.into(),
        loops: [
            encode("AAAA", true).unwrap(),
            encode("RRRR", true).unwrap(),
            encode("NNNN", true).unwrap(),
        ],
        cdr3: encode(cdr3, true).unwrap(),
    }
}
fn receptor(id: &str, seq: &str, paired: bool) -> Receptor {
    let mut chains = vec![chain("invented-v-one", seq)];
    if paired {
        chains.push(chain("invented-v-two", "CCCCAAAAAA"));
    }
    Receptor {
        id: id.into(),
        count: 1,
        metadata: BTreeMap::new(),
        chains,
    }
}
fn meta_record(id: &str, sequence: &str, phase: &str) -> Receptor {
    let mut r = receptor(id, sequence, false);
    r.metadata.insert("subject".into(), id.into());
    r.metadata.insert("cohort".into(), phase.into());
    r
}
fn meta_demo() -> Result<(), String> {
    let fg = vec![
        meta_record("d1", "AAAACAAAAA", "discovery"),
        meta_record("d2", "RAAACAAAAA", "discovery"),
        meta_record("d3", "NAAACAAAAA", "discovery"),
    ];
    let bg = vec![
        meta_record("db1", "AAWWWWWWAA", "discovery"),
        meta_record("db2", "RAWWWWWWAA", "discovery"),
    ];
    let held = vec![
        meta_record("h1", "DAAACAAAAA", "heldout"),
        meta_record("h2", "CAAACAAAAA", "heldout"),
    ];
    let held_bg = vec![
        meta_record("hb1", "NAWWWWWWAA", "heldout"),
        meta_record("hb2", "DAWWWWWWAA", "heldout"),
    ];
    let config = Config {
        threads: 1,
        tile: 1,
        memory: 64 * 1024 * 1024,
        optimized: true,
        radius: None,
    };
    let rules = Rules {
        radii: vec![0, 12],
        background_cap: 0.0,
        minimum_hits: 2,
        maximum_jaccard: 0.8,
        same_vj: false,
        motifs: vec![],
        resamples: 64,
        seed: 42,
    };
    let model = discover(
        &fg,
        &bg,
        Profile::pinned(1),
        rules,
        "human".into(),
        vec!["beta".into()],
        config.clone(),
        BTreeMap::new(),
    )?;
    let selected = model
        .candidates
        .iter()
        .filter(|c| c.reason == "discovery_selected")
        .count();
    let model = qualify(model, &held, &held_bg, config, BTreeMap::new())?;
    let qualified: Vec<_> = model.candidates.iter().filter(|c| c.qualified).collect();
    let application = vec![
        meta_record("a1", "EAAACAAAAA", "application"),
        meta_record("a2", "EAWWWWWWAA", "application"),
    ];
    let mut hits = 0;
    for c in &qualified {
        let scores = application
            .iter()
            .map(|r| model.profile.score(&c.center, r))
            .collect::<Result<Vec<_>, _>>()?;
        hits += matching(
            &c.center,
            &application,
            &scores,
            c.radius.unwrap(),
            &model.rules,
        )
        .len();
    }
    println!(
        "operational example: {selected} selected, {} qualified, {hits} application matches",
        qualified.len()
    );
    Ok(())
}
fn main() -> Result<(), String> {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "all".into());
    let allowed = [
        "all",
        "core",
        "extended",
        "paired",
        "gamma-delta",
        "neighborhood",
        "meta",
    ];
    if !allowed.contains(&mode.as_str()) {
        return Err(format!("mode must be one of {}", allowed.join(", ")));
    }
    let a = receptor("invented-1", "AAAACAAAAA", false);
    let b = receptor("invented-2", "AAAADAAAAA", false);
    let d = distance(&a, &b)?;
    if mode == "all" || mode == "core" {
        println!("core distance: {d}");
    }
    if mode == "all" || mode == "extended" {
        let p = Profile::pinned(1);
        p.validate(1, false)?;
        assert_eq!(p.score(&a, &b)?, d);
        println!("extended pinned-profile distance: {d}");
    }
    if mode == "all" || mode == "paired" {
        let a = receptor("invented-1", "AAAACAAAAA", true);
        let b = receptor("invented-2", "AAAADAAAAA", true);
        let p = Profile::pinned(2);
        p.validate(2, false)?;
        let n = p.score(&a, &b)?;
        assert_eq!(n, distance(&a, &b)?);
        println!("paired distance: {n}");
    }
    if mode == "all" || mode == "gamma-delta" {
        assert_eq!(chain_config("gamma,delta", "human")?, ["gamma", "delta"]);
        let a = receptor("invented-1", "AAAACAAAAA", true);
        let b = receptor("invented-2", "AAAADAAAAA", true);
        println!(
            "human gamma/delta synthetic score: {}",
            Profile::pinned(2).score(&a, &b)?
        );
    }
    if mode == "all" || mode == "neighborhood" {
        let targets = vec![a, b];
        let edges = ranked(&[0, d], &targets, d, None, None);
        let summary = summarize(&edges, &targets)?;
        println!(
            "neighborhood radius {d}: {} neighbors, supplied count {}",
            summary.neighbors, summary.weighted_count
        );
    }
    if mode == "all" || mode == "meta" {
        meta_demo()?;
    }
    Ok(())
}
