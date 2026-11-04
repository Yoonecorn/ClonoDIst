use clap::{Args, Parser, Subcommand};
use clonodist_core::{Receptor, Result};
use clonodist_distance::Config;
use clonodist_extended::{meta::*, *};
use clonodist_germline::Reference;
use clonodist_io::{hash, Transaction};
use std::{
    collections::BTreeMap,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};
#[derive(Parser)]
#[command(
    name = "clonodist-extended",
    version,
    about = "Profile-aware Extended workflows; Core behavior remains frozen"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}
#[derive(Subcommand)]
enum Command {
    Profile {
        #[arg(long, default_value = "pinned")]
        kind: String,
        #[arg(long, default_value_t = 1)]
        chains: usize,
        #[arg(long)]
        output: PathBuf,
    },
    Distance(Distance),
    Discover(Discover),
    Qualify(Evaluate),
    Apply(Apply),
}
#[derive(Args, Clone)]
struct Run {
    #[arg(long)]
    db: PathBuf,
    /// Use only the included invented demonstration reference; never use for real data.
    #[arg(long)]
    synthetic_demo: bool,
    #[arg(long)]
    mouse_trd_reference: Option<PathBuf>,
    #[arg(long)]
    input: PathBuf,
    #[arg(long)]
    output: PathBuf,
    #[arg(long, default_value = "human")]
    organism: String,
    #[arg(long, default_value = "beta")]
    chains: String,
    #[arg(long)]
    profile: Option<PathBuf>,
    #[arg(long)]
    profile_sha256: Option<String>,
    #[arg(long)]
    experimental: bool,
    #[arg(long, default_value_t = 1)]
    threads: usize,
    #[arg(long, default_value_t = 32)]
    tile: usize,
    #[arg(long, default_value_t = 1073741824)]
    memory: usize,
    #[arg(long)]
    scalar: bool,
}
type LoadedRun = (
    Reference,
    Vec<String>,
    Profile,
    Vec<Receptor>,
    BTreeMap<String, String>,
);

impl Run {
    fn read_input(&self, path: &Path, db: &Reference, chains: &[String]) -> Result<Vec<Receptor>> {
        if let Some(reference) = &self.mouse_trd_reference {
            let supplement = mouse_trd::MouseTrdReference::load(reference, db)?;
            read_with_mouse_trd(
                path,
                db,
                &self.organism,
                chains,
                self.memory / 2,
                &supplement,
            )
        } else {
            read(path, db, &self.organism, chains, self.memory / 2)
        }
    }
    fn config(&self) -> Config {
        Config {
            threads: self.threads,
            tile: self.tile,
            memory: self.memory,
            optimized: !self.scalar,
            radius: None,
        }
    }
    fn load(&self) -> Result<LoadedRun> {
        let chains = if self.mouse_trd_reference.is_some() {
            if self.organism != "mouse" || !["delta", "gamma,delta"].contains(&self.chains.as_str())
            {
                return Err("mouse TRD reference requires mouse delta or gamma,delta".into());
            }
            self.chains.split(',').map(str::to_owned).collect()
        } else {
            chain_config(&self.chains, &self.organism)?
        };
        let p = if let Some(path) = &self.profile {
            load_json::<Profile>(path, self.memory)?
        } else {
            Profile::pinned(chains.len())
        };
        p.validate(chains.len(), self.experimental)?;
        if self.profile.is_some() && self.profile_sha256.is_none() {
            return Err("external profile requires --profile-sha256".into());
        }
        if self
            .profile_sha256
            .as_ref()
            .is_some_and(|h| checksum(&p).ok().as_ref() != Some(h))
        {
            return Err("profile checksum mismatch".into());
        }
        let db = if self.synthetic_demo {
            if self.mouse_trd_reference.is_some() {
                return Err("synthetic demonstration cannot use a mouse TRD supplement".into());
            }
            eprintln!(
                "SYNTHETIC DEMONSTRATION ONLY: invented reference, not scientifically qualified"
            );
            Reference::load_synthetic_demo(&self.db)?
        } else {
            Reference::load(&self.db)?
        };
        let input_hash = hash(&self.input)?;
        let rs = self.read_input(&self.input, &db, &chains)?;
        let mut hashes = BTreeMap::from([
            ("input".into(), input_hash),
            ("db".into(), db.sha256.clone()),
        ]);
        if self.mouse_trd_reference.is_some() {
            hashes.insert(
                "mouse_trd_reference".into(),
                mouse_trd::REFERENCE_SHA256.into(),
            );
        }
        if let Some(path) = &self.profile {
            hashes.insert("profile_file".into(), hash(path)?);
        }
        Ok((db, chains, p, rs, hashes))
    }
}
#[derive(Args)]
struct Distance {
    #[command(flatten)]
    run: Run,
    #[arg(long)]
    target: Option<PathBuf>,
    #[arg(long, default_value = "dense")]
    format: String,
    #[arg(long)]
    radius: Option<u32>,
    #[arg(long)]
    top_k: Option<usize>,
    #[arg(long)]
    exclude_self: bool,
}
#[derive(Args)]
struct Discover {
    #[command(flatten)]
    run: Run,
    #[arg(long)]
    background: PathBuf,
    #[arg(long)]
    rules: PathBuf,
}
#[derive(Args)]
struct Evaluate {
    #[command(flatten)]
    run: Run,
    #[arg(long)]
    background: PathBuf,
    #[arg(long)]
    model: PathBuf,
}
#[derive(Args)]
struct Apply {
    #[command(flatten)]
    run: Run,
    #[arg(long)]
    model: PathBuf,
}
fn load_json<T: serde::de::DeserializeOwned>(p: &Path, memory: usize) -> Result<T> {
    use std::io::Read;
    let meta = std::fs::metadata(p).map_err(|e| e.to_string())?;
    if !meta.is_file() || meta.len() > memory as u64 / 64 {
        return Err("JSON exceeds admission".into());
    }
    let mut bytes = vec![];
    std::fs::File::open(p)
        .map_err(|e| e.to_string())?
        .take(meta.len() + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
    if bytes.len() as u64 != meta.len() {
        return Err("JSON changed during read".into());
    }
    serde_json::from_slice(&bytes).map_err(|e| e.to_string())
}
fn writer(tx: &Transaction, name: &str) -> Result<BufWriter<std::fs::File>> {
    Ok(BufWriter::with_capacity(
        65536,
        std::fs::File::create(tx.temporary.join(name)).map_err(|e| e.to_string())?,
    ))
}
fn finish(mut f: BufWriter<std::fs::File>) -> Result<()> {
    f.flush().map_err(|e| e.to_string())?;
    f.get_ref().sync_all().map_err(|e| e.to_string())
}
fn line<T: serde::Serialize>(f: &mut impl Write, v: &T) -> Result<()> {
    serde_json::to_writer(&mut *f, v).map_err(|e| e.to_string())?;
    f.write_all(b"\n").map_err(|e| e.to_string())
}
fn seal(
    tx: Transaction,
    p: &Profile,
    run: &Run,
    hashes: BTreeMap<String, String>,
    details: serde_json::Value,
) -> Result<()> {
    for (key, path) in [
        ("input", Some(&run.input)),
        ("db", Some(&run.db)),
        ("profile_file", run.profile.as_ref()),
        ("mouse_trd_reference", run.mouse_trd_reference.as_ref()),
    ] {
        if let Some(path) = path {
            if hashes.get(key) != Some(&hash(path)?) {
                return Err(format!("{key} changed during execution"));
            }
        }
    }
    tx.json("profile.json", p)?;
    let mut outputs = BTreeMap::new();
    for e in std::fs::read_dir(&tx.temporary).map_err(|e| e.to_string())? {
        let e = e.map_err(|e| e.to_string())?;
        outputs.insert(
            e.file_name().to_string_lossy().into_owned(),
            hash(&e.path())?,
        );
    }
    tx.json("manifest.json",&serde_json::json!({"namespace":"ClonoDist Extended","schema":"clonodist.extended.output.v1","profile_id":p.id,"profile_sha256":checksum(p)?,"profile":p,"organism":run.organism,"chains":run.chains,"reference_context":if run.synthetic_demo {"synthetic-demo-unqualified"} else {"pinned-reference"},"inputs":hashes,"outputs":outputs,"details":details,"memory_contract":"conservative managed allocation admission; no hard RSS guarantee"}))?;
    tx.commit()
}
fn model_compatible(m: &Model, p: &Profile, r: &Run, db: &Reference) -> Result<()> {
    m.verify()?;
    if m.profile_sha256 != checksum(p)?
        || m.organism != r.organism
        || m.chains.join(",") != r.chains
    {
        return Err("model/profile/organism/chain mismatch".into());
    }
    for candidate in &m.candidates {
        let center = &candidate.center;
        for (ch, c) in m.chains.iter().zip(&center.chains) {
            let l = &ch[..1];
            let v = center
                .metadata
                .get(&format!("v_{l}_gene"))
                .ok_or("model missing V")?;
            let j = center
                .metadata
                .get(&format!("j_{l}_gene"))
                .ok_or("model missing J")?;
            let aa = center
                .metadata
                .get(&format!("cdr3_{l}_aa"))
                .ok_or("model missing CDR3")?;
            let vg = db.gene(&m.organism, v)?;
            let jg = db.gene(&m.organism, j)?;
            let label = match ch.as_str() {
                "alpha" => "A",
                "beta" => "B",
                "gamma" => "G",
                "delta" => "D",
                _ => return Err("bad model chain".into()),
            };
            if vg.chain != label
                || jg.chain != label
                || vg.region != "V"
                || jg.region != "J"
                || vg.loops.len() < 3
                || c.v != *v
                || c.cdr3 != clonodist_core::encode(aa, true)?
            {
                return Err("model center reference mismatch".into());
            }
            for i in 0..3 {
                if c.loops[i] != clonodist_core::encode(&vg.loops[i], false)? {
                    return Err("model center germline mismatch".into());
                }
            }
        }
    }
    Ok(())
}
fn run() -> Result<()> {
    match Cli::parse().command {
        Command::Profile {
            kind,
            chains,
            output,
        } => {
            let p = match kind.as_str() {
                "pinned" => Profile::pinned(chains),
                "cdr3-focused" => Profile::cdr3(chains),
                _ => return Err("unknown fixed profile".into()),
            };
            p.validate(chains, true)?;
            let tx = Transaction::new(&output)?;
            tx.json("profile.json", &p)?;
            tx.json("identity.json",&serde_json::json!({"profile_id":p.id,"sha256":checksum(&p)?,"canonicalization":"sorted object keys, compact UTF-8 JSON, stable arrays, integer parameters"}))?;
            tx.commit()
        }
        Command::Distance(a) => {
            if !["dense", "edges", "csr"].contains(&a.format.as_str())
                || a.top_k == Some(0)
                || (a.format == "dense"
                    && (a.radius.is_some() || a.top_k.is_some() || a.exclude_self))
            {
                return Err("invalid distance mode/options".into());
            }
            let (db, chains, p, q, mut hashes) = a.run.load()?;
            let t = if let Some(path) = &a.target {
                hashes.insert("target".into(), hash(path)?);
                a.run.read_input(path, &db, &chains)?
            } else {
                q.clone()
            };
            let tx = Transaction::new(&a.run.output)?;
            tx.json("query.json", &q)?;
            tx.json("target.json", &t)?;
            let mut out = writer(
                &tx,
                if a.format == "dense" {
                    "dense.tsv"
                } else if a.format == "edges" {
                    "edges.tsv"
                } else {
                    "csr.indices.tsv"
                },
            )?;
            let mut data = if a.format == "csr" {
                Some(writer(&tx, "csr.data.tsv")?)
            } else {
                None
            };
            let mut ptr = if a.format == "csr" {
                Some(writer(&tx, "csr.indptr.tsv")?)
            } else {
                None
            };
            let mut summary = writer(&tx, "summaries.jsonl")?;
            if a.format == "dense" {
                write!(out, "query_id").map_err(|e| e.to_string())?;
                for r in &t {
                    write!(out, "\t{}", r.id).map_err(|e| e.to_string())?;
                }
                writeln!(out).map_err(|e| e.to_string())?;
            } else if a.format == "edges" {
                writeln!(out, "query_id\ttarget_id\tdistance").map_err(|e| e.to_string())?;
            }
            let mut nnz = 0u64;
            if let Some(f) = &mut ptr {
                writeln!(f, "0").map_err(|e| e.to_string())?;
            }
            let mut config = a.run.config();
            config.radius = a.radius;
            let decimals: Vec<Vec<u8>> = (0..1000u32)
                .map(|v| format!("\t{v}").into_bytes())
                .collect();
            let timing = rows(&q, &t, &p, config, |start, block| {
                for (offset, row) in block.iter().enumerate() {
                    let i = start + offset;
                    if a.format == "dense" {
                        write!(out, "{}", q[i].id).map_err(|e| e.to_string())?;
                        for d in row {
                            if let Some(bytes) = decimals.get(*d as usize) {
                                out.write_all(bytes).map_err(|e| e.to_string())?;
                            } else {
                                write!(out, "\t{d}").map_err(|e| e.to_string())?;
                            }
                        }
                        writeln!(out).map_err(|e| e.to_string())?;
                        continue;
                    }
                    let mut edges: Vec<(usize, u32)> = row
                        .iter()
                        .enumerate()
                        .filter(|(j, d)| {
                            a.radius.is_none_or(|r| **d <= r)
                                && (!a.exclude_self || q[i].id != t[*j].id)
                        })
                        .map(|(j, d)| (j, *d))
                        .collect();
                    edges.sort_by(|(a, x), (b, y)| x.cmp(y).then_with(|| t[*a].id.cmp(&t[*b].id)));
                    if let Some(k) = a.top_k {
                        edges.truncate(k);
                    }
                    let count = edges.iter().try_fold(0u64, |s, (j, _)| {
                        s.checked_add(t[*j].count).ok_or("summary count overflow")
                    })?;
                    let subjects: std::collections::BTreeSet<_> = edges
                        .iter()
                        .filter_map(|(j, _)| t[*j].metadata.get("subject"))
                        .collect();
                    let samples: std::collections::BTreeSet<_> = edges
                        .iter()
                        .filter_map(|(j, _)| t[*j].metadata.get("sample"))
                        .collect();
                    line(
                        &mut summary,
                        &serde_json::json!({"query_id":q[i].id,"neighbors":edges.len(),"weighted_count":count,"subjects":subjects.len(),"samples":samples.len()}),
                    )?;
                    if a.format == "csr" {
                        edges.sort_by_key(|e| e.0);
                    }
                    for (j, d) in edges {
                        if a.format == "edges" {
                            writeln!(out, "{}\t{}\t{d}", q[i].id, t[j].id)
                                .map_err(|e| e.to_string())?;
                        } else {
                            writeln!(out, "{j}").map_err(|e| e.to_string())?;
                            writeln!(data.as_mut().ok_or("CSR data missing")?, "{d}")
                                .map_err(|e| e.to_string())?;
                        }
                        nnz = nnz.checked_add(1).ok_or("nnz overflow")?;
                    }
                    if let Some(f) = &mut ptr {
                        writeln!(f, "{nnz}").map_err(|e| e.to_string())?;
                    }
                }
                Ok(())
            })?;
            finish(out)?;
            finish(summary)?;
            if let Some(f) = data {
                finish(f)?;
            }
            if let Some(f) = ptr {
                finish(f)?;
            }
            if let Some(path) = &a.target {
                if hashes["target"] != hash(path)? {
                    return Err("target changed".into());
                }
            }
            seal(
                tx,
                &p,
                &a.run,
                hashes,
                serde_json::json!({"workflow":"distance","shape":[q.len(),t.len()],"format":a.format,"radius":a.radius,"top_k":a.top_k,"exclude_self":a.exclude_self,"nnz":nnz,"compute_seconds":timing.0,"consume_seconds":timing.1}),
            )
        }
        Command::Discover(a) => {
            if a.run.mouse_trd_reference.is_some() {
                return Err("mouse delta meta inference is not qualified; supplemental reference is admitted for distance workflows only".into());
            }
            let (db, chains, p, fg, mut hashes) = a.run.load()?;
            hashes.insert("discovery_background".into(), hash(&a.background)?);
            hashes.insert("rules".into(), hash(&a.rules)?);
            let bg = read(
                &a.background,
                &db,
                &a.run.organism,
                &chains,
                a.run.memory / 2,
            )?;
            let rules: Rules = load_json(&a.rules, a.run.memory)?;
            let model = discover(
                &fg,
                &bg,
                p.clone(),
                rules,
                a.run.organism.clone(),
                chains,
                a.run.config(),
                hashes.clone(),
            )?;
            if hashes["discovery_background"] != hash(&a.background)?
                || hashes["rules"] != hash(&a.rules)?
            {
                return Err("discovery inputs changed".into());
            }
            let tx = Transaction::new(&a.run.output)?;
            tx.json("model.json", &model)?;
            seal(
                tx,
                &p,
                &a.run,
                hashes,
                serde_json::json!({"workflow":"discovery","candidates":model.candidates.len(),"selected":model.candidates.iter().filter(|c|c.reason=="discovery_selected").count()}),
            )
        }
        Command::Qualify(a) => {
            if a.run.mouse_trd_reference.is_some() {
                return Err("mouse delta meta inference is not qualified; supplemental reference is admitted for distance workflows only".into());
            }
            let (db, chains, p, fg, mut hashes) = a.run.load()?;
            let model: Model = load_json(&a.model, a.run.memory)?;
            model_compatible(&model, &p, &a.run, &db)?;
            hashes.insert("evaluation_background".into(), hash(&a.background)?);
            hashes.insert("discovery_model".into(), hash(&a.model)?);
            let bg = read(
                &a.background,
                &db,
                &a.run.organism,
                &chains,
                a.run.memory / 2,
            )?;
            let model = qualify(model, &fg, &bg, a.run.config(), hashes.clone())?;
            if hashes["evaluation_background"] != hash(&a.background)?
                || hashes["discovery_model"] != hash(&a.model)?
            {
                return Err("evaluation inputs changed".into());
            }
            let tx = Transaction::new(&a.run.output)?;
            tx.json("model.json", &model)?;
            seal(
                tx,
                &p,
                &a.run,
                hashes,
                serde_json::json!({"workflow":"qualification","qualified":model.candidates.iter().filter(|c|c.qualified).count(),"evaluation":model.qualification}),
            )
        }
        Command::Apply(a) => {
            if a.run.mouse_trd_reference.is_some() {
                return Err("mouse delta meta inference is not qualified; supplemental reference is admitted for distance workflows only".into());
            }
            let (db, _, p, rs, mut hashes) = a.run.load()?;
            let model: Model = load_json(&a.model, a.run.memory)?;
            model_compatible(&model, &p, &a.run, &db)?;
            if model.phase != "qualified" {
                return Err("application requires qualified workflow artifact".into());
            }
            hashes.insert("qualified_model".into(), hash(&a.model)?);
            let tx = Transaction::new(&a.run.output)?;
            let mut out = writer(&tx, "matches.jsonl")?;
            for c in &model.candidates {
                for r in &rs {
                    let d = p.score(&c.center, r)?;
                    let (g, m) = constraints(&c.center, r, &model.rules);
                    let neighbor = c.radius.is_some_and(|radius| d <= radius);
                    let direct =
                        signature(&c.center)? == signature(r)? && gene_identity(&c.center, r);
                    let mut reasons = vec![];
                    if !c.qualified {
                        reasons.push("candidate_not_qualified");
                    }
                    if !neighbor {
                        reasons.push("outside_radius");
                    }
                    if !g {
                        reasons.push("gene_constraint");
                    }
                    if !m {
                        reasons.push("motif_constraint");
                    }
                    line(
                        &mut out,
                        &serde_json::json!({"candidate":c.center.id,"record":r.id,"distance":d,"radius":c.radius,"direct_match":direct,"distance_neighbor":neighbor,"gene_compatible":g,"motif_compatible":m,"accepted":reasons.is_empty(),"reasons":reasons}),
                    )?;
                }
            }
            finish(out)?;
            if hashes["qualified_model"] != hash(&a.model)? {
                return Err("model changed".into());
            }
            seal(
                tx,
                &p,
                &a.run,
                hashes,
                serde_json::json!({"workflow":"application"}),
            )
        }
    }
}
fn main() {
    if let Err(e) = run() {
        eprintln!("ClonoDist Extended: {e}");
        std::process::exit(2);
    }
}
