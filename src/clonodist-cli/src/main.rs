use clap::{Args, Parser, Subcommand};
use clonodist_core::{self as core, Result};
use clonodist_distance::{Config, Engine};
use clonodist_germline::Reference;
use clonodist_io::{hash, read_tsv, Transaction};
use serde_json::{json, Value};
use std::{
    io::{BufWriter, Write},
    path::PathBuf,
    time::Instant,
};
#[derive(Parser)]
#[command(
    version,
    about = "Exact pinned TCR distances and descriptive neighborhoods"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}
#[derive(Subcommand)]
enum Command {
    Validate(Run),
    Distance(Run),
    Join(Run),
    Neighbors(Run),
    ExplainPair(Explain),
    /// JSON-lines scientific adapter for qualification (not strict repertoire input).
    #[command(hide = true)]
    Batch {
        #[arg(long)]
        db: PathBuf,
    },
}
#[derive(Args)]
struct Run {
    #[arg(long)]
    input: PathBuf,
    #[arg(long)]
    reference: Option<PathBuf>,
    #[arg(long)]
    db: PathBuf,
    /// Use only the included invented demonstration reference; never use for real data.
    #[arg(long)]
    synthetic_demo: bool,
    #[arg(long, default_value = "human")]
    organism: String,
    #[arg(long, value_delimiter = ',', default_value = "beta")]
    chains: Vec<String>,
    #[arg(long)]
    output: PathBuf,
    #[arg(long, default_value_t = 1)]
    threads: usize,
    #[arg(long, default_value_t = 32)]
    tile: usize,
    /// Managed memory budget in bytes; conservative allocation planning, not an OS RSS limit.
    #[arg(long, default_value_t = 536870912)]
    memory_limit: usize,
    #[arg(long)]
    radius: Option<u32>,
    #[arg(long)]
    top_k: Option<usize>,
    #[arg(long)]
    exclude_self: bool,
    #[arg(long)]
    deduplicate: bool,
    #[arg(long)]
    scalar: bool,
    #[arg(long,default_value="dense-tsv",value_parser=["dense-tsv","edges-tsv","csr"])]
    format: String,
}
#[derive(Args)]
struct Explain {
    #[command(flatten)]
    run: Run,
    #[arg(long)]
    id_a: String,
    #[arg(long)]
    id_b: String,
}
fn batch(db: PathBuf) -> Result<()> {
    use std::io::BufRead;
    let db = Reference::load(&db)?;
    let stdin = std::io::stdin();
    let mut out = BufWriter::new(std::io::stdout().lock());
    for line in stdin.lock().lines() {
        let v: Value =
            serde_json::from_str(&line.map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
        let answer = (|| -> Result<Value> {
            match v["op"].as_str().unwrap_or("region") {
                "region" => {
                    let a = core::encode(v["left"].as_str().ok_or("left")?, false)?;
                    let b = core::encode(v["right"].as_str().ok_or("right")?, false)?;
                    let p = if v["parameters"].is_null() {
                        core::CDR3
                    } else {
                        serde_json::from_value(v["parameters"].clone())
                            .map_err(|e| e.to_string())?
                    };
                    Ok(json!({"distance":core::region(&a,&b,p)?}))
                }
                "pair" => {
                    let species = v["organism"].as_str().ok_or("organism")?;
                    let mut breakdown = serde_json::Map::new();
                    let mut total = 0;
                    for c in v["chains"].as_array().ok_or("chains")? {
                        let c = c.as_str().ok_or("chain")?;
                        if !["alpha", "beta"].contains(&c) {
                            return Err("unsupported chain".into());
                        }
                        let l = &c[..1];
                        let resolve = |r: &Value| {
                            db.resolve(
                                species,
                                c,
                                r[format!("v_{l}_gene")].as_str().ok_or("v")?,
                                None,
                                r[format!("cdr3_{l}_aa")].as_str().ok_or("cdr3")?,
                                false,
                            )
                        };
                        let a = resolve(&v["rows"][0])?;
                        let b = resolve(&v["rows"][1])?;
                        let d = core::chain(&a, &b)?;
                        total = core::add(total, d.total)?;
                        breakdown.insert(c.into(), json!(d));
                    }
                    Ok(json!({"chains":breakdown,"total":total}))
                }
                "normalize" => clonodist_io::compat::normalize(&v["case"], &db),
                _ => Err("unknown batch op".into()),
            }
        })();
        let answer = match answer {
            Ok(a) => a,
            Err(e) => json!({"rust_error":e}),
        };
        serde_json::to_writer(&mut out, &answer).map_err(|e| e.to_string())?;
        out.write_all(b"\n").map_err(|e| e.to_string())?;
        out.flush().map_err(|e| e.to_string())?;
    }
    Ok(())
}
// Bounded immutable formatting table; full u32 fallback preserves every value.
const fn decimal_table() -> [([u8; 4], usize); 1000] {
    let mut table = [([0; 4], 0); 1000];
    let mut n = 0;
    while n < 1000 {
        let mut bytes = [b'\t', 0, 0, 0];
        let len = if n < 10 {
            bytes[1] = b'0' + n as u8;
            2
        } else if n < 100 {
            bytes[1] = b'0' + (n / 10) as u8;
            bytes[2] = b'0' + (n % 10) as u8;
            3
        } else {
            bytes[1] = b'0' + (n / 100) as u8;
            bytes[2] = b'0' + (n / 10 % 10) as u8;
            bytes[3] = b'0' + (n % 10) as u8;
            4
        };
        table[n] = (bytes, len);
        n += 1;
    }
    table
}
const DECIMALS: [([u8; 4], usize); 1000] = decimal_table();
fn run(args: Run, kind: &str, explain: Option<(&str, &str)>) -> Result<()> {
    let start = Instant::now();
    if args.memory_limit < 64 * 1024 * 1024 {
        return Err("memory-limit must be at least 64 MiB".into());
    }
    if args.chains.is_empty()
        || args.chains.len() > 2
        || args
            .chains
            .iter()
            .any(|s| !["alpha", "beta"].contains(&s.as_str()))
        || args.chains.windows(2).any(|w| w[0] >= w[1])
    {
        return Err("chains must be alpha, beta, or alpha,beta".into());
    }
    if kind == "neighbors" && args.radius.is_none() && args.top_k.is_none() {
        return Err("neighbors requires --radius or --top-k".into());
    }
    if kind == "join" && args.reference.is_none() {
        return Err("join requires --reference".into());
    }
    if args.top_k == Some(0) {
        return Err("top-k must be positive".into());
    }
    if args.exclude_self && args.reference.is_some() {
        return Err("exclude-self is only supported within a repertoire".into());
    }
    let input_hash = hash(&args.input)?;
    let reference_hash = args.reference.as_ref().map(|p| hash(p)).transpose()?;
    let db_clock = Instant::now();
    let db = if args.synthetic_demo {
        eprintln!("SYNTHETIC DEMONSTRATION ONLY: invented reference, not scientifically qualified");
        Reference::load_synthetic_demo(&args.db)?
    } else {
        Reference::load(&args.db)?
    };
    let reference_seconds = db_clock.elapsed().as_secs_f64();
    let (q, ledger) = read_tsv(
        &args.input,
        &db,
        &args.organism,
        &args.chains,
        args.deduplicate,
        args.memory_limit / 2,
    )?;
    let (target, target_ledger) = if let Some(p) = &args.reference {
        let (t, l) = read_tsv(
            p,
            &db,
            &args.organism,
            &args.chains,
            args.deduplicate,
            args.memory_limit / 2,
        )?;
        (Some(t), Some(l))
    } else {
        (None, None)
    };
    let t = target.as_deref().unwrap_or(&q);
    let input_seconds = start.elapsed().as_secs_f64();
    let tx = Transaction::new(&args.output)?;
    tx.json("records.json", &q)?;
    tx.json("ledger.json", &ledger)?;
    if let Some(l) = target_ledger {
        tx.json("target-ledger.json", &l)?;
        tx.json("target-records.json", &t)?;
    }
    let mut comparisons = 0u64;
    let mut edges = 0u64;
    let mut kernel_seconds = 0.;
    let mut compute_seconds = 0.;
    let mut serialization_seconds = 0.;
    let mut engine_setup_seconds = 0.;
    if let Some((id_a, id_b)) = explain {
        let a = q.iter().find(|r| r.id == id_a).ok_or("id-a not found")?;
        let b = t.iter().find(|r| r.id == id_b).ok_or("id-b not found")?;
        let breakdown: Result<Vec<_>> = a
            .chains
            .iter()
            .zip(&b.chains)
            .map(|(a, b)| core::chain(a, b))
            .collect();
        tx.json(
            "breakdown.json",
            &json!({"chains":breakdown?,"total":core::distance(a,b)?}),
        )?;
    } else if kind != "validate" {
        let config = Config {
            threads: args.threads,
            tile: args.tile,
            memory: args.memory_limit,
            optimized: !args.scalar,
            radius: args.radius,
        };
        let setup_clock = Instant::now();
        let engine = Engine::new(&q, t, config)?;
        engine_setup_seconds = setup_clock.elapsed().as_secs_f64();
        let filename = match args.format.as_str() {
            "dense-tsv" => "matrix.tsv",
            "csr" => "csr.tsv",
            _ => "edges.tsv",
        };
        if args.format == "dense-tsv"
            && (args.radius.is_some() || args.top_k.is_some() || args.exclude_self)
        {
            return Err("radius/top-k require edges-tsv or csr".into());
        }
        let mut out = BufWriter::new(
            std::fs::File::create(tx.temporary.join(filename)).map_err(|e| e.to_string())?,
        );
        let mut summary = BufWriter::new(
            std::fs::File::create(tx.temporary.join("neighborhoods.tsv"))
                .map_err(|e| e.to_string())?,
        );
        writeln!(
            summary,
            "query_id\tneighbors\tweighted_count\tsubjects\tsamples"
        )
        .map_err(|e| e.to_string())?;
        if args.format == "dense-tsv" {
            write!(out, "query_id").map_err(|e| e.to_string())?;
            for r in t {
                write!(out, "\t{}", r.id).map_err(|e| e.to_string())?;
            }
            writeln!(out).map_err(|e| e.to_string())?;
        } else if args.format == "csr" {
            writeln!(out, "row\tindptr_start\tindptr_end\tindices\tdistances")
                .map_err(|e| e.to_string())?;
        } else {
            writeln!(out, "query_id\ttarget_id\tdistance").map_err(|e| e.to_string())?;
        }
        let clock = Instant::now();
        let phase = engine.run_profile(|start, rows| {
            for (offset, row) in rows.iter().enumerate() {
                let i = start + offset;
                comparisons = comparisons
                    .checked_add(t.len() as u64)
                    .ok_or("comparison counter overflow")?;
                if args.format == "dense-tsv" {
                    write!(out, "{}", q[i].id).map_err(|e| e.to_string())?;
                    for d in row {
                        if let Some((bytes, len)) = DECIMALS.get(*d as usize) {
                            out.write_all(&bytes[..*len]).map_err(|e| e.to_string())?;
                        } else {
                            write!(out, "\t{d}").map_err(|e| e.to_string())?;
                        }
                    }
                    writeln!(out).map_err(|e| e.to_string())?;
                } else {
                    let mut e = clonodist_neighbors::ranked_unique(
                        row,
                        t,
                        args.radius.unwrap_or(u32::MAX - 1),
                        args.top_k,
                        if args.exclude_self {
                            Some(q[i].id.as_str())
                        } else {
                            None
                        },
                    );
                    let s = clonodist_neighbors::summarize(&e, t)?;
                    writeln!(
                        summary,
                        "{}\t{}\t{}\t{}\t{}",
                        q[i].id, s.neighbors, s.weighted_count, s.subjects, s.samples
                    )
                    .map_err(|e| e.to_string())?;
                    let edge_start = edges;
                    edges = edges
                        .checked_add(e.len() as u64)
                        .ok_or("edge counter overflow")?;
                    if args.format == "csr" {
                        e.sort_unstable_by_key(|x| x.0);
                        let indices = e
                            .iter()
                            .map(|x| x.0.to_string())
                            .collect::<Vec<_>>()
                            .join(",");
                        let distances = e
                            .iter()
                            .map(|x| x.1.to_string())
                            .collect::<Vec<_>>()
                            .join(",");
                        writeln!(out, "{i}\t{edge_start}\t{edges}\t{indices}\t{distances}")
                            .map_err(|e| e.to_string())?;
                    } else {
                        for (j, d) in e {
                            writeln!(out, "{}\t{}\t{d}", q[i].id, t[j].id)
                                .map_err(|e| e.to_string())?;
                        }
                    }
                }
            }
            Ok(())
        })?;
        compute_seconds = phase.0;
        serialization_seconds = phase.1;
        kernel_seconds = clock.elapsed().as_secs_f64();
        out.flush().map_err(|e| e.to_string())?;
        out.get_ref().sync_all().map_err(|e| e.to_string())?;
        summary.flush().map_err(|e| e.to_string())?;
        summary.get_ref().sync_all().map_err(|e| e.to_string())?;
    }
    let mut outputs = serde_json::Map::new();
    for entry in std::fs::read_dir(&tx.temporary).map_err(|e| e.to_string())? {
        let p = entry.map_err(|e| e.to_string())?.path();
        outputs.insert(p.file_name().unwrap().to_string_lossy().into(),json!({"sha256":hash(&p)?,"bytes":std::fs::metadata(&p).map_err(|e|e.to_string())?.len()}));
    }
    if hash(&args.input)? != input_hash
        || args.reference.as_ref().map(|p| hash(p)).transpose()? != reference_hash
    {
        return Err("input changed during run".into());
    }
    let manifest = json!({"schema_version":1,"version":env!("CARGO_PKG_VERSION"),"profile":core::PROFILE,"git_commit":option_env!("CLONODIST_GIT_SHA").unwrap_or("unrecorded-build"),"command":std::env::args().collect::<Vec<_>>(),"organism":args.organism,"chains":args.chains,"input_sha256":input_hash,"reference_input_sha256":reference_hash,"database_sha256":db.sha256,"reference_context":if args.synthetic_demo {"synthetic-demo-unqualified"} else {"pinned-reference"},"threads":args.threads,"tile":args.tile,"memory_budget_bytes":args.memory_limit,"optimized":!args.scalar,"input_seconds":input_seconds,"reference_seconds":reference_seconds,"engine_setup_seconds":engine_setup_seconds,"compute_seconds":compute_seconds,"serialization_seconds":serialization_seconds,"kernel_and_write_seconds":kernel_seconds,"elapsed_seconds":start.elapsed().as_secs_f64(),"comparisons":comparisons,"retained_edges":edges,"query_records":q.len(),"target_records":t.len(),"radius":args.radius,"top_k":args.top_k,"exclude_self":args.exclude_self,"outputs":outputs,"platform":std::env::consts::OS,"architecture":std::env::consts::ARCH});
    tx.json("manifest.json", &manifest)?;
    tx.commit()?;
    Ok(())
}
fn main() {
    let answer = match Cli::parse().command {
        Command::Batch { db } => batch(db),
        Command::Validate(a) => run(a, "validate", None),
        Command::Distance(a) => run(a, "distance", None),
        Command::Join(a) => run(a, "join", None),
        Command::Neighbors(a) => run(a, "neighbors", None),
        Command::ExplainPair(a) => run(a.run, "explain", Some((&a.id_a, &a.id_b))),
    };
    if let Err(e) = answer {
        eprintln!("clonodist: {e}");
        std::process::exit(2);
    }
}

#[cfg(test)]
mod formatting_tests {
    use super::DECIMALS;
    #[test]
    fn every_preformatted_value_is_canonical_decimal() {
        for (n, (bytes, len)) in DECIMALS.iter().enumerate() {
            assert_eq!(&bytes[..*len], format!("\t{n}").as_bytes());
        }
    }
}
