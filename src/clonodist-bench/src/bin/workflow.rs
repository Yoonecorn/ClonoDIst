//! Warm complete comparison traversal harness; output identity is checked externally.
use clonodist_core::{Receptor, Result};
use clonodist_distance::{Config, Engine};
use serde_json::json;
use std::{
    io::{BufWriter, Write},
    path::PathBuf,
    time::Instant,
};
fn execute() -> Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.len() != 8 {
        return Err("usage: workflow query.json target.json threads tile optimized radius-or-all repetitions output.bin".into());
    }
    let q: Vec<Receptor> =
        serde_json::from_slice(&std::fs::read(&args[0]).map_err(|e| e.to_string())?)
            .map_err(|e| e.to_string())?;
    let t: Vec<Receptor> =
        serde_json::from_slice(&std::fs::read(&args[1]).map_err(|e| e.to_string())?)
            .map_err(|e| e.to_string())?;
    let config = Config {
        threads: args[2].parse().map_err(|_| "threads")?,
        tile: args[3].parse().map_err(|_| "tile")?,
        optimized: args[4] == "true",
        radius: args[5].parse().ok(),
        memory: 1024 * 1024 * 1024,
    };
    let start = Instant::now();
    let engine = Engine::new(&q, &t, config)?;
    let setup = start.elapsed().as_secs_f64();
    let reps: usize = args[6].parse().map_err(|_| "repetitions")?;
    let mut times = vec![];
    let mut counts = vec![];
    for rep in 0..reps {
        let clock = Instant::now();
        let mut sum = 0u128;
        let mut retained = 0u64;
        engine.run(|start, rows| {
            for (i, row) in rows.iter().enumerate() {
                for (j, d) in row.iter().enumerate() {
                    if engine.config.radius.is_none_or(|r| *d <= r) {
                        retained += 1;
                        sum = sum.wrapping_add(
                            u128::from(*d)
                                + ((start + i) as u128) * 1000000007
                                + (j as u128) * 1000000009,
                        );
                    }
                }
            }
            Ok(())
        })?;
        times.push(clock.elapsed().as_secs_f64());
        counts.push(json!({"checksum":sum.to_string(),"retained":retained}));
        if rep == 0 && args[7] != "-" {
            let path = PathBuf::from(&args[7]);
            let mut out = BufWriter::new(std::fs::File::create(path).map_err(|e| e.to_string())?);
            engine.run(|_, rows| {
                for row in rows {
                    for d in row {
                        let v = if engine.config.radius.is_some_and(|r| *d > r) {
                            u32::MAX
                        } else {
                            *d
                        };
                        out.write_all(&v.to_le_bytes()).map_err(|e| e.to_string())?;
                    }
                }
                Ok(())
            })?;
            out.flush().map_err(|e| e.to_string())?;
        }
    }
    println!(
        "{}",
        json!({"setup_seconds":setup,"times_seconds":times,"checks":counts,"pairs":q.len()*t.len()})
    );
    Ok(())
}
fn main() {
    if let Err(e) = execute() {
        eprintln!("{e}");
        std::process::exit(2);
    }
}
