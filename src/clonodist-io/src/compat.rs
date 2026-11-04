//! Finite pinned DataFrame normalization adapter. Not the strict production reader.
use clonodist_core::{distance, Receptor, Result};
use clonodist_germline::Reference;
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, BTreeSet};
fn warning(s: String) -> Value {
    json!({"category":"UserWarning","message":s})
}
fn err(out: &mut Value, kind: &str, message: String) {
    out["error"] = json!({"type":kind,"message":message});
}
fn cmp(a: &Value, b: &Value) -> std::cmp::Ordering {
    match (a.as_f64(), b.as_f64()) {
        (Some(a), Some(b)) => a.total_cmp(&b),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        _ => a.as_str().unwrap_or("").cmp(b.as_str().unwrap_or("")),
    }
}
pub fn normalize(case: &Value, db: &Reference) -> Result<Value> {
    let mut out = json!({"warnings":[]});
    let k = &case["kwargs"];
    let organism = k["organism"].as_str().unwrap_or("human");
    if !["human", "mouse"].contains(&organism) {
        err(
            &mut out,
            "ValueError",
            "organism must be 'mouse', 'human', or 'rhesus'".into(),
        );
        return Ok(out);
    }
    let chains: Vec<String> = k["chains"]
        .as_array()
        .map(|a| {
            a.iter()
                .map(|v| v.as_str().unwrap_or("").to_owned())
                .collect()
        })
        .unwrap_or(vec!["beta".into()]);
    if chains
        .iter()
        .any(|s| !["alpha", "beta"].contains(&s.as_str()))
    {
        err(&mut out,"ValueError","TCRrep chains arg can be one or more of the following ['alpha', 'beta', 'gamma', 'delta'] case-sensitive".into());
        return Ok(out);
    }
    if !k["imgt_aligned"].is_null() && !k["imgt_aligned"].is_boolean() {
        err(
            &mut out,
            "ValueError",
            "TCRrep imgt_aligned argument must be a boolean".into(),
        );
        return Ok(out);
    }
    let mut rows: Vec<Map<String, Value>> = case["rows"]
        .as_array()
        .ok_or("rows must be array")?
        .iter()
        .map(|r| {
            r.as_object()
                .cloned()
                .ok_or("row must be object".to_string())
        })
        .collect::<Result<_>>()?;
    let mut columns: BTreeSet<String> = rows.iter().flat_map(|r| r.keys().cloned()).collect();
    let mut warnings = vec![];
    for ch in &chains {
        let l = &ch[..1];
        for (col,msg) in [(format!("cdr3_{l}_aa"),format!("cell_df needs a column called 'cdr3_{l}_aa' to track the CDR3 amino acid sequence\n")),(format!("v_{l}_gene"),format!("cell_df needs a column called 'v_{l}_gene' for default functions\n"))]{if !columns.contains(&col){warnings.push(warning(msg));}}
    }
    if !columns.contains("count") {
        warnings.push(warning(
            "cell_df needs a counts column to track clonal number of frequency\n".into(),
        ));
        warnings.push(warning(
            "No 'count' column provided; count column set to 1".into(),
        ));
        for r in &mut rows {
            r.insert("count".into(), json!(1));
        }
        columns.insert("count".into());
    }
    for row in &mut rows {
        for col in &columns {
            row.entry(col.clone()).or_insert(Value::Null);
        }
    }
    let mut inferred = vec![];
    for ch in &chains {
        let l = &ch[..1];
        let vcol = format!("v_{l}_gene");
        if !columns.contains(&vcol) {
            out["warnings"] = json!(warnings);
            err(
                &mut out,
                "AttributeError",
                format!("'DataFrame' object has no attribute '{vcol}'"),
            );
            return Ok(out);
        }
        for (idx, region) in ["cdr1", "cdr2", "pmhc"].iter().enumerate() {
            let col = format!("{region}_{l}_aa");
            inferred.push(col.clone());
            for row in &mut rows {
                let v = &row[&vcol];
                let g = v.as_str().and_then(|s| db.gene(organism, s).ok());
                let s = g.and_then(|g| g.loops.get(idx));
                if s.is_none() {
                    let name = v.as_str().unwrap_or("nan");
                    warnings.push(warning(format!("{name} gene was not recognized in reference db no cdr seq could be inferred")));
                }
                row.insert(col.clone(), s.map_or(Value::Null, |s| json!(s)));
            }
        }
    }
    let initial = columns
        .iter()
        .filter(|s| !["count", "cell_id", "clone_id"].contains(&s.as_str()))
        .cloned()
        .collect::<Vec<_>>();
    let mut index = initial.clone();
    index.extend(inferred);
    for ch in &chains {
        index.push(format!("cdr3_{}_aa", &ch[..1]));
    }
    let input_count: f64 = rows.iter().map(|r| r["count"].as_f64().unwrap_or(0.)).sum();
    if k["deduplicate"] != false {
        let mut groups: Vec<Map<String, Value>> = vec![];
        for r in rows {
            if index
                .iter()
                .filter(|s| r.contains_key(*s))
                .any(|s| r[s].is_null())
            {
                continue;
            }
            let keys: Vec<_> = index.iter().filter(|s| r.contains_key(*s)).collect();
            if let Some(g) = groups
                .iter_mut()
                .find(|g| keys.iter().all(|s| g.get(*s) == r.get(*s)))
            {
                let n = g["count"].as_f64().unwrap_or(0.) + r["count"].as_f64().unwrap_or(0.);
                g.insert("count".into(), number(n));
            } else {
                let mut g = r;
                g.remove("cell_id");
                g.remove("clone_id");
                g.insert("count".into(), number(g["count"].as_f64().unwrap_or(0.)));
                groups.push(g);
            }
        }
        groups.sort_by(|a, b| {
            for s in &index {
                let order = cmp(
                    a.get(s).unwrap_or(&Value::Null),
                    b.get(s).unwrap_or(&Value::Null),
                );
                if !order.is_eq() {
                    return order;
                }
            }
            std::cmp::Ordering::Equal
        });
        for (i, r) in groups.iter_mut().enumerate() {
            r.insert("clone_id".into(), json!(i + 1));
        }
        rows = groups;
    }
    let retained: f64 = rows.iter().map(|r| r["count"].as_f64().unwrap_or(0.)).sum();
    if input_count != retained {
        warnings.push(warning(format!("Not all cells/sequences could be grouped into clones.{} of {} were not captured. This occurs when any of the values in the index columns are null or missing for a given sequence. To see entries with missing values use: tcrdist.repertoire.TCRrep._show_incomplete()",input_count-retained,input_count)));
    }
    out["warnings"] = json!(warnings);
    out["index_cols"] = json!(index);
    out["input_count"] = json!(input_count);
    out["retained_count"] = json!(retained);
    out["normalized"] = json!(rows);
    if rows.is_empty() {
        err(
            &mut out,
            "ValueError",
            "zero-size array to reduction operation maximum which has no identity".into(),
        );
        return Ok(out);
    }
    let mut dense = Map::new();
    for ch in &chains {
        let l = &ch[..1];
        let ccol = format!("cdr3_{l}_aa");
        if !rows[0].contains_key(&ccol) {
            err(&mut out, "KeyError", format!("'{ccol}'"));
            return Ok(out);
        }
        if rows.len() > 1 && rows.iter().any(|r| !r[&ccol].is_string()) {
            err(
                &mut out,
                "TypeError",
                "'<' not supported between instances of 'str' and 'int'".into(),
            );
            return Ok(out);
        }
        let mut records = vec![];
        for r in &rows {
            let chain = db.resolve(
                organism,
                ch,
                r[&format!("v_{l}_gene")].as_str().ok_or("V type")?,
                None,
                r[&ccol].as_str().unwrap_or(""),
                false,
            )?;
            records.push(Receptor {
                id: String::new(),
                count: 1,
                metadata: BTreeMap::new(),
                chains: vec![chain],
            });
        }
        let mut matrix = vec![];
        for (i, a) in records.iter().enumerate() {
            let mut row = vec![];
            for (j, b) in records.iter().enumerate() {
                row.push(if i == j { 0 } else { distance(a, b)? });
            }
            matrix.push(row);
        }
        dense.insert(ch.clone(), json!(matrix));
    }
    out["dense"] = json!(dense);
    Ok(out)
}
fn number(n: f64) -> Value {
    if n.fract() == 0. && n >= i64::MIN as f64 && n <= i64::MAX as f64 {
        json!(n as i64)
    } else {
        json!(n)
    }
}
