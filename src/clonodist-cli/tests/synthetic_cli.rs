use std::path::PathBuf;
use std::process::Command;

#[test]
fn synthetic_distance_requires_explicit_opt_in_and_marks_output() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let out = std::env::temp_dir().join(format!("clonodist-cli-synthetic-{}", std::process::id()));
    assert!(!out.exists());
    let run = |synthetic: bool| {
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_clonodist"));
        cmd.current_dir(&root).args([
            "distance",
            "--db",
            "examples/synthetic-reference.txt",
            "--input",
            "examples/synthetic-receptors.tsv",
            "--chains",
            "beta",
            "--output",
            out.to_str().unwrap(),
        ]);
        if synthetic {
            cmd.arg("--synthetic-demo");
        }
        cmd.output().unwrap()
    };
    assert!(!run(false).status.success());
    assert!(!out.exists());
    let result = run(true);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(String::from_utf8_lossy(&result.stderr).contains("SYNTHETIC DEMONSTRATION ONLY"));
    assert_eq!(
        std::fs::read(out.join("matrix.tsv")).unwrap(),
        std::fs::read(root.join("examples/expected/single-chain-matrix.tsv")).unwrap()
    );
    let manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(out.join("manifest.json")).unwrap()).unwrap();
    assert_eq!(manifest["reference_context"], "synthetic-demo-unqualified");
    std::fs::remove_dir_all(out).unwrap();
}
