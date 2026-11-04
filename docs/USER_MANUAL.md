# ClonoDist user manual

Version 1.0.0 · Yoon Sing Yap · [yoonatworkz@gmail.com](mailto:yoonatworkz@gmail.com)

## 1. What ClonoDist does

ClonoDist calculates distances between T-cell receptors. It can also find receptors within a chosen distance of one another. It records the settings, input identities and reference context used in a run. This helps you check and repeat your work.

ClonoDist uses the established tcrdist3/pwseqdist scoring method. **It does not introduce a new receptor-distance formula.**

## 2. What ClonoDist does not do

ClonoDist does not predict an antigen or infer HLA. Distance similarity does not prove shared antigen specificity. Operational neighborhood matches are not sensitivity, specificity, HLA inference or clinical validation. They do not establish causal biology or population error rates.

## 3. Core and Extended

**Core** provides alpha and beta distances, paired alpha/beta distances, joins and neighborhoods. **Extended** adds defined gamma/delta contexts, explicit profiles, and separate discovery, qualification and application steps for operational neighborhood analysis. Both use the same scoring foundation.

The examples below use invented receptors. **Synthetic examples test software behavior only.** They do not reproduce the scientific qualification campaign or prove biological validity.

## Which workflow should I use?

| Your goal | Program | Workflow |
|---|---|---|
| Alpha or beta distance | `clonodist` | Core single-chain distance |
| Paired alpha/beta distance | `clonodist` | Core paired distance |
| Profile-aware distance | `clonodist-extended` | Extended distance |
| Gamma/delta distance | `clonodist-extended` | Extended gamma/delta distance |
| Radius neighborhood | `clonodist` | Core radius-neighborhood calculation |
| Operational discovery and matching | `clonodist-extended` | Discovery, qualification and application |

The v1.0 manual demonstrates radius-based neighborhood calculation. The operational discovery, qualification, and application workflow requires a prespecified study design and governed inputs. Read the command help and reproducibility guidance before using that workflow.

## 4. System requirements

You need a computer, a terminal, Rust and Cargo, and a C compiler with build tools. You need Git only if you choose to clone the repository. You need Python 3 for the optional mouse TRD preparation script. Allow disk space for the software and your outputs. Dense distance tables can grow with the square of the number of records.

The commands below use a macOS or Linux-style shell. The scientific qualification was performed in the documented local Apple Silicon environment. Windows is not scientifically qualified for v1.0, and no Windows workflow is validated here. A successful build on another system does not qualify that system scientifically.

## 5. Install Rust

1. Open the official [Rust installation page](https://www.rust-lang.org/tools/install/).
2. Follow its instructions for your operating system.
3. Install the C build tools for your system using the instructions below.
4. Open a new terminal.
5. Check Rust:

   ```bash
   rustc --version
   ```

6. Check Cargo:

   ```bash
   cargo --version
   ```

Both commands should print a version. If a command is not found, reopen the terminal. If it still fails, follow the installation page's PATH guidance.

### Install C build tools

On macOS, run:

```bash
xcode-select --install
```

On Ubuntu or Debian, run:

```bash
sudo apt update
sudo apt install build-essential
```

On Fedora, run:

```bash
sudo dnf group install development-tools
```

Follow the prompts. These commands install system software and may ask for your password. Windows is outside the qualified v1.0 environment; no Windows installation workflow is validated here.

## 6. Download and build ClonoDist

Repository: [https://github.com/Yoonecorn/ClonoDIst](https://github.com/Yoonecorn/ClonoDIst).

Pre-release note: The repository is private during release preparation. Authorized users must sign in to GitHub. After public release, no repository invitation is required.

1. Open a terminal in the folder where you want to place ClonoDist.
2. Clone the repository:

   ```bash
   git clone https://github.com/Yoonecorn/ClonoDIst.git
   ```

3. Enter the folder containing `Cargo.toml`:

   ```bash
   cd ClonoDIst
   ```

4. Build the workspace:

   ```bash
   cargo build --workspace --locked
   ```

The first build needs a network connection to fetch Cargo dependencies. `--locked` uses the versions in `Cargo.lock`. If dependencies are already cached and you have no network, add `--offline` before running the command. The examples below assume Cargo uses the default `target` folder. If you set a custom `CARGO_TARGET_DIR`, use that folder instead of `target` in the commands.

## 7. Check the installation

Run these commands from the folder containing `Cargo.toml`:

```bash
target/debug/clonodist --version
target/debug/clonodist-extended --version
```

Expected output:

```text
clonodist-cli 1.0.0
clonodist-extended 1.0.0
```

To see available commands, run:

```bash
target/debug/clonodist --help
target/debug/clonodist-extended --help
```

If a program is not found, repeat the build in Section 6 and check the target folder. For Windows without a Unix-style shell, executable names and commands may differ.

## 8. Prepare an input file

Use a UTF-8 tab-separated file. Each row represents one receptor record. The included [synthetic input](../examples/synthetic-receptors.tsv) contains two invented records. This excerpt shows its beta-chain fields; the complete file also has alpha, gamma and delta fields:

```tsv
id	count	v_b_gene	j_b_gene	cdr3_b_aa
demo-1	1	FAKE_B_V*01	FAKE_B_J*01	AAAACAAAAA
demo-2	1	FAKE_B_V*01	FAKE_B_J*01	AAAADAAAAA
```

The example above uses real tab separators. A real input file must use tabs, not spaces. Keep the header row. Give each record a unique `id`. A supplied `count` must be a positive whole number. Record what your counts mean; ClonoDist does not infer their units.

For a selected chain, provide its exact V allele, J allele and uppercase CDR3 amino-acid sequence. Use suffix `a` for alpha, `b` for beta, `g` for gamma and `d` for delta. For example, beta uses `v_b_gene`, `j_b_gene` and `cdr3_b_aa`. A paired input needs all three fields for both chains. CDR3 must contain 1–127 canonical amino acids. Optional `subject` and `sample` fields can identify the source of a record.

The `FAKE` allele names work only with the included demonstration reference. For real data, use exact alleles present in the separately acquired reference. ClonoDist does not guess or correct allele names.

## 9. Five-minute synthetic example

Start in a terminal outside the repository folder. This example uses only invented inputs. It tests software behavior, not biological validity. A first build can take longer than five minutes.

1. Clone the repository:

   ```bash
   git clone https://github.com/Yoonecorn/ClonoDIst.git
   ```

2. Enter the repository:

   ```bash
   cd ClonoDIst
   ```

3. Build the programs:

   ```bash
   cargo build --workspace --locked
   ```

4. Check the Core program version:

   ```bash
   target/debug/clonodist --version
   ```

   Expected terminal output:

   ```text
   clonodist-cli 1.0.0
   ```

5. Create an output folder:

   ```bash
   mkdir -p output
   ```

6. Calculate the synthetic beta-chain distances:

   ```bash
   target/debug/clonodist distance \
     --synthetic-demo \
     --db examples/synthetic-reference.txt \
     --input examples/synthetic-receptors.tsv \
     --chains beta \
     --output output/quick-start
   ```

7. Open `output/quick-start/matrix.tsv` in a text editor. The result is tab-separated:

   ```tsv
   query_id	demo-1	demo-2
   demo-1	0	12
   demo-2	12	0
   ```

The distance between the invented records `demo-1` and `demo-2` is 12 under this demonstration. It is not evidence of shared antigen specificity. The command prints a synthetic-only warning and marks its manifest `synthetic-demo-unqualified`. Never use `--synthetic-demo` for real receptors.

Sections 10–14 provide more examples. Each command writes a **new** output directory. If one exists, choose a different `--output` path.

## 10. Core single-chain distance

**Purpose:** compare the two invented beta chains in Section 8. The input is `examples/synthetic-receptors.tsv`; its beta columns are shown above. The example is synthetic and is not scientifically qualified.

Run:

```bash
target/debug/clonodist distance \
  --synthetic-demo \
  --db examples/synthetic-reference.txt \
  --input examples/synthetic-receptors.tsv \
  --chains beta \
  --output output/core
```

`distance` requests a dense matrix. `--db` identifies the demonstration reference. `--input` identifies the receptor file. `--chains beta` selects the beta fields. `--output` names the new result directory.

Open `output/core/matrix.tsv`. It should match [the expected single-chain matrix](../examples/expected/single-chain-matrix.tsv):

```tsv
query_id	demo-1	demo-2
demo-1	0	12
demo-2	12	0
```

`query_id` names each row. The other column headers name the target records. A cell is the distance between its row and column. The diagonal is zero for this example. The distance between the two invented records is 12.

**Trouble:** if `output/core` already exists, use a new directory name. If the demonstration reference is refused, check both `--synthetic-demo` and the `--db` path.

## 11. Extended single-chain distance

**Purpose:** run the profile-aware Extended program on the same two invented beta chains and beta fields shown in Section 8. This is a synthetic behavior check, not scientific qualification.

```bash
target/debug/clonodist-extended distance \
  --synthetic-demo \
  --db examples/synthetic-reference.txt \
  --input examples/synthetic-receptors.tsv \
  --chains beta \
  --output output/extended
```

The options have the same roles as in Section 10. The pinned profile is used by default. Open `output/extended/dense.tsv`. It has the same rows, columns and values as the [expected single-chain matrix](../examples/expected/single-chain-matrix.tsv): `demo-1` versus `demo-2` is 12. The row ID is the query and each other header is a target ID.

**Trouble:** use a new `--output` directory if the named directory exists. Check that you ran `clonodist-extended`, not the Core program, when looking for `dense.tsv`.

## 12. Paired alpha/beta distance

**Purpose:** compare both alpha and beta chains together. The included invented input has `v_a_gene`, `j_a_gene`, `cdr3_a_aa` and the corresponding beta fields. For both demo records, the alpha CDR3 is the same as the beta CDR3 shown in Section 8. This example is synthetic and is not scientifically qualified.

```bash
target/debug/clonodist distance \
  --synthetic-demo \
  --db examples/synthetic-reference.txt \
  --input examples/synthetic-receptors.tsv \
  --chains alpha,beta \
  --output output/paired
```

`--chains alpha,beta` selects both chains. The paired distance is their sum. The remaining options have the roles described in Section 10.

Open `output/paired/matrix.tsv`. It should match [the expected paired matrix](../examples/expected/paired-chain-matrix.tsv):

```tsv
query_id	demo-1	demo-2
demo-1	0	24
demo-2	24	0
```

`query_id` is the row record. Each other column is a target record. The cross-record distance is 24. **Trouble:** if a field is missing, check that both alpha and beta V, J and CDR3 columns are present.

## 13. Human gamma/delta distance

**Purpose:** exercise the Extended human gamma/delta path with invented receptors. The included input has `v_g_gene`, `j_g_gene`, `cdr3_g_aa` and the corresponding delta fields. The demo gamma and delta CDR3 values match the two invented values shown in Section 8. This example is synthetic and is not scientifically qualified.

```bash
target/debug/clonodist-extended distance \
  --synthetic-demo \
  --db examples/synthetic-reference.txt \
  --input examples/synthetic-receptors.tsv \
  --chains gamma,delta \
  --output output/gamma-delta
```

`--chains gamma,delta` selects both chains. The default organism is human. `--db`, `--input` and `--output` identify the demonstration reference, input and new output directory.

Open `output/gamma-delta/dense.tsv`. It should match [the expected paired matrix](../examples/expected/paired-chain-matrix.tsv). Its `query_id` rows and target columns are the two invented IDs; the cross-record distance is 24.

**Trouble:** check both gamma and delta fields. Do not use the invented reference for real data. Human delta qualification is limited to the documented dedicated-locus reference context. Mouse delta needs the separate setup in Section 15.

## 14. Radius neighborhood

**Purpose:** retain beta-chain pairs at distance 12 or less. This uses the two invented beta records shown in Section 8. It is a synthetic behavior check, not a biological assessment.

```bash
target/debug/clonodist neighbors \
  --synthetic-demo \
  --db examples/synthetic-reference.txt \
  --input examples/synthetic-receptors.tsv \
  --chains beta \
  --radius 12 \
  --format edges-tsv \
  --output output/neighbors
```

`neighbors` finds retained pairs. `--radius 12` includes distances equal to 12. `--format edges-tsv` asks for an edge table. The other options identify the invented reference, input, chain and new output directory. Self edges are included by default.

Open `output/neighbors/edges.tsv`. It should match [the expected radius edges](../examples/expected/radius-12-edges.tsv):

```tsv
query_id	target_id	distance
demo-1	demo-1	0
demo-1	demo-2	12
demo-2	demo-2	0
demo-2	demo-1	12
```

`query_id` is the starting record. `target_id` is a retained neighbor. `distance` is the pair's distance. Each record has two neighbors here, including itself. `output/neighbors/neighborhoods.tsv` summarizes neighbor counts and supplied counts.

**Trouble:** if you see fewer edges, check the radius and whether self edges were excluded. A missing sparse edge means no stored edge; it is different from a stored zero-distance edge.

## 15. Mouse TRD reference preparation

What is a checksum? A checksum is a long identifier calculated from a file. If a file changes, its checksum normally changes. ClonoDist uses checksums to confirm that the expected profile or reference was supplied. A matching checksum confirms file identity. It does not establish biological validity.

**Purpose:** prepare the exact supplemental reference for qualified mouse delta or paired mouse gamma/delta distance work. This is an optional governed-reference workflow, not a synthetic example. It uses source files that you obtain separately. The repository does not distribute those files.

1. Read [reference acquisition](../README.md#third-party-data-boundary), including source terms and expected hashes.
2. Obtain the exact historical base reference and IMGT TRDV and TRDJ FASTA snapshots from their original providers, subject to their terms. A current download may not match the archived snapshot.
3. Create a folder named `local-data`, then place the files in that folder. Name them `TRDV.fasta`, `TRDJ.fasta` and `combo_xcr_2024-03-05.tsv` for the following command.
4. Run:

   ```bash
   python3 scripts/prepare_mouse_trd.py \
     --trdv local-data/TRDV.fasta \
     --trdj local-data/TRDJ.fasta \
     --base local-data/combo_xcr_2024-03-05.tsv \
     --output local-data/mouse-trd.json
   ```

`--trdv` and `--trdj` identify the two user-acquired FASTA inputs. `--base` identifies the exact user-acquired base reference. `--output` names a new local supplement. The script checks file identities before it writes the result. On success it reports:

```text
Prepared pinned mouse TRD supplement; SHA-256: 482280a9dfacda0d45acf5243b68e088e7649819a05bc24b90a21f776536f649
```

The result is a JSON reference contract, not a distance table; it has no output columns. Supply it through Extended's `--mouse-trd-reference` option for a supported mouse delta distance run. Mouse delta meta inference is outside the qualified scope.

**Trouble:** a missing file or SHA-256 mismatch means the required exact input was not supplied. Check the filename and source identity. Do not bypass the check. If the exact historical snapshots cannot be obtained, these instructions cannot recreate the qualified supplement. Do not redistribute the generated reference without checking its underlying rights.

## 16. Choose a profile

**Purpose:** inspect the pinned one-chain profile. This command uses no receptor input or external reference. It creates a profile description; it is not a scientific test.

```bash
target/debug/clonodist-extended profile \
  --kind pinned \
  --chains 1 \
  --output output/pinned-profile
```

`--kind pinned` selects the pinned settings. `--chains 1` requests a one-chain profile. `--output` names a new directory. Open `output/pinned-profile/identity.json`. Its `profile_id` is `tcrdist3-pinned-0e6efa3`; its `sha256` is the canonical profile identity. The directory also contains `profile.json`. These JSON files have named fields, not table columns.

The pinned profile is the starting point for qualified work in its documented context. When using an external profile, Extended requires its canonical identity through `--profile-sha256`; a file-byte hash is not necessarily the canonical profile identity. Alternative profile mechanics have been tested, but their biological interpretation is experimental.

**Trouble:** choose a new output path if the directory exists. For a checksum mismatch, use the canonical identity for the exact profile and chain count.

## 17. Understand the output files

Each command creates a new directory. Files vary by command and format:

| File | Meaning |
|---|---|
| Core `matrix.tsv` or Extended `dense.tsv` | Dense distances. Rows are queries; columns after `query_id` are targets. |
| `edges.tsv` | Retained pairs. Columns are `query_id`, `target_id` and `distance`. |
| Core `neighborhoods.tsv` | Per-query summary. `neighbors` counts retained records; `weighted_count` sums their supplied counts. `subjects` and `samples` count distinct reported values. |
| Core `records.json` | Normalized input records and encoded chain information. |
| Core `ledger.json` | Input and output row/count totals and row-to-ID mapping. Use it to check conservation. |
| `manifest.json` | Run context and file identities. It includes settings and hashes; timing and other run metadata can vary. Synthetic runs identify their reference context. |
| Extended `query.json`, `target.json` | Normalized records used as queries and targets. |
| Extended `profile.json` | Effective profile used for the run. |
| Extended `summaries.jsonl` | Line-oriented summaries for applicable sparse output; it can be empty for a dense run. |
| Profile command `identity.json` | Profile ID and canonical SHA-256 identity. |
| CSR files | Sparse data, target indices and row pointers. A stored zero is distinct from an absent edge. |
| Operational `model.json`, `matches.jsonl` | Model or line-oriented matching results from the applicable discovery, qualification or application command. |

Keep the manifest with your input, reference and profile records. A hash identifies bytes; it does not by itself establish biological validity or authenticated provenance.

## 18. Interpret distances and neighborhoods

A smaller distance means the selected scoring method found the records closer under that reference and profile. It does not prove the receptors recognize the same antigen. A zero distance is a score, not proof of biological identity.

A radius is an inclusive cutoff. Changing it changes which edges are retained. Self edges can affect counts. Record whether you include them. `weighted_count` uses the counts you supplied; it is not automatically a cell or molecule count.

Extended operational analysis separates discovery, qualification and application. Fix the rules and background before evaluation. Keep the model and decision reasons. Study separation and controls matter. Operational matches do not provide sensitivity, specificity, HLA inference or clinical validation.

## 19. Common errors and repairs

- **`output ... exists` or refusal to overwrite:** choose a new `--output` directory. Keep the old result if you need it.
- **Reference size or SHA-256 mismatch:** acquire the exact reference identified in [reference acquisition](../README.md#third-party-data-boundary). Do not disable the check.
- **Synthetic marker rejected:** for the examples, use both `--synthetic-demo` and `examples/synthetic-reference.txt`. Do not use them for real data.
- **Unknown allele or chain mismatch:** check the exact V and J allele names, chain and organism against the selected reference. Do not silently substitute alleles.
- **Invalid input:** check for tabs, the required headers, unique IDs, canonical uppercase amino acids and positive counts.
- **Profile checksum mismatch:** use the canonical `sha256` from the identity for the exact profile. Check the chain count.
- **Memory planning refusal:** reduce the input, worker count or tile size, or choose an appropriate budget. The budget is not a hard process-memory limit.
- **Command not found:** build from the folder containing `Cargo.toml`. Check the target directory and executable name.

## 20. Reproducibility and reference rights

Record the software version, command, input identity, reference identity, profile identity, count units and output manifest. Keep records of how inputs, study groups and backgrounds were selected. Compare scientific outputs separately from timing or build metadata. See [reproducibility](../README.md#qualified-scope-and-limitations).

Third-party source records and reference datasets are not redistributed. Users must obtain governed inputs from their original providers. This includes IMGT, VDJdb and 10x Genomics materials where needed. Read [reference acquisition](../README.md#third-party-data-boundary) for source identities and terms. Do not assume that permission to use a source also permits redistribution of derived material.

## Checklist before real analysis

- Confirm the selected chain.
- Confirm the organism.
- Confirm exact V and J allele names.
- Confirm the reference source and version.
- Confirm the analysis profile.
- Confirm the meaning of `count`.
- Confirm whether self edges are included.
- Choose a new output directory.
- Keep the manifest with the results.
- Do not interpret distance as proof of antigen specificity.

## 21. Scientific limitations

Scientific qualification is limited to the documented profiles, references, inputs and environment. It does not imply qualification on every platform, species, receptor or setting. Alternative-profile biological interpretation is experimental. General gamma/delta antigen inference and mouse delta meta inference are outside scope. Local performance measurements do not establish universal speed or real-mouse scaling. Managed memory planning is not a hard RSS guarantee. Dense output can be quadratic.

## 22. Glossary

- **Allele:** a named version of a gene in the selected reference.
- **CDR3:** a receptor sequence region supplied in the input.
- **Chain:** one part of a receptor, such as alpha or beta.
- **CSR:** a sparse representation using data, indices and row pointers.
- **Distance:** a score calculated under a specified method, profile and reference.
- **Manifest:** a record of settings and identities for a run.
- **Profile:** the settings defining the scoring context.
- **Qualification:** evidence limited to stated conditions, not proof of biological or clinical validity.
- **Radius:** the inclusive maximum distance retained for a neighborhood.
- **Reference:** the external database used to resolve allele names and regions.

## 23. Get help

Read [the README](../README.md), [troubleshooting](#19-common-errors-and-repairs), and the manuscript links there. For command options, use the help commands in Section 7. Contact Yoon Sing Yap at [yoonatworkz@gmail.com](mailto:yoonatworkz@gmail.com). Do not send private or sensitive data without an appropriate agreement.
