# ClonoDist

ClonoDist calculates reference-qualified T-cell receptor distances and supports inspectable neighborhood analysis. It uses the established tcrdist3/pwseqdist scoring method; it does not introduce a new receptor-distance formula.

## Main capabilities

Core supports human and mouse alpha, beta and paired alpha/beta distances, joins and neighborhoods. Extended adds defined gamma/delta contexts, explicit profiles, and operational discovery, qualification and application. Outputs record the execution context and decisions.

## Qualified scope and limitations

Qualification is limited to the documented inputs, references, profiles, workflows and local Apple Silicon environment. Human delta uses a dedicated-locus context. Mouse delta requires a supplemental reference; mouse delta meta inference is outside scope. Alternative-profile biological interpretation is experimental. Synthetic examples test software behavior only. Distance similarity does not prove shared antigen specificity. Operational matches do not establish sensitivity, specificity, HLA inference, causal biology or clinical validity. Dense work and output can be quadratic; managed memory planning is not a hard RSS limit.

For reproducibility, retain input and reference identities, the profile, program version, command and output manifest. Record count units and how records and backgrounds were selected. Timing and build metadata can differ even when scientific output is equal.

## Installation

Install Rust, Cargo and a C build toolchain. Open a terminal in this repository. Build the workspace:

```bash
cargo build --workspace --locked
```

The first build needs network access for uncached dependencies.

## Synthetic example

This example uses invented input and reference records. It does not establish biological validity. Run it after the build:

```bash
mkdir -p output
target/debug/clonodist distance \
  --synthetic-demo \
  --db examples/synthetic-reference.txt \
  --input examples/synthetic-receptors.tsv \
  --chains beta \
  --output output/readme
```

Open `output/readme/matrix.tsv`. The distance between the two invented records is 12. The output directory must be new.

## User manual

Read the [User manual](docs/USER_MANUAL.md) for installation, examples, output interpretation, profiles and troubleshooting.

## Manuscript

Read the [Manuscript](docs/CLONODIST-MANUSCRIPT-v1.0.0.pdf) for the accepted methods, results and limitations.

## Third-party data boundary

Third-party source records and reference datasets are not redistributed. The manuscript reports project-generated aggregate analyses. Users must obtain governed inputs from their original providers and follow their terms. Permission to use a source does not automatically permit redistribution of derived material.

For Core and Extended reference-based work, obtain `tcrdist/db/combo_xcr_2024-03-05.tsv` from [tcrdist3 commit 0e6efa39092e2c64c9c0b52a72c2f5961c161ebb](https://github.com/kmayerb/tcrdist3/tree/0e6efa39092e2c64c9c0b52a72c2f5961c161ebb). Its expected SHA-256 is `1d7fcac37ad6d955b3a3f693526e36623186464aef2c528cf8e0322dc89706cd`. The software license does not establish all rights in that historical composite. Supply the file locally with `--db`; do not bypass identity checks or redistribute it without resolving its rights.

Mouse delta also requires the exact IMGT mouse TRDV and TRDJ FASTA snapshots retrieved on 2026-10-29. Obtain them from the [IMGT reference directory](https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/TR/) under the [IMGT terms](https://www.imgt.org/about/termsofuse.php). The required TRDV SHA-256 is `41956333b48e440e2adeaf7b9546a2fe1f60ba7853b2027fb4e1b5067729c49c`; TRDJ is `6276c6ea88f4f248fa5c37886a6194a5d6b1b47f1f0927828b4b497e8401ca34`. The manual describes the preparation command. Current files may differ; if the exact historical snapshots are unavailable, this qualified mouse-delta reference cannot be recreated. Do not distribute the generated supplement without checking underlying rights.

The scientific comparison used [pwseqdist commit 852e159ad12582973bfbf23dd33fef068723742e](https://github.com/agartland/pwseqdist/tree/852e159ad12582973bfbf23dd33fef068723742e) and the tcrdist3 source above. The optional observational study used [VDJdb release 2025-12-29](https://github.com/antigenomics/vdjdb-db/releases/tag/2025-12-29), whose recorded archive license is CC BY-ND 4.0, and the [10x Genomics healthy-donor dataset](https://www.10xgenomics.com/datasets/human-t-cells-from-a-healthy-donor-1-k-cells-multi-v-2-2-standard-5-0-0), identified by its provider as CC BY 4.0. Obtain these separately under their terms. The recorded VDJdb archive `vdjdb-2025-12-29.zip` has SHA-256 `36fac39bcc57817fd807e3de9a3ad90820ef8aa6db47578e2016be27dd5da25a`. The recorded 10x filtered contig file `sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_filtered_contig_annotations.csv` has SHA-256 `2356d5babeb29a0e5133e7b139ffa20016cd8a84fcb9d7c0b71431163be83d52`. The original selection materials are not supplied, so downloading them alone does not reproduce the published study. See [third-party notices](THIRD_PARTY_NOTICES.md).

## Citation

See [CITATION.cff](CITATION.cff) and cite the underlying methods described in the manuscript. Repository: https://github.com/Yoonecorn/ClonoDIst.

## License

Original software is licensed under [MIT](LICENSE). Original rights-cleared prose and figure presentation use [CC BY 4.0](CONTENT_LICENSE.md). Third-party materials retain their own terms.

## Contact

Yoon Sing Yap — yoonatworkz@gmail.com
