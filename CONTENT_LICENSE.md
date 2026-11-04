# Content license

Copyright © 2026 Yoon Sing Yap. Independently authored, rights-cleared ClonoDist documentation prose, manuscript prose, and original diagrams and figure presentations in this repository are licensed under the [Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/). Attribution must identify the creator, indicate changes, and link to the license. The full legal terms are at https://creativecommons.org/licenses/by/4.0/legalcode.en.

The software source code, project-authored executable examples, tests, and scripts are covered by [LICENSE](LICENSE), the MIT License. This content notice does not change the terms of code snippets incorporated into the software.

Neither license grants rights in third-party data, code, databases, source sequences, names, trademarks, or other material whose copyright belongs to others. In particular, IMGT, VDJdb, 10x Genomics, tcrdist3, pwseqdist, Parasail, historical composite references and Cargo dependencies retain their own terms. Their appearance as citations or method acknowledgments is not a grant of a license. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and [REFERENCE_ACQUISITION.md](README.md#third-party-data-boundary). Do not assume that reference-derived or VDJdb-derived artifacts can be redistributed under the ClonoDist licenses.

## Path boundaries

- `src/**`, executable code in `scripts/*.py`, and `Cargo.toml` files: original software, MIT, subject to the notices above.
- `*.md`, `docs/USER_MANUAL.md`, and `docs/CLONODIST-MANUSCRIPT-v1.0.0.pdf`: original prose and original figure presentations, CC BY 4.0 where owned and rights-cleared. Code excerpts remain MIT.
- `Cargo.lock` and `scripts/mouse_trd_provenance.json`: project-authored dependency or provenance metadata. Identifiers and source facts do not change the rights in the referenced third-party materials. No reference data or third-party software license is granted by their inclusion.
- `LICENSE`, `CITATION.cff`, and `.gitignore`: license or project metadata; the license text applies to original software as stated in `LICENSE`.

These boundaries do not authorize redistribution of third-party material. The license for original presentation does not extend to the underlying third-party sources or grant redistribution rights in their records or reference data.
