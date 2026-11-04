# Third-party notices and method attribution

ClonoDist uses the established tcrdist3/pwseqdist scoring method. Please cite the underlying research: Dash et al., *Nature* (2017), doi:10.1038/nature22383; Mayer-Blackwell et al., *eLife* (2021), doi:10.7554/eLife.68605; and the software projects [tcrdist3](https://github.com/kmayerb/tcrdist3) and [pwseqdist](https://github.com/agartland/pwseqdist). ClonoDist's independently authored Rust implementation is not a distribution of those Python packages. Their software licenses are MIT; their data, dependencies and other components retain separate rights.

The recorded comparison environment used [Parasail](https://github.com/jeffdaily/parasail) (Daily, *BMC Bioinformatics* 2016, doi:10.1186/s12859-016-0930-z). It is not linked into the production Rust executable by this source package.

The scientific manuscript acknowledges IMGT®, the international ImMunoGeneTics information system® (CNRS/Université de Montpellier; founder and former director Marie-Paule Lefranc), VDJdb developers and contributing studies, and 10x Genomics. IMGT data and metadata, VDJdb records, 10x data, and historical composite references are not licensed by ClonoDist. See [reference acquisition](README.md#third-party-data-boundary) for exact source identities and terms. This source package does not distribute downloaded IMGT, VDJdb or 10x data.

## Rust dependencies

The dependencies below are fetched separately by Cargo, not vendored in this source repository. Their versions are pinned in `Cargo.lock`. The license expressions and repository links below are the packages' declared Cargo metadata; consult each package's own license and notices when distributing compiled artifacts. Alternative choices in an expression apply according to the package's terms. The list includes dependencies for supported targets and build configurations, not necessarily every dependency linked on one platform.

| Package | Version | Declared license | Source |
|---|---|---|---|
| anstream | 0.6.18 | MIT OR Apache-2.0 | https://github.com/rust-cli/anstyle.git |
| anstyle | 1.0.14 | MIT OR Apache-2.0 | https://github.com/rust-cli/anstyle.git |
| anstyle-parse | 0.2.6 | MIT OR Apache-2.0 | https://github.com/rust-cli/anstyle.git |
| anstyle-query | 1.1.5 | MIT OR Apache-2.0 | https://github.com/rust-cli/anstyle.git |
| anstyle-wincon | 3.0.11 | MIT OR Apache-2.0 | https://github.com/rust-cli/anstyle.git |
| block-buffer | 0.10.4 | MIT OR Apache-2.0 | https://github.com/RustCrypto/utils |
| cc | 1.5.1 | MIT OR Apache-2.0 | https://github.com/rust-lang/cc-rs |
| cfg-if | 1.0.5 | MIT OR Apache-2.0 | https://github.com/rust-lang/cfg-if |
| clap | 4.5.26 | MIT OR Apache-2.0 | https://github.com/clap-rs/clap |
| clap_builder | 4.5.26 | MIT OR Apache-2.0 | https://github.com/clap-rs/clap |
| clap_derive | 4.5.24 | MIT OR Apache-2.0 | https://github.com/clap-rs/clap |
| clap_lex | 0.7.4 | MIT OR Apache-2.0 | https://github.com/clap-rs/clap |
| colorchoice | 1.0.5 | MIT OR Apache-2.0 | https://github.com/rust-cli/anstyle.git |
| cpufeatures | 0.2.17 | MIT OR Apache-2.0 | https://github.com/RustCrypto/utils |
| crossbeam-deque | 0.8.8 | MIT OR Apache-2.0 | https://github.com/crossbeam-rs/crossbeam |
| crossbeam-epoch | 0.9.21 | MIT OR Apache-2.0 | https://github.com/crossbeam-rs/crossbeam |
| crossbeam-utils | 0.8.23 | MIT OR Apache-2.0 | https://github.com/crossbeam-rs/crossbeam |
| crypto-common | 0.1.7 | MIT OR Apache-2.0 | https://github.com/RustCrypto/traits |
| csv | 1.3.1 | Unlicense/MIT | https://github.com/BurntSushi/rust-csv |
| csv-core | 0.1.13 | Unlicense/MIT | https://github.com/BurntSushi/rust-csv |
| digest | 0.10.7 | MIT OR Apache-2.0 | https://github.com/RustCrypto/traits |
| either | 1.18.0 | MIT OR Apache-2.0 | https://github.com/rayon-rs/either |
| find-msvc-tools | 0.1.14 | MIT OR Apache-2.0 | https://github.com/rust-lang/cc-rs |
| generic-array | 0.14.7 | MIT | https://github.com/fizyk20/generic-array.git |
| heck | 0.5.0 | MIT OR Apache-2.0 | https://github.com/withoutboats/heck |
| is_terminal_polyfill | 1.70.2 | MIT OR Apache-2.0 | https://github.com/polyfill-rs/is_terminal_polyfill |
| itoa | 1.0.18 | MIT OR Apache-2.0 | https://github.com/dtolnay/itoa |
| libc | 0.2.189 | MIT OR Apache-2.0 | https://github.com/rust-lang/libc |
| memchr | 2.8.3 | Unlicense OR MIT | https://github.com/BurntSushi/memchr |
| once_cell_polyfill | 1.70.2 | MIT OR Apache-2.0 | https://github.com/polyfill-rs/once_cell_polyfill |
| proc-macro2 | 1.0.107 | MIT OR Apache-2.0 | https://github.com/dtolnay/proc-macro2 |
| quote | 1.0.47 | MIT OR Apache-2.0 | https://github.com/dtolnay/quote |
| rayon | 1.12.0 | MIT OR Apache-2.0 | https://github.com/rayon-rs/rayon |
| rayon-core | 1.13.0 | MIT OR Apache-2.0 | https://github.com/rayon-rs/rayon |
| ryu | 1.0.23 | Apache-2.0 OR BSL-1.0 | https://github.com/dtolnay/ryu |
| serde | 1.0.228 | MIT OR Apache-2.0 | https://github.com/serde-rs/serde |
| serde_core | 1.0.228 | MIT OR Apache-2.0 | https://github.com/serde-rs/serde |
| serde_derive | 1.0.228 | MIT OR Apache-2.0 | https://github.com/serde-rs/serde |
| serde_json | 1.0.145 | MIT OR Apache-2.0 | https://github.com/serde-rs/json |
| sha2 | 0.10.9 | MIT OR Apache-2.0 | https://github.com/RustCrypto/hashes |
| sha2-asm | 0.6.4 | MIT | https://github.com/RustCrypto/asm-hashes |
| shlex | 2.0.1 | MIT OR Apache-2.0 | https://github.com/comex/rust-shlex |
| strsim | 0.11.1 | MIT | https://github.com/rapidfuzz/strsim-rs |
| syn | 2.0.119 | MIT OR Apache-2.0 | https://github.com/dtolnay/syn |
| typenum | 1.20.1 | MIT OR Apache-2.0 | https://github.com/paholg/typenum |
| unicode-ident | 1.0.27 | (MIT OR Apache-2.0) AND Unicode-3.0 | https://github.com/dtolnay/unicode-ident |
| utf8parse | 0.2.2 | Apache-2.0 OR MIT | https://github.com/alacritty/vte |
| version_check | 0.9.5 | MIT/Apache-2.0 | https://github.com/SergioBenitez/version_check |
| windows-link | 0.2.1 | MIT OR Apache-2.0 | https://github.com/microsoft/windows-rs |
| windows-sys | 0.61.2 | MIT OR Apache-2.0 | https://github.com/microsoft/windows-rs |
