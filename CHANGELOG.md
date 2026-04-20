# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pytantan/compare/v0.1.4...HEAD

## [v0.1.4] - 2026-04-20
[v0.1.4]: https://github.com/althonos/pytantan/compare/v0.1.3...v0.1.4

### Fixed
- Inadequate runtime dispatch of AVX2 code causing crashes on pre-Haswell x86-64 ([#2](https://github.com/althonos/pytantan/issues/2)).

### Changed
- Compile wheels for Python Limited API 3.11.


## [v0.1.3] - 2024-10-21
[v0.1.3]: https://github.com/althonos/pytantan/compare/v0.1.2...v0.1.3

### Changed
- Update `scoring-matrices` dependency to `v0.3.0`.


## [v0.1.2] - 2024-10-15
[v0.1.2]: https://github.com/althonos/pytantan/compare/v0.1.1...v0.1.2

### Changed
- Rewrite package build using `scikit-build-core`.
- Update documentation to use the PyData theme.

### Removed
- Support for Python 3.6.


## [v0.1.1] - 2024-05-15
[v0.1.1]: https://github.com/althonos/pytantan/compare/v0.1.0...v0.1.1

### Fixed
- Compilation on `armv7` platform.
- Patching of MacOS compiler flags for platform-specific support.


## [v0.1.0] - 2024-05-08
[v0.1.0]: https://github.com/althonos/pytantan/compare/324bdb80e...v0.1.0

Initial release.
