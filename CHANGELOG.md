# Changelog

## v1.1.0 - 2026-03-30 [v1.1.0 Release]

### Changed

- Refactored the main pipeline entrypoint in `QPPL.py` to use structured logging and cleaner task dispatch.
- Updated configuration handling in `modules/config.py` with centralized parameter schemas, typed parsing, and generated defaults.
- Updated default config (`qppl.conf`) keys and defaults:
  - Renamed filter keys to `cut_adapter`, `filter_quality`, `filter_length`.
- Refactored task modules (`modules/filter.py`, `modules/readquality.py`, `modules/assembly.py`, `modules/characterization.py`) to:
  - use reusable command generation helpers,
  - improve per-file and total runtime reporting,
  - standardize command execution via conda environments.
- Improved assembly workflow orchestration, including explicit step ordering and output directory constants.
- Improved characterization workflow orchestration and standardized intermediate/result table generation.
- Updated CLI parser in `modules/parser.py` with explicit `-h/--help` registration.
- Updated documentation in `README.md` with revised setup, task descriptions, configuration schema, and output layout.
- Bumped version metadata to `v1.1.0` in `modules/version.py` and made the logo version dynamic.

### Added

- Added colorized console logging and per-task log files (`readquality.log`, `filter.log`, `assembly.log`, `characterization.log`).
- Added new config keys in `qppl.conf`: `porechop_threads`, `shasta_read_length`, `polish_threads`, and `taxmyphage_db`.
- Expanded the filtering step to automatically fetch NanoLyse reference files when missing.

### Removed

- Removed legacy config field `prefix`.
- Removed legacy filter key names from defaults and parser logic.

### Fixed

- Improved FASTQ input detection and extension handling (including case-insensitive suffix stripping in task modules).
- Improved error visibility for failed external tool commands through structured logging.

## v1.0.0 - 2025-04-13 [ESCMID 2025 Release]

### Added

- Initial public release of QPPL.
