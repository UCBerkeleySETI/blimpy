This file is a version history of turbo_seti amendments, beginning with version 2.0.2.  Entries appear in version descending order (newest first, oldest last).
<br>
<br>
|    Date    | Version | Contents |
| :--: | :--: | :-- |
| 2021-08-07 | 2.0.21 | New signal_processing source file, "dedoppler.py" (discussion in PR #220).  |
| 2021-08-06 | 2.0.20 | New utility, "stix" (Issue #221).  |
| 2021-07-28 | 2.0.19 | Update fil2h5 to handle YUGE data matrixes.  |
| 2021-07-17 | 2.0.18 | Get rid of numpy "RuntimeWarning: Mean of empty slice" messages (Issue #212).  |
| 2021-07-13 | 2.0.17 | New utility: stax.  |
| 2021-07-08 | 2.0.16 | Increase test coverage of calc_n_coarse_chan().  |
| | | Improve messaging when calc_n_coarse_chan() emits warnings. |
| 2021-06-14 | 2.0.15 | Fix issue #210 - Guard against unusual Filterbank headers created by setigen apps.  |
| 2021-06-12 | 2.0.14 | Fix issue #208 - Miscalculated max_data_array_size when available RAM < 1 GB.  |
| 2021-06-10 | 2.0.13 | Fix issue #205 - Define MeerKAT in the list of observatories. |
| | | Fix issue #207 guppi.py generate_filterban_header(). |
| 2021-05-29 | 2.0.12 | Fix issue #203 - calc_n_coarse_chan default to 64. |
| 2021-04-14 | 2.0.11 | Fix issue #196 - automate memory requirements. |
| 2021-03-11 | 2.0.10 | Reopened enhancement #178 - calcload utility - added a verbose parameter. |
| 2021-03-08 | 2.0.9 | Implemented enhancement #182 - rawhdr utility (get header from raw files). |
| | | Amended setup.cfg to enable hdf5plugin to be installed optimized by installation. |
| 2021-03-05 | 2.0.8 | Yanked NUMBA from requirements.txt and waterfall.py due to observed instability in a large data array. |
| 2021-03-04 | 2.0.7 | Fix issue #177 - Amend waterfall.py by adding a \_\_del\_\_ function to ensure that HDF5 files are closed. |
| | | Fix issue #178 - Introduce a new utility (calcload) to calculate max_load for Waterfall. |
| 2021-03-01 | 2.0.6 | Fix issue #171 - grab_data() needed a clear error message when "heavy" data had not been loaded. |
| 2020-12-18 | 2.0.5 | Ignore documentation files in CI (PR #166). |
| | | Slice and dice by time as well as frequencies in the `dice` command (PR #167). |
| 2020-12-18 | 2.0.4 | Deprecate Travis CI in favor of Github Actions (PR #164). |
| 2020-12-16 | 2.0.3 | Numba acceleration in Waterfall.block_dc() (PR #162). |
| | | Removed references to nonexistent `filutil` command. |
| | | Removed generated pop-up plot windows while running `pytest`; the figures are still saved. |
| | | Removed outdated Docker files. |
| | | Updated setuptools build requirements. |
| | | Updated documentation. |
| 2020-12-15 | 2.0.2  | Current as of 2020-12-15. |
