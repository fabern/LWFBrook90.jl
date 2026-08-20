# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Guiding principles:

- Changelogs are for humans, not machines.
- There should be an entry for every single version.
- The same types of changes should be grouped.
- Target audience is package user, not package developer.

Types of changes:

- `Added` for new features.
- `Changed` for changes in existing functionality.
- `Deprecated` for soon-to-be removed features.
- `Removed` for now removed features.
- `Fixed` for any bug fixes.
- `Security` in case of vulnerabilities.

## [Unreleased]

### Added

- This `CHANGELOG.md` file and release instructions in `README-for-Developers.md`.

## [0.9.9] - 2024-04-05

### Added

- Added gamma-distributed root profiles.
- Added `get_forcing`, `get_states`, and `get_fluxes` post-processing functions.

### Changed

- Released the version used for the [LWFBrook90.jl description paper](https://doi.org/10.1029/2024MS004445).
- Streamlined the post-processing API.

### Fixed

- Corrected the Millington-Quirk tortuosity calculation.

## [0.9.8] - 2024-01-05

### Added

- Added daily, monthly, and yearly water-partitioning summaries.
- Added a helper for preparing simulations for LWFBrook90R.
- Allowed positive beta values in normalized beta root profiles.

### Changed

- Renamed `get_aboveground` to `get_amounts`.

### Fixed

- Corrected the isotope update time step for groundwater and soil water.
- Improved setup errors and Julia 1.10 compatibility.

## [0.9.7] - 2023-07-09

### Added

- Added `get_soil_` for extracting soil states and fluxes by depth.
- Added `saveSPAC` and a time column to above-ground results.

### Changed

- Forcing data now reports an error instead of extrapolating outside its period.

### Removed

- Removed `run_simulation` in favor of the current simulation API.

## [0.9.6] - 2023-06-21

### Added

- Added cumulative fluxes and water-balance errors to amount plots.
- Added plots that highlight saturated soil areas.

### Changed

- Improved simulation performance and plotting.
- Required Julia 1.8 or later.

### Fixed

- Avoided incorrect warnings for added soil layers and nonzero start times.

## [0.9.5] - 2023-04-05

### Added

- Added more soil-horizon options to `remakeSPAC`.
- Added an option to disable solver return-code checks in `simulate!`.

### Changed

- Reduced the default output from `simulate!`.

### Fixed

- Fixed `run_example`, plot labels, and LWFBrook90R input preparation.

## [0.9.4] - 2023-03-17

### Added

- Exposed solver options through `simulate!`.
- Added control over the tolerance used to select soil output depths.

## [0.9.3] - 2023-03-14

### Added

- Added support for a single soil layer.
- Extended `remakeSPAC` soil-layer updates.

### Fixed

- Fixed LAI remaking and root-water-uptake centroid calculations.

## [0.9.2] - 2023-03-06

### Added

- Allowed initial soil and scalar states to be remade.

### Fixed

- Handled zero LAI and zero interception depth.

## [0.9.1] - 2023-03-05

### Added

- Added `remakeSPAC` for copying a model with changed parameters.
- Added helpers for generating soil, initial-condition, stomatal, and canopy inputs.
- Added LAI interpolation from key canopy parameters.

### Changed

- Changed model setup to `model = loadSPAC(...)` followed by `setup(model)`.
- Integrated soil discretization into `loadSPAC`.

## [0.9.0] - 2023-02-17

### Added

- Added a new API for defining root-density profiles.
- Added compact display summaries for `SPAC` and `DiscretizedSPAC` objects.

## [0.8.0] - 2023-02-11

### Added

- Added `plotforcingandstates` for inspecting forcing and model states.

### Changed

- Renamed state fields to include clearer units, such as `CC.MJm2` and `RWU.mmday`.

## [0.7.5] - 2023-02-01

### Added

- Added validation for gaps in forcing data.

### Changed

- Clarified the calculation of changing root profiles.

### Fixed

- Fixed occasional soil-conductivity `NaN` values and simplified snowpack updates.

## [0.7.4] - 2022-12-01

### Fixed

- Fixed layouts, legends, and labels in amount and isotope plots.

## [0.7.3] - 2022-11-30

### Added

- Added the root-water-uptake centroid to heatmaps.

### Changed

- Refined the `plotamounts` layout.

### Fixed

- Fixed solver return-code handling with SciMLBase.

## [0.7.2] - 2022-11-21

### Changed

- Made `plotamounts` more flexible.

### Fixed

- Added checks for soil-discretization input files.

## [0.7.1] - 2022-11-11

### Changed

- Documented the post-processing functions.

## [0.7.0] - 2022-11-11

### Added

- Added `plotamounts` and `plotisotopes`.
- Included refined soil discretization in simulation results.

### Changed

- Simplified the exported API and improved default plots.

## [0.6.3] - 2022-11-08

### Fixed

- Fixed `LWFBrook90.ISO.plotisotopes`.

## [0.6.2] - 2022-11-08

### Changed

- Reduced isotope plots to daily time steps.

### Fixed

- Handled root-finding bracket warnings more safely.

## [0.6.1] - 2022-11-07

### Added

- Allowed functions to define initial soil state and root density.

### Fixed

- Fixed `run_example` and improved conversion between R and Julia inputs.

## [0.6.0] - 2022-11-05

### Added

- Included date and time values in discretized models.

### Changed

- Simplified model inputs and removed the unused `CS` parameter.
- Cleaned up the exported API.

## [0.5.1] - 2022-11-03

### Added

- Added continuous interpolation for time-varying model inputs.
- Allowed a custom vertical step in soil discretization.

### Changed

- Reworked model objects and state access around `SPAC`, `DiscretizedSPAC`, and component arrays.

### Fixed

- Fixed groundwater state handling and mixed-up groundwater isotope labels.

## [0.5.0] - 2022-07-26

### Added

- Added xylem storage and isotope tracking, including hydraulic lift.
- Added xylem volume as an input parameter.

### Changed

- Standardized isotope inputs as delta values.
- Normalized time-varying vegetation inputs.

## [0.4.0] - 2022-07-18

### Added

- Added isotope transport and isotope result plots.
- Added `solve_LWFB90` as a simpler simulation API.
- Added water-balance error output and daily root-water-uptake plots.

### Changed

- Made input validation stricter and improved simulation performance.

### Fixed

- Fixed irregular input interpolation, snow accumulation, and several isotope edge cases.

## [0.3.6] - 2022-03-09

### Changed

- Improved reliability checks against the published Hammel simulations.

## [0.3.5] - 2022-02-28

### Added

- Allowed manual soil and root discretization.
- Added a plot recipe for simulation results.

### Changed

- Renamed depth inputs to state their units explicitly.

### Fixed

- Fixed initial conditions, saved time steps, and soil discretization.

## [0.3.4] - 2021-09-26

### Added

- Allowed soil refinement at requested state-output depths.

### Fixed

- Fixed the example script.

## [0.3.3] - 2021-09-23

### Fixed

- Corrected a unit typo in the input data.

## [0.3.2] - 2021-08-24

### Changed

- Renamed and documented the input-reading interface.

## [0.3.1] - 2021-08-21

### Fixed

- Fixed unit handling in generated R input files, documentation, and examples.

## [0.3.0] - 2021-08-21

### Changed

- Split soil horizons, soil discretization, and initial conditions into separate inputs.
- Added units to input data and checked them when files are read.
- Simplified and renamed several input parameters and files.

### Removed

- Removed measured streamflow input and obsolete input files.

## [0.2.0] - 2021-03-16

### Changed

- Streamlined input preprocessing and time-varying input interpolation.
- Reorganized soil parameters and simplified the internal parameter set.

### Fixed

- Fixed bypass-flow activation, xylem flow, and initial root dimensions.

## [0.1.4+doc1] - 2021-03-10

### Changed

- Documentation-only tag; no package code changes.

## [0.1.4] - 2021-03-10

### Changed

- Improved documentation deployment and release tests.

## [0.1.3+doc1] - 2021-03-09

### Changed

- Documentation-only tag; no package code changes.

## [0.1.3] - 2021-03-09

### Changed

- Updated the user documentation.

## [0.1.2+doc1] - 2021-03-08

### Changed

- Documentation-only release; no package code changes.

## [0.1.2] - 2021-03-05

### Changed

- Reduced unnecessary package imports.

### Fixed

- Fixed the example-data path used by `run_example`.

## [0.1.1] - 2021-03-01

### Added

- Added `run_example` for running the included example data set.

### Fixed

- Fixed images and links in the README and documentation.

## [0.1.0] - 2021-03-01

### Added

- Initial release of the LWFBrook90 hydrological model for Julia.
- Included the Beatenberg example data set and simulation output plots.

[Unreleased]: https://github.com/fabern/LWFBrook90.jl/compare/v0.9.9...HEAD
[0.9.9]: https://github.com/fabern/LWFBrook90.jl/compare/v0.9.8...v0.9.9
[0.9.8]: https://github.com/fabern/LWFBrook90.jl/compare/v0.9.7...v0.9.8
[0.9.7]: https://github.com/fabern/LWFBrook90.jl/compare/v0.9.6...v0.9.7
[0.9.6]: https://github.com/fabern/LWFBrook90.jl/compare/v0.9.5...v0.9.6
[0.9.5]: https://github.com/fabern/LWFBrook90.jl/compare/v0.9.4...v0.9.5
[0.9.4]: https://github.com/fabern/LWFBrook90.jl/compare/v0.9.3...v0.9.4
[0.9.3]: https://github.com/fabern/LWFBrook90.jl/compare/v0.9.2...v0.9.3
[0.9.2]: https://github.com/fabern/LWFBrook90.jl/compare/v0.9.1...v0.9.2
[0.9.1]: https://github.com/fabern/LWFBrook90.jl/compare/v0.9.0...v0.9.1
[0.9.0]: https://github.com/fabern/LWFBrook90.jl/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/fabern/LWFBrook90.jl/compare/v0.7.5...v0.8.0
[0.7.5]: https://github.com/fabern/LWFBrook90.jl/compare/v0.7.4...v0.7.5
[0.7.4]: https://github.com/fabern/LWFBrook90.jl/compare/v0.7.3...v0.7.4
[0.7.3]: https://github.com/fabern/LWFBrook90.jl/compare/v0.7.2...v0.7.3
[0.7.2]: https://github.com/fabern/LWFBrook90.jl/compare/v0.7.1...v0.7.2
[0.7.1]: https://github.com/fabern/LWFBrook90.jl/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/fabern/LWFBrook90.jl/compare/v0.6.3...v0.7.0
[0.6.3]: https://github.com/fabern/LWFBrook90.jl/compare/v0.6.2...v0.6.3
[0.6.2]: https://github.com/fabern/LWFBrook90.jl/compare/v0.6.1...v0.6.2
[0.6.1]: https://github.com/fabern/LWFBrook90.jl/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/fabern/LWFBrook90.jl/compare/v0.5.1...v0.6.0
[0.5.1]: https://github.com/fabern/LWFBrook90.jl/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/fabern/LWFBrook90.jl/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/fabern/LWFBrook90.jl/compare/v0.3.6...v0.4.0
[0.3.6]: https://github.com/fabern/LWFBrook90.jl/compare/v0.3.5...v0.3.6
[0.3.5]: https://github.com/fabern/LWFBrook90.jl/compare/v0.3.4...v0.3.5
[0.3.4]: https://github.com/fabern/LWFBrook90.jl/compare/v0.3.3...v0.3.4
[0.3.3]: https://github.com/fabern/LWFBrook90.jl/compare/v0.3.2...v0.3.3
[0.3.2]: https://github.com/fabern/LWFBrook90.jl/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/fabern/LWFBrook90.jl/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/fabern/LWFBrook90.jl/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/fabern/LWFBrook90.jl/compare/v0.1.4...v0.2.0
[0.1.4+doc1]: https://github.com/fabern/LWFBrook90.jl/tree/v0.1.4%2Bdoc1
[0.1.4]: https://github.com/fabern/LWFBrook90.jl/compare/v0.1.3...v0.1.4
[0.1.3+doc1]: https://github.com/fabern/LWFBrook90.jl/tree/v0.1.3%2Bdoc1
[0.1.3]: https://github.com/fabern/LWFBrook90.jl/compare/v0.1.2...v0.1.3
[0.1.2+doc1]: https://github.com/fabern/LWFBrook90.jl/tree/v0.1.2%2Bdoc1
[0.1.2]: https://github.com/fabern/LWFBrook90.jl/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/fabern/LWFBrook90.jl/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/fabern/LWFBrook90.jl/releases/tag/v0.1.0
