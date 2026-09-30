# Multispecies design walkthrough

The revised proposal is in [the design document](../docs/multispecies-design.md).
The [Julia example](example-script-02-multispecies.jl) constructs the proposed input
data, with one complete configuration per species. It does not call the solver or
demonstrate an implemented `loadSPAC` API.

Each species owns its roots, LAI/phenology, xylem hydraulics, stomatal control,
interception and isotope initial conditions. Site forcing, soil properties,
soil water and the ground snowpack are shared by default, even with tiled canopies.
Ground snow representation is a separate choice from canopy mixing. A single `species` argument keeps related inputs together
and makes partial overrides and validation local to the selected component.

The canopy model must be explicit. Area tiles use local radiation intensity and
convert fluxes and stores to stand area once, using `cover_fraction`. Species in
an intermingled canopy need joint radiation/aerodynamic calculations and a different
area convention. With `ground_snowpack = :shared`, aggregate tile throughfall and
snow-surface forcing before updating one ground snowpack; keep intercepted snow
species-specific. The proposal documents these choices, together with the required
changes to loading, setup, solver state/caches, soil coupling, isotopes and output.

The `DAV2020-multispecies-full` and `DAV2020-multispecies-bare-minimum` CSVs retain
the original prototype format. They provide shared forcing and independent plant
parameter/root/phenology columns, but their existence does not establish solver
support or numerical validity. The design document describes the migration to a
strict ownership schema and the tests needed before claiming species independence
or single-species equivalence.

The existing uncommitted solver prototype and tests were preserved through the
rebase. Their known gaps are recorded in the design document. No multispecies
simulation results are asserted as validated by this revision.
