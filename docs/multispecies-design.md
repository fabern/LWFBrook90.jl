# Multispecies input and solver design

Status: design proposal, revised against `develop` at `4a7f361`. The existing
multispecies code is an unvalidated prototype. The interface below is proposed;
these keywords are not yet a supported `loadSPAC` interface.

## Recommendation

Supply one complete vegetation configuration for each species or cohort, alongside
one shared site and soil configuration. Each vegetation entry owns its parameters,
root distribution, canopy evolution, initial conditions, derived quantities, and
runtime caches. Species share soil water through a single soil solver. By default,
they also share the ground snowpack, while intercepted snow remains species-specific.

Use a single `species` argument instead of parallel species dictionaries in
`params`, `root_distribution`, `canopy_evolution`, and `cover_fractions`. For example,
`species.Beech` contains everything needed to construct the Beech plant component.
This also allows two cohorts of the same biological species to have different ages
or rooting profiles. The keys identify model components, not a fixed taxonomy.

Keep the current BROOK90 parameter names inside each entry's `params` field. A
loader can accept NamedTuples and CSVs, then normalize both into the same internal
configuration. Parameter ownership must be declared in a schema, rather than
inferred from whichever CSV cell is populated or whichever species is first.

The construction-only [Julia example](../examples/example-script-02-multispecies.jl)
shows the proposed shape. Its `site_params` contains shared parameter overrides;
its `species` entries each contain `params`, `root_distribution`,
`canopy_evolution`, and `initial_conditions`. `IC_soil`, the soil grid, weather,
precipitation isotopes, irrigation, storm durations, and solver options remain
shared inputs. `IC_scalar` supplies shared groundwater and ground snowpack initial
conditions. `ground_snowpack = :shared` is independent of `canopy_mixing = :tiles`.

## Ownership of parameters and states

The following is the proposed ownership for an area-tile canopy model. A shared
mixed-canopy model needs an additional stand canopy calculation, described below.
Equal defaults may be resolved independently for every species; equality of values
does not make the corresponding plant state shared.

| Owner | Inputs and derived quantities |
| --- | --- |
| Species roots | `root_distribution` or `Rootden_<species>`; `MXRTLN`, `RTRAD`, `INITRLEN`, `INITRDEP`, `RGRORATE`, `RGROPER`, `NOOUTF`; its own age-dependent `RELDEN(t, layer)` and root/rhizosphere resistances |
| Species canopy and phenology | `MAXLAI`, `DENSEF_baseline_`, `SAI_baseline_`, `HEIGHT_baseline_m`, `AGE_baseline_yrs`; its own LAI, SAI, height, density and age evolution |
| Species xylem hydraulics | `MXKPL`, `FXYLEM`; its own plant and xylem resistance and leaf water limitation |
| Species stomatal control | `GLMAX`, `GLMIN`, `R5`, `CVPD`, `RM`, `TL`, `T1`, `T2`, `TH`, `PSICR`; its own canopy conductance, potential transpiration and actual uptake |
| Species/tile canopy exchange | `CR`, `LWIDTH`, `RHOTP`, `NN`, `LPC`, `CZS`, `CZR`, `HS`, `HR`, `ZMINH`, `ALB`, `ALBSN`, `KSNVP`, `LAIMLT`, `SAIMLT`; its own exchange calculations in tile mode |
| Species interception | `FRINTLAI`, `FSINTLAI`, `FRINTSAI`, `FSINTSAI`, `CINTRL`, `CINTRS`, `CINTSL`, `CINTSS`; separate `INTR` and `INTS` amounts and isotope signatures |
| Species isotope reservoir | `VXYLEM_mm`, initial xylem isotope signatures; separate xylem composition and root-water-uptake diagnostics |
| Shared site and forcing | `LAT_DEG`, `ESLOPE_DEG`, `ASPECT_DEG`, `C1`, `C2`, `C3`, `WNDRAT`, `FETCH`, `Z0W`, `ZW`, `RSTEMP`, `Z0G`, `Z0S`; weather, precipitation, irrigation and their isotope inputs |
| Shared soil and groundwater | Soil horizons, hydraulic properties, grid, soil water/potential/isotopes, groundwater state; `DISPERSIVITY_mm`, `IDEPTH_m`, `QDEPTH_m`, `SLVPDEPTH_m`, `RSSA`, `RSSB`, `INFEXP`, `BYPAR`, `QFPAR`, `QFFC`, `IMPERV`, `DSLOPE`, `LENGTH_SLOPE`, `DRAIN`, `GSC`, `GSP` |
| Ground snowpack | Common coefficients `MELFAC`, `CCFAC`, `GRDMLT`, `MAXLQF`, `SNODEN`; shared `SNOW`, `CC`, `SNOWLQ` and snow isotope states by default, for either canopy mode. Optional separate ground snowpacks are a distinct spatial-resolution choice. |
| Solver | `DTIMAX`, `DSWMAX`, `DPSIMAX`, isotope/irrigation flags; one set per simulation |

`VXYLEM_mm` is currently a prescribed isotope mixing volume. It does not introduce
plant capacitance or a vulnerability curve. The existing RHS keeps xylem water
amount constant. A future capacitance model needs separate species water stores,
water potentials, capacitance/vulnerability parameters and a water balance that
distinguishes root uptake from transpiration. Do not invent those parameters in
the input schema before that hydraulic model is selected.

## Make the canopy assumption explicit

Species identity alone does not define how vegetation shares radiation. Require an
explicit `canopy_mixing` choice for multispecies input; do not silently choose it
from the number of columns.

| Choice | Meaning | Area convention |
| --- | --- | --- |
| `:tiles` | Non-overlapping canopy tiles above a common, horizontally mixed soil column. Each tile solves its own canopy exchange and interception. This approximates a mosaic with shared belowground water. | `cover_fraction = w_s` is the tile area divided by total stand area; all fractions are explicit and sum to one. LAI, root length, conductance and storage inputs are per tile ground area. |
| `:mixed` | Species coexist in a common canopy. Light allocation, aerodynamic coupling and ground energy balance are calculated jointly from canopy structure and species contributions. | Species LAI, root length, conductance and stores are expressed per total stand ground area. A tile `cover_fraction` must not be applied again. |

The current CSV examples describe the `:tiles` input convention. They do not
establish a radiation model for an intermingled forest. Supporting `:mixed` requires
a choice of canopy structure/light allocation and a coupled exchange calculation;
the implementation must reject that mode until those equations are implemented.

For tiles, each tile receives the same incident radiation and precipitation
**intensity** as the site. Evaluate the nonlinear canopy, stomatal and hydraulic
calculations with that intensity. Convert the resulting extensive flux or storage
to stand area exactly once. Shared ground stores already use stand area and are
not multiplied by individual tile fractions:

```text
TRANI_stand[i, s] = w_s * TRANI_tile[i, s]
TRANI_soil[i]     = sum(TRANI_stand[i, s] for s)
LAI_stand        = sum(w_s * LAI_tile[s] for s)
INTR_stand       = sum(w_s * INTR_tile[s] for s)
VXYLEM_stand[s]  = w_s * VXYLEM_tile[s]
```

Do not pass `w_s * radiation` to each tile's stomatal response. Do not multiply
already weighted LAI or plant conductivity by the fraction a second time. Xylem
turnover must use uptake and reservoir volume on the same area basis: either
`TRANI_tile / VXYLEM_tile` or `TRANI_stand / VXYLEM_stand`.

## Ground snowpack: independent of canopy tiling

Canopy tiling does not require separate ground snowpacks. Recommend
`ground_snowpack = :shared` for the first implementation in both canopy modes,
consistent with the shared soil column. Keep each species' intercepted rain `INTR`
and intercepted snow `INTS` separate. Store ground `SNOW` (SWE), cold content `CC`,
liquid water `SNOWLQ` and ground snow isotope composition once per stand area.
An optional `ground_snowpack = :tiles` would retain spatial snow histories under
each canopy tile; it is a separate, more detailed model choice, not implied by
`canopy_mixing = :tiles`.

This separation has a precedent in [CLM's subgrid hierarchy](https://escomp.github.io/CTSM/tech_note/Ecosystem/CLM50_Tech_Note_Ecosystem.html#surface-heterogeneity-and-data-structure):
vegetation patches can share one soil/snow column, with area-weighted patch fluxes
supplying its boundary conditions. That precedent supports the model structure;
it does not validate a particular BROOK90 aggregation or its site-level accuracy.

For tiled canopies above a shared ground snowpack, use this sequence each
precipitation interval:

1. Read the common ground snow state and temperature. All canopy tiles use that
   state for snow burial, roughness, snow-albedo switching and surface exchange.
2. Compute interception and wet-canopy corrections separately for each species.
   Compute each tile's rain/snow throughfall, potential snow evaporation/condensation
   (`PSNVP`) and snow energy flux (`SNOEN`) using its own canopy traits and the
   common ground state. These are local flux densities, before snow-availability
   and melt/refreezing constraints are applied.
3. Area-weight those boundary fluxes into stand fluxes. For example:

   ```text
   RTHR  = sum(w_s * RTHR_tile[s] for s)
   STHR  = sum(w_s * STHR_tile[s] for s)
   SNOEN = sum(w_s * SNOEN_tile[s] for s)
   PSNVP = sum(w_s * PSNVP_tile[s] for s)
   ```

4. Call `SNOWPACK` once with the shared `SNOW`, `CC`, `SNOWLQ` and aggregated
   forcing. Apply groundmelt, snow-availability limits, liquid retention and
   refreezing once. Its accepted `SNVP`, `SMLT` and `RSNO` enter the stand budgets;
   soil input is `RTHR - RSNO + SMLT`, plus the existing shared irrigation input.
5. Update the ground snow isotope balance once, using the accepted shared fluxes.
   Mix incoming isotope atom-fraction fluxes before the update. Apply the common
   snow-cover condition to soil evaporation for all canopy tiles.

Evaluate nonlinear canopy functions before averaging their fluxes. In particular,
BROOK90's `SNOENRGY` contains exponential LAI/SAI attenuation: energy calculated
from an average LAI/SAI is generally different from the weighted sum of tile energy
fluxes. Ground snow temperature is derived from the common cold content and SWE,
not from an arithmetic mean of hypothetical tile temperatures.

The model assumption is one effective ground snow state for the stand. Differences
in snowfall reaching the ground, shading and exchange still affect the boundary
fluxes, but their spatial covariance with snow depth, cold content and liquid water
is discarded. This represents a lumped column; it does not explicitly simulate
lateral snow or heat transport. Distinct persistent snow patches require separate
ground states (or a future snow-cover-fraction model).

Consequences include common snow depth and temperature, common snow presence for
albedo/soil-evaporation switching, and one meltwater isotope history. Differences
in snow disappearance below different canopies cannot be represented. Because
melt, refreezing, liquid retention and availability limits are nonlinear, stand
melt and recharge can differ from a model with separate snowpacks, even if both
models conserve water. There is no general fixed direction of that difference.

A simplified example illustrates the threshold effect. Start with 10 mm SWE on
each of two equal tiles, ignore cold content and liquid retention, and give them
potential melt of 20 and 0 mm during a step. Separate snowpacks produce 5 mm stand
melt after each local melt is capped at its snow inventory. A shared snowpack sees
10 mm potential stand melt and can melt all 10 mm. This is the change in model
assumption that must be accepted when choosing a shared ground state.

Averaging **updated** tile snowpacks and feeding their average back to all tiles
represents neither preserved local histories nor the shared model described here.
For a shared snowpack, aggregate boundary fluxes **before one snowpack update**.
For separate snowpacks, retain each state and aggregate outputs only for reporting
and input to the shared soil column. Ground snow representation and intercepted
snow representation remain independent decisions.

## Loading, overrides and backward compatibility

1. Declare the ordered species/cohort IDs once, from `species` or the parameter
   CSV. Match root, phenology and initial-condition entries against those IDs.
   Reordering species must only reorder the corresponding outputs.
2. Split the legacy default parameter list into shared, species and solver
   categories using the ownership schema. Apply shared defaults once and species
   defaults independently. Resolve each required `NaN` input for every species.
3. Use precedence: category defaults, file values, explicit Julia overrides.
   Merge partial overrides into the corresponding species entry; reject unknown
   fields and IDs. A shared parameter in a species entry, or a species parameter
   in `site_params`, is an error. No implicit inheritance from another species.
4. Keep wide CSVs as a file adapter: `param_id,site,Beech,Spruce`, species-suffixed
   phenology columns, `Rootden_Beech`, `Rootden_Spruce`, and species IC rows. Normalize
   them into the same species records used by programmatic input. In the proposed
   strict format, a species row has `NA` in `site`; a shared row has `NA` in all
   species columns. Identical species values are repeated, not placed in `site`.
5. Validate finite, nonnegative tile fractions and their unit sum; do not supply
   equal fractions automatically. Validate each root profile, rooting depth,
   phenology time coverage, conductance/temperature parameters, storage capacity,
   isotope initial values, and grid coverage independently. Inactive, zero-area
   tiles have zero stand contribution and must bypass divisions by their area.
6. Keep legacy single-species input working through an adapter that creates one
   species at unit area. Existing `params`, `root_distribution`, `canopy_evolution`
   and plant IC keywords retain their old meaning when `species` is absent. With
   `species` present, plant overrides belong under their species entry; reject
   ambiguous combinations. `remakeSPAC` must follow the same rules.
7. Keep the normalized configuration intact in `SPAC.pars`. Never discard the
   shared site values while extracting species records. Adapt `saveSPAC` and input
   export to round-trip the IDs, category values and area convention. A single
   species compatibility view must never be used as multispecies solver input.

When updating the CSV fixtures to the strict schema, move `RM`, canopy aerodynamic
coefficients and tile exchange/snow modifiers from `site` into every species column.
Keep solver settings and soil settings, including `SLVPDEPTH_m`, shared. With a
shared ground snowpack, retain the existing common snow/heat/liquid-water IC rows;
only plant and interception ICs belong under species. An optional separate-ground-
snowpack adapter would additionally need tile snow/heat/liquid-water IC rows. The
existing fixtures remain prototype-format examples pending the loader change.

## Required solver changes

| Location | Required change |
| --- | --- |
| `func_read_inputData.jl` | Normalize and validate the ownership schema and species records; retain shared parameters; implement override and legacy adapters. |
| `func_discretize_soil_domain.jl` and `setup` in `LWFBrook90.jl` | Refine one soil grid, preserve every species root profile on that grid, and build each species' transient roots/phenology using its own age and growth parameters. |
| `func_DiffEq_definition_p.jl` | Build typed shared parameters plus species parameters and preallocated caches. Use a layers-by-species uptake matrix and separate resistance/conductance caches. Preserve `develop`'s in-place calculations. |
| `func_MSB_functions.jl` / `module_PET.jl` / `module_EVP.jl` | Pass the complete species configuration to canopy, `SRSC`, `PLNTRES` and `TBYLAYER` calculations. For a mixed canopy, calculate shared light/aerodynamic constraints before species exchange. Split `MSBPREINT` so species interception/throughfall is computed before a single shared ground snow update. Area-weight tile `SNOEN`/`PSNVP` inputs before applying snow constraints; decide soil-evaporation snow cover from the shared state/input. |
| `func_DiffEq_definition_u0.jl`, callbacks and RHS | Give each species interception, xylem and uptake state/diagnostics. Keep one ground snow/heat/liquid-water/isotope state by default, independent of canopy mode; update it once per interval. Initialize all isotope states, set all derivative fields on every RHS call, and reset/save every species accumulator. Derive stand totals for reporting instead of integrating duplicate aggregate plant stores. |
| Soil coupling in callbacks / `func_DiffEq_definition_f.jl` | Compute all species demands against the same current soil state before advancing it. Sum uptake once into the shared soil balance; retain accepted species fluxes through each soil/transport step. Any supply limiting must update species fluxes and reported totals consistently, without giving earlier species priority. |
| Isotope callbacks / `module_ISO.jl` | Track each species' accepted uptake by layer and mix isotope atom fractions using those fluxes. Update its own xylem reservoir on the matching area basis. Define no-uptake behavior; account explicitly for signed root outflux when hydraulic redistribution is enabled. |
| `func_DiffEq_definition_ode.jl` and water balance helpers | Include every physical species/tile amount in solver norms and balance diagnostics without counting aggregate aliases twice. |
| `func_postprocess.jl` and input export | Expose stand totals and species series with IDs and area units. Keep legacy columns for one species. Reject unsupported multispecies export to the single-vegetation R interface. |

For the first hydraulic implementation, preserve BROOK90's existing no-capacitance,
daily demand/supply approximation per species. Share soil storage and potentials,
not plant resistances or stomatal parameters. Re-evaluating uptake within soil
substeps is a subsequent model decision; do not claim a fully coupled hydraulic
competition model from independently calculated daily potentials alone.

## Gaps in the preserved prototype

The worktree prototype is retained as implementation material. Review found:

- `define_LWFB90_p` uses the first species to supply the legacy parameter tuple.
  The loader copies site defaults into species parameters and subsequently drops
  the separate site record.
- The species loop selects `GLMAX`, `GLMIN`, `PSICR`, `MXKPL`, `MXRTLN` and some
  canopy inputs, but still passes common `R5`, `CVPD`, `RM`, `TL`, `T1`, `T2`, `TH`,
  `RTRAD`, `NOOUTF`, `RHOTP` and interception coefficients. Accepting separate CSV
  values therefore does not ensure separate species behavior.
- Species xylem initialization hardcodes its amount and leaves species isotope
  initialization incomplete. Stand interception initialization and updates use
  inconsistent weighting, and initial aggregate isotope values use the first
  species rather than a mass-weighted mixture.
- Species uptake is weighted to stand area before xylem turnover, while
  `VXYLEM_mm` is left on the tile basis. The prototype runs separate snowpack
  updates from the same old state and averages their updated results. For the
  proposed shared ground model, replace that with input aggregation followed by
  one snowpack update. The existing RHS/norm do not cover all added plant fields.
- The one-species equivalence test runs the same legacy loading path twice. The
  two-species tests use identical traits and the first 20 days of winter. These
  tests cannot establish independent stomatal or hydraulic control.

The rebase conflict resolution retains the prototype alongside the new
`SLVPDEPTH_m`/`aux_du_SLVPI` interface and isotope/cache improvements from `develop`.
It is not a validation of the prototype's numerical results.

## Acceptance criteria for implementation

- Legacy single-species regression results and a genuine named one-species CSV
  input agree within solver tolerances, including nonzero interception and isotope
  initial states.
- Splitting one tile into identical 50/50 and 30/70 species conserves stand water,
  interception, snow and isotope trajectories, including species xylem turnover,
  during an active growing season.
- Distinct root profiles use their intended layers. Changing one species' LAI,
  `R5`, `CVPD`, temperature thresholds, `RTRAD`, `MXKPL`, `FXYLEM`, `PSICR`,
  `NOOUTF`, interception coefficients or `VXYLEM_mm` changes that component's
  calculation. Effects on other species occur through declared shared coupling.
- Species permutation, zero-area species, no uptake and dry soil do not change
  stand results spuriously or introduce undefined state. Combined withdrawal
  respects available soil water and all reported sinks match the water balance.
- Shared ground snow is updated once with conserved area-weighted water/energy
  inputs. Snowfall, rain retention, sublimation, melt and recharge budgets close;
  no ground store or flux is counted once per species. Heterogeneous-canopy tests
  cover cold-content removal, liquid retention and snow disappearance, and distinguish
  the shared update from averaging locally constrained snowpack updates.
- Soil, ground snow and plant isotope mass balances close; totals are mixed in
  atom-fraction space, and signed hydraulic redistribution follows the selected
  transport model.
- Load/save/remake round trips preserve independent species inputs and ownership;
  validation catches misspelled IDs, wrong categories and missing required values.
- Unit, integration, regression and allocation checks run on the revised solver.
  None of these numerical criteria is asserted as passed by this design revision.
