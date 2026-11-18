# landisR.bda

Base BDA 4.0.1 port for landisR 0.1.0.9011 and
landisR.succession.biomass 0.1.0.9002 for full coupled parity.

The algorithm follows the LANDIS-II Foundation's
[Base BDA 4.0.1 source](https://github.com/LANDIS-II-Foundation/Extension-Base-BDA/tree/f92a067eca53df7122e96c361c8c1638608310fd).
See the upstream notices and Apache 2.0 license in that repository.

Read the Base BDA file with `bda_read()`. For a coupled scenario, initialize
the simulation with `recruitment_rng = "landis_v7"` and pass it to
`bda_prepare_coupling()`. Include the prepared BDA extension in the desired
position of the disturbance list passed to `landisR::landis_sim()`.
The port supports fixed and variable
pulses, uniform and normal intervals, start and end years, and the background
ROS. It requires a Base BDA timestep equal to the canonical Biomass Succession
timestep.

bda_log() returns the event summary. Each run event keeps the damaged sites
and the removed cohorts. For diagnostics, select cell IDs through
trace_cells. With trace_cells = "all", bda_maps() reconstructs the four
map matrices. Use this option on validation grids when the full diagnostics
are needed.

Dispersal supports MaxRadius and the 4N, 8N, 12N, and 24N percolation templates.
For spatial diagnostics, use trace_dispersal = TRUE to retain the selected
epicenters, candidate counts, and RNG positions. With trace_cells enabled,
spatial events also report zone membership for the requested cells.

The resource neighborhood supports uniform, linear, and Gaussian weights.
`NeighborSpeedUp` controls which positions contribute after the full
LANDIS-II shuffle. Sets no larger than the stride retain all their neighbors.
Linear edges at the specified radius retain zero-weight neighbors in the RNG
sequence. Gaussian weights use the integer half-radius and single-precision
rounding of Base BDA 4.0.1.

Set trace_calendar = TRUE to collect interval and RNG diagnostics without
any per-cell data. Retrieve the table with bda_calendar(). The first BDA
step runs the initial scheduling draw before that step's succession
when using bda_extension() directly. For LANDIS-II parity with other
disturbances, call bda_prepare_coupling() on the year-zero simulation,
then use its returned simulation and extension together. This draws the
initial interval before any disturbance runs. Resume using
readRDS() and landisR::landis_sim().
The extension checkpoint records the pending interval, elapsed time,
announced next year, and the normal generator's cached Boolean bits. Spatial
checkpoints also retain the previous zone, the sites damaged in the latest
epidemic, and the epicenter count.
The MT19937 stream itself lives in the shared simulation.

DisturbanceModifiers supports effects from fire severity and harvest
prescriptions. Wrap the Base Fire extension with
`bda_disturbance_source(fire, "fire")` so the fire report is visible in
the configured disturbance order. Biomass Harvest 4.4 lacks the
native site variable used by Base BDA 4.0.1 for harvest modifiers. The
optional `bda_disturbance_source(harvest, "harvest")` bridge supplies
complete cohort removals when requested.

## Multiple agents

Put the first agent file on the `BDAInputFiles` line and each additional
file on a separate line. Agents run in file order and share the forest CSR
and random stream. Each keeps its own schedule and dispersion state.
The Boolean cache used for normal intervals is shared across agents.

A modifier selector matches `BDA` or the case-sensitive name of the
last agent that killed a cohort at the site. The multiplier uses the elapsed
time since that event. A mortality event in the current step replaces the
previous record for a later agent even though its elapsed time is zero.
There is one sparse, checkpointed history of the latest BDA mortality
when a rule uses it.

Use `bda_log(run)` for the ordered event summary. In a multi-agent run an
individual event can be retrieved from `run$events$bda[["10"]][["budworm"]]`.
This event can be passed to `bda_maps()` when the run traces all cells.
