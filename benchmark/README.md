# Performance benchmarks

The benchmark suite measures the ODE solve for the `DAV2020-full` example over 30 days,
both with and without isotope transport. Loading input files and constructing the
discretized model happen in BenchmarkTools' `setup` phase and are not timed.

Run the suite from the repository root:

```sh
julia --project=. benchmark/runbenchmarks.jl
```

Optionally save the raw BenchmarkTools results for later analysis:

```sh
BENCHMARK_OUTPUT=benchmark-results.json julia --project=. benchmark/runbenchmarks.jl
```

`benchmark/benchmarks.jl` follows the conventional `SUITE::BenchmarkGroup` layout used
by PkgBenchmark.jl and AirspeedVelocity.jl. The benchmark GitHub workflow compares pull
requests with the target branch and publishes timing and allocation changes in the job
summary. It reports regressions instead of failing correctness tests, because timings on
shared CI runners are inherently noisy.
