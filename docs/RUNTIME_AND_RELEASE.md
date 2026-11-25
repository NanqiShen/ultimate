# Runtime and release checks

Use the installable `ultimate-bioinfo` package for the CLI and the source checkout
for cluster launch templates. This doesn't replace the existing `preflight`,
production approval, Slurm, or delivery gates.

1. `ultimate doctor --config path/to/project.yaml` — metadata-only installation
   check (Python, dependencies, packaged templates, configuration, output
   parent, and declared sample sheet readability). Returns JSON and a nonzero
   exit on a failed required check. It doesn't inspect sample contents or
   prove that backends or production approval are ready.
2. `ultimate debug --config path/to/project.yaml --threads 8` — preview the
   budget without creating a run. `ultimate debug --trace ./trace.jsonl` —
   find the first failed stage from one explicit trace. Optional run tracing
   is enabled by `ultimate run --config path/to/project.yaml --debug-log ./trace.jsonl`.
   Use a new destination; traces contain stage/status/timing and exception
   types, not raw inputs or exception messages.
3. CPU budget precedence: `--threads`, `ULTIMATE_THREADS`, YAML
   `runtime.threads`, then Slurm's per-task allocation or a local default of
   at most four. The result is capped by the detected CPU affinity and the
   Slurm per-task limit. A Slurm job without a declared count defaults to one.
   The effective budget is recorded as `thread_budget` in the run manifest.
   Supported native thread pools and subprocess environment variables are
   limited per run. Module order remains sequential and backends with their
   own explicit worker settings retain those settings. This controls thread
   oversubscription; it is not a measured throughput improvement.
4. From the package checkout: `uv build` produces the wheel and sdist. Install
   the wheel into an isolated Python 3.11+ environment and invoke
   `ultimate doctor` from outside the checkout, to verify packaged resources.
   Publish only the intended release branch and do not include environments,
   logs, patient inputs, or analysis outputs in a release.
