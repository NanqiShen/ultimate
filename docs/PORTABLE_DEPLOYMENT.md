# Portable Ultimate checkout

The package root is the directory containing its own `pyproject.toml`,
`src/ultimate`, `config`, and `slurm` directories (the `ultimate/` subdirectory
in the agent_hpc Git repository). This directory can be placed on any filesystem.

- Start with Python 3.11+ and install the package from this directory
  (`python -m pip install -e .`). Recreate environments instead of moving conda binaries.
  For 10x Matrix Market and scRNA backends, install the optional `scrna` extra.
- Project YAML input and output paths are relative to that YAML's directory.
  An explicit `project.server_root: ..` is relative to the same location and
  overrides the environment. A missing root uses `ULTIMATE_ROOT` or checkout
  discovery. CLI `--root` takes precedence for commands accepting it.
- Set `ULTIMATE_RESOURCES` for reference files (the YAML shorthand
  `@resources/`). `@project/` is relative to the normalized project root.
  Set `ULTIMATE_EXTERNAL_ROOT` for optional, separately managed legacy inputs.
- For Slurm, export `ULTIMATE_ROOT` or submit from the package root. Examples
  are templates, not automatically submitted jobs. Scheduler stdout/stderr
  default to the submission directory; pass explicit `--output`/`--error`
  to put them elsewhere. Set `CONDA_SH` and `ULTIMATE_ENV` for site-specific
  environments. Rendered handoff examples are relative to the project root;
  run their commands from that root after review and substitute real inputs.
- A wheel install outside a checkout uses the caller's working directory unless
  `ULTIMATE_ROOT` is set. For a moved checkout re-install editable entrypoints.
- Previously generated job files and run reports retain their storage provenance.
  Recreate job envelopes for the new deployment; do not rewrite archived outputs.
