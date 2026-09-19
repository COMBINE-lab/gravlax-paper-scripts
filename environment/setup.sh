# Source this before any project command: `source env/setup.sh`
# Project-local toolchain homes below the project root (home NFS is nearly full;
# ~/.cargo and ~/.rustup are dangling symlinks; the default TMPDIR is on a nearly full volume).
export PROJ_ROOT=${GRAVLAX_PROJECT_ROOT:?set GRAVLAX_PROJECT_ROOT to the analysis project root}
export CARGO_HOME=$PROJ_ROOT/env/cargo
export RUSTUP_HOME=$PROJ_ROOT/env/rustup
export TMPDIR=$PROJ_ROOT/env/tmp
export CONDA_PKGS_DIRS=$PROJ_ROOT/env/conda-pkgs
export UV_CACHE_DIR=$PROJ_ROOT/env/uv-cache
mkdir -p "$TMPDIR" "$CONDA_PKGS_DIRS"

# Toolchains: rust 1.98 is shared with the coresets project (same host, same rustup home layout).
export PATH=$CARGO_HOME/bin:$PROJ_ROOT/env/sc/bin:$PATH

# Shared 256-core box, no scheduler: cap parallelism (be a considerate neighbour).
export CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-32}
export AIE_THREADS=${AIE_THREADS:-32}

# Fast-iteration default: restrict the whole pipeline to one gene-dense chromosome.
# Unset (`unset AIE_CHROM`) for whole-genome runs.
export AIE_CHROM=${AIE_CHROM:-chr19}
