# pipeline/config.sh - paths and threads for the clean pipeline. Source it before any step:
#   source pipeline/config.sh
# Override any variable in the environment before sourcing (e.g. CH06_WORK=/data2/... source pipeline/config.sh).
P=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)                       # pipeline scripts
export CH06_WORK=${CH06_WORK:-$(cd "$P/.." && pwd)/run}                # all outputs (not tracked by git)
export CH06_INPUTS=${CH06_INPUTS:-$CH06_WORK/inputs/chapter06-inputs}  # Zenodo package 10.5281/zenodo.23690055
export GTDB_SSU=${GTDB_SSU:-$CH06_WORK/inputs/ssu_all_r226.fna.gz}     # GTDB r226 SSU
export GTDB_TREE=${GTDB_TREE:-$CH06_WORK/inputs/bac120_r226.tree}      # GTDB r226 bac120 tree
export SILVA_NR99=${SILVA_NR99:-$CH06_WORK/inputs/SILVA_138.2_NR99.fasta}  # SILVA 138.2 NR99 (headers: acc Bacteria;...;Genus;organism)
export THREADS=${THREADS:-16}
export CH06_PIPELINE=$P
mkdir -p "$CH06_WORK"/{inputs,refs,work,results,logs,tree,place,figures}
