#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10G
#SBATCH -c 5
#SBATCH -o vast_combine_%A.out

function usage() {
  printf "
    Utility script to submit a SLURM job. It merges vast-tools quantification output into one table. If merge.tab with
    SAMPLE_NAME\t\GROUP_NAME format exists in the project directory, vast-tools first merges samples before combining the output.

    Arguments:

      -s [species]: VastDB species code, for more information check: https://github.com/vastgroup/vast-tools#vastdb-libraries
      -d [project_dir]: path to a project directory

    Usage: sbatch vast_combine.sh -s species_code -d /path/to/project_directory/
    " 1>&2
}

module load VAST-TOOLS/2.5.1 # vast-tools installed as env module: https://github.com/vastgroup/vast-tools#installation

while getopts "h?s:d:" option; do
  case "${option}" in

  h)
    usage
    exit 1
    ;;
  s)
    species=$OPTARG
    ;;
  d)
    project_dir=$OPTARG
    ;;
  esac
done

VASTDB="/path/to/VAST-TOOLS/2.5.1/VASTDB/" # path to VastDB

if [ -f "${project_dir}/merge.tab" ]; then
  # first, merge vast-tools runs
  vast-tools merge \
    --outDir "$project_dir" \
    --groups "${project_dir}/merge.tab" \
    --sp "$species" \
    --move_to_PARTS &&
    # then, run vast-tools combine with merged samples
  vast-tools combine \
      -o "${project_dir}/vast_out/" \
      -sp "$species" \
      --dbDir "$VASTDB" \
      --cores 5
else
  # combine without merging
  vast-tools combine \
    -o "${project_dir}/vast_out" \
    -sp "$species" \
    --dbDir "$VASTDB" \
    --cores 5
fi
