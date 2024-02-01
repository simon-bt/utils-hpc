#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10G
#SBATCH -c 5
#SBATCH -o vast_align_%A_%a.out

function usage() {
  printf "
    Utility script to submit a SLURM job array of N jobs. It quantifies splicing levels with vast-tools align.

    Arguments:

      -t [table]: path to a job array table formatted as follows: JOB_NUMBER\tSRA\tSAMPLE_NAME;
      if using in-house generated data, keep the 2nd column empty
      -s [species]: VastDB species code, for more information check: https://github.com/vastgroup/vast-tools#vastdb-libraries
      -d [project_dir]: path to a project directory with fastq subdirectory

    Usage: sbatch --array 1-N run_vast_tools.sh -t /path/to/table.tab -s species_code -d /path/to/project_directory/
    " 1>&2
}

module load VAST-TOOLS/2.5.1 # vast-tools installed as env module

while getopts "h?t:s:d:" option; do
  case "${option}" in

  h)
    usage
    exit 1
    ;;
  t)
    table=${OPTARG}
    ;;
  s)
    species=${OPTARG}
    ;;
  d)
    project_dir=${OPTARG}
    ;;
  esac
done

VASTDB="/path/to/VAST-TOOLS/2.5.1/VASTDB/" # path to VastDB
sample_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $table)

# project_dir must contain fastq/ subdirectory
if [ -f "${project_dir}/fastq/${sample_name}_2.fastq.gz" ]; then # check for the paired-end reads
  vast-tools align \
    "${project_dir}/fastq/${sample_name}_1.fastq.gz" "${project_dir}/fastq/${sample_name}_2.fastq.gz" \
    --sp "$species" \
    --name "$sample_name" \
    --dbDir "$VASTDB" \
    --cores 5 \
    --resume
else
  vast-tools align \
    "${project_dir}/fastq/${sample_name}_1.fastq.gz" \
    --sp "$species" \
    --name "$sample_name" \
    --dbDir "$VASTDB" \
    --cores 5 \
    --resume
fi
