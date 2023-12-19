#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10G
#SBATCH -c 5
#SBATCH -o salmon_quant_%A_%a.out

: '
  Utility script to submit a SLURM job array of N jobs with salmon quant.

  Arguments:

    -t [table]: path to a job array table formatted as follows: JOB_NUMBER\tSRA\tSAMPLE_NAME;
    if using in-house generated data, keep the 2nd column empty
    -i [index]: path to a Salmon index directory
    -b [bootstraps]: number of bootstraps to perform
    -d [project_dir]: path to a project directory with fastq subdirectory

  Usage: sbatch --array 1-N salmon_quant.sh -t /path/to/table.tab -i /path/to/salmon_index/ -b n_bootstraps -d /path/to/project_directory/
'

module load Python                  # Python installed as env module
conda activate /path/to/conda/envs/ # must have Salmon installed

while getopts t:i:b:d: option; do
  case "${option}" in
  t)
    table=${OPTARG}
    ;;
  i)
    index=${OPTARG}
    ;;
  b)
    bootstraps=${OPTARG}
    ;;
  d)
    project_dir=${OPTARG}
    ;;
  esac
done

# create a salmon_out/ subdirectory if does not exist
if [ ! -d "${project_dir}/salmon_out" ]; then
  mkdir "${project_dir}/salmon_out"
fi

sample_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $table)

# project_dir must contain fastq/ subdirectory
if [ -f "${project_dir}/fastq/${sample_name}_2.fastq.gz" ]; then # check for paired-end reads
  salmon quant \
    --index "$index" \
    --libType A \
    --mates1 "${project_dir}/fastq/${sample_name}_1.fastq.gz" \
    --mates2 "${project_dir}/fastq/${sample_name}_2.fastq.gz" \
    --validateMappings \
    --numBootstraps $bootstraps \
    --threads 5 \
    --output "${project_dir}/salmon_out/${sample_name}"
else
  salmon quant \
    --index "$index" \
    --libType A \
    -r "${project_dir}/fastq/${sample_name}_1.fastq.gz" \
    --validateMappings \
    --numBootstraps $bootstraps \
    --threads 5 \
    --output "${project_dir}/salmon_out/${sample_name}"
fi
