#!/usr/bin/env bash
#SBATCH --qos=regular
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH -c 2
#SBATCH -o download_fastqs_%A_%a.out

function usage(){
  printf "
    Utility script to submit a SLURM job array of N jobs. It downloads FASTQ files with SRA toolkit and performs quality
    control with fastqc.

    Arguments:

      -t [table]: path to a job array table formatted as follows: JOB_NUMBER\tSRA\tSAMPLE_NAME;
      -d [project_dir]: path to a project directory

    Usage: sbatch --array 1-N download_fastqs.sh -t /path/to/table.tab -d /path/to/project_directory/
    " 1>&2
}

module load Python                           # Python installed as env module
conda activate /path/to/conda/envs/           # must have fastqc installed
module load SRA-Toolkit/3.0.0-centos_linux64 # SRA-Toolkit installed as env module

while getopts "h?t:d:" option; do
  case "${option}" in

  h)
    usage
    exit 1
    ;;
  t)
    table=${OPTARG}
    ;;
  d)
    project_dir=${OPTARG}
    ;;
  esac
done

# create a fastq/ subdirectory if does not exist
if [ ! -d "${project_dir}/fastq" ]; then
  mkdir "${project_dir}/fastq"
fi

# create a fastqc/ subdirectory if does not exist
if [ ! -d "${project_dir}/fastq/fastqc" ]; then
  mkdir "${project_dir}/fastq/fastqc"
fi

sra=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $table)
sample_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $table)

# download FASTQ file/s from SRA repository
prefetch "$sra" \
  --output-directory "${project_dir}/fastq/" \
  --output-file "${sample_name}" &&
fastq-dump --outdir "${project_dir}/fastq/" \
    --accession "${sample_name}" \
    --gzip \
    --readids \
    --dumpbase \
    --skip-technical \
    --split-files &&
rm "${project_dir}/${sample_name}" &&

# run fastqc
if [ -f "${project_dir}/fastq/${sample_name}_2.fastq.gz" ]; then # check for the paired-end reads
  fastqc "${project_dir}/fastq/${sample_name}_1.fastq.gz" "${project_dir}/fastq/${sample_name}_2.fastq.gz" \
    --outdir "${project_dir}/fastq/fastqc" \
    --threads 2
else
  fastqc "${project_dir}/fastq/${sample_name}_1.fastq.gz" \
    --outdir "${project_dir}/fastq/fastqc" \
    --threads 2
fi
