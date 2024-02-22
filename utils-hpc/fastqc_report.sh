#!/usr/bin/env bash
#SBATCH --qos=regular
#SBATCH --time=12:00:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10G
#SBATCH -c 2
#SBATCH -o fastqc_report_%A.out

function usage() {
  printf "
    Utility script to submit a SLURM job. It performs quality control of FASTQ files with fastqc and generates a report
    with multiqc.

    Arguments:
      -d [fastq_dir]: path to a directory with FASTQ files in fastq.gz format

    Usage: sbatch fastqc.sh -d /path/to/fastq_directory/
    " 1>&2
}

module load Python                  # Python installed as env module
conda activate /path/to/conda/envs/ # must have fastqc and multiqc installed

while getopts "h?d:" option; do
  case "${option}" in

  h)
    usage
    exit 1
    ;;
  d)
    fastq_dir=${OPTARG}
    ;;
  esac
done

if [ ! -d "${fastq_dir}/fastqc" ]; then
  mkdir "${fastq_dir}/fastqc"
fi

cd "${fastq_dir}/"
fastqc *fastq.gz -t 2 --outdir "${fastq_dir}/fastqc" &&
multiqc "${fastq_dir}/fastqc" --verbose --module fastqc --outdir "${fastq_dir}/fastqc"
