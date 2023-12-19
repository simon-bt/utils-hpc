#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10G
#SBATCH -c 5
#SBATCH -o salmon_index_%A.out

: '
  Utility script to submit a SLURM job with salmon index.

  Arguments:

    -k [kmer]: kmer length for Salmon index
    -d [project_dir]: path to a project directory

  Usage: sbatch salmon_index.sh -k kmer -d /path/to/project_directory/
'

module load Python                  # Python installed as env module
conda activate /path/to/conda/envs/ # must have Salmon installed

while getopts k:d: option; do
  case "${option}" in
  k)
    kmer=${OPTARG}
    ;;
  d)
    project_dir=${OPTARG}
    ;;
  esac
done

RELEASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"
TRANSCRIPTS="gencode.v38.transcripts.fa.gz"

# download transcripts files
curl "${RELEASE}/${TRANSCRIPTS}" -o "${project_dir}/transcripts.fa.gz" &&
# generate salmon index
salmon index \
    --transcripts "${project_dir}/transcripts.fa.gz" \
    --kmerLen $kmer \
    --gencode \
    --index "${project_dir}/salmon_index" \
    --threads 5 &&
rm "${project_dir}/transcripts.fa.gz"
