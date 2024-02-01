#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH -c 5
#SBATCH --mem-per-cpu=15G
#SBATCH -o star_index_%A.out

function usage() {
  printf "
    Utility script to submit a SLURM job. It generates STAR index.

    Arguments:

      -l [read_len]: read length - 1 for STAR index
      -d [project_dir]: path to a project directory

    Usage: sbatch star_index.sh -l read_len -d /path/to/project_directory/
    " 1>&2
}

module load Python                  # Python installed as env module
conda activate /path/to/conda/envs/ # must have STAR and samtools installed

while getopts "h?l:d:" option; do
  case "${option}" in

  h)
    usage
    exit 1
    ;;
  l)
    read_len=${OPTARG}
    ;;
  d)
    project_dir=${OPTARG}
    ;;
  esac
done

# create a star_out/ subdirectory if does not exist
if [ ! -d "${project_dir}/star_out" ]; then
  mkdir "${project_dir}/star_out"
fi

# create a star_index/ subdirectory if does not exist
if [ ! -d "${project_dir}/star_out/star_index" ]; then
  mkdir "${project_dir}/star_out/star_index"
fi

RELEASE="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"
GTF="gencode.v38.primary_assembly.annotation.gtf.gz"
GENOME="GRCh38.primary_assembly.genome.fa.gz"

# download annotation and genome files
curl "${RELEASE}/${GTF}" -o "${project_dir}/annotation.gtf.gz" &&
curl "${RELEASE}/${GENOME}" -o "${project_dir}/genome.fa.gz" &&
gunzip "${project_dir}/annotation.gtf.gz" "${project_dir}/genome.fa.gz" &&
STAR \
    --runThreadN 5 \
    --runMode genomeGenerate \
    --genomeDir "${project_dir}/star_out/star_index" \
    --sjdbGTFfile "${project_dir}/annotation.gtf" \
    --genomeFastaFiles "${project_dir}/genome.fa" \
    --sjdbOverhang $read_len &&
rm "${project_dir}/annotation.gtf" "${project_dir}/genome.fa"
