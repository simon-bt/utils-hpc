#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10G
#SBATCH -c 5
#SBATCH -o salmon_quant_multi_%A_%a.out

function usage(){
  printf "
    Utility script to submit a SLURM job array of N jobs. It quantifies transcript abundances with salmon quant using
    multiple FASTQ files per mate.

    Arguments:

      -t [table]: path to a job array table formatted as follows: JOB_NUMBER\tSAMPLE_NAME\tSAMPLE_1,SAMPLE_2,SAMPLE_n
      -i [index]: path to a Salmon index directory
      -b [bootstraps]: number of bootstraps to perform
      -d [project_dir]: path to a project directory with fastq subdirectory

    Usage: sbatch --array 1-N salmon_quant_multi.sh -t /path/to/table.tab -i /path/to/salmon_index/ -b n_bootstraps -d /path/to/project_directory/
    " 1>&2
}

module load Python                  # Python installed as env module
conda activate /path/to/conda/envs/ # must have Salmon installed

while getopts "h?t:i:b:d:" option; do
  case "${option}" in

  h)
    usage
    exit 1
    ;;
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

sample_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $table)

# get names and paths of files for merging
file_names=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $table)
file_paths=$(awk '$0="'${project_dir}'/fastq/"$0' <<<"${file_names//,/ ${project_dir}/fastq/}")
# get first file for merging
IFS=',' read -r file _ <<<"$file_names"

# project_dir must contain fastq/ subdirectory
if [ -f "${project_dir}/fastq/${file}_2.fastq.gz" ]; then # check for paired-end reads
  # get all forward mate files
  mates1_files="${file_paths// /_1.fastq.gz }_1.fastq.gz"
  # get all reverse mate files
  mates2_files="${file_paths// /_2.fastq.gz }_2.fastq.gz"
  salmon quant \
    --index "$index" \
    --libType A \
    --mates1 $mates1_files \
    --mates2 $mates2_files \
    --validateMappings \
    --numBootstraps $bootstraps \
    --threads 5 \
    --output "${project_dir}/salmon_out/${sample_name}"
else
  # get all single mate files
  mates1_files="${file_paths// /_1.fastq.gz }_1.fastq.gz"
  salmon quant \
    --index "$index" \
    --libType A \
    -r $mates1_files \
    --validateMappings \
    --numBootstraps $bootstraps \
    --threads 5 \
    --output "${project_dir}/salmon_out/${sample_name}"
fi
