#!/usr/bin/env bash
#SBATCH --qos=regular
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH -c 2
#SBATCH -o parse_bams_%A_%a.out

function usage(){
  printf "
    Utility script to submit a SLURM job array of N jobs. It parses BAM files for selected chromosomes.

    Arguments:

    -t [table]: path to a job array table formatted as follows: JOB_NUMBER\tSRA\tSAMPLE_NAME\tORIENTATION;
    -c [chromosomes]: chromosomes, comma-separated
    -d [project_dir]: path to a directory with BAM files

    Usage: sbatch --array 1-N parse_bams.sh -t /path/to/table.tab -c chromosomes -d /path/to/project_dir/
    " 1>&2
}

module load Python                           # Python installed as env module
conda activate /path/to/conda/envs/        # must have fastqc installed

while getopts "h?t:c:d:" option; do
  case "${option}" in

  h)
    usage
    exit 1
    ;;
  t)
    table=${OPTARG}
    ;;
  c)
    chromosomes=${OPTARG}
    ;;
  d)
    project_dir=${OPTARG}
    ;;
  esac
done

sample_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $table)

IFS=',' read -r -a array <<< "$chromosomes"
for chromosome in "${array[@]}"
  do
    dir="${project_dir}/${chromosome}"
    if [ ! -d "$dir" ]; then
      mkdir "$dir"
    fi
    samtools view -@ 2 -b "${project_dir}/${sample_name}.Aligned.sortedByCoord.out.bam" "$chromosome" > \
      "${dir}/${sample_name}_${chromosome}.bam" &&
    samtools index -@ 2 "${dir}/${sample_name}_${chromosome}.bam"
done
