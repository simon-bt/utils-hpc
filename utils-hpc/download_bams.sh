#!/usr/bin/env bash
#SBATCH --qos=regular
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH -c 2
#SBATCH -o download_bams_%A_%a.out

function usage(){
  printf "
    Utility script to submit a SLURM job array of N jobs. It downloads BAM files with SRA toolkit and converts them to
    FASTQ files.

    Arguments:

    -t [table]: path to a job array table formatted as follows: JOB_NUMBER\tSRA\tSAMPLE_NAME\tORIENTATION;
    -d [project_dir]: path to a project directory

    Usage: sbatch --array 1-N download_bams.sh -t /path/to/table.tab -d /path/to/project_directory/
    " 1>&2
}

module load Python                           # Python installed as env module
conda activate /path/to/conda/envs/        # must have fastqc installed
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

# create a bams/ subdirectory if does not exist
if [ ! -d "${project_dir}/bams" ]; then
  mkdir "${project_dir}/bams"
fi

# create a fastq/ subdirectory if does not exist
if [ ! -d "${project_dir}/fastq" ]; then
  mkdir "${project_dir}/fastq"
fi

sra=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $table)
sample_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $table)
orientation=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $table)

# download BAM file from SRA repository
prefetch "$sra" \
  --output-directory "${project_dir}/bams/" \
  --output-file "$sample_name" &&
sam-dump "$sample_name" \
    --seqid \
    --header \
    --prefix "$sra" \
    --output-file "${project_dir}/bams/${sample_name}.bam" &&
rm "${project_dir}/${sample_name}" &&

# convert BAM file into FASTQ file/s
if [ "$orientation" == "paired" ]; then
  # first, sort BAM file
  samtools sort -n -o "${project_dir}/bams/${sample_name}.sorted.bam" "${project_dir}/bams/${sample_name}.bam" &&
  # then, split BAM file into paired end files
  bedtools bamtofastq \
      -i "${project_dir}/bams/${sample_name}.sorted.bam" \
      -fq "${project_dir}/fastq/${sample_name}_1.fastq" \
      -fq2 "${project_dir}/fastq/${sample_name}_2.fastq" &&
  gzip "${project_dir}/fastq/${sample_name}_1.fastq" "${project_dir}/fastq/${sample_name}_2.fastq"
else
  # convert BAM file into single end file
  bedtools bamtofastq \
    -i "${project_dir}/bams/${sample_name}.bam" \
    -fq "${project_dir}/fastq/${sample_name}_1.fastq" &&
  gzip "${project_dir}/fastq/${sample_name}_1.fastq"
fi &&
samtools index -@ 2 "${project_dir}/bams/${sample_name}*.bam"
