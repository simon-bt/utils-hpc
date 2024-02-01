#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH -c 5
#SBATCH --mem-per-cpu=15G
#SBATCH -o star_align_multi_%A_%a.out

function usage() {
  printf "
    Utility script to submit a SLURM job array of N jobs. It aligns RNAseq data with STAR and generates BAM and
    bedGraph files using multiple FASTQ files per mate.

    Arguments:

      -t [table]: path to a job array table formatted as follows: JOB_NUMBER\tSAMPLE_NAME\tSAMPLE_1,SAMPLE_2,SAMPLE_n;
      if using in-house generated data, keep the 2nd column empty
      -i [index]: path to a STAR index directory
      -d [project_dir]: path to a project directory with fastq subdirectory
      -s [strandness]: strandedness of wiggle/bedGraph output, Stranded/Unstranded

    Usage: sbatch --array 1-N star_align_multi.sh -t /path/to/table.tab -i /path/to/star_index/ -d /path/to/project_directory/ -s strandness
    " 1>&2
}

module load Python                  # Python installed as env module
conda activate /path/to/conda/envs/ # must have STAR and samtools installed

while getopts "h?t:i:d:s:" option; do
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
  d)
    project_dir=${OPTARG}
    ;;
  s)
    strandness=${OPTARG}
    ;;
  esac
done

# create a star_out/ subdirectory if does not exist
if [ ! -d "${project_dir}/star_out" ]; then
  mkdir "${project_dir}/star_out"
fi

sample_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $table)

# get names and paths of files for merging
file_names=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $table)
file_paths=$(awk '$0="'${project_dir}'/fastq/"$0' <<<"${file_names//,/,${project_dir}/fastq/}")
# get first file for merging
IFS=',' read -r file _ <<<"$file_names"

# create a sample_name/ output subdirectory if does not exist
if [ ! -d "${project_dir}/star_out/${sample_name}" ]; then
  mkdir "${project_dir}/star_out/${sample_name}"
fi

# project_dir must contain fastq/ subdirectory
if [ -f "${project_dir}/fastq/${file}_2.fastq.gz" ]; then # check for paired-end reads
  # get all forward mate files
  mates1_files="${file_paths//,/_1.fastq.gz,}_1.fastq.gz"
  # get all reverse mate files
  mates2_files="${file_paths//,/_2.fastq.gz,}_2.fastq.gz"
  STAR \
    --runThreadN 5 \
    --readFilesIn $mates1_files $mates2_files \
    --genomeDir "$index" \
    --genomeLoad NoSharedMemory \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --outFilterType BySJout \
    --outSAMattributes Standard \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${project_dir}/star_out/${sample_name}/${sample_name}."
else
  mates1_files="${file_paths//,/_1.fastq.gz,}_1.fastq.gz"
  STAR \
    --runThreadN 5 \
    --readFilesIn $mates1_files \
    --genomeDir "$index" \
    --genomeLoad NoSharedMemory \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --outFilterType BySJout \
    --outSAMattributes Standard \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${project_dir}/star_out/${sample_name}/${sample_name}."
fi &&
samtools index -@ 5 "${project_dir}/star_out/${sample_name}/${sample_name}.Aligned.sortedByCoord.out.bam" &&
STAR \
    --runThreadN 5 \
    --runMode inputAlignmentsFromBAM \
    --inputBAMfile "${project_dir}/star_out/${sample_name}/${sample_name}.Aligned.sortedByCoord.out.bam" \
    --outWigType bedGraph \
    --outWigNorm RPM \
    --outFileNamePrefix "${project_dir}/star_out/${sample_name}/${sample_name}." \
    --outWigStrand "$strandness"
