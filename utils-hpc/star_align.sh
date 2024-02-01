#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH -c 5
#SBATCH --mem-per-cpu=15G
#SBATCH -o star_align_%A_%a.out

function usage(){
  printf "
    Utility script to submit a SLURM job array of N jobs. It aligns RNAseq data with STAR and generates BAM and
    bedGraph files.

    Arguments:

      -t [table]: path to a job array table formatted as follows: JOB_NUMBER\tSRA\tSAMPLE_NAME;
      if using in-house generated data, keep the 2nd column empty
      -i [index]: path to a STAR index directory
      -d [project_dir]: path to a project directory with fastq subdirectory
      -s [strandness]: strandedness of wiggle/bedGraph output, Stranded/Unstranded

    Usage: sbatch --array 1-N star_align.sh -t /path/to/table.tab -i /path/to/star_index/ -d /path/to/project_directory/ -s strandness
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

sample_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $table)

# create a sample_name/ output subdirectory if does not exist
if [ ! -d "${project_dir}/star_out/${sample_name}" ]; then # check for paired-end reads
  mkdir "${project_dir}/star_out/${sample_name}"
fi

# project_dir must contain fastq/ subdirectory
if [ -f "${project_dir}/fastq/${sample_name}_2.fastq.gz" ]; then
  STAR \
    --runThreadN 5 \
    --readFilesIn "${project_dir}/fastq/${sample_name}_1.fastq.gz" "${project_dir}/fastq/${sample_name}_2.fastq.gz" \
    --genomeDir "$index" \
    --genomeLoad NoSharedMemory \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --outFilterType BySJout \
    --outSAMattributes Standard \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${project_dir}/star_out/${sample_name}/${sample_name}."
else
  STAR \
    --runThreadN 5 \
    --readFilesIn "${project_dir}/fastq/${sample_name}_1.fastq.gz" \
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
