#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH -c 5
#SBATCH --mem-per-cpu=15G
#SBATCH -o generate_bigwigs_%A_%a.out

function usage(){
  printf "
    Utility script to submit a SLURM job array of N jobs. It sorts bedGraph files generated by STAR and converts them
    to bigWig files with bedGraphToBigWig.

    Arguments:

      -t [table]: path to a job array table formatted as follows: JOB_NUMBER\tSRA\tSAMPLE_NAME;
      if using in-house generated data, keep the 2nd column empty
      -f [index]: path to a genome.fa.fai file
      -d [project_dir]: path to a project directory with star_out subdirectory
      -s [strandness]: strandedness of wiggle/bedGraph output, Stranded/Unstranded
      -c [cleanup]: bool flag to cleanup directory, true/false

    Usage: sbatch --array 1-N generate_bigwigs.sh -t /path/to/table.tab -f /path/to/fa.fai/ -d /path/to/project_directory/ -s strandness -c cleanup
    " 1>&2
}

module load Python                  # Python installed as env module
conda activate /path/to/conda/envs/ # must have bedGraph and bedGraphToBigWig
module load OpenSSL/.1.1

while getopts "h?t:f:d:s:c:" option; do
  case "${option}" in

  h)
    usage
    exit 1
    ;;
  t)
    table=${OPTARG}
    ;;
  f)
    fai=${OPTARG}
    ;;
  d)
    project_dir=${OPTARG}
    ;;
  s)
    strandness=${OPTARG}
    ;;
  c)
    cleanup=${OPTARG}
    ;;
  esac
done

sample_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $table)
sample_dir="${project_dir}/star_out/${sample_name}"
bg="${sample_name}.Signal.Unique.str1.out.bg"

if [ "$strandness" == "Stranded" ]; then # sort and convert stranded bedGraphs to bigwigs
  str2_bg="${sample_name}.Signal.Unique.str2.out.bg"
  bedSort "${sample_dir}/${bg}" "${sample_dir}/${bg}" &&
  bedGraphToBigWig "${sample_dir}/${bg}" "${fai}" "${sample_dir}/${sample_name}_plus.bigwig"
  bedSort "${sample_dir}/${str2_bg}" "${sample_dir}/${str2_bg}" &&
  bedGraphToBigWig "${sample_dir}/${str2_bg}" "${fai}" "${sample_dir}/${sample_name}_minus.bigwig"
else
  bedSort "${sample_dir}/${bg}" "${sample_dir}/${bg}" &&
  bedGraphToBigWig "${sample_dir}/${bg}" "${fai}" "${sample_dir}/${sample_name}.bigwig"
fi
# cleanup
if [ "$cleanup" = true ]; then
  cd "${sample_dir}" || exit
  rm *.out.bg *.bam* # remove bedGraph and bam files
fi
