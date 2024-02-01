# utils-hpc

Utility scripts for SLURM jobs.

# Content

List of stand-alone utilities:

* [download_fastqs.sh](utils-hpc/download_fastqs.sh) - downloads FASTQ files from SRA repository with [SRA Toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/) and perform quality control with [FastQC](https://github.com/s-andrews/FastQC)
* [download_bams.sh](utils-hpc/download_bam.sh) - downloads BAM files from SRA repository with [SRA Toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/) and convert them to FASTQ files with [bedtools](https://github.com/arq5x/bedtools2)
* [star_index.sh](utils-hpc/star_index.sh) - generates [STAR](https://github.com/alexdobin/STAR) index
* [star_align.sh](utils-hpc/star_align.sh) - performs RNA-sequencing data alignment to [STAR](https://github.com/alexdobin/STAR) index and generates gene counts, BAM files and bigwig files
* [star_align_multi.sh](utils-hpc/star_align_multi.sh) - performs RNA-sequencing data alignment to [STAR](https://github.com/alexdobin/STAR) index using multiple FASTQ files per sample, and generates gene counts, BAM files and bigwig files
* [generate_bigwigs.sh](utils-hpc/generate_bigwigs.sh) - sorts STAR-generated bedGraph files and converts them to bigWig files
* [salmon_index.sh](utils-hpc/salmon_index.sh) - generates [Salmon](https://github.com/COMBINE-lab/salmon) index
* [salmon_quant.sh](utils-hpc/salmon_quant.sh) - quantifies transcript abundances from RNA-sequencing data using [Salmon](https://github.com/COMBINE-lab/salmon)
* [salmon_quant_multi.sh](utils-hpc/salmon_quant_multi.sh) - quantifies transcript abundances from RNA-sequencing data using [Salmon](https://github.com/COMBINE-lab/salmon) for merged FASTQ files
* [vast_align.sh](utils-hpc/vast_align.sh) - quantifies splicing from RNA-sequencing data using [VASTdb database](https://github.com/vastgroup/vast-tools)
* [vast_combine.sh](utils-hpc/vast_combine.sh) - combines vast-tools aling output into one table
* [parse_bams.sh](utils-hpc/parse_bams.sh) - parses BAM files for selected chromosomes

List of example config files:

* [config.tab](configs/config.tab) - config to download FASTQ or BAM files
* [config_multi.tab](configs/config_multi.tab) - config to execute Salmon quantification and STAR alignment with multiple 
FASTQ files per sample
* [merge.tab](configs/merge.tab) - config to merge vast-tools quantification files before combining into inclusion table

# License

This project is under MIT License.