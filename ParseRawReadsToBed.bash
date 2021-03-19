#!/bin/bash
# This script takes an sra file as input and runs the file through the SRA toolkit, trimmomatic, 
# bowtie, samtools, and bedtools to create a bed file.
# The first argument should be the input sra or gzipped fastq data, the second input should be the path to the fasta file of sequences for trimmomatic,
#   and the third argument should be the path to the basename for the bowtie2 index files.  (e.g. for files basename.1.bt2 in dir, give: dir/basename)

# Get trimmomatic's jar path.
trimmomaticPath=$(dpkg -L trimmomatic | grep .jar$ | head -1)

# Get the main directory that data is being stored in.
dataDirectory=${1%/*}

echo

# Determine if the input is sra or gzipped fastq.
if [[ $1 == *\.sr ]]
then
    # If given an sra file, generate the rawFastq file from it.
    echo "sra format given."
    dataName=${1%.sr}
    rawFastq="$dataName.sr.fastq.gz"
    echo "Converting to fastq format..."
    fastq-dump --gzip -O $dataDirectory $1

elif [[ $1 == *\.fastq\.gz ]]
then
    # If given a gzipped fastq file, set the dataName and rawFastq variables accordingly.
    echo "gzipped fastq given."
    dataName=${1%.fastq.gz}
    rawFastq=$1

else
    echo "Error: given file: $1 is not an sra file or a gzipped fastq file."
    exit 1

fi

echo "Working with $1"

# Create the names of all other intermediate and output files.
trimmedFastq="${dataName}_trimmed.fastq.gz"
bowtieSAMOutput="$dataName.sam"
BAMOutput="$dataName.bam.gz"
finalBedOutput="$dataName.bed"

# Trim the data (Single End)
echo "Trimming adaptors..."
java -jar $trimmomaticPath SE $rawFastq $trimmedFastq "ILLUMINACLIP:$2:2:30:10"

# Align the reads to the genome.
echo "Aligning reads with bowtie2..."
bowtie2 -x $3 -U $trimmedFastq -S $bowtieSAMOutput

# Convert from sam to bam.
echo "Converting from sam to bam..."
samtools view -b -o $BAMOutput $bowtieSAMOutput

# Gzip the sam file.  (Can't find a way to have bowtie do this to the output by default...)
echo "Gzipping sam file..."
gzip $bowtieSAMOutput

# Convert to final bed output.
echo "Converting to bed..."
bedtools bamtobed -i $BAMOutput > $finalBedOutput
