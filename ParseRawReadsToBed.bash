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
inputData=$1; shift
adaptorFile=$1; shift
bt2IndexBasename=$1; shift

if [[ $# > 0 ]]
then
    bowtieBinary=$1
else
    bowtieBinary="bowtie2"
fi    

# Determine if the input is sra or gzipped fastq.
if [[ $inputData == *\.sr ]]
then
    # If given an sra file, generate the rawFastq file from it.
    echo "sra format given."
    dataName=${inputData%.sr}
    rawFastq="$dataName.sr.fastq.gz"
    echo "Converting to fastq format..."
    fastq-dump --gzip -O $dataDirectory $inputData

elif [[ $inputData == *\.fastq\.gz ]]
then
    # If given a gzipped fastq file, set the dataName and rawFastq variables accordingly.
    echo "gzipped fastq given."
    dataName=${inputData%.fastq.gz}
    rawFastq=$inputData

else
    echo "Error: given file: $inputData is not an sra file or a gzipped fastq file."
    exit 1

fi

echo "Working with $inputData"

# Create the names of all other intermediate and output files.
trimmedFastq="${dataName}_trimmed.fastq.gz"
bowtieSAMOutput="$dataName.sam"
BAMOutput="$dataName.bam.gz"
finalBedOutput="$(dirname $(dirname $dataName))/$(basename $dataName).bed"

# Trim the data (Single End)
if [[ $adaptorFile != "NONE" ]]
then
    echo "Trimming adaptors..."
    java -jar $trimmomaticPath SE $rawFastq $trimmedFastq "ILLUMINACLIP:$adaptorFile:2:30:10"
else
    echo "Skipping adaptor trimming."
    trimmedFastq=$rawFastq
fi

# Align the reads to the genome.
echo "Aligning reads with bowtie2..."
$bowtieBinary -x $bt2IndexBasename -U $trimmedFastq -S $bowtieSAMOutput

# Convert from sam to bam.
echo "Converting from sam to bam..."
samtools view -b -o $BAMOutput $bowtieSAMOutput

# Gzip the sam file.  (Can't find a way to have bowtie do this to the output by default...)
echo "Gzipping sam file..."
gzip -f $bowtieSAMOutput

# Convert to final bed output.
echo "Converting to bed..."
bedtools bamtobed -i $BAMOutput > $finalBedOutput
