#!/bin/bash
# This script takes a gzipped fastq file and a file of adapater sequences to trim and uses them to run trimmomatic
# and trim the adaptor sequences off of the reads.

# Get trimmomatic's jar path.
trimmomaticPath=$(dpkg -L trimmomatic | grep .jar$ | head -1)

# Get the main directory that data is being stored in.
dataDirectory=${1%/*}

echo
fastqFile=$1; shift
adaptorFile=$1; shift

# Ensure the input is a gzipped fastq file.
if [[ $fastqFile == *\.fastq\.gz ]]
then
    # If given a gzipped fastq file, set the dataName and rawFastq variables accordingly.
    dataName=${fastqFile%.fastq.gz}

else
    echo "Error: given file: $fastqFile is not a gzipped fastq file."
    exit 1

fi

echo "Working with $fastqFile"

# Create the name of the trimmed output file.
trimmedFastq="${dataName}_trimmed.fastq.gz"

# Trim the data (Single End)
echo "Trimming adaptors..."
java -jar $trimmomaticPath SE $fastqFile $trimmedFastq "ILLUMINACLIP:$adaptorFile:2:30:10"
