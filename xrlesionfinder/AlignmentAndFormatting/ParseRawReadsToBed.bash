#!/bin/bash
# This script takes an sra file as input and runs the file through the SRA toolkit, trimmomatic, 
# bowtie, samtools, and bedtools to create a bed file.
# The first argument should be the input sra or gzipped fastq data, the second input should be the path to the fasta file of sequences for trimmomatic,
# and the third argument should be the path to the basename for the bowtie2 index files.  (e.g. for files basename.1.bt2 in dir, give: dir/basename)

# Get trimmomatic's jar path.
trimmomaticPath=$(dpkg -L trimmomatic | grep .jar$ | head -1)

inputReads2=""

while [[ $# > 0 ]]; do
  case $1 in
    -1|--read-file-1)
      inputReads1="$2"
      shift
      shift
      ;;
    -2|--read-file-2)
      inputReads2="$2"
      shift
      shift
      ;;
    -a|--adapterFile)
      adapterFile="$2"
      shift
      shift
      ;;
    -i|--bt2-index)
      bt2IndexBasename="$2"
      shift
      shift
      ;;
    -t|--threads)
      threads="$2"
      shift
      shift
      ;;
    -c|--custom-bowtie2-arguments)
      customBowtie2Arguments="$2"
      shift
      shift
      ;;
    -b|--bowtie2-binary)
      bowtie2Binary="$2"
      shift
      shift
      ;;
    -*|--*)
      echo "Unknown option: $1"
      exit 1
      ;;
    *)
      echo "Unexpected positional argument: $1"
      shift
      ;;
  esac
done

# Get the main directory and .tmp directories.
dataDirectory=${inputReads1%/*}
tmpDirectory="$dataDirectory/.tmp"

# If no custom bowtie2 binary was given, use the default.
if [[ -z $bowtie2Binary ]]
then
    bowtie2Binary="bowtie2"
fi

# Determine whether the file is gzipped or not and set the dataName variable accordingly.
dataName=${inputReads1##*/}
if [[ $inputReads1 == *\.fastq ]]
then
    echo "fastq given."
    dataName=${dataName%.fastq}
elif [[ $inputReads1 == *\.fastq\.gz ]]
then
    echo "gzipped fastq given."
    dataName=${dataName%.fastq.gz}
else
    echo "Error: given file: $inputReads1 is not a fastq file or a gzipped fastq file."
    exit 1
fi

# Don't forget to trim the read identifier off of the data name if the data is paired-end.
if [[ ! -z "$inputReads2" ]]
then
    dataName=${dataName%_R1}
fi

if [[ -z "$inputReads2" ]]
then
    echo "Working with $inputData"
else
    echo "Working with $inputReads1 and $inputReads2"
fi

# Create the names of all other intermediate and output files.
trimmedFastq="$tmpDirectory/${dataName}_trimmed.fastq.gz"
bowtieSAMOutput="$tmpDirectory/$dataName.sam"
bowtieStatsOutput="$tmpDirectory/${dataName}_bowtie2_stats.txt"
BAMOutput="$tmpDirectory/$dataName.bam.gz"
finalBedOutput="$dataDirectory/$dataName.bed"

if [[ ! -z "$inputReads2" ]]
then
    trimmedFastqP1="$tmpDirectory/${dataName}_trimmed_1P.fastq.gz"
    trimmedFastqP2="$tmpDirectory/${dataName}_trimmed_2P.fastq.gz"
fi

# Trim the data
if [[ $adapterFile != "NONE" ]]
then
    if [[ -z "$inputReads2" ]]
    then
        echo "Trimming adapters in single-end mode..."
        java -jar $trimmomaticPath SE -threads $threads $inputReads1 $trimmedFastq "ILLUMINACLIP:$adapterFile:2:30:10"
    else
        echo "Trimming adapters in paired-end mode..."
        java -jar $trimmomaticPath PE -threads $threads $inputReads1 $inputReads2 \
        -baseout $trimmedFastq "ILLUMINACLIP:$adapterFile:2:30:10"
    fi

else
    echo "Skipping adapter trimming."
    if [[ -z "$inputReads2" ]]
    then
        trimmedFastq=$inputReads1
    else
        trimmedFastqP1=$inputReads1
        trimmedFastqP2=$inputReads2
    fi
fi

# Align the reads to the genome.
echo "Aligning reads with bowtie2..."
if [[ -z "$inputReads2" ]]
then
    if [[ -z "$customBowtie2Arguments" ]]
    then
        $bowtie2Binary -x $bt2IndexBasename -U $trimmedFastq -S $bowtieSAMOutput -p $threads \
        |& tail -6 | tee $bowtieStatsOutput
    else
        $bowtie2Binary -x $bt2IndexBasename -1 $trimmedFastqP1 -2 $trimmedFastqP2 -S $bowtieSAMOutput -p $threads \
        $customBowtie2Arguments |& tail -6 | tee $bowtieStatsOutput
    fi
else
    if [[ -z "$customBowtie2Arguments" ]]
    then
        $bowtie2Binary -x $bt2IndexBasename -1 $trimmedFastqP1 -2 $trimmedFastqP2 -S $bowtieSAMOutput \
        -p $threads |& tail -20 | tee $bowtieStatsOutput
    else
        $bowtie2Binary -x $bt2IndexBasename -1 $trimmedFastqP1 -2 $trimmedFastqP2 -S $bowtieSAMOutput \
        -p $threads $customBowtie2Arguments |& tail -20 | tee $bowtieStatsOutput
    fi
fi

# Convert from sam to bam.
echo "Converting from sam to bam..."
samtools view -b --threads $((threads-1)) -o $BAMOutput $bowtieSAMOutput

# Gzip the sam file.  (Can't find a way to have bowtie do this to the output by default...)
echo "Gzipping sam file..."
pigz -f -p $threads $bowtieSAMOutput

# Convert to final bed output.
echo "Converting to bed..."
bedtools bamtobed -i $BAMOutput > $finalBedOutput
