#!/usr/bin/bash

#
# Wonderdump is a workaround to download SRA files directly
# when fastq-dump's internet connection does not work.
# Which can happen surprisingly frequently.
#
# Usage:
#   wonderdump -X 10000 SRR1553500

set -ue

# This is where we will store the file.
SRA_DIR=~/ncbi/public/sra

# Make the directory if it does not exist.
mkdir -p $SRA_DIR

# All other parameters up to last.
PARAMS=${@:1:$(($# - 1))}

# The last parameter must be the SRR number.
SRR=${@:$#}

if [ -z "$SRR" ]
then
      echo "*** Please specify an SRR number"
      exit;
fi

echo "*** Getting SRR run: $SRR"

# Create the full path to the file.
SRA_FILE="$SRA_DIR/$SRR.sra"
TMP_FILE="$SRA_DIR/$SRR.tmp"

# Download only if it does not exist.
if [ ! -f $SRA_FILE ];
then
    PATH1=${SRR:0:6}
    PATH2=${SRR:0:10}
    URL="ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${PATH1}/${PATH2}/${SRR}.sra"
    echo "*** Downloading: $URL"
    echo "*** Saving to: $SRA_FILE"
    curl $URL > $TMP_FILE

    # Move to local file only if successful.
    mv $TMP_FILE $SRA_FILE
else
    echo "*** SRA file found: $SRA_FILE"
fi

# Run the fastq-dump.
fastq-dump $PARAMS ${SRA_FILE}