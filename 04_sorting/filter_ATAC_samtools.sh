#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --partition=mack
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --mail-user=e959b751@ku.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Script to perform filtering of blacklist regions
# paired end reads
# usage: sbatch filter_ATAC_samtools.sh
#Submitted batch job 33262643

# Required modules for servers
module load bedtools
module load samtools

# set sorting method flags from input
methodTag="coordinate"

# retrieve file with black list regions
blacklistPath=$(grep "blacklistRegions:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/blacklistRegions://g")

# retrieve analysis outputs absolute path
outputsPath=$(grep "outputs:" ../"inputData/inputPaths.txt" | tr -d " " | sed "s/outputs://g")

# setup the inputs path
inputsPath=$outputsPath"/filterd_sorted_"$methodTag"_ATAC"

# create outputs directory
outputFolder=$inputsPath

# move to outputs directory
cd "$outputFolder"

# name output file of inputs
inputOutFile="$outputFolder"/"$outputFolder"_summary.txt
# add software version to output summary file
bedtools --version > $inputOutFile
samtools --version >> $inputOutFile

# loop through all reads and sort sam/bam files for input to samtools
for f1 in $inputsPath"/"*/; do
	# trim file path from current folder name
	curSampleNoPath=$(basename "$f1")
	# create directory for current sample outputs
	mkdir "$outputFolder"/"$curSampleNoPath"
	# status message
	echo "Sample $curSampleNoPath is being processed..."
	# remove blacklist regions
	bedtools intersect -abam "$outputFolder"/"$curSampleNoPath"/noDups.bam -b $blacklistPath -v > "$outputFolder"/"$curSampleNoPath"/noBlacklist.bam
	samtools sort -@ 8 -o "$outputFolder"/"$curSampleNoPath"/sorted.bam -T /tmp/"$curSampleNoPath".sorted.bam "$outputFolder"/"$curSampleNoPath"/noBlacklist.bam
	samtools index -@ 8 "$outputFolder"/"$curSampleNoPath"/sorted.bam
	# status message
	echo "Sample $curSampleNoPath has been processed!"
done
