#!/bin/bash

## add flags for numThreads and output

sam=$1
threads=$2

echo CONVERTING :: $sam

samtools view -Sb $sam | samtools sort -@ $threads > $sam".bam"
