#!/bin/bash -euo pipefail
mkdir BismarkIndex
cp genome.fa BismarkIndex/
bismark_genome_preparation --bowtie2  BismarkIndex
