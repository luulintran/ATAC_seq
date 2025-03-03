#!/bin/bash
#BSUB -J homer_annP
#BSUB -o logs/homer.%J.out
#BSUB -e logs/homer.%J.err
#BSUB -n 12
#BSUB -R rusage[mem=50]

# Define motif name
MOTIF_NAME="SOX9" # Change this as needed. This adds the motif name to the output file name

# Create directories
mkdir -p logs
mkdir -p annotatePeaks_motif_inst

# Load modules to use homer
. /usr/share/Modules/init/bash
module load modules modules-init
module load homer

# Define paths to BED files, motif file, and output directory
BED1="diffbind_bed_files/CTRL_e13_enriched.bed"
BED2="diffbind_bed_files/NICD_e13_enriched.bed"
MOTIF="NICD_motifs/homerResults/motif8.motif"
output_dir="annotatePeaks_motif_inst"

# Run HOMER annotatePeaks.pl for both BED files using the motif file -m
echo "Running HOMER annotatePeaks.pl with motif file for: $MOTIF_NAME"
annotatePeaks.pl "$BED1" mm10 -m "$MOTIF" > "$output_dir/${MOTIF_NAME}_CTRL_homerMotifs_annP.txt"

annotatePeaks.pl "$BED2" mm10 -m "$MOTIF" > "$output_dir/${MOTIF_NAME}_NICD_homerMotifs_annP.txt"

echo "Annotation complete. Results saved to $output_dir"