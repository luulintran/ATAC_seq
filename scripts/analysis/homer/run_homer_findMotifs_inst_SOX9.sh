#!/bin/bash  
#BSUB -J homer
#BSUB -o logs/homer.%J.out  
#BSUB -e logs/homer.%J.err  
#BSUB -n 12
#BSUB -R rusage[mem=50]  

# Create directories
mkdir -p logs
mkdir -p findMotifsGenome_motif_inst

# Load modules
. /usr/share/Modules/init/bash
module load modules modules-init
module load homer

# Define motif name
MOTIF_NAME="SOX9" # Change this as needed. This adds the motif name to the output file name.

# Define paths to BED files, motif file, and genome FASTA file
BED1="diffbind_bed_files/CTRL_e13_enriched.bed"
BED2="diffbind_bed_files/NICD_e13_enriched.bed"
MOTIF="SOX9motif.motif" # Motif file corresponding to $MOTIF_NAME. Change as needed. SOX9 motif was downloaded from HOMER
GENOME="../genome_file/mm10.fa"

# Define output directory
OUTPUT_DIR="findMotifsGenome_motif_inst"

# Make sure genome file and motif file exist
if [[ ! -f "$GENOME" ]]; then
  echo "Error: Genome file not found at $GENOME"
  exit 1
fi

if [[ ! -f "$MOTIF" ]]; then
  echo "Error: Motif file not found at $MOTIF"
  exit 1
fi

# Run HOMER findMotifsGenome.pl for both BED files using the motif file (-find)
echo "Running HOMER findMotifsGenome.pl with motif file for: $MOTIF_NAME"

findMotifsGenome.pl "$BED1" "$GENOME" "$OUTPUT_DIR" -find "$MOTIF" > "$OUTPUT_DIR/${MOTIF_NAME}_CTRL_findMotifs_inst.txt"

findMotifsGenome.pl "$BED2" "$GENOME" "$OUTPUT_DIR" -find "$MOTIF" > "$OUTPUT_DIR/${MOTIF_NAME}_NICD_findMotifs_inst.txt"

echo "findMotifsGenome complete. Results saved to $OUTPUT_DIR"