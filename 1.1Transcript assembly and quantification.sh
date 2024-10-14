## Transcript assembly and quantification from single cell long read RNA-seq data.

#!/bin/bash
# FLAMES version=0.1

sample=$1

## match cell barcode

match_cell_barcode='FLAMES/src/bin/match_cell_barcode'
${match_cell_barcode} Data/fastq/${sample}.fq Result/match_cell_barcode/${sample}/${sample}.stat.txt Result/match_cell_barcode/${sample}/${sample}.matched.fastq.gz Data/nanopore/barcode/${sample}/${sample}_barcodes.tsv 2


## assemble and quantify based on FLAMES sc long pipeline
### merge all fastq

mer_fq='Data/merged.fastq.gz'
find Data/fastq -name *.fastq.gz | xargs cat > ${mer_fq}
file_gff='reference/gencode.v38.primary_assembly.annotation.gff3'
file_fa='reference/GRCh38.primary_assembly.genome.fa'
dir_out='Result/FLAMES'
python FLAMES/python/sc_long_pipeline.py -a ${file_gff}  -i ${mer_fq} -o ${dir_out} -f ${file_fa} 

## transcript classification
assembled_gtf='Result/isoform_annotated.filtered.gff3'
python ~/software/SQANTI3-4.3/sqanti3_qc.py ${assembled_gtf} ${file_gff} ${file_fa} -o ${sample} -d Result/SQANTI_output --cage_peak Reference/human.refTSS_v3.1.hg38.bed
