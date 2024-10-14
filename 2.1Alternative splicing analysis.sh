#!/bin/bash
# suppa version=2.3
## step 1
assembled_gtf='Result/isoform_annotated.filtered.gtf'
suppa.py generateEvents -i ${assembled_gtf} -o Result/isoform_annotated.events -e SE SS MX RI FL -f ioe
cd Result
awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' *.ioe > isoform_annotated.all.events.ioe
## step 2
ioe_merge_file='Result/isoform_annotated.all.events.ioe'
transcript_exp='Result/transcript_exp.txt'
suppa.py psiPerEvent -i ${ioe_merge_file} -e ${transcript_exp} -o Result/total_project_events
