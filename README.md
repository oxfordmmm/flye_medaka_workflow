# flye_medaka_workflow

Assembles nanopore data with Flye, polish with medaka. Then runs MLST and BAKTA and puts all GFF files in ROARY. 

```bash
nextflow run  https://github.com/oxfordmmm/flye_medaka_workflow \
        -r v0.1 \
        --inputFastq data/ \
        -resume \
        --inputmodel r1041_e82_400bps_hac_g615_model.tar.gz \
        -with-trace -profile conda,bmrc

```

## Outputs

CSVs - seqkit contig summaries
MLST - output from MLST
assemblies - contigs from flye or spades
bakta - output from bakta annotation
roary - output from ROARY
