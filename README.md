# ManeSelectBed

This repo contains a script that periodically generates BED files from MANE Select transcripts. The script produces several types of BED files, including exon, CDS, intron, and UTR BED files. Additionally, it generates a corrected BED file based on the exon BED file, with specific adjustments.



## Features

**Exon BED**: Contains the coordinates of exons from MANE Select transcripts.

**CDS BED**: Contains the coordinates of coding sequences (CDS) from MANE Select transcripts.

**Intron BED**: Contains the coordinates of introns from MANE Select transcripts.

**UTR BED**: Contains the coordinates of untranslated regions (UTR) from MANE Select transcripts.



### Corrected BED

A refined version of the exon BED file with the following adjustments:

**UTR Correction**: Exons that include UTR regions are corrected using the CDS BED file.

**Polymorphic Genes**: For genes with multiple transcripts, only the transcript with the smallest RefSeq ID is retained.

**Intersection Handling**: If genes have overlapping regions, the coordinates are adjusted to prioritize the region that appears first in the genome sequence. The corrected BED file ensures that there are no overlapping entries.


## Link

[Gencode](https://www.gencodegenes.org/)

[NCBI MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/)


