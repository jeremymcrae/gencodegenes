
### GENCODEGenes

This code loads genes from GENCODE GTF files into, groups transcripts by gene, 
and provides methods for transcripts, so you can find CDS distances and sequences.

### Install
```sh
pip install gencodegenes
```

### Usage

```py
from gencodegenes import Gencode

gencode = Gencode('PATH_TO_GTF')

# get gene by HGNC symbol
gene = gencode['OR5A1]
transcripts = gene.transcripts
canonical = gene.canonical

```