
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

gencode = Gencode(GTF_PATH)  # or Gencode(GTF, FASTA_PATH) to give transcripts DNA sequence

# get gene by HGNC symbol
gene = gencode['OR5A1']
transcripts = gene.transcripts
canonical = gene.canonical  # picks MANE transcript if available, if none named
                            # as MANE, picks the one tagged as appris_principal
                            # (or longest CDS if multiple), if none tagged, picks
                            # the longest protein coding, if none protein coding,
                            # picks the longest cDNA 
gene.start, gene.end, gene.chrom, gene.strand, gene.symbol # other attributes available


# find gene nearest a genomic position, or overlapping a genomic region
gencode.nearest('chr1', 1000000)
gencode.in_region('chr1', 1000000, 2000000)

# and the transcript has a bunch of methods
tx = gene.canonical
tx.in_exons(pos)                         # check if pos in exons
tx.in_coding_region(pos)                 # check if pos in CDS
tx.get_coding_distance(pos)              # get distance in CDS to CDS start
tx.get_closest_exon(pos)                 # find exon closest to position
tx.get_position_on_chrom(cds_pos)        # convert CDS pos to genomic pos
tx.get_codon_info(pos)                   # get info about codon for a site
tx.get_codon_number_for_cds_pos(cds_pos) # convert CDS pos to codon number
tx.translate(seq)                        # translate DNA to AA (if opened with Fasta)

# the transcript also has associated data fields
tx.name         # transcript ID
tx.chrom        # transcript chromosome
tx.start        # transcript start (TSS)
tx.end          # transcript end
tx.cds_start    # CDS start position
tx.cds_end      # CDS end position 
tx.type         # transcript type e.g. protein_coding
tx.strand       # strand (+ or -)
tx.exons        # list of exon coordinates
tx.cds          # list of CDS coordinates
tx.cds_sequence # get cDNA sequence (if Gencode was opened with fasta)

```