
### GENCODEGenes

This package loads genes from GENCODE GTF files, groups transcripts by gene, 
and provides methods for transcripts, so you can find exon coordinates, CDS 
distances and sequences.

### Install
```sh
pip install gencodegenes
```

### Usage

```py
from gencodegenes import Gencode

gencode = Gencode(GTF_PATH)
# full function arguments are Gencode(gencode, fasta=None, coding_only=True)
#  - gencode: path to GTF file (plain or gzipped)
#  - fasta: pass in path to fasta file to get gene transcripts with sequence
#  - coding_only: pass in False to include all transcripts, not just protein coding

# use as a context manager to close the fasta when done
with Gencode(GTF_PATH, fasta=FASTA_PATH) as gencode:
    ...

# get gene by HGNC symbol
gene = gencode['OR5A1']
transcripts = gene.transcripts
canonical = gene.canonical  # picks the transcript tagged Ensembl_canonical (the
                            # MANE Select transcript in human, where one exists), if
                            # none tagged, picks from those tagged appris_principal,
                            # if none tagged, picks from all transcripts. Within
                            # these, picks the longest CDS, or the longest cDNA if
                            # none are protein coding
gene.start, gene.end, gene.chrom, gene.strand, gene.symbol # other attributes available
gene.alternate_ids  # gene_id and hgnc_id from the GTF


# find gene nearest a genomic position, or overlapping a genomic region. 
# Chromosomes match with or without the 'chr' prefix
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
tx.get_codon_number_for_cds_position(cds_pos) # convert CDS pos to codon number
tx.translate(seq)                        # translate DNA to AA
tx.consequence(pos, ref, alt)            # get variant consequence (if opened with fasta)

# the transcript also has associated data fields
tx.name         # transcript ID
tx.chrom        # transcript chromosome
tx.start        # transcript start (lowest position on chromosome)
tx.end          # transcript end (highest position on chromosome)
tx.cds_start    # CDS start position
tx.cds_end      # CDS end position 
tx.type         # transcript type e.g. protein_coding
tx.strand       # strand (+ or -)
tx.exons        # list of exon coordinates
tx.cds          # list of CDS coordinates
tx.cds_sequence # get cDNA sequence (if Gencode was opened with fasta)
tx.attributes   # dict-like view of the GTF attributes e.g. tx.attributes['gene_id']

```