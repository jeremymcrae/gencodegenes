
import gzip
from pathlib import Path
import unittest
import tempfile

from gencodegenes.gencode import Gencode, Gene, _parse_gtfline, _open_gencode
from gencodegenes.transcript import Transcript

def write_gtf(path, lines):
    with open(path, 'wt') as output:
        for line in lines:
            output.write(line)

def make_fasta(path, chroms):
    with open(path, 'wt') as output:
        for chrom in chroms:
            output.write(f'>{chrom}\n')
            lines = ['A' * 50 + '\n'] * 50
            output.writelines(lines)

class TestGencode(unittest.TestCase):
    
    def setUp(self):
        ''' set path to folder with test data
        '''
        self.folder = Path(__file__).parent /  "data"
        self.gtf_path = self.folder / 'example.grch38.gtf'
        self.fasta_path = self.folder / 'example.grch38.fa'
        temp_gtf = tempfile.NamedTemporaryFile(delete=False)
        temp_fasta = tempfile.NamedTemporaryFile(delete=False)
        self.temp_gtf_path = temp_gtf.name
        self.temp_fasta_path = temp_fasta.name
        temp_gtf.close()
        temp_fasta.close()
        self.maxDiff = None
    
    def tearDown(self):
        try:
            Path(self.temp_gtf_path).unlink()
            Path(self.temp_fasta_path).unlink()
            Path(self.temp_fasta_path + '.fai').unlink(missing_ok=True)
        except:
            pass

    def test_gencode_opens(self):
        ''' test we can open a gencode object
        '''
        gencode = Gencode(self.gtf_path, self.fasta_path)
        self.assertEqual(len(gencode), 1)
        gene = gencode['OR4F5']
        with self.assertRaises(KeyError):
            gencode['ZZZZZZZ']
    
    def test_gencode_genome(self):
        ''' test each Gencode uses its own genome fasta, and closes it on exit
        '''
        with Gencode(self.gtf_path, self.fasta_path) as gencode:
            seq = gencode['OR4F5'].canonical.genomic_sequence
        self.assertEqual(len(seq), 927)
        
        # exiting without a fasta is fine
        with Gencode(self.gtf_path) as gencode:
            pass
        gencode.__exit__()
        
        # genomes aren't shared between Gencode objects
        with open(self.temp_fasta_path, 'wt') as output:
            output.write('>chr1\n' + ('A' * 60 + '\n') * 1200)
        real = Gencode(self.gtf_path, self.fasta_path)
        poly_a = Gencode(self.gtf_path, self.temp_fasta_path)
        no_fasta = Gencode(self.gtf_path)
        self.assertEqual(real['OR4F5'].canonical.genomic_sequence, seq)
        self.assertEqual(poly_a['OR4F5'].canonical.genomic_sequence, 'A' * 927)
        self.assertEqual(no_fasta['OR4F5'].canonical.genomic_sequence, '')
    
    def test_gencode_genome_chrom_edges(self):
        ''' test transcripts near the ends of a chromosome get their sequence
        '''
        seq = 'ACGT' * 20
        with open(self.temp_fasta_path, 'wt') as output:
            output.write(f'>chr1\n{seq}\n')
        lines = []
        for tx_id, strand, start, end in [('START', '+', 2, 20), ('END', '-', 62, 80)]:
            # 15 bp CDS, so the CDS isn't padded from the flanking sequence
            for feature, x, y in [('transcript', start, end), ('exon', start, end),
                    ('CDS', start + 2, start + 16)]:
                lines.append(f'chr1\tHAVANA\t{feature}\t{x}\t{y}\t.\t{strand}\t.\t'
                    f'transcript_id "{tx_id}"; gene_name "{tx_id}"; transcript_type "protein_coding";\n')
        write_gtf(self.temp_gtf_path, lines)
        
        with Gencode(self.temp_gtf_path, self.temp_fasta_path) as gencode:
            start = gencode['START'].transcripts[0]
            end = gencode['END'].transcripts[0]
        
        self.assertEqual(start.genomic_offset, 1)
        self.assertEqual(start.genomic_sequence, seq[0:20])
        self.assertEqual(start.cds_sequence, seq[3:18])
        
        self.assertEqual(end.genomic_offset, 1)
        self.assertEqual(end.genomic_sequence, seq[60:80])
        self.assertEqual(end.cds_sequence, end.reverse_complement(seq[63:78]))
    
    def test_gene_alternate_ids(self):
        ''' test Gene accepts alternate IDs as str, list of str, or bytes
        '''
        self.assertEqual(Gene('TEST').alternate_ids, [])
        self.assertEqual(Gene('TEST', 'ENSG1').alternate_ids, ['ENSG1'])
        self.assertEqual(Gene('TEST', ['ENSG1', 'HGNC:1']).alternate_ids, ['ENSG1', 'HGNC:1'])
        self.assertEqual(Gene('TEST', [b'ENSG1']).alternate_ids, ['ENSG1'])
    
    def test_gencode_add_gene(self):
        ''' test adding genes doesn't duplicate genes in the region index
        '''
        gencode = Gencode(self.gtf_path)
        gene = Gene('NEW')
        gene.add_transcript(Transcript('ENST_NEW', 'chr1', 500000, 501000, '+',
            exons=[(500000, 501000)], cds=[(500000, 501000)]))
        gencode.add_gene(gene)
        gencode.add_gene(gene)
        
        self.assertEqual([x.symbol for x in gencode.in_region('chr1', 69000, 70100)], ['OR4F5'])
        self.assertEqual([x.symbol for x in gencode.in_region('chr1', 500500, 500600)], ['NEW'])
        
        # genes without transcripts lack coordinates, so can't be indexed
        with self.assertRaises(ValueError):
            gencode.add_gene(Gene('EMPTY'))
    
    def test_gencode_in_region(self):
        ''' test that in_region pulls out the correct genes
        '''
        lines = '##format: gtf\n' \
                'chr1\tHAVANA\tgene\t10\t20\t.\t-\t.\tgene_name "TEST1";\n' \
                'chr1\tHAVANA\ttranscript\t10\t20\t.\t-\t.\ttranscript_id "ENST_A";gene_name "TEST1"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\tCDS\t15\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\tgene\t100\t110\t.\t-\t.\tgene_name "TEST2";\n' \
                'chr1\tHAVANA\ttranscript\t100\t110\t.\t-\t.\ttranscript_id "ENST_B";gene_name "TEST2"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\texon\t100\t110\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\tCDS\t105\t110\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding;"\n'\
                'chr2\tHAVANA\tgene\t100\t110\t.\t-\t.\tgene_name "TEST3";\n' \
                'chr2\tHAVANA\ttranscript\t100\t110\t.\t-\t.\ttranscript_id "ENST_C";gene_name "TEST3"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr2\tHAVANA\texon\t100\t110\t.\t-\t.\ttranscript_id "ENST_C" gene_name "TEST3"; transcript_type "protein_coding"\n' \
                'chr2\tHAVANA\tCDS\t105\t110\t.\t-\t.\ttranscript_id "ENST_C" gene_name "TEST3"; transcript_type "protein_coding"\n'
        
        write_gtf(self.temp_gtf_path, lines)
        make_fasta(self.temp_fasta_path, ['chr1', 'chr2'])
        gencode = Gencode(self.temp_gtf_path, self.temp_fasta_path)
        
        genes = gencode.in_region('chr1', 5, 6)  # shouldn't get anything
        self.assertEqual(genes, [])
        genes = gencode.in_region('chr1', 5, 15)  # should get TEST1 (via start)
        self.assertEqual(set(x.symbol for x in genes), set(["TEST1"]))
        genes = gencode.in_region('chr1', 5, 25)  # should get TEST1 (via start, end)
        self.assertEqual(set(x.symbol for x in genes), set(["TEST1"]))
        genes = gencode.in_region('chr1', 15, 25)  # should get TEST1 (via end)
        self.assertEqual(set(x.symbol for x in genes), set(["TEST1"]))
        genes = gencode.in_region('chr1', 12, 13) # should get TEST1 (via encapsulate)
        self.assertEqual(set(x.symbol for x in genes), set(["TEST1"]))
        genes = gencode.in_region('chr1', 5, 105) # should get TEST1 (via start, end), TEST2 (via start)
        self.assertEqual(set(x.symbol for x in genes), set(["TEST1", "TEST2"]))
        genes = gencode.in_region('chr1', 5, 115) # should get TEST1 (via start, end), TEST2 (via start, end)
        self.assertEqual(set(x.symbol for x in genes), set(["TEST1", "TEST2"]))
        genes = gencode.in_region('chr1', 15, 105) # should get TEST1 (via end), TEST2 (via start)
        self.assertEqual(set(x.symbol for x in genes), set(["TEST1", "TEST2"]))
        genes = gencode.in_region('chr1', 25, 105)  # should get TEST2 (via start)
        self.assertEqual(set(x.symbol for x in genes), set(["TEST2"]))
        genes = gencode.in_region('chr1', 102, 105) # should get TEST2 (via encapsulate)
        self.assertEqual(set(x.symbol for x in genes), set(["TEST2"]))
        genes = gencode.in_region('chr1', 120, 150) # shouldn't get anything
        self.assertEqual(genes, [])
        
        # no overlapping genes right up to the border of a gene
        in_region = gencode.in_region('chr1', 0, 9)
        self.assertEqual(len(in_region), 0)
        in_region = gencode.in_region('chr1', 9, 9) # check a single base before
        self.assertEqual(len(in_region), 0)
        in_region = gencode.in_region('chr1', 21, 30)
        self.assertEqual(len(in_region), 0)
        
        # we find matches if even a single base of the gene is in the region
        in_region = gencode.in_region('chr1', 0, 11)
        self.assertEqual(len(in_region), 1)
        del gencode
    
    def test_gencode_nearest(self):
        ''' test that we can find the nearest gene from Gencode
        '''
        lines = ['##format: gtf\n',
                'chr1\tHAVANA\tgene\t10\t20\t.\t-\t.\tgene_name "TEST1";\n',
                'chr1\tHAVANA\ttranscript\t10\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding"; tag "appris_principal_1";\n',
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding;"\n',
                'chr1\tHAVANA\tCDS\t15\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding;"\n',
                'chr1\tHAVANA\tgene\t100\t110\t.\t-\t.\tgene_name "TEST2";\n',
                'chr1\tHAVANA\ttranscript\t100\t110\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding"; tag "appris_principal_1";\n',
                'chr1\tHAVANA\texon\t100\t110\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding;"\n',
                'chr1\tHAVANA\tCDS\t105\t110\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding;"\n',
                'chr2\tHAVANA\tgene\t100\t110\t.\t-\t.\tgene_name "TEST3";\n',
                'chr2\tHAVANA\ttranscript\t100\t110\t.\t-\t.\ttranscript_id "ENST_C"; gene_name "TEST3"; transcript_type "protein_coding"; tag "appris_principal_1";\n',
                'chr2\tHAVANA\texon\t100\t110\t.\t-\t.\ttranscript_id "ENST_C"; gene_name "TEST3"; transcript_type "protein_coding;"\n',
                'chr2\tHAVANA\tCDS\t105\t110\t.\t-\t.\ttranscript_id "ENST_C"; gene_name "TEST3"; transcript_type "protein_coding;"\n']
        
        write_gtf(self.temp_gtf_path, lines)
        make_fasta(self.temp_fasta_path, ['chr1', 'chr2'])
        gencode = Gencode(self.temp_gtf_path, self.temp_fasta_path)
        
        self.assertEqual(gencode.nearest('chr1', 5).symbol, 'TEST1')
        self.assertEqual(gencode.nearest('chr1', 10).symbol, 'TEST1')
        self.assertEqual(gencode.nearest('chr1', 11).symbol, 'TEST1')
        self.assertEqual(gencode.nearest('chr1', 21).symbol, 'TEST1')
        self.assertEqual(gencode.nearest('chr1', 80).symbol, 'TEST2')
        self.assertEqual(gencode.nearest('chr1', 105).symbol, 'TEST2')
        self.assertEqual(gencode.nearest('chr1', 2000).symbol, 'TEST2')
    
        # and look for genes on the final chrom, both before and after the final gene
        self.assertEqual(gencode.nearest('chr2', 0).symbol, 'TEST3')
        self.assertEqual(gencode.nearest('chr2', 2000).symbol, 'TEST3')
        
        with self.assertRaises(ValueError):
            gencode.nearest('chrZZZ', 2000)
        
        del gencode
    
    def test_gencode_chrom_prefix(self):
        ''' test region lookups work whether or not chromosomes have a chr prefix
        '''
        lines = ['##format: gtf\n',
                '1\tHAVANA\ttranscript\t10\t20\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding";\n',
                '1\tHAVANA\texon\t10\t20\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding";\n',
                '1\tHAVANA\tCDS\t10\t20\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding";\n',
                'chr2\tHAVANA\ttranscript\t10\t20\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding";\n',
                'chr2\tHAVANA\texon\t10\t20\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding";\n',
                'chr2\tHAVANA\tCDS\t10\t20\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding";\n']
        
        write_gtf(self.temp_gtf_path, lines)
        gencode = Gencode(self.temp_gtf_path)
        
        for chrom in ['1', 'chr1']:
            self.assertEqual([x.symbol for x in gencode.in_region(chrom, 5, 15)], ['TEST1'])
            self.assertEqual(gencode.nearest(chrom, 100).symbol, 'TEST1')
        for chrom in ['2', 'chr2']:
            self.assertEqual([x.symbol for x in gencode.in_region(chrom, 5, 15)], ['TEST2'])
            self.assertEqual(gencode.nearest(chrom, 100).symbol, 'TEST2')
        
        with self.assertRaises(ValueError):
            gencode.in_region('3', 5, 15)
        with self.assertRaises(ValueError):
            gencode.nearest('chr3', 100)
    
    def test_gencode_symbol_at_multiple_loci(self):
        ''' test genes sharing a symbol at different loci are kept apart
        '''
        attrs = 'gene_name "SHOX"; transcript_type "protein_coding";'
        lines = []
        for chrom, gene_id, start, end in [('chrX', 'ENSG1', 100, 200),
                ('chrY', 'ENSG1_PAR_Y', 5000, 6000), ('chrY', 'ENSG2', 8000, 9000)]:
            for feature in ['transcript', 'exon', 'CDS']:
                lines.append(f'{chrom}\tHAVANA\t{feature}\t{start}\t{end}\t.\t+\t.\t'
                    f'gene_id "{gene_id}"; transcript_id "{gene_id}_T"; {attrs}\n')
        
        write_gtf(self.temp_gtf_path, lines)
        gencode = Gencode(self.temp_gtf_path)
        
        self.assertEqual(len(gencode), 1)
        self.assertEqual(list(gencode), ['SHOX'])
        gene = gencode['SHOX']
        self.assertEqual((gene.chrom, gene.start, gene.end), ('chrX', 100, 200))
        self.assertEqual([x.name for x in gene.transcripts], ['ENSG1_T'])
        
        self.assertEqual(gencode.in_region('chrX', 1, 10000), [gene])
        self.assertEqual(gencode.nearest('chrX', 150), gene)
        
        genes = gencode.in_region('chrY', 1, 10000)
        self.assertEqual(sorted((x.start, x.end) for x in genes), [(5000, 6000), (8000, 9000)])
        self.assertEqual(gencode.nearest('chrY', 8500).start, 8000)
        
        # genes with a symbol already present are only added on other chromosomes
        gencode.add_gene(gene)
        new = Gene('SHOX')
        new.add_transcript(Transcript('ENST3', 'chr1', 10, 20, '+'))
        gencode.add_gene(new)
        self.assertEqual(gencode.in_region('chrX', 1, 10000), [gene])
        self.assertEqual(gencode.in_region('chr1', 1, 100), [new])
    
    def test_gencode_canonical(self):
        ''' test we find the correct canonical transcript
        '''
        lines = ['##format: gtf\n',
                'chr1\tHAVANA\tgene\t20\t30\t.\t-\t.\tgene_name "TEST1";\n',
                'chr1\tHAVANA\ttranscript\t20\t30\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding"; tag "appris_principal_1";\n',
                'chr1\tHAVANA\texon\t20\t30\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding;"\n',
                'chr1\tHAVANA\tCDS\t25\t30\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding;"\n',
                'chr1\tHAVANA\tgene\t100\t110\t.\t-\t.\tgene_name "TEST1";\n',
                'chr1\tHAVANA\ttranscript\t100\t110\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST1"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\texon\t100\t110\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST1"; transcript_type "protein_coding;"\n',
                'chr1\tHAVANA\tCDS\t110\t100\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST1"; transcript_type "protein_coding;"\n',
        ]
        
        write_gtf(self.temp_gtf_path, lines)
        make_fasta(self.temp_fasta_path, ['chr1', 'chr2'])
        gencode = Gencode(self.temp_gtf_path, self.temp_fasta_path)
        
        gene = gencode['TEST1']
        canonical = gene.canonical 
        self.assertEqual(canonical.name, 'ENST_A')
        del gencode
        
        # give the second transcript the appris_principal tag as well, which
        # given that it has the longer CDS, would make it the canonical now
        lines[-3] = lines[-3].strip() + ' tag "appris_principal_1";\n'
        write_gtf(self.temp_gtf_path, lines)
        make_fasta(self.temp_fasta_path, ['chr1', 'chr2'])
        gencode = Gencode(self.temp_gtf_path, self.temp_fasta_path)
        
        gene = gencode['TEST1']
        canonical = gene.canonical
        self.assertEqual(canonical.name, 'ENST_B')
        del gencode
        
        # Ensembl_canonical takes priority over the longer appris_principal CDS,
        # even when the Ensembl_canonical transcript is also appris_principal
        lines[2] = lines[2].strip() + ' tag "Ensembl_canonical";\n'
        write_gtf(self.temp_gtf_path, lines)
        gencode = Gencode(self.temp_gtf_path, self.temp_fasta_path)
        
        gene = gencode['TEST1']
        canonical = gene.canonical
        self.assertEqual(canonical.name, 'ENST_A')
        del gencode
    
    def test_gencode_canonical_noncoding(self):
        ''' test a non-coding transcript isn't picked over a coding transcript
        '''
        lines = ['##format: gtf\n',
                'chr1\tHAVANA\ttranscript\t20\t30\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\texon\t20\t30\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\tCDS\t25\t30\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\ttranscript\t100\t500\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "TEST1"; transcript_type "retained_intron";\n',
                'chr1\tHAVANA\texon\t100\t500\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "TEST1"; transcript_type "retained_intron";\n',
                'chr1\tHAVANA\ttranscript\t1000\t1100\t.\t+\t.\ttranscript_id "ENST_C"; gene_name "TEST2"; transcript_type "lncRNA";\n',
                'chr1\tHAVANA\texon\t1000\t1020\t.\t+\t.\ttranscript_id "ENST_C"; gene_name "TEST2"; transcript_type "lncRNA";\n',
                'chr1\tHAVANA\ttranscript\t1000\t1100\t.\t+\t.\ttranscript_id "ENST_D"; gene_name "TEST2"; transcript_type "lncRNA";\n',
                'chr1\tHAVANA\texon\t1000\t1100\t.\t+\t.\ttranscript_id "ENST_D"; gene_name "TEST2"; transcript_type "lncRNA";\n',
        ]
        
        write_gtf(self.temp_gtf_path, lines)
        gencode = Gencode(self.temp_gtf_path, coding_only=False)
        
        self.assertEqual(gencode['TEST1'].canonical.name, 'ENST_A')
        
        # without any coding transcripts, fall back to the longest exonic length
        self.assertEqual(gencode['TEST2'].canonical.name, 'ENST_D')
        del gencode
    
    def test_gencode_nearest_enveloping_gene(self):
        ''' test nearest finds a gene that envelops another gene
        '''
        lines = ['##format: gtf\n',
                'chr1\tHAVANA\ttranscript\t1000\t500000\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "BIG"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\texon\t1000\t500000\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "BIG"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\tCDS\t1000\t500000\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "BIG"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\ttranscript\t10000\t11000\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "SMALL"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\texon\t10000\t11000\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "SMALL"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\tCDS\t10000\t11000\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "SMALL"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\ttranscript\t600000\t601000\t.\t+\t.\ttranscript_id "ENST_C"; gene_name "LATE"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\texon\t600000\t601000\t.\t+\t.\ttranscript_id "ENST_C"; gene_name "LATE"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\tCDS\t600000\t601000\t.\t+\t.\ttranscript_id "ENST_C"; gene_name "LATE"; transcript_type "protein_coding";\n']
        
        write_gtf(self.temp_gtf_path, lines)
        gencode = Gencode(self.temp_gtf_path)
        
        # BIG ends 40 kb before the site, LATE starts 60 kb after
        self.assertEqual(gencode.nearest('chr1', 540000).symbol, 'BIG')
        self.assertEqual(gencode.nearest('chr1', 570000).symbol, 'LATE')
        self.assertEqual(gencode.nearest('chr1', 1).symbol, 'BIG')
        self.assertEqual(gencode.nearest('chr1', 700000).symbol, 'LATE')
        del gencode
    
    def test_gencode_nearest_adjacent_gene(self):
        ''' test nearest picks the gene containing a site over an adjacent gene
        '''
        lines = ['##format: gtf\n',
                'chr1\tHAVANA\ttranscript\t1000\t5000\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "LONG"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\texon\t1000\t5000\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "LONG"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\tCDS\t1000\t5000\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "LONG"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\ttranscript\t5001\t6000\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "SHORT"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\texon\t5001\t6000\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "SHORT"; transcript_type "protein_coding";\n',
                'chr1\tHAVANA\tCDS\t5500\t5600\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "SHORT"; transcript_type "protein_coding";\n']
        
        write_gtf(self.temp_gtf_path, lines)
        gencode = Gencode(self.temp_gtf_path)
        
        # the site is in the UTR of SHORT, next to the end of LONG, which has
        # a longer CDS
        self.assertEqual(gencode.nearest('chr1', 5001).symbol, 'SHORT')
        self.assertEqual(gencode.nearest('chr1', 5000).symbol, 'LONG')
    
    def test_parse_gtf_gene_line(self):
        ''' test we can parse a GTF line for a gene feature
        '''
        line = 'chr1\tHAVANA\tgene\t69091\t70008\t.\t+\t.\tgene_id "ENSG00000186092.4";' \
            'gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5";' \
            'level 2; havana_gene "OTTHUMG00000001094.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'gene', 
            'start': 69091,
            'end': 70008,
            'strand': b'+',
            'symbol': b'',
            'alternate_ids': [],
            'tx_id': b'',
            'transcript_type': b'',
            'is_canonical': 0,
            'attributes': {},
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_transcript_line(self):
        ''' test we can parse a GTF line for a transcript feature
        '''
        line = 'chr1\tHAVANA\ttranscript\t69091\t70008\t.\t+\t.\tgene_id "ENSG00000186092.4"; '\
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; ' \
            'gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; ' \
            'transcript_status "KNOWN"; transcript_name "OR4F5-001"; level 2; ' \
            'protein_id "ENSP00000334393.3"; tag "basic"; transcript_support_level "NA"; ' \
            'tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS30547.1"; ' \
            'havana_gene "OTTHUMG00000001094.2"; havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'transcript', 
            'start': 69091,
            'end': 70008,
            'strand': b'+',
            'symbol': b'OR4F5',
            'alternate_ids': [b'ENSG00000186092.4'],
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_canonical': 5,
            'attributes': {
                b'gene_id': b'ENSG00000186092.4',
                b'transcript_id': b'ENST00000335137.3',
                b'gene_type': b'protein_coding',
                b'gene_status': b'KNOWN',
                b'gene_name': b'OR4F5',
                b'transcript_type': b'protein_coding',
                b'transcript_status': b'KNOWN',
                b'transcript_name': b'OR4F5-001',
                b'level': b'2',
                b'protein_id': b'ENSP00000334393.3',
                b'tag': b'basic,appris_principal_1,CCDS',
                b'transcript_support_level': b'NA',
                b'ccdsid': b'CCDS30547.1',
                b'havana_gene': b'OTTHUMG00000001094.2',
                b'havana_transcript': b'OTTHUMT00000003223.2',
                },
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_transcript_line_alternate_ids(self):
        ''' test we can parse a GTF line which includes both alternate IDs
        '''
        line = 'chr1\tHAVANA\ttranscript\t69091\t70008\t.\t+\t.\tgene_id "ENSG00000186092.4"; '\
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; ' \
            'gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; ' \
            'transcript_status "KNOWN"; transcript_name "OR4F5-001"; level 2; ' \
            'protein_id "ENSP00000334393.3"; tag "basic"; transcript_support_level "NA"; ' \
            'hgnc_id "HGNC:14825"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS30547.1"; ' \
            'havana_gene "OTTHUMG00000001094.2"; havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'transcript', 
            'start': 69091,
            'end': 70008,
            'strand': b'+',
            'symbol': b'OR4F5',
            'alternate_ids': [b'ENSG00000186092.4', b'HGNC:14825'],
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_canonical': 5,
            'attributes': {
                b'gene_id': b'ENSG00000186092.4',
                b'transcript_id': b'ENST00000335137.3',
                b'gene_type': b'protein_coding',
                b'gene_status': b'KNOWN',
                b'gene_name': b'OR4F5',
                b'transcript_type': b'protein_coding',
                b'transcript_status': b'KNOWN',
                b'transcript_name': b'OR4F5-001',
                b'level': b'2',
                b'protein_id': b'ENSP00000334393.3',
                b'tag': b'basic,appris_principal_1,CCDS',
                b'transcript_support_level': b'NA',
                b'hgnc_id': b'HGNC:14825',
                b'ccdsid': b'CCDS30547.1',
                b'havana_gene': b'OTTHUMG00000001094.2',
                b'havana_transcript': b'OTTHUMT00000003223.2',
                },
            }
        self.assertEqual(obj, expected)
        
        # try without alternate gene ID from ensembl gene ID
        line = 'chr1\tHAVANA\ttranscript\t69091\t70008\t.\t+\t.\t '\
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; ' \
            'gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding";' \
            'hgnc_id "HGNC:14825"; tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS30547.1"; '
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'transcript', 
            'start': 69091,
            'end': 70008,
            'strand': b'+',
            'symbol': b'OR4F5',
            'alternate_ids': [b'HGNC:14825'],
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_canonical': 5,
            'attributes': {
                b'transcript_id': b'ENST00000335137.3',
                b'gene_type': b'protein_coding',
                b'gene_status': b'KNOWN',
                b'gene_name': b'OR4F5',
                b'transcript_type': b'protein_coding',
                b'hgnc_id': b'HGNC:14825',
                b'tag': b'appris_principal_1,CCDS',
                b'ccdsid': b'CCDS30547.1',
                },
            }
        self.assertEqual(obj, expected)
        
        # try without any alternate ID
        line = 'chr1\tHAVANA\ttranscript\t69091\t70008\t.\t+\t.\t '\
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; ' \
            'gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding";' \
            'tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS30547.1"; '
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'transcript', 
            'start': 69091,
            'end': 70008,
            'strand': b'+',
            'symbol': b'OR4F5',
            'alternate_ids': [],
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_canonical': 5,
            'attributes': {
                b'transcript_id': b'ENST00000335137.3',
                b'gene_type': b'protein_coding',
                b'gene_status': b'KNOWN',
                b'gene_name': b'OR4F5',
                b'transcript_type': b'protein_coding',
                b'tag': b'appris_principal_1,CCDS',
                b'ccdsid': b'CCDS30547.1',
                },
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_exon_line(self):
        '''test we can parse a GTF line for an exon feature
        '''
        line = 'chr1\tHAVANA\texon\t69091\t70008\t.\t+\t.\tgene_id "ENSG00000186092.4"; ' \
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; ' \
            'gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; ' \
            'transcript_status "KNOWN"; transcript_name "OR4F5-001"; exon_number 1; ' \
            'exon_id "ENSE00002319515.1"; level 2; protein_id "ENSP00000334393.3"; ' \
            'tag "basic"; transcript_support_level "NA"; tag "appris_principal_1"; ' \
            'tag "CCDS"; ccdsid "CCDS30547.1"; havana_gene "OTTHUMG00000001094.2"; ' \
            'havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'exon', 
            'start': 69091,
            'end': 70008,
            'strand': b'+',
            'symbol': b'',
            'alternate_ids': [],
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_canonical': 0,  ## exons don't get checked for principal tag
            'attributes': {},
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_cds_line(self):
        '''test we can parse a GTF line for a CDS feature
        '''
        line = 'chr1\tHAVANA\tCDS\t69091\t70005\t.\t+\t0\tgene_id "ENSG00000186092.4"; ' \
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; ' \
            'gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; ' \
            'transcript_status "KNOWN"; transcript_name "OR4F5-001"; exon_number 1; ' \
            'exon_id "ENSE00002319515.1"; level 2; protein_id "ENSP00000334393.3"; ' \
            'tag "basic"; transcript_support_level "NA"; tag "appris_principal_1"; ' \
            'tag "CCDS"; ccdsid "CCDS30547.1"; havana_gene "OTTHUMG00000001094.2"; ' \
            'havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'CDS', 
            'start': 69091,
            'end': 70005,
            'strand': b'+',
            'symbol': b'',
            'alternate_ids': [],
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_canonical': 0,  ## CDS don't get checked for principal tag
            'attributes': {},
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_start_codon_line(self):
        '''test we can parse a GTF line for a start codon feature
        '''
        line = 'chr1\tHAVANA\tstart_codon\t69091\t69093\t.\t+\t0\tgene_id "ENSG00000186092.4"; ' \
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status "KNOWN"; ' \
            'gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; ' \
            'transcript_name "OR4F5-001"; exon_number 1; exon_id "ENSE00002319515.1"; level 2; ' \
            'protein_id "ENSP00000334393.3"; tag "basic"; transcript_support_level "NA"; ' \
            'tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS30547.1"; ' \
            'havana_gene "OTTHUMG00000001094.2"; havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'start_codon', 
            'start': 69091,
            'end': 69093,
            'strand': b'+',
            'symbol': b'',
            'alternate_ids': [],
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_canonical': 0,
            'attributes': {},
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_stop_codon_line(self):
        '''test we can parse a GTF line for a stop codon feature
        '''
        line = 'chr1\tHAVANA\tstop_codon\t70006\t70008\t.\t+\t0\tgene_id "ENSG00000186092.4"; ' \
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status "KNOWN";' \
            'gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; ' \
            'transcript_name "OR4F5-001"; exon_number 1; exon_id "ENSE00002319515.1"; level 2; ' \
            'protein_id "ENSP00000334393.3"; tag "basic"; transcript_support_level "NA"; ' \
            'tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS30547.1"; ' \
            'havana_gene "OTTHUMG00000001094.2"; havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'stop_codon', 
            'start': 70006,
            'end': 70008,
            'strand': b'+',
            'symbol': b'',
            'alternate_ids': [],
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_canonical': 0,
            'attributes': {},
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_UTR_line(self):
        '''test we can parse a GTF line for a UTR feature
        '''
        line = 'chr1\tHAVANA\tUTR\t70006\t70008\t.\t+\t.\tgene_id "ENSG00000186092.4"; ' \
            'transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status "KNOWN"; ' \
            'gene_name "OR4F5"; transcript_type "protein_coding"; transcript_status "KNOWN"; ' \
            'transcript_name "OR4F5-001"; exon_number 1; exon_id "ENSE00002319515.1"; level 2; ' \
            'protein_id "ENSP00000334393.3"; tag "basic"; transcript_support_level "NA"; ' \
            'tag "appris_principal_1"; tag "CCDS"; ccdsid "CCDS30547.1"; ' \
            'havana_gene "OTTHUMG00000001094.2"; havana_transcript "OTTHUMT00000003223.2";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'UTR', 
            'start': 70006,
            'end': 70008,
            'strand': b'+',
            'symbol': b'',
            'alternate_ids': [],
            'tx_id': b'ENST00000335137.3',
            'transcript_type': b'protein_coding',
            'is_canonical': 0,
            'attributes': {},
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_minus_strand(self):
        '''test we can parse a GTF line for a feature on the minus strand
        '''
        line = 'chr1\tHAVANA\tCDS\t70006\t70008\t.\t-\t.\tgene_name "TEST";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1', 
            'feature': b'CDS', 
            'start': 70006,
            'end': 70008,
            'strand': b'-',
            'symbol': b'',
            'alternate_ids': [],
            'tx_id': b'',
            'transcript_type': b'',
            'is_canonical': 0,
            'attributes': {},
            }
        self.assertEqual(obj, expected)
    
    def test_parse_gtf_field_widths(self):
        '''test GTF fields are found regardless of the width of each field
        '''
        line = 'chr1\t.\ttranscript\t70006\t70008\t1000\t-\t2\ttranscript_id "ENST_A"; gene_name "TEST";\n'
        obj = _parse_gtfline(line.encode('utf8'))
        expected = {'chrom': b'chr1',
            'feature': b'transcript',
            'start': 70006,
            'end': 70008,
            'strand': b'-',
            'symbol': b'TEST',
            'alternate_ids': [],
            'tx_id': b'ENST_A',
            'transcript_type': b'',
            'is_canonical': 0,
            'attributes': {b'transcript_id': b'ENST_A', b'gene_name': b'TEST'},
            }
        self.assertEqual(obj, expected)
        
        with self.assertRaises(ValueError):
            _parse_gtfline(b'chr1\tHAVANA\tCDS\t70006\t70008\n')
    
    def test__open_gencode_multi_gene(self):
        '''test we can parse a GTF with multiple genes
        '''
        lines = '##format: gtf\n' \
                'chr1\tHAVANA\tgene\t10\t20\t.\t-\t.\tgene_name "TEST1";\n' \
                'chr1\tHAVANA\ttranscript\t10\t20\t.\t-\t.\ttranscript_id "ENST_A";gene_name "TEST1"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\tCDS\t15\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\tgene\t10\t30\t.\t-\t.\tgene_name "TEST2";\n' \
                'chr1\tHAVANA\ttranscript\t10\t30\t.\t-\t.\ttranscript_id "ENST_B";gene_name "TEST2"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\texon\t10\t30\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\tCDS\t15\t30\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding;"\n'
        
        write_gtf(self.temp_gtf_path, lines)
        data = _open_gencode(self.temp_gtf_path)
        
        self.assertEqual(len(data), 2)
        symbol1, tx1, is_principal = data[0]
        symbol2, tx2, is_principal = data[1]
        self.assertNotEqual(symbol1, symbol2)
    
    def test__open_gencode_blank_lines(self):
        '''test blank and comment lines mid-file are skipped, in plain and gzipped GTFs
        '''
        lines = '##format: gtf\n' \
                'chr1\tHAVANA\ttranscript\t10\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding";\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding";\n' \
                '\n' \
                '# a comment\n' \
                '\r\n' \
                '\r\r\n' \
                'chr1\tHAVANA\ttranscript\t10\t30\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding";\n' \
                'chr1\tHAVANA\texon\t10\t30\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding";\n' \
                '\n'
        
        write_gtf(self.temp_gtf_path, lines)
        data = _open_gencode(self.temp_gtf_path)
        self.assertEqual([x[0] for x in data], ['TEST1', 'TEST2'])
        
        gz_path = self.temp_gtf_path + '.gz'
        with gzip.open(gz_path, 'wt') as output:
            output.write(lines)
        try:
            data = _open_gencode(gz_path)
        finally:
            Path(gz_path).unlink()
        self.assertEqual([x[0] for x in data], ['TEST1', 'TEST2'])
    
    def test__open_gencode_without_transcript_lines(self):
        '''test GTFs without transcript lines (e.g. from UCSC) get gene fields and spans
        '''
        attrs = 'transcript_type "protein_coding"; tag "Ensembl_canonical";'
        lines = []
        for gene_id, symbol, exons, cds in [
                ('ENSG1', 'GENEA', [(100, 200), (300, 400)], [(150, 200), (300, 350)]),
                ('ENSG2', 'GENEB', [(1000, 2000)], [(1100, 1900)])]:
            for i, (start, end) in enumerate(exons):
                lines.append(f'chr1\tucsc\texon\t{start}\t{end}\t.\t+\t.\tgene_id "{gene_id}"; '
                    f'transcript_id "{gene_id}_T"; gene_name "{symbol}"; exon_number "{i + 1}"; {attrs}\n')
            for start, end in cds:
                lines.append(f'chr1\tucsc\tCDS\t{start}\t{end}\t.\t+\t0\tgene_id "{gene_id}"; '
                    f'transcript_id "{gene_id}_T"; gene_name "{symbol}"; {attrs}\n')
        
        write_gtf(self.temp_gtf_path, lines)
        data = _open_gencode(self.temp_gtf_path)
        self.assertEqual([x[0] for x in data], ['GENEA', 'GENEB'])
        self.assertEqual([(x[1].start, x[1].end) for x in data], [(100, 400), (1000, 2000)])
        self.assertEqual([x[2] for x in data], [10, 10])
        self.assertEqual(data[0][1].attributes['gene_id'], 'ENSG1')
        self.assertNotIn('exon_number', data[0][1].attributes)
        
        gencode = Gencode(self.temp_gtf_path)
        self.assertEqual(gencode['GENEA'].alternate_ids, ['ENSG1'])
        self.assertEqual([x.symbol for x in gencode.in_region('chr1', 50, 150)], ['GENEA'])
    
    def test__open_gencode_position_sorted(self):
        '''test transcripts are combined when their lines are interleaved
        '''
        def line(chrom, feature, start, end, tx_id):
            return f'{chrom}\tHAVANA\t{feature}\t{start}\t{end}\t.\t+\t.\tgene_id "G1"; ' \
                f'transcript_id "{tx_id}"; gene_name "A"; transcript_type "protein_coding";\n'
        
        lines = [line('chr1', 'transcript', 10, 40, 'T1'),
                line('chr1', 'exon', 10, 20, 'T1'),
                line('chr1', 'transcript', 15, 25, 'T2'),
                line('chr1', 'exon', 15, 25, 'T2'),
                line('chr1', 'exon', 30, 40, 'T1'),
                line('chr1', 'transcript', 50, 60, 'T3'),
                line('chr1', 'exon', 50, 60, 'T3'),
                # the same transcript ID on another chromosome (e.g. PAR genes)
                line('chrX', 'exon', 100, 200, 'T1'),
                line('chrX', 'exon', 300, 400, 'T1'),
                line('chrY', 'exon', 500, 600, 'T1')]
        
        write_gtf(self.temp_gtf_path, lines)
        data = _open_gencode(self.temp_gtf_path)
        self.assertEqual([(x[1].name, x[1].chrom, x[1].start, x[1].end) for x in data],
            [('T1', 'chr1', 10, 40), ('T2', 'chr1', 15, 25), ('T3', 'chr1', 50, 60),
             ('T1', 'chrX', 100, 400), ('T1', 'chrY', 500, 600)])
        self.assertEqual(data[0][1].exons, [{'start': 10, 'end': 20}, {'start': 30, 'end': 40}])
        self.assertEqual(data[3][1].exons, [{'start': 100, 'end': 200}, {'start': 300, 'end': 400}])
    
    def test__open_gencode_line_endings(self):
        '''test CRLF line endings and a missing final newline, in plain and gzipped GTFs
        '''
        lines = 'chr1\tHAVANA\ttranscript\t10\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding";\r\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding"'
        
        write_gtf(self.temp_gtf_path, lines)
        data = _open_gencode(self.temp_gtf_path)
        self.assertEqual([x[0] for x in data], ['TEST1'])
        self.assertEqual(data[0][1].exons, [{'start': 10, 'end': 20}])
        
        # gzipped files are detected by content, not by the file extension
        with gzip.open(self.temp_gtf_path, 'wt') as output:
            output.write(lines)
        data = _open_gencode(self.temp_gtf_path)
        self.assertEqual([x[0] for x in data], ['TEST1'])
    
    def test__open_gencode_skips_bad_transcripts(self):
        '''test malformed transcripts are skipped, without stopping the GTF load
        '''
        lines = '##format: gtf\n' \
                'chr1\tHAVANA\ttranscript\t10\t20\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding";\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding";\n' \
                'chr1\tHAVANA\tCDS\t12\t18\t.\t+\t.\ttranscript_id "ENST_A"; gene_name "TEST1"; transcript_type "protein_coding";\n' \
                'chr1\tHAVANA\ttranscript\t30\t40\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding";\n' \
                'chr1\tHAVANA\texon\t30\t40\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding";\n' \
                'chr1\tHAVANA\tCDS\t32\t38\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding";\n' \
                'chr1\tHAVANA\tstop_codon\t50\t52\t.\t+\t.\ttranscript_id "ENST_B"; gene_name "TEST2"; transcript_type "protein_coding";\n' \
                'chr1\tHAVANA\ttranscript\t60\t70\t.\t.\t.\ttranscript_id "ENST_C"; gene_name "TEST3"; transcript_type "protein_coding";\n' \
                'chr1\tHAVANA\texon\t60\t70\t.\t.\t.\ttranscript_id "ENST_C"; gene_name "TEST3"; transcript_type "protein_coding";\n'
        
        # ENST_B has a stop codon outside its exons, ENST_C lacks a valid strand
        write_gtf(self.temp_gtf_path, lines)
        data = _open_gencode(self.temp_gtf_path)
        self.assertEqual([x[1].name for x in data], ['ENST_A'])
    
    def test__open_gencode_missing_file(self):
        '''test a GTF which can't be opened raises an error
        '''
        with self.assertRaises(ValueError):
            _open_gencode(self.temp_gtf_path + '.missing')
    
    def test__open_gencode_not_coding(self):
        '''test we can parse GTFs without protein coding transcripts
        '''
        lines = '##format: gtf\n' \
                'chr1\tHAVANA\tgene\t10\t20\t.\t-\t.\tgene_name "TEST";\n' \
                'chr1\tHAVANA\ttranscript\t10\t20\t.\t-\t.\ttranscript_id "ENST_A";gene_name "TEST"; transcript_type "processed_transcript"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST"; transcript_type "processed_transcript;"\n' \
                'chr1\tHAVANA\tCDS\t15\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST"; transcript_type "processed_transcript;"\n' \
                        
        write_gtf(self.temp_gtf_path, lines)
        data = _open_gencode(self.temp_gtf_path)
        
        self.assertEqual(len(data), 0)
        
        # but if we allow all transcript types (not just protein_coding, we)
        write_gtf(self.temp_gtf_path, lines)
        data = _open_gencode(self.temp_gtf_path, coding_only=False)
        
        self.assertEqual(len(data), 1)
        
    def test__open_gencode_multi_transcript(self):
        '''test we can parse a GTF with multiple transcripts for the same gene
        '''
        lines = '##format: gtf\n' \
                'chr1\tHAVANA\tgene\t10\t20\t.\t-\t.\tgene_name "TEST";\n' \
                'chr1\tHAVANA\ttranscript\t10\t20\t.\t-\t.\ttranscript_id "ENST_A";gene_name "TEST"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\tCDS\t15\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\tgene\t10\t30\t.\t-\t.\tgene_name "TEST";\n' \
                'chr1\tHAVANA\ttranscript\t10\t30\t.\t-\t.\ttranscript_id "ENST_B";gene_name "TEST"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\texon\t10\t30\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\tCDS\t15\t30\t.\t-\t.\ttranscript_id "ENST_B"; gene_name "TEST"; transcript_type "protein_coding;"\n'
        
        write_gtf(self.temp_gtf_path, lines)
        data = _open_gencode(self.temp_gtf_path)
        
        self.assertEqual(len(data), 2)
        symbol1, tx1, is_principal = data[0]
        symbol2, tx2, is_principal = data[1]
        self.assertEqual(symbol1, symbol2)
        self.assertEqual(tx1.name, 'ENST_A')
        self.assertEqual(tx1.cds, [{'start': 15, 'end': 20}])
        self.assertEqual(tx2.name, 'ENST_B')
        self.assertEqual(tx2.cds, [{'start': 15, 'end': 30}])
        
    def test__open_gencode_multi_exon(self):
        '''test we can parse a GTF into transcripts
        '''
        lines = '##format: gtf\n' \
                'chr1\tHAVANA\tgene\t10\t100\t.\t-\t.\tgene_name "TEST";\n' \
                'chr1\tHAVANA\ttranscript\t10\t100\t.\t-\t.\ttranscript_id "ENST_A";gene_name "TEST"; transcript_type "protein_coding"; tag "appris_principal_1";\n' \
                'chr1\tHAVANA\tUTR\t10\t15\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST"; transcript_type "protein_coding";\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\tCDS\t15\t20\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\texon\t30\t40\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\tCDS\t30\t40\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\texon\t90\t100\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST"; transcript_type "protein_coding;"\n' \
                'chr1\tHAVANA\tUTR\t90\t100\t.\t-\t.\ttranscript_id "ENST_A"; gene_name "TEST"; transcript_type "protein_coding;"\n'
        
        write_gtf(self.temp_gtf_path, lines)
        data = _open_gencode(self.temp_gtf_path)
        
        self.assertEqual(len(data), 1)
        symbol, tx, is_principal = data[0]
        self.assertEqual(symbol, 'TEST')
        self.assertEqual(tx.name, 'ENST_A')
        self.assertEqual(tx.strand, '-')
        self.assertEqual(tx.exons, [{'start': 10, 'end': 20}, {'start': 30, 'end': 40}, {'start': 90, 'end': 100}])
        self.assertEqual(tx.cds, [{'start': 15, 'end': 20}, {'start': 30, 'end': 40}])
    
    def test__open_gencode_unquoted(self):
        '''test we can parse a GTF without quoted attributes
        '''
        lines = '##format: gtf\n' \
                'chr1\tHAVANA\tgene\t10\t100\t.\t-\t.\tgene_name "TEST";\n' \
                'chr1\tHAVANA\ttranscript\t10\t100\t.\t-\t.\ttranscript_id ENST_A;gene_name TEST; transcript_type protein_coding; tag appris_principal_1;\n' \
                'chr1\tHAVANA\tUTR\t10\t15\t.\t-\t.\ttranscript_id ENST_A; gene_name TEST; transcript_type protein_coding;\n' \
                'chr1\tHAVANA\texon\t10\t20\t.\t-\t.\ttranscript_id ENST_A; gene_name TEST; transcript_type protein_coding;\n' \
                'chr1\tHAVANA\tCDS\t15\t20\t.\t-\t.\ttranscript_id ENST_A; gene_name TEST; transcript_type protein_coding;\n' \
                'chr1\tHAVANA\texon\t30\t40\t.\t-\t.\ttranscript_id ENST_A; gene_name TEST; transcript_type protein_coding;\n' \
                'chr1\tHAVANA\tCDS\t30\t40\t.\t-\t.\ttranscript_id ENST_A; gene_name TEST; transcript_type protein_coding;\n' \
                'chr1\tHAVANA\texon\t90\t100\t.\t-\t.\ttranscript_id ENST_A; gene_name TEST; transcript_type protein_coding;\n' \
                'chr1\tHAVANA\tUTR\t90\t100\t.\t-\t.\ttranscript_id ENST_A; gene_name TEST; transcript_type protein_coding;\n'
        
        write_gtf(self.temp_gtf_path, lines)
        data = _open_gencode(self.temp_gtf_path)
        
        self.assertEqual(len(data), 1)
        symbol, tx, is_principal = data[0]
        self.assertEqual(symbol, 'TEST')
        self.assertEqual(tx.name, 'ENST_A')
        self.assertEqual(tx.strand, '-')
        self.assertEqual(tx.exons, [{'start': 10, 'end': 20}, {'start': 30, 'end': 40}, {'start': 90, 'end': 100}])
        self.assertEqual(tx.cds, [{'start': 15, 'end': 20}, {'start': 30, 'end': 40}])
    

