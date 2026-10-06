""" class to test the Transcript class
"""

import unittest

from gencodegenes.transcript import Transcript

class TestTranscriptPy(unittest.TestCase):
    """ unit test the Transcript class
    """
    
    def setUp(self):
        """ construct a Transcript object for unit tests
        """
        
        self.gene = self.construct_gene()
    
    def construct_gene(self, name='TEST', chrom='1', start=1000, end=2000,
            strand='+', exons=[(1000, 1200), (1800, 2000)],
            cds=[(1100, 1200), (1800, 1900)]):
        
        tx = Transcript(name, chrom, start, end, strand)
        tx.exons = exons
        tx.cds = cds
        
        return tx
    
    def test_set_exons(self):
        """ test that set_exons() works correctly
        """
        
        exons = [(0, 200), (800, 1000)]
        cds = [(1100, 1200), (1800, 1900)]
        self.gene.exons = exons
        
        self.assertEqual(self.gene.exons, [{'start': 0, 'end': 200}, {'start': 800, 'end': 1000}])
        
        self.gene = self.construct_gene(strand='-')
        self.gene.exons = exons
        self.assertEqual(self.gene.exons, [{'start': 0, 'end': 200}, {'start': 800, 'end': 1000}])
    
    def test_set_exons_missing_exon(self):
        """ test that set_exons() works correctly when we lack coordinates
        """
        
        # also test that we can determine the exons when we don't have any given
        # for a transcript with a single CDS region
        exons = []
        cds = [(1100, 1200)]
        self.gene.exons = exons
        self.assertEqual(self.gene.exons, [])
        self.gene.cds = cds
        self.assertEqual(self.gene.exons, [{'start': 1000, 'end': 2000}])
        
        # check that missing exons, but 2+ CDS regions raises an error.
        cds = [(1100, 1200), (1300, 1400)]
        self.gene.exons = []
        with self.assertRaises(ValueError):
            self.gene.cds = cds
    
    def test_set_cds(self):
        """ test that set_cds() works correctly
        """
        
        exons = [(0, 200), (800, 1000)]
        cds = [(100, 200), (800, 900)]
        
        # make sure we raise an error if we try to set the CDS before the exons
        with self.assertRaises(ValueError):
            tx = Transcript('test', '1', 0, 1000, '+')
            tx.cds = cds
        
        # check CDS positions
        self.gene.exons = exons
        self.gene.cds = cds
        self.assertEqual(self.gene.cds, [{'start': 100, 'end': 200}, {'start': 800, 'end': 900}])
        
        # check that CDS ends outside an exon are corrected
        exons = [(0, 200), (300, 400), (800, 1000)]
        cds = [(100, 200), (300, 402)]
        self.gene.exons = exons
        self.gene.cds = cds
        self.assertEqual(self.gene.cds, [{'start': 100, 'end': 200},
            {'start': 300, 'end': 400}, {'start': 800, 'end': 802}])
        
        cds = [(298, 400), (800, 1000)]
        self.gene.exons = exons
        self.gene.cds = cds
        self.assertEqual(self.gene.cds, [{'start': 198, 'end': 200},
            {'start': 300, 'end': 400}, {'start': 800, 'end': 1000}])
    
    def test_fix_cds_boundary(self):
        """ test that _fix_out_of_exon_cds_boundary() works correctly
        """
        
        exons = [(1100, 1200), (1300, 1400), (1800, 1900)]
        cds = [(1300, 1400)]
        
        tx = Transcript('test', '1', 0, 1000, '+')
        
        tx.exons = exons
        tx.cds = cds
        
        self.assertEqual(tx._fix_cds_boundary(1295), {'start': 1195, 'end': 1200})
        self.assertEqual(tx._fix_cds_boundary(1205), {'start': 1300, 'end': 1305})
        
        self.assertEqual(tx._fix_cds_boundary(1402), {'start': 1800, 'end': 1802})
        self.assertEqual(tx._fix_cds_boundary(1798), {'start': 1398, 'end': 1400})
        
        # raise an error if the position is within the exons
        with self.assertRaises(ValueError):
            self.gene._fix_cds_boundary(1105)
    
    def test_cdsange(self):
        """ unit test checking the CDS end points by strand
        """
        
        # self.gene.cds_min = 100
        # self.gene.cds_max = 200
        
        self.gene = self.construct_gene(exons=[(0, 300)], cds=[(100, 200)])
        
        self.assertEqual(self.gene.cds_start, 100)
        self.assertEqual(self.gene.cds_end, 200)
        
        self.gene = self.construct_gene(exons=[(0, 300)], cds=[(100, 200)], strand='-')
        self.assertEqual(self.gene.cds_start, 200)
        self.assertEqual(self.gene.cds_end, 100)
        
        with self.assertRaises(ValueError):
            self.construct_gene(strand='x')
    
    def test___add__(self):
        """ test that __add__() works correctly
        """
        
        exons = [(10, 20), (50, 60), (90, 100)]
        cds_2 = [(50, 60), (90, 95)]
        
        a = Transcript("a", "1", 10, 100, "+")
        b = Transcript("b", "1", 10, 100, "+")
        c = Transcript("c", "1", 10, 100, "+")
        d = Transcript("d", "1", 10, 100, "+")
        
        a.exons = exons
        a.cds = [(55, 60), (90, 100)]
        
        b.exons = exons
        b.cds = [(50, 60), (90, 95)]
        
        c.exons = [(45, 65)]
        c.cds = [(45, 65)]
        
        d.exons = [(30, 40)]
        d.cds = [(30, 40)]
        
        # check that adding two Transcripts gives the union of CDS regions
        self.assertEqual((a + b).cds, [{'start': 50, 'end': 60}, {'start': 90, 'end': 100}])
        self.assertEqual((a + c).cds, [{'start': 45, 'end': 65}, {'start': 90, 'end': 100}])
        
        # check that addition is reversible
        self.assertEqual((c + a).cds, [{'start': 45, 'end': 65}, {'start': 90, 'end': 100}])
        
        # check that adding previously unknown exons works
        self.assertEqual((a + d).cds, [{'start': 30, 'end': 40}, {'start': 55, 'end': 60}, {'start': 90, 'end': 100}])
        
        # check that we can add transcript + None correctly
        self.assertEqual(a + None, a)
        self.assertEqual(None + a, a)
    
    def test___add__not_overlapping(self):
        ''' test that __add__() works correctly when transcripts do not overlap
        '''
        
        a = Transcript("a", "1", 10, 50, "+")
        b = Transcript("b", "1", 60, 80, "+")
        
        a.exons = [(10, 50)]
        a.cds = [(10, 50)]
        a.genomic_sequence = 'N' * 40
        
        b.exons = [(60, 80)]
        b.cds = [(60, 80)]
        b.genomic_sequence = 'N' * 20
        
        self.assertEqual(len((a + b).genomic_sequence), 70)
    
    def test___add__cds_length_fixed(self):
        """ check that we can merge transcripts, even with fixed CDS coords
        """
        
        a = Transcript("a", "1", 10, 20, "+")
        a.exons = [(10, 20)]
        a.cds = [(10, 20)]
        
        a.cds_sequence = 'ACTGTACGCAT'
        a.genomic_offset = 5
        a.genomic_sequence = 'CGTAGACTGTACGCATCGATT'
        
        b = Transcript("b", "1", 0, 10, "+")
        b.exons = [(0, 10)]
        b.cds = [(0, 10)]
        
        b.cds_sequence = 'ACTGTACGCAT'
        b.genomic_offset = 5
        b.genomic_sequence = 'CGTAGACTGTACGCATCGTAG'
        
        # without a fix to tx.cpp to adjust an exon coordinate simultaneously,
        # the line below would give an error.
        c = a + b
    
    def test_merge_coordinates(self):
        """ test that we can merge transcripts with odd overlaps
        """
        
        a = Transcript("a", "1", 10, 20, "+")
        
        exons1 = [{'start': 10, 'end': 20}, {'start': 25, 'end': 40}]
        exons2 = [{'start': 10, 'end': 30}]
        
        self.assertEqual(a._merge_coordinates(exons1, exons2),
            a._merge_coordinates(exons2, exons1))
    
    def test_in_exons(self):
        """ test that in_exons() works correctly
        """
        
        # self.gene.exons = [(1000, 1200), (1800, 2000)]
        
        # check for positions inside the exon ranges
        self.assertTrue(self.gene.in_exons(1000))
        self.assertTrue(self.gene.in_exons(1001))
        self.assertTrue(self.gene.in_exons(1200))
        self.assertTrue(self.gene.in_exons(1800))
        self.assertTrue(self.gene.in_exons(1801))
        self.assertTrue(self.gene.in_exons(1999))
        self.assertTrue(self.gene.in_exons(2000))
        
        # check positions outside the exon ranges
        self.assertFalse(self.gene.in_exons(999))
        self.assertFalse(self.gene.in_exons(1201))
        self.assertFalse(self.gene.in_exons(1799))
        self.assertFalse(self.gene.in_exons(2001))
        self.assertFalse(self.gene.in_exons(-1100))
    
    def test_get_closest_exon(self):
        """ test that get_closest_exon() works correctly
        """
        #
        exon_1 = {'start': 1000, 'end': 1200}
        exon_2 = {'start': 1800, 'end': 2000}
        
        # find for positions closer to the first exon
        self.assertEqual(self.gene.get_closest_exon(0), exon_1)
        self.assertEqual(self.gene.get_closest_exon(999), exon_1)
        self.assertEqual(self.gene.get_closest_exon(1000), exon_1)
        self.assertEqual(self.gene.get_closest_exon(1100), exon_1)
        self.assertEqual(self.gene.get_closest_exon(1200), exon_1)
        self.assertEqual(self.gene.get_closest_exon(1201), exon_1)
        
        # a site equidistant from the exons will pick the later exon
        self.assertEqual(self.gene.get_closest_exon(1500), exon_2)
        
        # find for positions closer to the second exon
        self.assertEqual(self.gene.get_closest_exon(1501), exon_2)
        self.assertEqual(self.gene.get_closest_exon(1799), exon_2)
        self.assertEqual(self.gene.get_closest_exon(1800), exon_2)
        self.assertEqual(self.gene.get_closest_exon(1900), exon_2)
        self.assertEqual(self.gene.get_closest_exon(2000), exon_2)
        self.assertEqual(self.gene.get_closest_exon(2001), exon_2)
        self.assertEqual(self.gene.get_closest_exon(10000), exon_2)
    
    def test_construct_without_cds(self):
        """ check exons, offset and sequence are kept when constructed without CDS
        """
        tx = Transcript('TEST', '1', 1, 10, '+', exons=[(1, 10)],
            sequence='AACCGGTTAACCGG', offset=2)
        self.assertEqual(tx.exons, [{'start': 1, 'end': 10}])
        self.assertEqual(tx.cds, [])
        self.assertEqual(tx.genomic_offset, 2)
        self.assertEqual(tx.genomic_sequence, 'AACCGGTTAACCGG')
        
        # a lone CDS still needs exons, unless the CDS fits within the transcript
        with self.assertRaises(ValueError):
            Transcript('TEST', '1', 1, 10, '+', cds=[(1, 5), (7, 10)])
    
    def test_exon_lookups_without_exons(self):
        """ check exon lookups raise ValueError if the transcript lacks exons
        """
        tx = Transcript('TEST', '1', 1000, 2000, '+')
        with self.assertRaises(ValueError):
            tx.in_exons(1100)
        with self.assertRaises(ValueError):
            tx.get_closest_exon(1100)
    
    def test_in_coding_region(self):
        """ test that in_coding_region() works correctly
        """
        
        # self.gene.cds = [(1100, 1200), (1800, 1900)]
        
        # check for positions inside the exon ranges
        self.assertTrue(self.gene.in_coding_region(1100))
        self.assertTrue(self.gene.in_coding_region(1101))
        self.assertTrue(self.gene.in_coding_region(1200))
        self.assertTrue(self.gene.in_coding_region(1800))
        self.assertTrue(self.gene.in_coding_region(1801))
        self.assertTrue(self.gene.in_coding_region(1899))
        self.assertTrue(self.gene.in_coding_region(1900))
        
        # check positions outside the exon ranges
        self.assertFalse(self.gene.in_coding_region(1099))
        self.assertFalse(self.gene.in_coding_region(1201))
        self.assertFalse(self.gene.in_coding_region(1799))
        self.assertFalse(self.gene.in_coding_region(1901))
        self.assertFalse(self.gene.in_coding_region(-1100))
    
    # def test_get_exon_containing_position(self):
    #     """ test that get_exon_containing_position() works correctly
    #     """
    #
    #     exons = [(1000, 1200), (1800, 2000)]
    #
    #     self.assertEqual(self.gene.get_exon_containing_position(1000, exons), 0)
    #     self.assertEqual(self.gene.get_exon_containing_position(1200, exons), 0)
    #     self.assertEqual(self.gene.get_exon_containing_position(1800, exons), 1)
    #     self.assertEqual(self.gene.get_exon_containing_position(2000, exons), 1)
    #
    #     # raise an error if the position isn't within the exons
    #     with self.assertRaises(RuntimeError):
    #         self.gene.get_exon_containing_position(2100, exons)
    
    def test_get_coding_distance(self):
        """ test that get_coding_distance() works correctly
        """
        
        # self.gene.cds = [(1100, 1200), (1800, 1900)]
        
        # raise an error for positions outside the CDS
        self.assertEqual(self.gene.get_coding_distance(900), {'pos': -100, 'offset': -100})
        self.assertEqual(self.gene.get_coding_distance(1000), {'pos': -100, 'offset': 0})
        self.assertEqual(self.gene.get_coding_distance(1051), {'pos': -49, 'offset': 0})
        self.assertEqual(self.gene.get_coding_distance(1300), {'pos': 100, 'offset': 100})
        self.assertEqual(self.gene.get_coding_distance(1700), {'pos': 101, 'offset': -100})
        self.assertEqual(self.gene.get_coding_distance(2000), {'pos': 301, 'offset': 0})
        self.assertEqual(self.gene.get_coding_distance(2100), {'pos': 301, 'offset': 100})
        
        # zero distance between a site and itself
        self.assertEqual(self.gene.get_coding_distance(1100), {'pos': 0, 'offset': 0})
        
        # within a single exon, the distance is between the start and end
        self.assertEqual(self.gene.get_coding_distance(1200), {'pos': 100, 'offset': 0})
        
        # if we traverse exons, the distance bumps up at exon boundaries
        self.assertEqual(self.gene.get_coding_distance(1800), {'pos': 101, 'offset': 0})
        
        # check full distance across gene
        self.assertEqual(self.gene.get_coding_distance(1900), {'pos': 201, 'offset': 0})
        
        # check that the distance bumps up for each exon boundary crossed
        cds = [(1100, 1200), (1300, 1400), (1800, 1900)]
        exons = [(1100, 1200), (1300, 1400), (1800, 1900)]
        self.gene = self.construct_gene(exons=exons, cds=cds)
        self.assertEqual(self.gene.get_coding_distance(1900), {'pos': 302, 'offset': 0})
        
        # now try a gene where the site is in an upstream exon
        self.gene = self.construct_gene(exons=[(10, 20), (30, 40), (90, 100)],
            cds= [(30, 40), (90, 95)])
        self.assertEqual(self.gene.get_coding_distance(15), {'pos': -6, 'offset': 0})
    
    def test_chrom_pos_to_cds(self):
        """ test that chrom_pos_to_cds() works correctly
        """
        # self.gene.cds = [(1100, 1200), (1800, 1900)]
        
        # note that all of these chr positions are 0-based (ie pos - 1)
        self.assertEqual(self.gene.get_coding_distance(1100), {'pos': 0, 'offset': 0})
        self.assertEqual(self.gene.get_coding_distance(1101), {'pos': 1, 'offset': 0})
        self.assertEqual(self.gene.get_coding_distance(1199), {'pos': 99, 'offset': 0})
        
        # check that outside exon boundaries gets the closest exon position, if
        # the variant is close enough
        self.assertEqual(self.gene.get_coding_distance(1200), {'pos': 100, 'offset': 0})
        self.assertEqual(self.gene.get_coding_distance(1201), {'pos': 100, 'offset': 1})
        self.assertEqual(self.gene.get_coding_distance(1798), {'pos': 101, 'offset': -2})
        self.assertEqual(self.gene.get_coding_distance(1799), {'pos': 101, 'offset': -1})
        
        # # check that sites sufficiently distant from an exon raise an error, or
        # # sites upstream of a gene, just outside the CDS, but within an exon
        # with self.assertRaises(RuntimeError):
        #     self.gene.get_coding_distance(1215)
        # with self.assertRaises(RuntimeError):
        #     self.gene.get_coding_distance(1098)
        self.assertEqual(self.gene.get_coding_distance(1215), {'pos': 100, 'offset': 15})
        self.assertEqual(self.gene.get_coding_distance(1098), {'pos': -2, 'offset': 0})
        
        # check that sites in a different exon are counted correctly
        self.assertEqual(self.gene.get_coding_distance(1799), {'pos': 101, 'offset': -1})
        
        # check that sites on the reverse strand still give the correct CDS
        self.gene = self.construct_gene(strand="-")
        self.assertEqual(self.gene.get_coding_distance(1900), {'pos': 0, 'offset': 0})
        self.assertEqual(self.gene.get_coding_distance(1890), {'pos': 10, 'offset': 0})
        self.assertEqual(self.gene.get_coding_distance(1799), {'pos': 100, 'offset': 1})
        self.assertEqual(self.gene.get_coding_distance(1792), {'pos': 100, 'offset': 8})
        self.assertEqual(self.gene.get_coding_distance(1792), {'pos': 100, 'offset': 8})
        
        self.assertEqual(self.gene.get_coding_distance(1205), {'pos': 101, 'offset': -5})
        self.assertEqual(self.gene.get_coding_distance(1200), {'pos': 101, 'offset': 0})
    
    def test_get_boundary_distance(self):
        """ check the function to get distances to the nearest intron/exon boundary
        """
        print(self.gene)
        # check a site upstream of the gene
        self.assertEqual(self.gene.get_boundary_distance(50), 950)
        
        # check a site at the start of a gene
        self.assertEqual(self.gene.get_boundary_distance(1000), 0)
        
        # check some sites within the first exon
        self.assertEqual(self.gene.get_boundary_distance(1100), 101)
        self.assertEqual(self.gene.get_boundary_distance(1150), 51)
        
        # check sites in the first intron
        self.assertEqual(self.gene.get_boundary_distance(1250), 50)
        self.assertEqual(self.gene.get_boundary_distance(1400), 200)
        
        # check a site in the first exon, as it becomes closer to the next intron
        self.assertEqual(self.gene.get_boundary_distance(1101), 100)
        
        # check a site downstream of the gene
        self.assertEqual(self.gene.get_boundary_distance(2200), 200)
    
    def test_get_codon_info(self):
        """ check the function that checks the codon information for a position
        """
        
        self.gene.cds_sequence = 'ATGTCCATGTTGATGTTG'

        # make sure a site well outside the gene raises an error
        with self.assertRaises(ValueError):
            self.gene.get_codon_info(50)
        
        # a position near the start site, but upstream of the CDS will raise a
        # different error
        with self.assertRaises(ValueError):
            self.gene.get_codon_info(1050)
        
        # the first base after the CDS end (CDS position 202 of 202) is outside
        # the CDS too
        with self.assertRaises(ValueError):
            self.gene.get_codon_info(1901)
        
        # check the first base of the CDS
        self.assertEqual(self.gene.get_codon_info(1100),
            {'cds_pos': 0, 'codon_seq': 'ATG', 'intra_codon': 0,
                "codon_number": 0, 'initial_aa': 'M', 'offset': 0})
        
        # check the second base of the CDS
        self.assertEqual(self.gene.get_codon_info(1101),
            {'cds_pos': 1, 'codon_seq': 'ATG', 'intra_codon': 1,
                "codon_number": 0, 'initial_aa': 'M', 'offset': 0})
        
        # check the third base of the CDS
        self.assertEqual(self.gene.get_codon_info(1102),
            {'cds_pos': 2, 'codon_seq': 'ATG', 'intra_codon': 2,
                "codon_number": 0, 'initial_aa': 'M', 'offset': 0})
        
        # check the fourth base of the CDS
        self.assertEqual(self.gene.get_codon_info(1103),
            {'cds_pos': 3, 'codon_seq': 'TCC', 'intra_codon': 0,
                "codon_number": 1, 'initial_aa': 'S', 'offset': 0})
        
        # check a site 2 bp into the first intron. We assign this as the
        # position of the closest exon boundary, but without any codon info
        self.assertEqual(self.gene.get_codon_info(1202),
            {'cds_pos': 100, 'codon_seq': None, 'intra_codon': None,
                "codon_number": None, 'initial_aa': None, 'offset': 2})
    
    def test_noncoding_transcript(self):
        """ check transcripts without a CDS give consistent values
        """
        tx = self.construct_gene(cds=[])
        self.assertEqual((tx.cds_start, tx.cds_end), (0, 0))
        self.assertFalse(tx.in_coding_region(1100))
        
        self.assertEqual(tx.consequence(1100, 'A', 'G'), 'non_coding_transcript_exon_variant')
        self.assertEqual(tx.consequence(1500, 'A', 'G'), 'intron_variant')
        self.assertEqual(tx.consequence(1201, 'A', 'G'), 'splice_donor_variant')
        self.assertEqual(tx.consequence(500, 'A', 'G'), 'upstream_gene_variant')
    
    def test_consequence_mnv(self):
        """ check multi-base substitutions apply every base within the CDS
        """
        # CDS is ATG GCC TGG TAA (M A W *) at positions 4-15
        tx = Transcript('TEST', '1', 1, 20, '+', exons=[(1, 20)], cds=[(4, 15)],
            sequence='CCCATGGCCTGGTAAGGGGG')
        self.assertEqual(tx.consequence(9, 'C', 'T'), 'synonymous_variant')
        self.assertEqual(tx.consequence(10, 'TG', 'CA'), 'missense_variant')
        
        # GCC>GCT is synonymous, but TGG>TAG gains a stop in the next codon
        self.assertEqual(tx.consequence(9, 'CTG', 'TTA'), 'stop_gained')
        
        # the MNV starts upstream of the CDS, but alters the start codon
        self.assertEqual(tx.consequence(2, 'CCA', 'GGG'), 'start_lost')
        
        # same CDS on the - strand, at + strand positions 6-17. Position 8 is the
        # first base of the stop codon (TAA>CAA), position 9 is the last base of
        # the TGG codon (TGG>TGA)
        tx = Transcript('TEST', '1', 1, 20, '-', exons=[(1, 20)], cds=[(6, 17)],
            sequence='GGGATGGCCTGGTAACCCCC')
        self.assertEqual(tx.consequence(8, 'AC', 'GT'), 'stop_gained')
        self.assertEqual(tx.consequence(8, 'A', 'G'), 'stop_lost')
    
    def test_consequence_indel_frame(self):
        """ check indels are classed as inframe or frameshift by length change
        """
        self.assertEqual(self.gene.consequence(1150, 'G', 'GCCC'), 'inframe_insertion')
        self.assertEqual(self.gene.consequence(1150, 'G', 'GC'), 'frameshift_variant')
        self.assertEqual(self.gene.consequence(1150, 'G', 'GCC'), 'frameshift_variant')
        self.assertEqual(self.gene.consequence(1150, 'GCCC', 'G'), 'inframe_deletion')
        self.assertEqual(self.gene.consequence(1150, 'GC', 'G'), 'frameshift_variant')
