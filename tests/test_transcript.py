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
        
        # raise an error if there is no exon to shift the position into
        with self.assertRaises(ValueError):
            tx._fix_cds_boundary(1095)
        with self.assertRaises(ValueError):
            tx._fix_cds_boundary(1905)
        
        # and if the CDS is set to start before the first exon, or end after
        # the last exon
        with self.assertRaises(ValueError):
            tx.cds = [(1095, 1150)]
        with self.assertRaises(ValueError):
            tx.cds = [(1850, 1905)]
    
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
        a.genomic_sequence = 'N' * 41
        
        b.exons = [(60, 80)]
        b.cds = [(60, 80)]
        b.genomic_sequence = 'N' * 21
        
        self.assertEqual(len((a + b).genomic_sequence), 71)
    
    def test___add__one_without_sequence(self):
        ''' check transcripts combine when only one has a genomic sequence, without
        merging sequences
        '''
        a = Transcript("a", "1", 10, 20, "+", exons=[(10, 20)], cds=[(10, 15)],
            sequence='A' * 11)
        b = Transcript("b", "1", 30, 40, "+", exons=[(30, 40)], cds=[(30, 35)])
        expected = [{'start': 10, 'end': 15}, {'start': 30, 'end': 35}]
        for c in (a + b, b + a):
            self.assertEqual(c.cds, expected)
            self.assertEqual(c.genomic_sequence, '')
    
    def test___add__not_overlapping_with_offset(self):
        ''' check the gap between non-overlapping transcripts is filled with Ns,
        so the second transcript's sequence stays at the right positions
        '''
        # a spans 8-22 with the 2 bp flanks, b spans 38-52
        a = Transcript("a", "1", 10, 20, "+", exons=[(10, 20)], cds=[(12, 17)],
            offset=2, sequence='CCACGTACGTACGGG')
        b = Transcript("b", "1", 40, 50, "+", exons=[(40, 50)], cds=[(42, 47)],
            offset=2, sequence='TTGATTACAGATTAA')
        
        # a's sequence (8-22), Ns for 23-39, then b's sequence from 40
        c = a + b
        self.assertEqual(c.genomic_sequence, 'CCACGTACGTACGGG' + 'N' * 17 + 'GATTACAGATTAA')
        self.assertEqual(c.cds_sequence, 'GTACGT' + 'TTACAG')
    
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
    
    def test_equality(self):
        """ check transcripts are only equal with the same ID, position and structure
        """
        tx = self.construct_gene()
        self.assertEqual(tx, self.gene)
        self.assertEqual(hash(tx), hash(self.gene))
        self.assertEqual(len({tx, self.gene}), 1)
        
        # same span, but different ID, strand, exons or CDS
        self.assertNotEqual(self.construct_gene(name='OTHER'), self.gene)
        self.assertNotEqual(self.construct_gene(strand='-'), self.gene)
        self.assertNotEqual(self.construct_gene(exons=[(1000, 1300), (1800, 2000)]), self.gene)
        self.assertNotEqual(self.construct_gene(cds=[(1150, 1200), (1800, 1900)]), self.gene)
        
        self.assertNotEqual(self.gene, None)
        self.assertNotEqual(self.gene, 'TEST')
    
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
        
        # on the - strand, upstream exons are at higher coordinates. The CDS
        # starts at 35, so 36-40 are -1 to -5, then 90 is -6 and 100 is -16
        self.gene = self.construct_gene(start=10, end=100, strand='-',
            exons=[(10, 20), (30, 40), (90, 100)], cds=[(10, 20), (30, 35)])
        self.assertEqual(self.gene.get_coding_distance(40), {'pos': -5, 'offset': 0})
        self.assertEqual(self.gene.get_coding_distance(90), {'pos': -6, 'offset': 0})
        self.assertEqual(self.gene.get_coding_distance(100), {'pos': -16, 'offset': 0})
    
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
    
    def construct_cq_gene(self, strand='+'):
        ''' build a small two exon transcript for checking variant consequences
        
        Layout (+ strand coordinates):
            exon 1: 10001-10020 (5' UTR at 10001-10010, CDS at 10011-10020)
            intron: 10021-10040 (GT at 10021-10022, AG at 10039-10040)
            exon 2: 10041-10060 (CDS at 10041-10051, 3' UTR at 10052-10060)
        The 21 bp CDS is ATG GCC AAA G|GC TGG CTG TAA (M A K G W L *), with
        codon 3 split across the intron.
        '''
        seq = 'C' * 10 + 'ATGGCCAAAG' + 'GTAAGT' + 'T' * 12 + 'AG' + 'GCTGGCTGTAA' + 'G' * 9
        tx = Transcript('TEST', '1', 10001, 10060, strand,
            exons=[(10001, 10020), (10041, 10060)],
            cds=[(10011, 10020), (10041, 10051)])
        if strand == '+':
            tx.genomic_sequence = seq
        return tx
    
    def test_consequence_outside_transcript(self):
        ''' check variants outside the transcript, up to 5 kb away
        '''
        tx = self.construct_cq_gene()
        self.assertEqual(tx.consequence(10000, 'A', 'G'), 'upstream_gene_variant')
        self.assertEqual(tx.consequence(5001, 'A', 'G'), 'upstream_gene_variant')
        self.assertEqual(tx.consequence(5000, 'A', 'G'), 'intergenic_variant')
        self.assertEqual(tx.consequence(10061, 'A', 'G'), 'downstream_gene_variant')
        self.assertEqual(tx.consequence(15060, 'A', 'G'), 'downstream_gene_variant')
        self.assertEqual(tx.consequence(15061, 'A', 'G'), 'intergenic_variant')
        
        # upstream and downstream swap for transcripts on the - strand
        tx = self.construct_cq_gene(strand='-')
        self.assertEqual(tx.consequence(10000, 'A', 'G'), 'downstream_gene_variant')
        self.assertEqual(tx.consequence(10061, 'A', 'G'), 'upstream_gene_variant')
    
    def test_consequence_utr(self):
        ''' check variants in exons outside the CDS are UTR variants
        '''
        tx = self.construct_cq_gene()
        self.assertEqual(tx.consequence(10001, 'C', 'T'), '5_prime_UTR_variant')
        self.assertEqual(tx.consequence(10010, 'C', 'T'), '5_prime_UTR_variant')
        self.assertEqual(tx.consequence(10052, 'G', 'A'), '3_prime_UTR_variant')
        self.assertEqual(tx.consequence(10060, 'G', 'A'), '3_prime_UTR_variant')
        
        # on the - strand, the 5' UTR is at the higher coordinates
        tx = self.construct_cq_gene(strand='-')
        self.assertEqual(tx.consequence(10052, 'G', 'A'), '5_prime_UTR_variant')
        self.assertEqual(tx.consequence(10060, 'G', 'A'), '5_prime_UTR_variant')
        self.assertEqual(tx.consequence(10001, 'C', 'T'), '3_prime_UTR_variant')
        self.assertEqual(tx.consequence(10010, 'C', 'T'), '3_prime_UTR_variant')
    
    def test_consequence_splice_and_intron(self):
        ''' check intronic variants by distance to the nearest exon
        
        Splice donor and acceptor variants are in the first or last 2 bp of the
        intron. Donor bases 3-6 are splice donor region (with base 5 classed
        separately), 7-8 are splice region, acceptor bases 3-8 are splice
        region, and acceptor bases 3-17 are the polypyrimidine tract.
        '''
        # distance into the 20 bp intron from the donor end
        expected = {1: 'splice_donor_variant', 2: 'splice_donor_variant',
            3: 'splice_donor_region_variant', 4: 'splice_donor_region_variant',
            5: 'splice_donor_5th_base_variant', 6: 'splice_donor_region_variant',
            7: 'splice_region_variant', 8: 'splice_region_variant',
            # 9-12 are 12-9 bp from the acceptor end
            9: 'splice_polypyrimidine_tract_variant',
            12: 'splice_polypyrimidine_tract_variant',
            13: 'splice_region_variant', 18: 'splice_region_variant',
            19: 'splice_acceptor_variant', 20: 'splice_acceptor_variant'}
        tx = self.construct_cq_gene()
        for dist, cq in expected.items():
            pos = 10020 + dist
            ref = tx.genomic_sequence[pos - tx.start]
            self.assertEqual(tx.consequence(pos, ref, 'C'), cq, msg=f'pos={pos}')
        
        # on the - strand, the donor is at the 3' end of the intron (in + strand
        # coordinates), next to the exon which is first in transcript order
        tx = self.construct_cq_gene(strand='-')
        for dist, cq in expected.items():
            pos = 10041 - dist
            self.assertEqual(tx.consequence(pos, 'T', 'C'), cq, msg=f'pos={pos}')
    
    def test_consequence_coding_snvs(self):
        ''' check SNVs in the CDS of ATG GCC AAA G|GC TGG CTG TAA (M A K G W L *)
        '''
        tx = self.construct_cq_gene()
        self.assertEqual(tx.cds_sequence, 'ATGGCCAAAGGCTGGCTGTAA')
        self.assertEqual(tx.consequence(10012, 'T', 'C'), 'start_lost')         # ATG>ACG
        self.assertEqual(tx.consequence(10016, 'C', 'T'), 'synonymous_variant') # GCC>GCT
        self.assertEqual(tx.consequence(10017, 'A', 'G'), 'missense_variant')   # AAA>GAA
        
        # the codon split across the intron, GGC (G)
        self.assertEqual(tx.consequence(10020, 'G', 'T'), 'missense_variant')   # GGC>TGC
        # GGC>GGA is synonymous, but within 3 bp of the exon start, so it's a
        # splice region variant, which is more severe
        self.assertEqual(tx.consequence(10042, 'C', 'A'), 'splice_region_variant')
        
        self.assertEqual(tx.consequence(10044, 'G', 'A'), 'stop_gained')        # TGG>TAG
        self.assertEqual(tx.consequence(10045, 'G', 'A'), 'stop_gained')        # TGG>TGA
        self.assertEqual(tx.consequence(10049, 'T', 'C'), 'stop_lost')          # TAA>CAA
        self.assertEqual(tx.consequence(10050, 'A', 'G'), 'stop_retained_variant') # TAA>TGA
        self.assertEqual(tx.consequence(10051, 'A', 'G'), 'stop_retained_variant') # TAA>TAG
    
    def test_consequence_empty_ref(self):
        ''' check an empty reference allele raises an error
        '''
        tx = self.construct_cq_gene()
        with self.assertRaisesRegex(ValueError, 'reference allele is empty'):
            tx.consequence(10015, '', 'A')
    
    def test_consequence_without_sequence(self):
        ''' check coding SNVs need the genomic sequence, but other sites don't
        '''
        tx = self.construct_cq_gene('-')
        with self.assertRaisesRegex(ValueError, 'TEST lacks genomic sequence'):
            tx.consequence(10015, 'A', 'G')
        # 10 bp from the acceptor, which is at the lower end on the - strand
        self.assertEqual(tx.consequence(10030, 'A', 'G'), 'splice_polypyrimidine_tract_variant')
        self.assertEqual(tx.consequence(10005, 'A', 'G'), '3_prime_UTR_variant')
    
    def test_consequence_mnv_into_intron(self):
        ''' check multi-base substitutions running from the CDS into an intron
        disrupt the splice site, as SNVs and deletions there do
        '''
        # + strand, exon 1 ends at 10020 (last two bases GG) and the intron
        # starts with GT, so the donor is 10021-10022. Exon 2 starts at 10041,
        # and the intron ends with AG, so the acceptor is 10039-10040.
        tx = self.construct_cq_gene()
        self.assertEqual(tx.consequence(10020, 'GG', 'CC'), 'splice_donor_variant')
        self.assertEqual(tx.consequence(10019, 'AGG', 'CCC'), 'splice_donor_variant')
        self.assertEqual(tx.consequence(10040, 'GG', 'CC'), 'splice_acceptor_variant')
        self.assertEqual(tx.consequence(10040, 'GGC', 'CCC'), 'splice_acceptor_variant')
        # entirely within the exon, the splice region is less severe than missense
        self.assertEqual(tx.consequence(10018, 'AA', 'GG'), 'missense_variant')
        
        # on the - strand the donor and acceptor swap, as the donor is at the
        # higher coordinates
        plus = self.construct_cq_gene()
        tx = self.construct_cq_gene(strand='-')
        tx.genomic_sequence = tx.reverse_complement(plus.genomic_sequence)
        self.assertEqual(tx.consequence(10020, 'GG', 'CC'), 'splice_acceptor_variant')
        self.assertEqual(tx.consequence(10040, 'GG', 'CC'), 'splice_donor_variant')
        self.assertEqual(tx.consequence(10040, 'GGC', 'CCC'), 'splice_donor_variant')
    
    def test_consequence_ref_mismatch(self):
        ''' check ref alleles must match the genomic sequence, where available
        '''
        tx = self.construct_cq_gene()
        # the CDS starts with ATG at 10011-10013
        with self.assertRaisesRegex(ValueError, "doesn't match the genomic sequence"):
            tx.consequence(10011, 'G', 'A')
        with self.assertRaisesRegex(ValueError, "doesn't match the genomic sequence"):
            tx.consequence(10011, 'ATC', 'A')
        self.assertEqual(tx.consequence(10012, 't', 'c'), 'start_lost')
        
        # only the part overlapping the sequence is checked (10001 is C)
        self.assertEqual(tx.consequence(10000, 'AC', 'A'), '5_prime_UTR_variant')
        with self.assertRaisesRegex(ValueError, "doesn't match the genomic sequence"):
            tx.consequence(10000, 'AG', 'A')
        
        # without a sequence, the ref allele isn't checked
        tx = self.construct_cq_gene(strand='-')
        self.assertEqual(tx.consequence(10001, 'G', 'A'), '3_prime_UTR_variant')
    
    def test_consequence_lowercase(self):
        ''' check lowercase (e.g. soft-masked) sequence and alleles are handled
        '''
        tx = Transcript('TEST', '1', 1, 9, '+', exons=[(1, 9)], cds=[(1, 9)],
            sequence='atggcctaa')
        self.assertEqual(tx.genomic_sequence, 'ATGGCCTAA')
        self.assertEqual(tx.cds_sequence, 'ATGGCCTAA')
        # GCC (Ala) -> ACC (Thr), and GCC -> GCT (Ala)
        self.assertEqual(tx.consequence(4, 'G', 'A'), 'missense_variant')
        self.assertEqual(tx.consequence(4, 'g', 'a'), 'missense_variant')
        self.assertEqual(tx.consequence(6, 'c', 't'), 'synonymous_variant')
    
    def test_without_exons(self):
        ''' check methods needing exons raise an error, rather than crashing
        '''
        tx = Transcript('TEST', '1', 100, 200, '+')
        for func in (lambda: tx.get_coding_distance(150),
                     lambda: tx.get_codon_info(150),
                     lambda: tx.consequence(150, 'A', 'G'),
                     lambda: tx.consequence(150, 'AT', 'A')):
            with self.assertRaisesRegex(ValueError, 'no exons assigned'):
                func()
    
    def test_consequence_last_codon_not_stop(self):
        ''' check synonymous changes in the last codon are only stop_retained
        if that codon is a stop codon, e.g. not for a 3' incomplete CDS
        '''
        tx = Transcript('TEST', '1', 1, 20, '+', exons=[(1, 20)], cds=[(4, 15)])
        tx.genomic_sequence = 'CCCATGGCCAAAGCCGGGGG'  # ATG GCC AAA GCC, no stop
        self.assertEqual(tx.consequence(15, 'C', 'T'), 'synonymous_variant')
    
    def test_consequence_start_retained(self):
        ''' check synonymous changes in a (non-AUG) start codon retain the start
        '''
        # CDS is CTG GCC TGG TAA (L A W *) at positions 4-15
        tx = Transcript('TEST', '1', 1, 20, '+', exons=[(1, 20)], cds=[(4, 15)],
            sequence='CCCCTGGCCTGGTAAGGGGG')
        self.assertEqual(tx.consequence(6, 'G', 'A'), 'start_retained_variant') # CTG>CTA
        self.assertEqual(tx.consequence(5, 'T', 'A'), 'start_lost')             # CTG>CAG
    
    def test_consequence_short_intron(self):
        ''' check variants in very short introns (< 10 bp) count as coding
        
        Introns this short are likely annotation artefacts, rather than spliced
        out, so variants in them aren't classed as intronic or splice variants.
        '''
        tx = Transcript('TEST', '1', 1, 40, '+', exons=[(1, 20), (26, 40)],
            cds=[(5, 20), (26, 30)])
        for pos in range(21, 26):
            self.assertEqual(tx.consequence(pos, 'A', 'G'), 'coding_sequence_variant',
                msg=f'pos={pos}')
    
    def test_consequence_deletion_across_splice_site(self):
        ''' check deletions spanning an exon boundary are classed by the splice
        site they disrupt, whether or not they remove an inframe length of exon
        '''
        # deletes 10018-10020 (end of exon 1) and 10021-10022 (intron)
        exon_end = (10017, 'AAAGGT', 'A')
        # deletes 10039-10040 (intron) and 10041-10043 (start of exon 2)
        exon_start = (10038, 'TAGGCT', 'T')
        
        tx = self.construct_cq_gene('+')
        self.assertEqual(tx.consequence(*exon_end), 'splice_donor_variant')
        self.assertEqual(tx.consequence(*exon_start), 'splice_acceptor_variant')
        
        tx = self.construct_cq_gene('-')
        self.assertEqual(tx.consequence(*exon_end), 'splice_acceptor_variant')
        self.assertEqual(tx.consequence(*exon_start), 'splice_donor_variant')
        
        # removing a non-inframe length of exon is still a splice variant, since
        # that is more severe than a frameshift
        tx = self.construct_cq_gene('+')
        # deletes 10039-10040 (intron) and 10041-10042 (start of exon 2)
        self.assertEqual(tx.consequence(10038, 'TAGGC', 'T'), 'splice_acceptor_variant')
        # deletes 10019-10020 (end of exon 1) and 10021 (intron)
        self.assertEqual(tx.consequence(10018, 'AAGG', 'A'), 'splice_donor_variant')
    
    def test_consequence_inframe_indel_codons(self):
        ''' check inframe indels which gain a stop, or lose the start or stop
        codon, using the + strand CDS ATG GCC AAA G|GC TGG CTG TAA, with
        3' UTR GGGGGGGGG.
        '''
        tx = self.construct_cq_gene('+')
        # deletes TGG (10012-10014): ATG GCC -> ACC
        self.assertEqual(tx.consequence(10011, 'ATGG', 'A'), 'start_lost')
        # deletes CAT (10010-10012), running from the 5' UTR into the start codon
        self.assertEqual(tx.consequence(10009, 'CCAT', 'C'), 'start_lost')
        # deletes CA (10010-10011), which isn't a frameshift, as it isn't
        # wholly within the CDS
        self.assertEqual(tx.consequence(10009, 'CCA', 'C'), 'start_lost')
        # deletes the whole GCC codon (10014-10016)
        self.assertEqual(tx.consequence(10013, 'GGCC', 'G'), 'inframe_deletion')
        # inserts TAA after GCC: GCC -> GCC TAA
        self.assertEqual(tx.consequence(10016, 'C', 'CTAA'), 'stop_gained')
        # inserts CCC after GCC
        self.assertEqual(tx.consequence(10016, 'C', 'CCCC'), 'inframe_insertion')
        # deletes GTA (10048-10050): CTG TAA -> CTA
        self.assertEqual(tx.consequence(10047, 'TGTA', 'T'), 'stop_lost')
        # deletes the whole stop codon (10049-10051)
        self.assertEqual(tx.consequence(10048, 'GTAA', 'G'), 'stop_lost')
        # deletes AAG (10050-10052), into the 3' UTR: TAA GGG -> TGG
        self.assertEqual(tx.consequence(10049, 'TAAG', 'T'), 'stop_lost')
        # deletes AGG (10051-10053), into the 3' UTR: TAA GGG -> TAG, still a stop
        self.assertEqual(tx.consequence(10050, 'AAGG', 'A'), 'stop_retained_variant')
        
        # without sequence, inframe indels can't be checked further
        tx = self.construct_cq_gene('-')
        self.assertEqual(tx.consequence(10048, 'GTAA', 'G'), 'inframe_deletion')
        
        # - strand transcript at 1-30, CDS at 4-18. In transcript orientation
        # this is 12 bp of 5' UTR, ATG GCC AAA GCC TAA, then GGG of 3' UTR.
        tx = Transcript('TEST', '1', 1, 30, '-', exons=[(1, 30)], cds=[(4, 18)])
        tx.genomic_sequence = 'C' * 12 + 'ATGGCCAAAGCCTAA' + 'GGG'
        # + strand: CCC TTAGGCTTTGGCCAT GGGGGGGGGGGG
        self.assertEqual(tx.cds_sequence, 'ATGGCCAAAGCCTAA')
        # inserts TTA after genomic 12, i.e. TAA before AAA in transcript
        # orientation: AAA -> TAA AAA
        self.assertEqual(tx.consequence(12, 'T', 'TTTA'), 'stop_gained')
        # deletes the stop codon (genomic 4-6, TTA)
        self.assertEqual(tx.consequence(3, 'CTTA', 'C'), 'stop_lost')
        # deletes genomic 3-5 (CTT), i.e. AAG in transcript orientation, into
        # the 3' UTR: TAA GGG -> TGG
        self.assertEqual(tx.consequence(2, 'CCTT', 'C'), 'stop_lost')
        # deletes genomic 15-17 (CCA), i.e. TGG in transcript orientation
        self.assertEqual(tx.consequence(14, 'GCCA', 'G'), 'start_lost')
        # deletes genomic 18-19 (TG), i.e. the A of ATG and a 5' UTR base
        self.assertEqual(tx.consequence(17, 'ATG', 'A'), 'start_lost')
        # inserts CCC before the start codon, which leaves it intact
        self.assertEqual(tx.consequence(18, 'T', 'TGGG'), '5_prime_UTR_variant')
    
    def test_consequence_start_codon_indels(self):
        ''' check indels in the start codon are start_lost if the altered
        peptide doesn't start or end with the reference peptide, or (for
        deletions) the CDS no longer begins with ATG, as in VEP
        '''
        # CDS ATG ATG AAA GCC TAA (M M K A *) at 4-18
        tx = Transcript('TEST', '1', 1, 20, '+', exons=[(1, 20)], cds=[(4, 18)])
        tx.genomic_sequence = 'CCC' + 'ATGATGAAAGCCTAA' + 'GG'
        # deletes TGA (5-7): ATG ATG -> ATG, i.e. MM -> M
        self.assertEqual(tx.consequence(4, 'ATGA', 'A'), 'start_lost')
        # deletes the whole start codon (4-6), so the next ATG takes its place
        self.assertEqual(tx.consequence(3, 'CATG', 'C'), 'inframe_deletion')
        # deleting the whole start codon of ATG GCC leaves GCC at the start
        tx = self.construct_cq_gene('+')
        self.assertEqual(tx.consequence(10010, 'CATG', 'C'), 'start_lost')
        
        tx = Transcript('TEST', '1', 1, 20, '+', exons=[(1, 20)], cds=[(4, 18)])
        tx.genomic_sequence = 'CCC' + 'ATGATGAAAGCCTAA' + 'GG'
        # inserts CGA within the start codon: ATG -> ACG ATG, i.e. M -> TM
        self.assertEqual(tx.consequence(4, 'A', 'ACGA'), 'inframe_insertion')
    
    def test_consequence_insertion_codons(self):
        ''' check insertions in the + strand CDS ATG GCC AAA G|GC TGG CTG TAA
        '''
        tx = self.construct_cq_gene('+')
        # inserts TAC within GCC: GCC -> GTA CCC (V P), which loses the A
        self.assertEqual(tx.consequence(10014, 'G', 'GTAC'), 'protein_altering_variant')
        # inserts CCG within GCC: GCC -> GCC GCC (A A)
        self.assertEqual(tx.consequence(10014, 'G', 'GCCG'), 'inframe_insertion')
        # inserts A within TGG: TGG -> TGA G, a frameshift creating a stop
        self.assertEqual(tx.consequence(10044, 'G', 'GA'), 'stop_gained')
        # inserts T between GCC and AAA. The neighbouring codons are left
        # intact, so this doesn't create a stop, despite T AAA -> TAA A
        self.assertEqual(tx.consequence(10016, 'C', 'CT'), 'frameshift_variant')
        # deletes the first T of TAA (10049): TAA G -> AAG
        self.assertEqual(tx.consequence(10048, 'GT', 'G'), 'stop_lost')
    
    def test_consequence_insertion_splice_sites(self):
        ''' check insertions near exon/intron boundaries
        
        Layout: exon 1 at 1001-1020 (all 5' UTR on + strand), intron at
        1021-1040, exon 2 at 1041-1100 with CDS at 1050-1070.
        '''
        tx = Transcript('TEST', '1', 1001, 1100, '+',
            exons=[(1001, 1020), (1041, 1100)], cds=[(1050, 1070)])
        # insertions are classed by the base before the inserted bases
        expected = {1017: '5_prime_UTR_variant',
            # within the 3 exonic bases next to the intron
            1018: 'splice_region_variant',
            # at the exon/intron junction, which counts as exonic
            1020: 'splice_region_variant',
            # 1-7 bp into the intron, by the closer flanking base
            1021: 'splice_donor_variant', 1022: 'splice_region_variant',
            1023: 'splice_donor_region_variant', 1025: 'splice_donor_region_variant',
            1026: 'splice_region_variant', 1027: 'splice_region_variant',
            # 8+ bp from the donor, and 8-16 bp from the acceptor
            1028: 'splice_polypyrimidine_tract_variant',
            1033: 'splice_region_variant', 1038: 'splice_region_variant',
            1039: 'splice_acceptor_variant',
            1040: 'splice_region_variant', 1043: '5_prime_UTR_variant'}
        for pos, cq in expected.items():
            self.assertEqual(tx.consequence(pos, 'A', 'AT'), cq, msg=f'pos={pos}')
        
        # insertions at the junction of a coding exon follow the CDS change. The
        # CDS is ATG GCC AAA G|GC TGG CTG TAA, and both insert A before GC,
        # i.e. GGC -> GAG C
        tx = self.construct_cq_gene('+')
        self.assertEqual(tx.consequence(10020, 'G', 'GA'), 'frameshift_variant')
        self.assertEqual(tx.consequence(10040, 'G', 'GA'), 'frameshift_variant')
    
    def test_consequence_stop_flank_across_exons(self):
        ''' check frameshifts in a stop codon at the end of an exon use the
        3' UTR from the next exon
        '''
        # CDS ATG GCC AAA TAA at 4-15, then a 14 bp intron, and exon 2
        for exon2, cq in (('G' * 11, 'splice_region_variant'), ('C' * 11, 'stop_lost')):
            tx = Transcript('TEST', '1', 1, 40, '+', exons=[(1, 15), (30, 40)],
                cds=[(4, 15)])
            tx.genomic_sequence = 'CCC' + 'ATGGCCAAATAA' + 'GTAAGTTTTTTTAG' + exon2
            # deletes the A at 14: TAA -> TA + G, i.e. stop retained (but in
            # the splice region), or TA + C, which loses the stop
            self.assertEqual(tx.consequence(13, 'TA', 'T'), cq, msg=f'exon2={exon2}')
    
    def test_consequence_splice_site_outside_cds(self):
        ''' check variants spanning an exon/intron boundary outside the CDS are
        classed by the splice site they disrupt
        
        Layout: exon 1 at 1001-1020 (all 5' UTR on + strand), intron at
        1021-1040, exon 2 at 1041-1100 with CDS at 1050-1070.
        '''
        # deletes 1019-1020 (end of exon 1) and 1021-1022 (intron)
        exon_end = (1018, 'AAAGT', 'A')
        # deletes 1038-1040 (intron) and 1041 (start of exon 2, before the CDS)
        exon_start = (1037, 'TTAGG', 'T')
        
        for cds in ([(1050, 1070)], []):
            tx = Transcript('TEST', '1', 1001, 1100, '+',
                exons=[(1001, 1020), (1041, 1100)], cds=cds)
            self.assertEqual(tx.consequence(*exon_end), 'splice_donor_variant')
            self.assertEqual(tx.consequence(*exon_start), 'splice_acceptor_variant')
            
            tx = Transcript('TEST', '1', 1001, 1100, '-',
                exons=[(1001, 1020), (1041, 1100)], cds=cds)
            self.assertEqual(tx.consequence(*exon_end), 'splice_acceptor_variant')
            self.assertEqual(tx.consequence(*exon_start), 'splice_donor_variant')
        
        # deletions crossing the transcript start have no splice site to disrupt
        tx = Transcript('TEST', '1', 1001, 1100, '+',
            exons=[(1001, 1020), (1041, 1100)], cds=[(1050, 1070)])
        self.assertEqual(tx.consequence(998, 'AAAAC', 'A'), '5_prime_UTR_variant')
        
        # deletions within the exon, away from its ends, are unaffected
        self.assertEqual(tx.consequence(1005, 'AAAA', 'A'), '5_prime_UTR_variant')
    
    def test_consequence_exonic_splice_region(self):
        ''' check low impact variants in the 3 exonic bases next to an intron
        are splice region variants
        
        Layout: exon 1 at 1001-1020 (all 5' UTR on + strand), intron at
        1021-1040, exon 2 at 1041-1100 with CDS at 1050-1070.
        '''
        for strand in '+-':
            tx = Transcript('TEST', '1', 1001, 1100, strand,
                exons=[(1001, 1020), (1041, 1100)], cds=[(1050, 1070)])
            utr = '5_prime_UTR_variant' if strand == '+' else '3_prime_UTR_variant'
            self.assertEqual(tx.consequence(1017, 'A', 'G'), utr)
            for pos in (1018, 1020, 1041, 1043):
                self.assertEqual(tx.consequence(pos, 'A', 'G'), 'splice_region_variant',
                    msg=f'strand={strand}, pos={pos}')
            self.assertEqual(tx.consequence(1044, 'A', 'G'), utr)
            # the transcript start isn't next to an intron
            self.assertEqual(tx.consequence(1001, 'A', 'G'), utr)
        
        tx = Transcript('TEST', '1', 1001, 1100, '+',
            exons=[(1001, 1020), (1041, 1100)])
        self.assertEqual(tx.consequence(1010, 'A', 'G'), 'non_coding_transcript_exon_variant')
        self.assertEqual(tx.consequence(1020, 'A', 'G'), 'splice_region_variant')
        
        # exonic bases next to introns under 10 bp aren't splice region
        tx = Transcript('TEST', '1', 1, 40, '+', exons=[(1, 20), (26, 40)])
        self.assertEqual(tx.consequence(19, 'A', 'G'), 'non_coding_transcript_exon_variant')
    
    def test_consequence_short_intron_outside_cds(self):
        ''' check variants in very short introns follow the surrounding exons,
        when the intron is in the UTR, or the transcript is non-coding
        '''
        for strand, utr in (('+', '5_prime_UTR_variant'), ('-', '3_prime_UTR_variant')):
            tx = Transcript('TEST', '1', 1, 40, strand, exons=[(1, 20), (26, 40)],
                cds=[(30, 38)])
            for pos in range(21, 26):
                self.assertEqual(tx.consequence(pos, 'A', 'G'), utr, msg=f'pos={pos}')
        
        tx = Transcript('TEST', '1', 1, 40, '+', exons=[(1, 20), (26, 40)])
        for pos in range(21, 26):
            self.assertEqual(tx.consequence(pos, 'A', 'G'),
                'non_coding_transcript_exon_variant', msg=f'pos={pos}')
    
    def test_consequence_indel_frame(self):
        """ check indels are classed as inframe or frameshift by length change
        """
        self.assertEqual(self.gene.consequence(1150, 'G', 'GCCC'), 'inframe_insertion')
        self.assertEqual(self.gene.consequence(1150, 'G', 'GC'), 'frameshift_variant')
        self.assertEqual(self.gene.consequence(1150, 'G', 'GCC'), 'frameshift_variant')
        self.assertEqual(self.gene.consequence(1150, 'GCCC', 'G'), 'inframe_deletion')
        self.assertEqual(self.gene.consequence(1150, 'GC', 'G'), 'frameshift_variant')
