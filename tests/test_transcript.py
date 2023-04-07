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
        
        self.assertEqual(a.merge_coordinates(exons1, exons2),
            a.merge_coordinates(exons2, exons1))
    
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
    
