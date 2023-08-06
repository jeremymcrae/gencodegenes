# cython: language_level=3, boundscheck=False
'''
Copyright (c) 2015 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

from itertools import combinations

cdef class Transcript:
    def __cinit__(self, name, chrom, start, end, strand, 
            transcript_type='protein_coding', exons=None, cds=None, sequence=None, 
            offset=0):
        ''' construct a Transcript object
        
        Args:
            name: ID of the transcript
            start: position in bp at 5' edge of transcript (on + strand)
            end: position in bp at 3' edge of transcript (on + strand)
            exons: list of tuples defining start and end positions of exons
            cds: list of tuples defining start and end positions of CDS regions
            sequence: DNA sequence of genome region of the transcript.
            offset: how many base pairs the DNA sequence extends outwards
        '''
        
        name = name.encode('utf8')
        chrom = chrom.encode('utf8')
        transcript_type = transcript_type.encode('utf8')
        self.thisptr = new Tx(name, chrom, start, end, ord(strand), transcript_type)
        
        if exons is not None and cds is not None:
            self.exons = exons
            self.cds = cds
            
            self.genomic_offset = offset
            if sequence is not None:
                self.genomic_sequence = sequence
    
    def __dealloc__(self):
        del self.thisptr
    
    def __repr__(self):
        
        exons = [ (x['start'], x['end']) for x in self.exons ]
        cds = [ (x['start'], x['end']) for x in self.cds ]
        seq = self.genomic_sequence
        
        if len(seq) > 40:
            seq = seq[:20] + '...[{} bp]...'.format(len(seq) - 40) + seq[-20:]
        
        if exons == []:
            exons = None
        
        if cds == []:
            cds = None
        
        if seq == '':
            seq = None
        else:
            seq = '"' + seq + '"'
        
        return f'Transcript(name="{self.name}", chrom="{self.chrom}", ' \
            f'start={self.start}, end={self.end}, strand="{self.strand}", ' \
            f'transcript_type="{self.type}", exons={exons}, cds={cds}, ' \
            f'sequence={seq}, offset={self.genomic_offset})'
    
    def __str__(self):
        return self.__repr__()
    
    def __hash__(self):
        return hash((self.chrom, self.start, self.end))
    
    def __eq__(self, other):
        ''' check if transcripts occupy exact same genomic region '''
        return self.__hash__() == other.__hash__()
    
    def _get_overlaps(self, exon, regions):
        ''' find all regions which overlap a given region
        '''
        return [ i for i, x in enumerate(regions) if
            exon['start'] <= x['end'] and exon['end'] >= x['start'] ]
    
    def _insert_region(self, coords, region):
        ''' include a region into a list of regions
        
        To include a region, we have to check which pre-existing regions the new
        region overlaps, so any overlaps can be merged into a single region.
        
        Args:
            coords: list of {start: X, end: Y} dictionaries
            region: dict of {'start': X, 'end': Y} positions
        '''
        indices = self._get_overlaps(region, coords)
        overlaps = [ coords[i] for i in indices ]
        start = min( x['start'] for x in overlaps + [region] )
        end = max( x['end'] for x in overlaps + [region] )
        
        for i in sorted(indices, reverse=True):
            del coords[i]
        
        return coords + [{'start': start, 'end': end}]
    
    def _merge_coordinates(self, first, second):
        ''' merge two sets of coordinates, to get the union of regions
        
        This uses an inefficient approach, looping over and over, but we won't
        need to perform this often.
        
        Args:
            first: list of {'start': x, 'end': y} dictionaries for first transcript
            second: list of {'start': x, 'end': y} dictionaries for second transcript
        
        Returns:
            list of [start, end] lists, sorted by position.
        '''
        coords = []
        for a, b in combinations(first + second, 2):
            if a['start'] <= b['end'] and a['end'] >= b['start']:
                region = {'start': min(a['start'], b['start']),
                    'end': max(a['end'], b['end'])}
                a, b = region, region
            
            coords = self._insert_region(coords, a)
            coords = self._insert_region(coords, b)
        
        return [ (x['start'], x['end']) for x in sorted(coords, key=lambda x: x['start']) ]
    
    def _merge_genomic_seq(self, other):
        ''' merge the genomic sequence from two transcripts
        
        We have two transcripts A, and B. We need to get the contiguous sequence
        from the start of the first transcript on the chromosome to the end of
        the second transcript. The transcripts may or may not overlap. There are
        three scenarios we need to account for:
        
        
        overlap without   A   =================
        enveloping                           ================= B
        
        
            overlap       A   ==============================
        and envelop               ==================  B
        
        
        no overlap        A  ===============
                                               ===============  B
        
        I've called the transcript whose sequence is first along the chromosome
        as 'lead', and the transcript whose sequence is last as 'lag', and the
        converse as 'not_lead', and 'not_lag'. Note that in the envelope case,
        the lead transcript can also be the lag transcript.
        '''
        
        # make sure that the surrounding sequence is the same length in both
        # transcripts.
        # TODO: this could be worked around, by figuring the minimum offset length,
        # TODO: then trimming the respective DNA offset sequences to that length.
        assert self.genomic_offset == other.genomic_offset
        
        # figure out which transcripts hold the leading and lagging sections
        lead, not_lead = self, other
        if self.start > other.start:
            lead, not_lead = other, self
        
        lag, not_lag = self, other
        if self.end < other.end:
            lag, not_lag = other, self
        
        lead_offset = lead.genomic_offset
        lead_gdna = lead.genomic_sequence
        initial = lead_gdna[:not_lead.start - lead.start + lead_offset]
        
        if self.start <= other.end and self.end >= other.start:
            intersect_start = not_lead.start - lead.start + lead_offset
            intersect_end = not_lag.end - lead.start + lead_offset
            intersect = lead_gdna[intersect_start:intersect_end]
        else:
            intersect = 'N' * (lag.start - lead.end - lead_offset * 2)
        
        lag_offset = lag.genomic_offset
        lag_gdna = lag.genomic_sequence
        
        # some transcripts overlap, but some do not. We need to find the position
        # where the lagging transcript takes over, which is either at the end of
        # not lagging transcript, or the start of the lagging transcript,
        # whichever is higher
        pos = max(not_lag.end, lag.start)
        final = lag_gdna[pos - lag.start + lead_offset:]
        
        return initial + intersect + final
    
    def __add__(self, other):
        """ combine the coding sequences of two Transcript objects
        
        When we determine the sites for sampling, occasioally we want to
        use sites from multiple alternative transcripts. We determine the sites
        for each transcript in turn, but mask the sites that have been collected
        in the preceeding transcripts. In order to be able to mask all previous
        trabnscripts, we need to combine the coding sequence of the transcripts
        as we move through them. This function performs the union of coding
        sequence regions between different transcripts.
        
        We do this outside of the c++ class, so as to be able to set up a
        Transcript object correctly.
        
        Args:
            other: a transcript to be combined with the current object.
        
        Returns:
            an altered instance of the class, where the coding sequence regions
            are the union of the coding regions of two Transcript objects. This
            disrupts the ability to get meaningingful sequence from the object,
            so don't try to extract sequence from the returned object.
        """
        
        # if we try transcript + None or None + transcript, return the original
        # transcript, rather than raising an error.
        if other is None:
            return self
        if self is None:
            return other
        
        altered = Transcript('{}:{}'.format(self.name, other.name),
            self.chrom, min(self.start, other.start),
            max(self.end, other.end), self.strand, self.type)
        
        exons = self._merge_coordinates(self.exons, other.exons)
        cds = self._merge_coordinates(self.cds, other.cds)
        
        altered.exons = exons
        altered.cds = cds
        
        if self.genomic_sequence != "":
            altered.genomic_offset = self.genomic_offset
            altered.genomic_sequence = self._merge_genomic_seq(other)
        
        return altered
    
    def __radd__(self, other):
        return self.__add__(other)
    
    @property
    def name(self):
        '''transcript ID'''
        return self.thisptr.get_name().decode('utf8')
    @property
    def chrom(self):
        '''chromosome the gene is on'''
        return self.thisptr.get_chrom().decode('utf8')
    @property
    def type(self):
        '''transcript functional type (protein_coding etc.)'''
        return self.thisptr.get_type().decode('utf8')
    @property
    def start(self):
        '''transcript start position (TSS) on chromosome'''
        return self.thisptr.get_start()
    @property
    def end(self):
        '''transcript end position on chromosome'''
        return self.thisptr.get_end()
    @property
    def strand(self):
        '''strand the transcript is on (+ or -)'''
        return chr(self.thisptr.get_strand())
    @property
    def cds_start(self):
        '''CDS start position on chromosome'''
        return self.thisptr.get_cds_start()
    @property
    def cds_end(self):
        '''CDS end position on chromosome'''
        return self.thisptr.get_cds_end()
    
    @property
    def exons(self):
        ''' list of exon coords as {'start': int, 'end': int} dicts '''
        return self.thisptr.get_exons()
    @exons.setter
    def exons(self, exon_ranges):
        self.thisptr.set_exons(exon_ranges)

    @property
    def cds(self):
        '''list of CDS coords as {'start': int, 'end': int} dicts '''
        return self.thisptr.get_cds()
    @cds.setter
    def cds(self, cds_ranges):
        self.thisptr.set_cds(cds_ranges)

    @property
    def cds_sequence(self):
        ''' CDS sequence for the transcript
        '''
        return self.thisptr.get_cds_sequence().decode('utf8')
    
    @cds_sequence.setter
    def cds_sequence(self, text):
        self.thisptr.add_cds_sequence(text.encode('utf8'))

    @property
    def genomic_offset(self):
        '''distance the DNA sequence extends symmetrically beyond gene boundaries'''
        return self.thisptr.get_genomic_offset()
    
    @genomic_offset.setter
    def genomic_offset(self, int offset):
        self.thisptr.set_genomic_offset(offset)

    @property
    def genomic_sequence(self):
        ''' genomic sequence spanning the transcript 
        '''
        return self.thisptr.get_genomic_sequence().decode('utf8')
    
    @genomic_sequence.setter
    def genomic_sequence(self, str text):
        self.thisptr.add_genomic_sequence(text.encode('utf8'))

    def _fix_cds_boundary(self, pos):
        ''' adjust CDS boundary

        Args:
            pos: nucleotide position on chromosome
        '''
        return self.thisptr.fix_cds_boundary(pos)
    
    def in_exons(self, position):
        ''' check if a site lies within the exon ranges
        
        Args:
            position: an integer-based chromosome position.
        '''
        
        return self.thisptr.is_exonic(position)
    
    def get_closest_exon(self, pos):
        ''' finds the the exon closest to a site

        Args:
            pos: nucleotide position on chromosome
        '''
        return self.thisptr.get_closest_exon(pos)
    
    def in_coding_region(self, pos):
        ''' determine if a genomic position lies within the coding region

        Args:
            pos: nucleotide position on chromosome
        '''
        return self.thisptr.in_coding_region(pos)
    
    def get_coding_distance(self, pos):
        ''' get distance to CDS start (and intronic offset)

        Args:
            pos: nucleotide position on chromosome
        '''
        coords = self.thisptr.get_coding_distance(pos)
        
        return {'pos': coords.position, 'offset': coords.offset}
    
    def get_position_on_chrom(self, pos, offset=0):
        ''' convert CDS coordinate to position on chromosome

        Args:
            pos: number of bases to CDS position
            offset: if intronic, number of bases into the intron
        '''
        return self.thisptr.get_position_on_chrom(pos, offset)
    
    def get_codon_number_for_cds_position(self, cds_pos):
        ''' find which codon number a CDS position is at

        Args:
            cds_pos: CDS distance in bp to CDS start
        '''
        return self.thisptr.get_codon_number_for_cds_position(cds_pos)
    
    def get_position_within_codon(self, cds_pos):
        ''' find the position within 

        Args:
            cds_pos: CDS distance in bp to CDS start
        '''
        return self.thisptr.get_position_within_codon(cds_pos)
    

    def reverse_complement(self, text):
        ''' reverse complement a sequence

        Args:
            text: DNA sequence
        '''
        return self.thisptr.reverse_complement(text).decode('utf8')
    
    def get_centered_sequence(self, pos, length=3):
        ''' get DNA sequence around a chromosome position

        Args:
            pos: chromosome position
            length: number of bases to include
        '''
        return self.thisptr.get_centered_sequence(pos, length).decode('utf8')
    
    def get_codon_sequence(self, codon_num):
        ''' get the cDNA sequence for a codon

        Args:
            codon_num: number of codon to extract
        '''
        return self.thisptr.get_codon_sequence(codon_num).decode('utf8')
    
    def translate(self, text):
        ''' translate DNA sequence to amino acid sequence

        Args:
            text: DNA/mRNA sequence
        '''
        return self.thisptr.translate(text.encode('utf8')).decode('utf8')
    
    def get_codon_info(self, pos):
        ''' find information about the codon that a position occurs at

        Args:
            pos: nucleotide position on chromosome
        '''
        codon = dict(self.thisptr.get_codon_info(pos))
        
        if codon['codon_number'] == -9999999:
            codon['codon_number'] = None
            codon['intra_codon'] = None
            codon['codon_seq'] = None
            codon['initial_aa'] = None
        
        if codon['codon_seq'] is not None:
            codon['codon_seq'] = codon['codon_seq'].decode('utf8')
        
        if codon['initial_aa'] is not None:
            codon['initial_aa'] = chr(codon['initial_aa'])
        
        return codon
    
    def get_boundary_distance(self, pos):
        ''' get the distance in bp between a site and the nearest exon boundary

        Args:
            pos: nucleotide position on chromosome
        '''
        return self.thisptr.get_boundary_distance(pos)
    
    def consequence(self, pos, ref, alt):
        ''' get consequence of variant on transcript

        Args:
            pos: nucleotide position of variant on chromosome
            ref: reference allele DNA sequence
            alt: alternate allele DNA sequence
        '''
        cq = self.thisptr.consequence(pos, ref.encode('utf8'), alt.encode('utf8'))
        return cq.decode('utf8')
