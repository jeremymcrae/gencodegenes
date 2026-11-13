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

from libcpp.map cimport map
from libcpp.string cimport string

# passed as the name by _wrap_tx, which assigns thisptr itself
cdef object _WRAP = object()

cdef class Transcript:
    def __cinit__(self, name=None, chrom=None, start=0, end=0, strand='+', 
            transcript_type='protein_coding', exons=None, cds=None, sequence=None, 
            offset=0, attributes=None):
        ''' construct a Transcript object
        
        Args:
            name: ID of the transcript
            start: position in bp at 5' edge of transcript (on + strand)
            end: position in bp at 3' edge of transcript (on + strand)
            exons: list of tuples defining start and end positions of exons
            cds: list of tuples defining start and end positions of CDS regions
            sequence: DNA sequence of genome region of the transcript, in the
                direction of transcription (so reverse complemented on the - strand)
            offset: how many base pairs the DNA sequence extends outwards
            attributes: dict of key/value pairs from the GTF attributes field
        '''
        
        if name is _WRAP:
            return
        if name is None:
            raise ValueError('transcript needs a name')
        if chrom is None:
            raise ValueError(f'transcript {name} needs a chromosome')

        name = name.encode('utf8')
        chrom = chrom.encode('utf8')
        transcript_type = transcript_type.encode('utf8')
        
        cdef map[string, string] _attributes
        if attributes is not None:
            for key, value in attributes.items():
                _attributes[key.encode('utf8')] = value.encode('utf8')
        self.thisptr = new Tx(name, chrom, start, end, ord(strand), transcript_type, _attributes)
        
        if exons is not None:
            self.exons = exons
        if cds is not None:
            self.cds = cds
        
        self.genomic_offset = offset
        if sequence is not None:
            self.genomic_sequence = sequence
    
    def __dealloc__(self):
        del self.thisptr
    
    def __repr__(self):
        exons = [(x['start'], x['end']) for x in self.exons] or None
        cds = [(x['start'], x['end']) for x in self.cds] or None
        seq = self.genomic_sequence
        if len(seq) > 40:
            seq = f'{seq[:20]}...[{len(seq) - 40} bp]...{seq[-20:]}'
        seq = f'"{seq}"' if seq else None
        
        return f'Transcript(name="{self.name}", chrom="{self.chrom}", ' \
            f'start={self.start}, end={self.end}, strand="{self.strand}", ' \
            f'transcript_type="{self.type}", exons={exons}, cds={cds}, ' \
            f'sequence={seq}, offset={self.genomic_offset})'
    
    def _key(self):
        return (self.name, self.chrom, self.start, self.end, self.strand)
    
    def __hash__(self):
        return hash(self._key())
    
    def __eq__(self, other):
        ''' check if transcripts have the same ID, position, exons and CDS '''
        if not isinstance(other, Transcript):
            return NotImplemented
        return self._key() == other._key() and self.exons == other.exons \
            and self.cds == other.cds
    
    @staticmethod
    def _merge_coordinates(first, second):
        ''' merge two sets of coordinates, to get the union of regions
        
        Args:
            first: list of {'start': x, 'end': y} dictionaries for first transcript
            second: list of {'start': x, 'end': y} dictionaries for second transcript
        
        Returns:
            list of (start, end) tuples, sorted by position.
        '''
        merged = []
        for x in sorted(first + second, key=lambda x: x['start']):
            if merged and x['start'] <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], x['end']))
            else:
                merged.append((x['start'], x['end']))
        return merged
    
    def _merge_genomic_seq(self, other):
        ''' merge the genomic sequences of two transcripts (on the + strand), to
        cover both transcripts and their flanks. Any gap between the sequences is
        filled with Ns. Where the sequences overlap, the bases come from the
        transcript starting first.
        '''
        if self.genomic_offset != other.genomic_offset:
            raise ValueError(f'cannot merge sequences with different offsets: '
                f'{self.name} ({self.genomic_offset}) and {other.name} ({other.genomic_offset})')
        offset = self.genomic_offset
        first, second = sorted([self, other], key=lambda x: x.start)
        
        # bases between the end of the first sequence and the start of the
        # second, which is negative if the sequences overlap
        gap = (second.start - offset) - (first.end + offset) - 1
        if gap >= 0:
            return first.genomic_sequence + 'N' * gap + second.genomic_sequence
        return first.genomic_sequence + second.genomic_sequence[-gap:]
    
    def __add__(self, other):
        """ combine the coding sequences of two Transcript objects
        
        When we determine the sites for sampling, occasionally we want to
        use sites from multiple alternative transcripts. We determine the sites
        for each transcript in turn, but mask the sites that have been collected
        in the preceding transcripts. In order to be able to mask all previous
        transcripts, we need to combine the coding sequence of the transcripts
        as we move through them. This function performs the union of coding
        sequence regions between different transcripts.
        
        We do this outside of the c++ class, so as to be able to set up a
        Transcript object correctly.
        
        Args:
            other: a transcript to be combined with the current object.
        
        Returns:
            an altered instance of the class, where the coding sequence regions
            are the union of the coding regions of two Transcript objects. This
            disrupts the ability to get meaningful sequence from the object,
            so don't try to extract sequence from the returned object.
        """
        
        # if we try transcript + None or None + transcript, return the original
        # transcript, rather than raising an error.
        if other is None:
            return self
        if not isinstance(other, Transcript):
            return NotImplemented
        if self.chrom != other.chrom or self.strand != other.strand:
            raise ValueError(f'cannot combine transcripts on different chromosomes or '
                f'strands: {self.name} ({self.chrom}{self.strand}) and '
                f'{other.name} ({other.chrom}{other.strand})')
        
        altered = Transcript('{}:{}'.format(self.name, other.name),
            self.chrom, min(self.start, other.start),
            max(self.end, other.end), self.strand, self.type)
        
        exons = self._merge_coordinates(self.exons, other.exons)
        cds = self._merge_coordinates(self.cds, other.cds)
        
        altered.exons = exons
        altered.cds = cds
        
        if self.genomic_sequence != "" and other.genomic_sequence != "":
            altered.genomic_offset = self.genomic_offset
            # the merged sequence is already on the + strand
            seq = self._merge_genomic_seq(other)
            altered.thisptr.add_genomic_sequence(seq.encode('utf8'), True)
        
        return altered
    
    def __radd__(self, other):
        # handles None + transcript, e.g. when building up a merged transcript
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
        '''lowest position of the transcript on the chromosome'''
        return self.thisptr.get_start()
    @property
    def end(self):
        '''highest position of the transcript on the chromosome'''
        return self.thisptr.get_end()
    @property
    def tss(self):
        '''transcription start site on the chromosome (end on the - strand)'''
        if self.thisptr.get_strand() == b'-':
            return self.thisptr.get_end()
        return self.thisptr.get_start()
    @property
    def strand(self):
        '''strand the transcript is on (+ or -)'''
        return chr(self.thisptr.get_strand())
    @property
    def cds_start(self):
        '''first CDS position in the direction of transcription, so the
        highest CDS position on the - strand'''
        return self.thisptr.get_cds_start()
    @property
    def cds_end(self):
        '''last CDS position in the direction of transcription, so the
        lowest CDS position on the - strand'''
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
        ''' genomic sequence spanning the transcript (plus genomic_offset bases on
        each side), on the + strand. This is set in the direction of
        transcription, so is reverse complemented when set on the - strand.
        '''
        return self.thisptr.get_genomic_sequence().decode('utf8')
    
    @genomic_sequence.setter
    def genomic_sequence(self, str text):
        self.thisptr.add_genomic_sequence(text.encode('utf8'))

    @property
    def attributes(self):
        ''' dict-like view of the raw key/value pairs from the GTF attributes field

        The view reads each field lazily from the underlying Tx object, so
        indexing a single key (e.g. ``tx.attributes['gene_id']``) does not
        materialize the whole map. Values for keys that were repeated in the GTF
        (e.g. "tag") are joined into a single comma-separated string.
        '''
        return _Attributes(self)

    def _fix_cds_boundary(self, pos):
        ''' find where a CDS end lying in an intron is moved to, in the adjacent
        exon (for testing)

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
        ''' finds the exon closest to a site

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
        if self.thisptr.get_cds().size() == 0:
            raise ValueError(f'{self.name} has no CDS')
        
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
        ''' find the position within a codon (0, 1 or 2) of a CDS position

        Args:
            cds_pos: CDS distance in bp to CDS start
        '''
        return self.thisptr.get_position_within_codon(cds_pos)
    

    def reverse_complement(self, text):
        ''' reverse complement a sequence

        Args:
            text: DNA sequence
        '''
        if isinstance(text, str):
            text = text.encode('utf8')
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
        if codon.pop('in_coding'):
            codon['codon_seq'] = codon['codon_seq'].decode('utf8')
            codon['initial_aa'] = chr(codon['initial_aa'])
        else:
            codon.update(codon_number=None, intra_codon=None, codon_seq=None,
                initial_aa=None)
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


cdef Transcript _wrap_tx(Tx & tx):
    ''' wrap an existing C++ Tx object in a Transcript, without re-parsing it

    This copy-constructs the given Tx (preserving its exons, CDS, type,
    positions, genomic sequence and attributes) and hands it to a Transcript,
    avoiding the need to disassemble the Tx into python values and rebuild an
    equivalent object through Transcript.__init__.
    '''
    cdef Transcript transcript = Transcript.__new__(Transcript, _WRAP)
    transcript.thisptr = new Tx(tx)
    return transcript


cdef class _Attributes:
    ''' a lazy, dict-like view over a Transcript's GTF attributes

    Each lookup reads directly from the underlying Tx object, so indexing a
    single key does not build the full attribute map. Assigning to a key
    (``tx.attributes[key] = value``) writes straight through to the Tx. Repeated
    GTF keys (e.g. "tag") are exposed as a single comma-separated value.
    '''
    cdef Transcript _tx

    def __cinit__(self, Transcript tx):
        self._tx = tx

    def __getitem__(self, key):
        if not self._tx.thisptr.has_attribute(key.encode('utf8')):
            raise KeyError(key)
        return self._tx.thisptr.get_attribute(key.encode('utf8')).decode('utf8')

    def __setitem__(self, key, value):
        if not isinstance(value, str):
            raise ValueError(f'value ({value}) must be str, got {type(value)}')
        self._tx.thisptr.set_attribute(key.encode('utf8'), value.encode('utf8'))

    def __contains__(self, key):
        return self._tx.thisptr.has_attribute(key.encode('utf8'))

    def keys(self):
        return [k for k, v in self.items()]

    def values(self):
        return [v for k, v in self.items()]

    def items(self):
        cdef dict attributes = self._tx.thisptr.get_attributes()
        return [(k.decode('utf8'), v.decode('utf8')) for k, v in attributes.items()]

    def __iter__(self):
        return iter(self.keys())

    def __len__(self):
        return self._tx.thisptr.get_attributes().size()

    def __eq__(self, other):
        return dict(self.items()) == other

    def __repr__(self):
        return repr(dict(self.items()))
