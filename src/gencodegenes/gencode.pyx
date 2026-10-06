# cython: language_level=3, boundscheck=False, emit_linenums=True

import bisect
import logging
from pathlib import Path

from cython.operator cimport dereference as deref
from libcpp.algorithm cimport lower_bound, upper_bound
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.map cimport map

from pyfaidx import Fasta

from gencodegenes.transcript cimport (
    Tx,
    Region,
    CDS_coords,
    Transcript,
    _wrap_tx,
    )

cdef extern from "gtf.h" namespace "gencode":
    cdef struct GTFLine:
        string chrom
        string feature
        int start
        int end
        string strand
        string symbol
        vector[string] alternate_ids
        string tx_id
        string transcript_type
        int is_canonical
        map[string, string] attributes
        
    GTFLine parse_gtfline(string line) except +

cdef extern from "gencode.h" namespace "gencode":
    cdef struct NamedTx:
        string symbol
        vector[string] alternate_ids
        Tx tx
        int is_canonical
    
    cdef struct GenePoint:
        int pos
        string symbol
    
    vector[NamedTx] open_gencode(string, bool) except +
    bool CompFunc(const GenePoint &l, const GenePoint &r)
    vector[string] _in_region(string chrom, int start, int end, 
        map[string, vector[GenePoint]] & starts, map[string, vector[GenePoint]] & ends,
        int max_window) except+

cpdef _parse_gtfline(string line):
    ''' python function for unit testing GTF parsing
    '''
    return parse_gtfline(line)

cdef _convert_exons(vector[Region] exons):
    ''' convert vector of exon Regions to list of lists
    
    We need exons and CDS as lists of lists for constructing the python 
    Transcript object.
    '''
    return [[y.start, y.end] for y in exons]

cpdef _open_gencode(gtf_path, coding_only=True):
    ''' python function for unit testing loading transcripts from GTF
    '''
    cdef vector[NamedTx] _transcripts = open_gencode(gtf_path.encode('utf8'), coding_only)
    cdef Tx tx
    
    transcripts = []
    for x in _transcripts:
        tx = x.tx
        transcripts.append((x.symbol.decode('utf8'), _wrap_tx(tx), x.is_canonical))
    return transcripts

cdef class Gene:
    cdef string _symbol
    cdef vector[Tx] _transcripts
    cdef vector[int] _canonical
    cdef str _chrom
    cdef int _start, _end
    cdef vector[string] _alternate_ids
    cdef object _genome  # pyfaidx Fasta shared with the parent Gencode, or None
    def __cinit__(self, symbol, alt_ids=None):
        if isinstance(symbol, str):
            symbol = symbol.encode('utf8')
        self._symbol = symbol
        
        if alt_ids is None:
            alt_ids = []
        elif isinstance(alt_ids, str):
            alt_ids = [alt_ids]
        self.alternate_ids = [x.decode('utf8') if isinstance(x, bytes) else x for x in alt_ids]
        self.start = 999999999
        self.end = -999999999
    
    cdef add_tx(self, Tx tx, int is_canonical):
        self._transcripts.push_back(tx)
        self._canonical.push_back(is_canonical)
        self.chrom = tx.get_chrom().decode('utf8')
        self.start = min(self.start, tx.get_start())
        self.end = max(self.end, tx.get_end())
    
    def add_transcript(self, _tx):
        ''' add a Transcript to the gene object
        '''
        assert isinstance(_tx, Transcript)
        cdef Transcript txn = _tx          # typed handle so Cython sees thisptr
        self.add_tx(deref(txn.thisptr), False)
    
    def __repr__(self):
        chrom = self.chrom
        return f'Gene("{self.symbol}", {chrom}:{self.start}-{self.end})'
    
    @property
    def symbol(self):
        return self._symbol.decode('utf8')
    
    @property
    def alternate_ids(self):
        return [x.decode('utf8') for x in self._alternate_ids]
    
    @alternate_ids.setter
    def alternate_ids(self, ids):
        self._alternate_ids.clear()
        if isinstance(ids, str):
            ids = [ids]
        for x in ids:
            self._alternate_ids.push_back(x.encode('utf8'))
    
    @property
    def chrom(self):
        return self._chrom
    @chrom.setter
    def chrom(self, value):
        self._chrom = value
    
    @property
    def start(self):
        return self._start
    @start.setter
    def start(self, value):
        self._start = value
    
    @property
    def end(self):
        return self._end
    @end.setter
    def end(self, value):
        self._end = value
    
    @property
    def strand(self):
        if self._transcripts.size() > 0:
            return chr(self._transcripts[0].get_strand())
        raise IndexError('no transcripts in gene yet')
    
    cdef _to_Transcript(self, Tx tx):
        ''' construct Transcript (python object) from Tx (c++ object)
        '''
        cdef Transcript transcript = _wrap_tx(tx)
        
        # if the transcript lacks a genomic sequence, pull one from the genome
        # fasta (if available), matching the transcript's strand orientation
        if tx.get_genomic_sequence().size() == 0 and self._genome is not None:
            offset = 5 if tx.get_genomic_offset() == 0 else tx.get_genomic_offset()
            chrom = tx.get_chrom().decode('utf8')
            start = tx.get_start()
            end = tx.get_end()
            # shrink the flanks for transcripts near the ends of the chromosome
            offset = max(0, min(offset, start - 1, len(self._genome[chrom]) - end + 1))
            seq = self._genome[chrom][start-1-offset:end-1+offset].seq.upper()
            if chr(tx.get_strand()) == '-':
                seq = tx.reverse_complement(seq.encode('utf8')).decode('utf8')
            transcript.genomic_offset = offset
            transcript.genomic_sequence = seq
        
        return transcript
    
    @property
    def transcripts(self):
        ''' get list of Transcripts for gene, with genomic DNA included
        '''
        return [self._to_Transcript(x) for x in self._transcripts]
    
    cdef int _cds_len(self, Tx tx):
        ''' get length of coding sequence for a Tx object based transcript
        '''
        if tx.get_cds().size() == 0:
            return 0
        cdef CDS_coords coords = tx.get_coding_distance(tx.get_cds_end())
        return coords.position + 1
    
    cdef Tx _max_by_cds(self, vector[Tx] transcripts) except *:
        ''' get longest transcript by CDS length
        '''
        cdef Tx max_tx
        length = 0
        for tx in transcripts:
            curr_len = self._cds_len(tx)
            if curr_len > length:
                length = curr_len
                max_tx = tx
        if length == 0:
            raise ValueError('no coding transcripts')
        return max_tx
    
    cdef int _exonic_len(self, Tx tx):
        ''' get length of exonic sequence for a Tx object based transcript
        '''
        length = 0
        for start, end in _convert_exons(tx.get_exons()):
            length += abs(end - start) + 1
        return length
        # return sum(abs(x.end - x.start) + 1 for x in tx.get_exons())
    
    cdef Tx _max_by_exonic(self, vector[Tx] transcripts) except *:
        ''' get longest transcript by CDS length
        '''
        cdef Tx max_tx
        length = 0
        for tx in transcripts:
            curr_len = self._exonic_len(tx)
            if curr_len > length:
                length = curr_len
                max_tx = tx
        if length == 0:
            raise ValueError('no exonic transcripts')
        return max_tx
    
    @property
    def canonical(self):
        ''' find the canonical transcript for a gene.
        
        Canonical is defined as:
            - transcript with Ensembl_canonical tag (peferred)
            - transcript with longest CDS tagged with appris_principal in the GTF
            - if no appris_principal tags for any tx, use the tx with longest CDS
            # TODO: for the last case, maybe check the protein coding subset first
        
        Occasionally there are multiple transcripts tagged as appris_principal
        and with the same longest CDS, we use the first one of those.
        '''
        cdef vector[Tx] canonical
        for i in range(self._transcripts.size()):
            max_score = max(self._canonical)
            if self._canonical[i] == max_score:
                canonical.push_back(self._transcripts[i])
        
        if canonical.size() == 0:
            canonical = self._transcripts
        
        cdef Tx max_tx
        try:
            max_tx = self._max_by_cds(canonical)
        except ValueError:
            max_tx = self._max_by_exonic(canonical)
        
        return self._to_Transcript(max_tx)
    
    def in_any_tx_cds(self, pos):
        ''' find if a pos is in coding region of any transcript of a gene
        '''
        return any(tx.in_coding_region(pos) for tx in self._transcripts)
    
    def distance(self, chrom, pos):
        ''' get distance to nearest boundary of a gene
        '''
        # sanatize the chromosome first
        if self.chrom.startswith('chr') and not chrom.startswith('chr'):
            chrom = f'chr{chrom}'
        elif not self.chrom.startswith('chr') and chrom.startswith('chr'):
            chrom = chrom[3:]
        
        if self.chrom != chrom:
            return None
        if self.start <= pos <= self.end:
            return 0
        return min(abs(self.start - pos), abs(self.end - pos))

cdef class Gencode:
    # maps symbol to a list of Genes, as some symbols are used at multiple loci,
    # e.g. PAR genes on chrX and chrY
    cdef dict genes
    cdef map[string, vector[GenePoint]] starts, ends
    cdef object _genome
    def __cinit__(self, gencode=None, fasta=None, coding_only=True):
        ''' initialise Gencode
        
        Args:
            gencode: path to gencode annotations file
            fasta: path to fasta for genome matching annotations build
            coding: restrict to protein_coding only by default
        '''
        if gencode is not None and not Path(gencode).exists():
            raise ValueError(f'cannot find gencode at: {gencode}')
        if fasta is not None and not Path(fasta).exists():
            raise ValueError(f'cannot find fasta at: {fasta}')
        self.genes = {}
        self._genome = None
        if fasta:
            logging.info(f'opening genome fasta: {fasta}')
            self._genome = Fasta(str(fasta))
        logging.info(f'opening gencode annotations: {gencode}')
        cdef vector[NamedTx] transcripts
        cdef Gene curr
        loci = {}
        if gencode is not None:
            transcripts = open_gencode(str(gencode).encode('utf8'), coding_only)
            for x in transcripts:
                # group transcripts into genes by chrom and gene_id, so genes
                # sharing a symbol at different loci are kept apart
                gene_id = x.symbol
                if x.tx.has_attribute(b'gene_id'):
                    gene_id = x.tx.get_attribute(b'gene_id')
                key = (x.tx.get_chrom(), gene_id)
                if key not in loci:
                    curr = Gene(x.symbol, x.alternate_ids)
                    curr._genome = self._genome
                    loci[key] = curr
                    self.genes.setdefault(curr.symbol, []).append(curr)
                curr = loci[key]
                curr.add_tx(x.tx, x.is_canonical)
        self._sort()
    
    def _sort(self):
        ''' index by starts and ends, to speed finding genes in a region
        '''
        self.starts.clear()
        self.ends.clear()
        for symbol, genes in self.genes.items():
            for i, gene in enumerate(genes):
                chrom = gene.chrom.encode('utf8')
                key = f'{symbol}\t{i}'.encode('utf8')
                
                # ensure the chromosome is present
                if self.starts.count(chrom) == 0:
                    self.starts[chrom] = []
                if self.ends.count(chrom) == 0:
                    self.ends[chrom] = []
                
                self.starts[chrom].push_back(GenePoint(gene.start, key))
                self.ends[chrom].push_back(GenePoint(gene.end, key))
        
        # sort start and end coords by position
        for x, values in self.starts:
            self.starts[x] = sorted(values, key=lambda x: x['pos'])
        for x, values in self.ends:
            self.ends[x] = sorted(values, key=lambda x: x['pos'])
    
    def __repr__(self):
        return f'Gencode(n_genes={len(self)})'
    def __len__(self):
        return len(self.genes)
    def __getitem__(self, symbol):
        ''' get the gene for a symbol (the first loaded, if at multiple loci)
        '''
        return self.genes[symbol][0]
    def __iter__(self):
        for x in self.genes:
            yield x
    
    cdef Gene _gene_at(self, string key):
        ''' get the Gene for a key from the starts/ends index
        '''
        symbol, idx = key.decode('utf8').rsplit('\t', 1)
        return self.genes[symbol][int(idx)]
    
    def add_gene(self, Gene gene):
        ''' add another gene to the Gencode object
        
        The gene is skipped if a gene with the same symbol is already on the
        same chromosome.
        '''
        if gene.chrom is None:
            raise ValueError(f'cannot add gene without transcripts: {gene.symbol}')
        genes = self.genes.setdefault(gene.symbol, [])
        if all(x.chrom != gene.chrom for x in genes):
            if gene._genome is None:
                gene._genome = self._genome
            genes.append(gene)
        self._sort()
    
    def nearest(self, str chrom, int pos):
        ''' find the nearest gene to a genomic chrom, pos coordinate
        '''
        _chrom = self._match_chrom(chrom)
        chrom = _chrom.decode('utf8')
        
        # first, account for any overlapping genes
        overlaps = self.in_region(chrom, pos, pos)
        if len(overlaps) > 0:
            # if we have > 0 prioritise if the position is in the CDS
            cds_overlaps = [x for x in overlaps if x.in_any_tx_cds(pos)]
            if len(cds_overlaps) > 0:
                overlaps = cds_overlaps
            # prioritise the gene with longest CDS (in the canonical tx)
            txs = [x.canonical for x in overlaps]
            lengths = [x.get_coding_distance(x.cds_end)['pos'] if x.cds else 0 for x in txs]
            idx = lengths.index(max(lengths))
            return overlaps[idx]
        
        # no overlaps observed, look for the nearest upstream or downstream gene
        cdef GenePoint site = GenePoint(pos, b'A');
        cdef int i = lower_bound(self.ends[_chrom].begin(), self.ends[_chrom].end(), site, &CompFunc) - self.ends[_chrom].begin()
        cdef int j = lower_bound(self.starts[_chrom].begin(), self.starts[_chrom].end(), site, &CompFunc) - self.starts[_chrom].begin()
        
        # upstream is the gene ending closest before pos, downstream is the gene
        # starting closest after pos. starts and ends are sorted independently,
        # so each index is only valid for its own vector
        i = max(i - 1, 0)
        j = min(j, <int>self.starts[_chrom].size() - 1)
        
        upstream = self._gene_at(self.ends[_chrom][i].symbol)
        downstream = self._gene_at(self.starts[_chrom][j].symbol)
        
        if upstream.distance(chrom, pos) <= downstream.distance(chrom, pos):
            return upstream
        else:
            return downstream
    
    def in_region(self, str _chrom, int start, int end, int max_window=2500000):
        ''' find genes within a genomic region
        
        Args:
            chrom: chromosome to search on
            start: start position of region
            end: end position of region
            max_window: some genes encapsulate the region, which means we have 
                to account for gene lengths of up to 2.3 Mb in the human genome.
                This permits extra search space in other organisms.
        
        Returns:
            list of Gene objects
        '''
        symbols = _in_region(self._match_chrom(_chrom), start, end, self.starts,
            self.ends, max_window)
        return [self._gene_at(x) for x in symbols]
    
    cdef bytes _match_chrom(self, str chrom):
        ''' find the chromosome name used in the annotations, allowing for
        differences in the 'chr' prefix (e.g. 'chr1' vs '1')
        '''
        alternate = chrom[3:] if chrom.startswith('chr') else f'chr{chrom}'
        for name in [chrom, alternate]:
            key = name.encode('utf8')
            if self.starts.count(key) > 0:
                return key
        raise ValueError(f'unknown_chrom: {chrom}')
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type=None, exc_value=None, traceback=None):
        ''' close the genome fasta (if one was opened)
        '''
        if self._genome is not None:
            self._genome.close()
