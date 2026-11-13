# cython: language_level=3, boundscheck=False, emit_linenums=True

import logging
from pathlib import Path
import weakref

from cython.operator cimport dereference as deref, preincrement as inc
from libc.limits cimport INT_MAX, INT_MIN
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.map cimport map
from libcpp.utility cimport move, pair

from pyfaidx import Fasta

from gencodegenes.transcript cimport (
    Tx,
    Region,
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
        int canonical_score
        map[string, string] attributes
        
    void parse_gtfline(string & line, GTFLine & info, bool all_fields) except +

cdef extern from "gencode.h" namespace "gencode":
    cdef struct NamedTx:
        string symbol
        vector[string] alternate_ids
        Tx tx
        int canonical_score
    
    cdef struct GeneInterval:
        int start
        int end
        int max_end
        int idx
        int max_end_idx
    
    vector[NamedTx] open_gencode(string, bool) except +
    void index_intervals(vector[GeneInterval] & intervals)
    vector[int] genes_in_region(const vector[GeneInterval] & intervals, int start,
        int end) except +
    pair[int, int] flanking_genes(const vector[GeneInterval] & intervals, int pos)

cpdef _parse_gtfline(string line):
    ''' python function for unit testing GTF parsing
    '''
    cdef GTFLine info
    parse_gtfline(line, info, False)
    return info

cpdef _open_gencode(gtf_path, coding_only=True):
    ''' python function for unit testing loading transcripts from GTF
    '''
    cdef vector[NamedTx] _transcripts = open_gencode(gtf_path.encode('utf8'), coding_only)
    cdef Tx tx
    
    transcripts = []
    for x in _transcripts:
        tx = x.tx
        transcripts.append((x.symbol.decode('utf8'), _wrap_tx(tx), x.canonical_score))
    return transcripts

# bases of genomic sequence loaded beyond each end of a transcript, unless the
# transcript sets its own genomic offset
cdef int DEFAULT_FLANK = 5

cdef list _chrom_names(str chrom):
    ''' get the names a chromosome could have, allowing for differences in the
    'chr' prefix (e.g. 'chr1' vs '1'), and in naming the mitochondrial genome
    '''
    bare = chrom[3:] if chrom.startswith('chr') else chrom
    names = [chrom, bare, f'chr{bare}']
    if bare in ('M', 'MT'):
        names += ['MT', 'chrM', 'M', 'chrMT']
    return list(dict.fromkeys(names))

class _Genome:
    ''' genome fasta, which keeps the sequences of the most recently read
    regions, so a gene's transcripts only need one read from the fasta
    '''
    def __init__(self, path, max_cached=5):
        self._fasta = Fasta(str(path))
        self._max_cached = max_cached
        self._cache = []  # (chrom, start, end, seq) tuples, most recent last
        self._lengths = {}
    
    def close(self):
        self._fasta.close()
    
    def chrom_name(self, chrom):
        ''' find the name the fasta uses for a chromosome
        '''
        if chrom not in self._fasta:
            for name in _chrom_names(chrom):
                if name in self._fasta:
                    return name
        return chrom
    
    def chrom_length(self, chrom):
        if chrom not in self._lengths:
            self._lengths[chrom] = len(self._fasta[chrom])
        return self._lengths[chrom]
    
    def region(self, chrom, start, end, region_start, region_end):
        ''' get a cached sequence which covers start to end (0-based, end
        exclusive), reading the surrounding region from the fasta if needed
        
        Returns:
            tuple of the sequence (as bytes), and the index of start within it
        '''
        for i, (cached_chrom, cached_start, cached_end, seq) in enumerate(self._cache):
            if cached_chrom == chrom and cached_start <= start and end <= cached_end:
                if i != len(self._cache) - 1:
                    self._cache.append(self._cache.pop(i))
                return seq, start - cached_start
        
        region_start = max(0, min(start, region_start))
        region_end = min(self.chrom_length(chrom), max(end, region_end))
        seq = self._fasta[chrom][region_start:region_end].seq.encode('utf8')
        self._cache.append((chrom, region_start, region_end, seq))
        if len(self._cache) > self._max_cached:
            del self._cache[0]
        return seq, start - region_start
    
    def sequence(self, chrom, start, end, region_start, region_end):
        ''' get the sequence from start to end (0-based, end exclusive) as bytes
        '''
        seq, idx = self.region(chrom, start, end, region_start, region_end)
        return seq[idx:idx + end - start]
    
    def add_sequence(self, Transcript transcript, int gene_start, int gene_end):
        ''' add the genomic sequence (with flanks) to a transcript, matching
        the transcript's strand orientation
        
        Args:
            transcript: Transcript lacking a genomic sequence
            gene_start: start of the gene containing the transcript. The whole
                gene (plus default flanks) is read, for the other transcripts.
            gene_end: end of the gene containing the transcript
        '''
        cdef Tx * tx = transcript.thisptr
        cdef const char * data
        cdef bytes seq
        cdef int idx, start = tx.get_start(), end = tx.get_end()
        cdef size_t length
        
        offset = DEFAULT_FLANK if tx.get_genomic_offset() == 0 else tx.get_genomic_offset()
        chrom = self.chrom_name(tx.get_chrom().decode('utf8'))
        # shrink the flanks for transcripts near the ends of the chromosome
        offset = max(0, min(offset, start - 1, self.chrom_length(chrom) - end))
        # region() takes 0-based coordinates, with exclusive ends
        seq, idx = self.region(chrom, start - 1 - offset, end + offset,
            gene_start - 1 - DEFAULT_FLANK, gene_end + DEFAULT_FLANK)
        transcript.genomic_offset = offset
        # copy straight from the cached bytes, without slicing them first.
        # Sequence missing from the fasta ends up too short, which raises
        # an error in add_genomic_sequence
        data = seq
        idx = min(idx, len(seq))
        length = min(end + offset - (start - 1 - offset), len(seq) - idx)
        tx.add_genomic_sequence(move(string(data + idx, length)), True)

cdef class Gencode

cdef class Gene:
    cdef string _symbol
    cdef vector[Tx] _transcripts
    cdef vector[int] _canonical_scores
    cdef str _chrom
    cdef int _start, _end
    cdef vector[string] _alternate_ids
    cdef object _genome  # _Genome shared with the parent Gencode, or None
    cdef list _parents   # weakrefs to the Gencode objects indexing this gene
    def __cinit__(self, symbol, alt_ids=None):
        if isinstance(symbol, str):
            symbol = symbol.encode('utf8')
        self._symbol = symbol
        self._parents = []
        
        if alt_ids is None:
            alt_ids = []
        elif isinstance(alt_ids, str):
            alt_ids = [alt_ids]
        self.alternate_ids = [x.decode('utf8') if isinstance(x, bytes) else x for x in alt_ids]
        self.start = INT_MAX
        self.end = INT_MIN
    
    cdef add_tx(self, Tx tx, int canonical_score):
        self.chrom = tx.get_chrom().decode('utf8')
        self.start = min(self.start, tx.get_start())
        self.end = max(self.end, tx.get_end())
        self._transcripts.push_back(move(tx))
        self._canonical_scores.push_back(canonical_score)
    
    def add_transcript(self, _tx):
        ''' add a Transcript to the gene object
        '''
        assert isinstance(_tx, Transcript)
        cdef Transcript txn = _tx          # typed handle so Cython sees thisptr
        if self._chrom is not None and txn.chrom != self._chrom:
            raise ValueError(f'cannot add {txn.name} on {txn.chrom} to {self.symbol} on {self._chrom}')
        self.add_tx(deref(txn.thisptr), False)
    
    cdef _mark_stale(self):
        ''' flag that the Gencode objects indexing this gene need to reindex it
        '''
        cdef Gencode parent
        for ref in self._parents:
            parent = ref()
            if parent is not None:
                parent._stale = True
    
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
        self._mark_stale()
    
    @property
    def start(self):
        return self._start
    @start.setter
    def start(self, value):
        self._start = value
        self._mark_stale()
    
    @property
    def end(self):
        return self._end
    @end.setter
    def end(self, value):
        self._end = value
        self._mark_stale()
    
    @property
    def strand(self):
        if self._transcripts.size() > 0:
            return chr(self._transcripts[0].get_strand())
        raise IndexError('no transcripts in gene yet')
    
    cdef _to_Transcript(self, Tx & tx):
        ''' construct Transcript (python object) from Tx (c++ object)
        '''
        cdef Transcript transcript = _wrap_tx(tx)
        if tx.get_genomic_sequence().size() == 0 and self._genome is not None:
            self._genome.add_sequence(transcript, self._start, self._end)
        return transcript
    
    @property
    def transcripts(self):
        ''' get list of Transcripts for gene, with genomic DNA included
        '''
        cdef size_t i
        return [self._to_Transcript(self._transcripts[i]) for i in range(self._transcripts.size())]
    
    cdef int _exonic_len(self, Tx & tx):
        ''' get length of exonic sequence for a Tx object based transcript
        '''
        cdef int length = 0
        cdef size_t i
        for i in range(tx.get_exons().size()):
            length += (tx.get_exons()[i].end - tx.get_exons()[i].start) + 1
        return length
    
    cdef size_t _longest(self, vector[size_t] & indices, bint coding) except *:
        ''' get index of the longest transcript, by CDS or exonic length
        '''
        cdef size_t i, max_i = 0
        cdef int curr_len, length = 0
        for i in indices:
            if coding:
                curr_len = self._transcripts[i].get_cds_length()
            else:
                curr_len = self._exonic_len(self._transcripts[i])
            if curr_len > length:
                length = curr_len
                max_i = i
        if length == 0:
            raise ValueError(f'no {"coding" if coding else "exonic"} transcripts')
        return max_i
    
    @property
    def canonical(self):
        ''' find the canonical transcript for a gene.
        
        Canonical is defined as:
            - transcript with Ensembl_canonical tag (preferred)
            - transcript with longest CDS tagged with appris_principal in the GTF
            - if no appris_principal tags for any tx, use the tx with longest CDS
            # TODO: for the last case, maybe check the protein coding subset first
        
        Occasionally there are multiple transcripts tagged as appris_principal
        and with the same longest CDS, we use the first one of those.
        '''
        return self._to_Transcript(self._transcripts[self._canonical_idx()])
    
    cdef size_t _canonical_idx(self) except *:
        ''' get the index of the canonical transcript (see canonical)
        '''
        cdef size_t i
        cdef int score, max_score = 0
        for score in self._canonical_scores:
            max_score = max(max_score, score)
        
        cdef vector[size_t] canonical
        for i in range(self._transcripts.size()):
            if self._canonical_scores[i] == max_score:
                canonical.push_back(i)
        
        try:
            return self._longest(canonical, True)
        except ValueError:
            return self._longest(canonical, False)
    
    def in_any_tx_cds(self, int pos):
        ''' find if a pos is in coding region of any transcript of a gene
        '''
        cdef size_t i
        for i in range(self._transcripts.size()):
            if self._transcripts[i].in_coding_region(pos):
                return True
        return False
    
    def distance(self, chrom, pos):
        ''' get distance to nearest boundary of a gene
        '''
        if chrom != self.chrom and chrom not in _chrom_names(self.chrom):
            return None
        if self.start <= pos <= self.end:
            return 0
        return min(abs(self.start - pos), abs(self.end - pos))

cdef class Gencode:
    # maps symbol to a list of Genes, as some symbols are used at multiple loci,
    # e.g. PAR genes on chrX and chrY
    cdef dict genes
    cdef list _indexed  # genes in the order they are indexed in intervals
    cdef map[string, vector[GeneInterval]] intervals
    cdef object _genome
    cdef bint _stale  # set when a gene's span may have changed since indexing
    cdef object __weakref__
    def __cinit__(self, gencode=None, fasta=None, coding_only=True):
        ''' initialise Gencode
        
        Args:
            gencode: path to gencode annotations file
            fasta: path to fasta for genome matching annotations build
            coding_only: restrict to protein_coding only by default
        '''
        if gencode is not None and not Path(gencode).exists():
            raise ValueError(f'cannot find gencode at: {gencode}')
        if fasta is not None and not Path(fasta).exists():
            raise ValueError(f'cannot find fasta at: {fasta}')
        self.genes = {}
        self._genome = None
        if fasta:
            logging.info(f'opening genome fasta: {fasta}')
            self._genome = _Genome(fasta)
        logging.info(f'opening gencode annotations: {gencode}')
        cdef vector[NamedTx] transcripts
        cdef NamedTx * x
        cdef size_t i
        cdef Gene curr
        cdef bytes chrom, gene_id
        loci = {}
        parent = weakref.ref(self)
        if gencode is not None:
            transcripts = open_gencode(str(gencode).encode('utf8'), coding_only)
            for i in range(transcripts.size()):
                x = &transcripts[i]
                # group transcripts into genes by chrom and gene_id, so genes
                # sharing a symbol at different loci are kept apart
                gene_id = x.symbol
                if x.tx.has_attribute(b'gene_id'):
                    gene_id = x.tx.get_attribute(b'gene_id')
                chrom = x.tx.get_chrom()
                key = (chrom, gene_id)
                if key not in loci:
                    curr = Gene(x.symbol, x.alternate_ids)
                    curr._genome = self._genome
                    curr._parents.append(parent)
                    loci[key] = curr
                    self.genes.setdefault(curr.symbol, []).append(curr)
                curr = loci[key]
                curr.add_tx(move(x.tx), x.canonical_score)
        self._index()
    
    def _index(self):
        ''' index gene coordinates per chromosome, to speed finding genes by
        position
        '''
        self._stale = False
        self._indexed = []
        self.intervals.clear()
        cdef Gene gene
        cdef string chrom
        for genes in self.genes.values():
            for gene in genes:
                chrom = gene._chrom.encode('utf8')
                self.intervals[chrom].push_back(GeneInterval(gene._start, gene._end,
                    0, len(self._indexed), 0))
                self._indexed.append(gene)
        
        cdef map[string, vector[GeneInterval]].iterator it = self.intervals.begin()
        while it != self.intervals.end():
            index_intervals(deref(it).second)
            inc(it)
    
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
            gene._parents.append(weakref.ref(self))
            genes.append(gene)
            self._stale = True
    
    def nearest(self, str chrom, int pos):
        ''' find the nearest gene to a genomic chrom, pos coordinate
        '''
        if self._stale:
            self._index()
        cdef vector[GeneInterval] * intervals = &self.intervals[self._match_chrom(chrom)]
        
        # first, account for any overlapping genes
        cdef Gene gene
        overlaps = [self._indexed[i] for i in genes_in_region(deref(intervals), pos, pos)]
        if len(overlaps) > 0:
            # if we have > 0 prioritise if the position is in the CDS
            cds_overlaps = [x for x in overlaps if x.in_any_tx_cds(pos)]
            if len(cds_overlaps) > 0:
                overlaps = cds_overlaps
            # prioritise the gene with longest CDS (in the canonical tx), read
            # from the Tx directly, to avoid loading sequence from the genome
            lengths = []
            for gene in overlaps:
                lengths.append(gene._transcripts[gene._canonical_idx()].get_cds_length())
            idx = lengths.index(max(lengths))
            return overlaps[idx]
        
        # no overlaps observed, so pick the closer of the gene ending closest
        # before pos, and the gene starting closest after pos
        cdef pair[int, int] flanks = flanking_genes(deref(intervals), pos)
        if flanks.first < 0:
            return self._indexed[flanks.second]
        if flanks.second < 0:
            return self._indexed[flanks.first]
        cdef Gene upstream = self._indexed[flanks.first]
        cdef Gene downstream = self._indexed[flanks.second]
        if pos - upstream._end <= downstream._start - pos:
            return upstream
        return downstream
    
    def in_region(self, str chrom, int start, int end, int max_window=2500000):
        ''' find genes within a genomic region
        
        Args:
            chrom: chromosome to search on
            start: start position of region
            end: end position of region
            max_window: unused, only kept for backwards compatibility
        
        Returns:
            list of Gene objects, ordered by gene start
        '''
        if self._stale:
            self._index()
        cdef vector[GeneInterval] * intervals = &self.intervals[self._match_chrom(chrom)]
        return [self._indexed[i] for i in genes_in_region(deref(intervals), start, end)]
    
    cdef bytes _match_chrom(self, str chrom):
        ''' find the chromosome name used in the annotations, allowing for
        differences in the 'chr' prefix (e.g. 'chr1' vs '1'), and in naming
        the mitochondrial genome
        '''
        key = chrom.encode('utf8')
        if self.intervals.count(key) > 0:
            return key
        for name in _chrom_names(chrom):
            key = name.encode('utf8')
            if self.intervals.count(key) > 0:
                return key
        raise ValueError(f'unknown_chrom: {chrom}')
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type=None, exc_value=None, traceback=None):
        ''' close the genome fasta (if one was opened)
        '''
        if self._genome is not None:
            self._genome.close()
