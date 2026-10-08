from typing import Iterator, TypedDict

class _Region(TypedDict):
    start: int
    end: int

class _CodonInfo(TypedDict):
    cds_pos: int
    codon_seq: str | None
    intra_codon: int | None
    codon_number: int | None
    initial_aa: str | None
    offset: int

class _CodingDistance(TypedDict):
    pos: int
    offset: int

class _Attributes:
    """ a lazy, dict-like view over a Transcript's GTF attributes """
    def __getitem__(self, key: str) -> str: ...
    def __setitem__(self, key: str, value: str) -> None: ...
    def __contains__(self, key: str) -> bool: ...
    def keys(self) -> list[str]: ...
    def values(self) -> list[str]: ...
    def items(self) -> list[tuple[str, str]]: ...
    def __iter__(self) -> Iterator[str]: ...
    def __len__(self) -> int: ...
    def __eq__(self, other: object) -> bool: ...
    def __repr__(self) -> str: ...

class Transcript:
    def __init__(
        self,
        name: str | None = ...,
        chrom: str | None = ...,
        start: int = ...,
        end: int = ...,
        strand: str = ...,
        transcript_type: str = ...,
        exons: list[tuple[int, int]] | None = ...,
        cds: list[tuple[int, int]] | None = ...,
        sequence: str | None = ...,
        offset: int = ...,
        attributes: dict[str, str] | None = ...,
    ) -> None: ...

    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...
    def __hash__(self) -> int: ...
    def __eq__(self, other: object) -> bool: ...
    def __add__(self, other: Transcript | None) -> Transcript: ...
    def __radd__(self, other: Transcript | None) -> Transcript: ...

    @property
    def name(self) -> str:
        """transcript ID"""
        ...
    @property
    def chrom(self) -> str:
        """chromosome the gene is on"""
        ...
    @property
    def type(self) -> str:
        """transcript functional type (protein_coding etc.)"""
        ...
    @property
    def start(self) -> int:
        """transcript start position (TSS) on chromosome"""
        ...
    @property
    def end(self) -> int:
        """transcript end position on chromosome"""
        ...
    @property
    def strand(self) -> str:
        """strand the transcript is on (+ or -)"""
        ...
    @property
    def cds_start(self) -> int:
        """CDS start position on chromosome"""
        ...
    @property
    def cds_end(self) -> int:
        """CDS end position on chromosome"""
        ...

    @property
    def exons(self) -> list[_Region]:
        """list of exon coords as {'start': int, 'end': int} dicts"""
        ...
    @exons.setter
    def exons(self, exon_ranges: list[tuple[int, int]]) -> None: ...

    @property
    def cds(self) -> list[_Region]:
        """list of CDS coords as {'start': int, 'end': int} dicts"""
        ...
    @cds.setter
    def cds(self, cds_ranges: list[tuple[int, int]]) -> None: ...

    @property
    def cds_sequence(self) -> str:
        """coding sequence, in the direction of transcription"""
        ...
    @cds_sequence.setter
    def cds_sequence(self, text: str) -> None: ...

    @property
    def genomic_offset(self) -> int:
        """distance the DNA sequence extends symmetrically beyond gene boundaries"""
        ...
    @genomic_offset.setter
    def genomic_offset(self, offset: int) -> None: ...

    @property
    def genomic_sequence(self) -> str:
        """DNA sequence of the transcript region, plus genomic_offset bp on each
        side. This is set in the direction of transcription, but returned on
        the + strand.
        """
        ...
    @genomic_sequence.setter
    def genomic_sequence(self, text: str) -> None: ...

    @property
    def attributes(self) -> _Attributes:
        """dict-like view of the raw key/value pairs from the GTF attributes field"""
        ...

    def in_exons(self, position: int) -> bool:
        """check if a chromosome position lies within any exon (coding or UTR)

        Args:
            position: chromosome position (1-based)

        Returns:
            True if the position is inside an exon, False otherwise

        Raises:
            ValueError: if the transcript has no exons
        """
        ...
    def get_closest_exon(self, pos: int) -> _Region:
        """find the exon containing, or nearest to, a chromosome position

        For positions outside exons, this is the exon with the closest boundary.

        Args:
            pos: chromosome position (1-based)

        Returns:
            exon coords as a {'start': int, 'end': int} dict

        Raises:
            ValueError: if the transcript has no exons
        """
        ...
    def in_coding_region(self, pos: int) -> bool:
        """check if a chromosome position lies within the coding sequence

        Only bases within CDS regions count, so intronic and UTR positions
        return False, as do all positions in non-coding transcripts.

        Args:
            pos: chromosome position (1-based)

        Returns:
            True if the position is a coding base, False otherwise
        """
        ...
    def get_coding_distance(self, pos: int) -> _CodingDistance:
        """get the distance from the CDS start to a chromosome position

        The distance counts exonic bases only (skipping introns), in the
        direction of transcription. The CDS start itself is 0, positions
        upstream (5') of the CDS start are negative, and positions
        downstream (3') are positive.

        Intronic positions are first moved to the nearest exon boundary, and
        the distance into the intron is reported separately as the offset.

        Args:
            pos: chromosome position (1-based)

        Returns:
            dict with:
                pos: exonic distance in bp from the CDS start
                offset: bp into the intron from the nearest exon boundary (0
                    for exonic positions). Positive offsets are 3' of an exon
                    (splice donor side), negative offsets are 5' of an exon
                    (splice acceptor side).
        """
        ...
    def get_position_on_chrom(self, pos: int, offset: int = ...) -> int:
        """convert a CDS position back to a chromosome position

        This is the inverse of get_coding_distance(), for positions within
        the CDS.

        Args:
            pos: distance in coding bp from the CDS start (0-based)
            offset: bp into the adjacent intron, in the direction of
                transcription (default 0, for coding positions)

        Returns:
            chromosome position (1-based)

        Raises:
            ValueError: if the position isn't within the CDS
        """
        ...
    def get_codon_number_for_cds_position(self, cds_pos: int) -> int:
        """find which codon a CDS position falls in

        Args:
            cds_pos: distance in coding bp from the CDS start (0-based)

        Returns:
            codon number (0-based, so the start codon is 0)
        """
        ...
    def get_position_within_codon(self, cds_pos: int) -> int:
        """find the position of a CDS base within its codon

        Args:
            cds_pos: distance in coding bp from the CDS start (0-based)

        Returns:
            position within the codon: 0, 1 or 2
        """
        ...
    def reverse_complement(self, text: str) -> str:
        """reverse complement a DNA or RNA sequence

        Case is preserved, IUPAC ambiguity codes are complemented, U pairs with
        A, and any other character becomes N.

        Args:
            text: DNA or RNA sequence

        Returns:
            reverse complemented DNA sequence
        """
        ...
    def get_centered_sequence(self, pos: int, length: int = ...) -> str:
        """get the genomic sequence centered on a chromosome position

        The sequence is from the + strand, in upper case, and requires the
        genomic sequence to be loaded (e.g. via Gencode(..., fasta=path)).

        Args:
            pos: chromosome position (1-based) at the center of the sequence
            length: number of bases to return, which must be odd (default 3)

        Returns:
            DNA sequence e.g. 'TCC' for a trinucleotide around pos

        Raises:
            ValueError: if length is even, pos is negative, or the sequence
                extends beyond the loaded genomic sequence
        """
        ...
    def get_codon_sequence(self, codon_num: int) -> str:
        """get the coding sequence for a codon

        Args:
            codon_num: codon number (0-based, so the start codon is 0)

        Returns:
            3 bp codon sequence, oriented in the direction of transcription

        Raises:
            ValueError: if the codon number is outside the coding sequence
        """
        ...
    def translate(self, text: str) -> str:
        """translate a DNA sequence to amino acids, using the standard genetic code

        Args:
            text: coding DNA sequence (case-insensitive), with a length that's
                a multiple of 3

        Returns:
            protein sequence as single letter amino acid codes, with '*' for
            stop codons

        Raises:
            ValueError: if the length isn't a multiple of 3, or a codon can't be
                translated (e.g. it contains an N)
        """
        ...
    def get_codon_info(self, pos: int) -> _CodonInfo:
        """get details of the codon at a chromosome position

        Requires the coding sequence to be available (from the genomic
        sequence, or set via cds_sequence).

        Args:
            pos: chromosome position (1-based)

        Returns:
            dict with:
                cds_pos: exonic distance in bp from the CDS start (see
                    get_coding_distance)
                codon_seq: 3 bp codon sequence, in the direction of transcription
                intra_codon: position within the codon (0, 1 or 2)
                codon_number: codon number (0-based)
                initial_aa: single letter amino acid that the codon encodes
                offset: bp into the intron from the nearest exon boundary (0 for
                    exonic positions)
            For intronic positions, codon_seq, intra_codon, codon_number and
            initial_aa are None.

        Raises:
            ValueError: if the position lies before the CDS start or after the
                CDS end
        """
        ...
    def get_boundary_distance(self, pos: int) -> int:
        """get the distance in bp from a position to the nearest exon boundary

        Uses the closest exon. The first intronic base next to an exon is 1 bp
        away. Coding positions get 1 bp added, so the first or last base of a
        coding exon is also 1 bp away, whereas a UTR or non-coding exon's
        boundary base is 0 bp away.

        Args:
            pos: chromosome position (1-based)

        Returns:
            distance in bp to the nearest exon boundary

        Raises:
            ValueError: if the transcript has no exons
        """
        ...
    def consequence(self, pos: int, ref: str, alt: str) -> str:
        """predict the consequence of a variant on this transcript

        Alleles are given on the + strand. Shared leading and trailing bases
        are trimmed before classification, so VCF-style alleles with a padding
        base are fine. SNVs and MNVs in the CDS need the coding sequence
        to be available; MNVs return the most severe consequence across the
        affected codons.

        Args:
            pos: chromosome position (1-based) of the first reference base
            ref: reference allele DNA sequence
            alt: alternate allele DNA sequence

        Returns:
            a Sequence Ontology term, one of:
                outside the transcript: upstream_gene_variant,
                    downstream_gene_variant (within 5 kb), intergenic_variant
                outside the CDS: 5_prime_UTR_variant, 3_prime_UTR_variant,
                    non_coding_transcript_exon_variant, intron_variant,
                    splice_acceptor_variant, splice_donor_variant (within 2 bp
                    of an exon), splice_region_variant (within 8 bp of an exon)
                indels in the CDS: frameshift_variant, inframe_insertion,
                    inframe_deletion, splice_donor_variant
                SNVs and MNVs in the CDS: stop_gained, stop_lost, start_lost,
                    missense_variant, start_retained_variant,
                    stop_retained_variant, synonymous_variant
                coding_sequence_variant, for variants in introns shorter than
                    10 bp, or that overlap the CDS but change no coding bases

        Raises:
            ValueError: if the variant can't be placed on the transcript (e.g.
                an indel without an overlapping exon), or the ref allele doesn't
                match the genomic sequence
        """
        ...
