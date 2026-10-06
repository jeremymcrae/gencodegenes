from os import PathLike
from types import TracebackType
from typing import Iterator

from gencodegenes.transcript import Transcript

StrPath = str | PathLike[str]

class Gene:
    def __init__(self, symbol: str | bytes, alt_ids: list[str] | str | None = ...) -> None: ...
    def __repr__(self) -> str: ...

    @property
    def symbol(self) -> str: ...

    @property
    def alternate_ids(self) -> list[str]: ...
    @alternate_ids.setter
    def alternate_ids(self, ids: list[str] | str) -> None: ...

    @property
    def chrom(self) -> str: ...
    @chrom.setter
    def chrom(self, value: str) -> None: ...

    @property
    def start(self) -> int: ...
    @start.setter
    def start(self, value: int) -> None: ...

    @property
    def end(self) -> int: ...
    @end.setter
    def end(self, value: int) -> None: ...

    @property
    def strand(self) -> str: ...

    @property
    def transcripts(self) -> list[Transcript]:
        """get list of Transcripts for gene, with genomic DNA included"""
        ...

    @property
    def canonical(self) -> Transcript:
        """find the canonical transcript for a gene"""
        ...

    def add_transcript(self, _tx: Transcript) -> None:
        """add a Transcript to the gene object"""
        ...

    def in_any_tx_cds(self, pos: int) -> bool:
        """find if a pos is in coding region of any transcript of a gene"""
        ...

    def distance(self, chrom: str, pos: int) -> int | None:
        """get distance to nearest boundary of a gene"""
        ...

class Gencode:
    def __init__(
        self,
        gencode: StrPath | None = ...,
        fasta: StrPath | None = ...,
        coding_only: bool = ...,
    ) -> None: ...

    def __repr__(self) -> str: ...
    def __len__(self) -> int: ...
    def __getitem__(self, symbol: str) -> Gene: ...
    def __iter__(self) -> Iterator[str]: ...

    def add_gene(self, gene: Gene) -> None:
        """add another gene to the Gencode object"""
        ...

    def nearest(self, chrom: str, pos: int) -> Gene:
        """find the nearest gene to a genomic chrom, pos coordinate"""
        ...

    def in_region(
        self, _chrom: str, start: int, end: int, max_window: int = ...
    ) -> list[Gene]:
        """find genes within a genomic region"""
        ...

    def __enter__(self) -> Gencode: ...
    def __exit__(
        self,
        exc_type: type[BaseException] | None = ...,
        exc_value: BaseException | None = ...,
        traceback: TracebackType | None = ...,
    ) -> None:
        """close the genome fasta (if one was opened)"""
        ...
