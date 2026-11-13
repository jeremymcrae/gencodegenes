from importlib.metadata import version

__version__ = version('gencodegenes')

from gencodegenes.gencode import Gencode, Gene
from gencodegenes.transcript import Transcript

__all__ = ['Gencode', 'Gene', 'Transcript']
