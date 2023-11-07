from importlib.metadata import version

__name__ = 'gencodegenes'
__version__ = version('gencodegenes')

from gencodegenes.gencode import Gencode, Gene
from gencodegenes.transcript import Transcript
