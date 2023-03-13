from pkg_resources import get_distribution

__version__ = get_distribution('gencodegenes').version

from gencodegenes.gencode import Gencode, Gene
from gencodegenes.transcript import Transcript
