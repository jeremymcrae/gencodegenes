from pkg_resources import get_distribution

__version__ = get_distribution('bgen').version

from gencodegenes.gencode import Gencode
from gencodegenes.transcript import Transcript
