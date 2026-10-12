#lets people write "from clean_dbs import read_tdt" instead of "from clean_dbs.read_tdt import read_tdt"
from .read_tdt import read_tdt
from .read_bv import read_bv
from .find_blocks import find_tdt_blocks
