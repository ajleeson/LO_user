import trapsfun
from importlib import reload
reload(trapsfun)

# run testtraps.py -g cas6 -r backfill -s continuation -d 2021.01.01 -f TRAPS_newest -test True

trapsfun.traps_placement('riv')
trapsfun.traps_placement('wwtp')