import trapsfun
from importlib import reload
reload(trapsfun)

# run testtraps.py -g cas6 -r backfill -s continuation -d 2020.01.01 -f traps13 -test True

trapsfun.traps_placement('riv')
trapsfun.traps_placement('wwtp')