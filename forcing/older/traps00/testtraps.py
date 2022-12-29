import trapsfun
from importlib import reload
reload(trapsfun)

# run testtraps.py -g cas6 -r backfill -s continuation -d 2020.01.01 -f traps00 -test True

trapsfun.traps_placement('riv') # riv or wwtp
trapsfun.traps_placement('wwtp')