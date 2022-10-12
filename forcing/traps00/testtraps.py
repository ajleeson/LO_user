import trapsfun
from importlib import reload
reload(trapsfun)

# run testtraps.py -g cas6 -t v0 -r backfill -s continuation -d 2010.01.01 -f traps00 -test True

result = trapsfun.traps_placement('riv') # riv or wwtp
# print(result)