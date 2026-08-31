"""
Combine two continuous box extractions into one file,
it works if both extractions contain the same variables.
"""

import xarray as xr
from pathlib import Path
from lo_tools import Lfun
Ldir = Lfun.Lstart()

#--------------------------------------
# Perturbation Run

# get files to concatenate
file1_pert = Path('../../LO_output/extract/cas7_t1dgeWB_x11abd3monthscont/box/cresst3_2020.06.01_2020.06.30.nc')
file2_pert = Path('../../LO_output/extract/cas7_t1dgeWB_x11abd/box/cresst3_2020.07.01_2020.10.31.nc')
# where to put output figures
out_dir = Ldir['LOo'] / 'extract' / 'cresst3_box_extractions'
Lfun.make_dir(out_dir)
out_file_pert = 'perturbation_testalkdose_2020.06.01_2020.10.31.nc'                                                                       

# One-line combine + write:
with xr.open_mfdataset([file1_pert, file2_pert], combine='by_coords') as ds:
    ds.sortby('ocean_time').to_netcdf(out_file_pert)

#--------------------------------------
# Baseline Run

# get files to concatenate
file1_base = Path('../../LO_output/extract/cas7_t1_x11ab/box/cresst3_2020.06.01_2020.08.31.nc')
file2_base = Path('../../LO_output/extract/cas7_t1_x11ab/box/cresst3_2020.09.01_2020.10.31.nc')
# where to put output figures
out_file_pert = 'baseline_testalkdose_2020.06.01_2020.10.31.nc'                                                                       

# One-line combine + write:
with xr.open_mfdataset([file1_pert, file2_pert], combine='by_coords') as ds:
    ds.sortby('ocean_time').to_netcdf(out_file_pert)