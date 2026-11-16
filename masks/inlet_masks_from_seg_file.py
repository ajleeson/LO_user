"""Create subdomain inlet masks from a TEF segment-index pickle."""

from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np
import xarray as xr


DEFAULT_SEG_FILE = Path(
    "/Users/aleeson/Documents/Research/LO_output/extract/tef2/"
    "seg_info_dict_cas7_c21_traps00.p"
)
DEFAULT_GRID_FILE = Path("/Users/aleeson/Documents/Research/LO_data/grids/cas7/grid.nc")
DEFAULT_REFERENCE_MASK_FILE = Path(
    "/Users/aleeson/Documents/Research/LO_output/hypVol_for_intrmdl_cmprsn/"
    "basin_masks_from_pugetsoundDObox.nc"
)
DEFAULT_OUTPUT_DIR = Path(
    "/Users/aleeson/Documents/Research/LO_output/hypVol_for_intrmdl_cmprsn"
)
DEFAULT_OUTPUT_FILE = DEFAULT_OUTPUT_DIR / "inlet_masks_from_seg_file.nc"


def find_subdomain_start(
    full_grid: xr.Dataset, reference_mask: xr.Dataset
) -> tuple[int, int]:
    """Find the reference mask's upper-left index in the full grid."""
    full_lat = full_grid["lat_rho"].values
    full_lon = full_grid["lon_rho"].values
    ref_lat = reference_mask["lat_rho"].values
    ref_lon = reference_mask["lon_rho"].values

    lat_matches = np.flatnonzero(
        np.isclose(full_lat[:, 0], ref_lat[0, 0], rtol=0, atol=1e-10)
    )
    lon_matches = np.flatnonzero(
        np.isclose(full_lon[0, :], ref_lon[0, 0], rtol=0, atol=1e-10)
    )
    if len(lat_matches) != 1 or len(lon_matches) != 1:
        raise ValueError("Could not uniquely locate the reference mask in the full grid.")

    j0, i0 = int(lat_matches[0]), int(lon_matches[0])
    j1 = j0 + ref_lat.shape[0]
    i1 = i0 + ref_lon.shape[1]
    if j1 > full_lat.shape[0] or i1 > full_lon.shape[1]:
        raise ValueError("Reference mask extends beyond the full grid.")

    if not np.allclose(
        full_lat[j0:j1, i0:i1], ref_lat, rtol=0, atol=1e-10
    ) or not np.allclose(
        full_lon[j0:j1, i0:i1], ref_lon, rtol=0, atol=1e-10
    ):
        raise ValueError("Reference mask coordinates do not match the full grid.")

    return j0, i0


def make_inlet_masks(
    seg_file: Path,
    grid_file: Path,
    reference_mask_file: Path,
    output_file: Path,
) -> None:
    """Create one binary RHO-grid mask for each inlet in ``seg_file``."""
    with seg_file.open("rb") as file:
        segment_info = pickle.load(file)

    with xr.open_dataset(grid_file) as full_grid, xr.open_dataset(
        reference_mask_file
    ) as reference_mask:
        j0, i0 = find_subdomain_start(full_grid, reference_mask)
        eta_size, xi_size = reference_mask.sizes["eta_rho"], reference_mask.sizes["xi_rho"]

        output = reference_mask[
            ["lat_rho", "lon_rho", "h", "mask_rho"]
        ].load()
        for inlet_name, inlet_info in segment_info.items():
            mask = np.zeros((eta_size, xi_size), dtype=np.int8)
            for j_full, i_full in inlet_info["ji_list"]:
                j_sub, i_sub = j_full - j0, i_full - i0
                if not (0 <= j_sub < eta_size and 0 <= i_sub < xi_size):
                    raise IndexError(
                        f"{inlet_name} point {(j_full, i_full)} is outside "
                        "the reference-mask subdomain."
                    )
                mask[j_sub, i_sub] = 1

            variable_name = inlet_name.removesuffix("_p")
            output[variable_name] = xr.DataArray(
                mask,
                dims=("eta_rho", "xi_rho"),
                attrs={
                    "long_name": f"{inlet_name} inlet mask",
                    "flag_values": np.array([0, 1], dtype=np.int8),
                    "flag_meanings": "not_inlet inlet",
                    "grid": "cas7",
                },
            )

        output.attrs["source_segment_file"] = str(seg_file)
        output.attrs["subdomain_origin_in_full_grid"] = f"j={j0}, i={i0}"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output.to_netcdf(output_file)


make_inlet_masks(
    DEFAULT_SEG_FILE,
    DEFAULT_GRID_FILE,
    DEFAULT_REFERENCE_MASK_FILE,
    DEFAULT_OUTPUT_FILE,
)
