"""
Read in all of the desired fits files that will be fed into lightcurve in the very near future, and produce other versions of them that include the timestamps converted to be in bmjd_tdb in stead of mjd_(tt, utc)?
"""
from astropy.time import Time
from astropy.io import fits

def convert_times(fits_file):
    """
    function that should read in the fits files and output new ones somewhere with BMJD_tdb times in a place that the new local_lightcurve.cos.py file will find them within run_step_grid.py
    """
    pass
