# to run this script, put the followind command in the terminal:

# python get_stars_at_oc.py Gulliver_27 146.0801 -54.1154 0.12 19 

# where: Gulliver_27 is the name of the object (in our case OC), 146.0801 is the RA, -54.1154 is the DEC, 0.12 is the search radius in degrees, and 19 is the magnitude limit

import warnings
warnings.simplefilter('ignore')
import sys
import numpy as np
from matplotlib import pyplot as plt
from astropy import coordinates as coord
from astropy import units as u
from astroquery.gaia import Gaia
from astropy.table import Table
from zero_point import zpt
zpt.load_tables()

oc_identifier = str(sys.argv[1])
oc_ra = float(sys.argv[2])
oc_dec = float(sys.argv[3])
oc_search_radius = float(sys.argv[4])
oc_mag_limit = float(sys.argv[5])


def get_dr3_data(ra, dec, search_radius, mag_limit):  
    query = \
    """
    select
        g.source_id, g.ra, g.ra_error, g.dec, g.dec_error, g.parallax, g.parallax_error, g.parallax_over_error,
       
        g.ruwe, g.nu_eff_used_in_astrometry, g.pseudocolour, g.ecl_lat,
        g.astrometric_excess_noise, g.astrometric_excess_noise_sig, g.visibility_periods_used,
        g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag, g.bp_rp, g.phot_bp_rp_excess_factor, phot_g_mean_flux, g.duplicated_source,
        g.l, g.b
    from
        gaiadr3.gaia_source as g
    where
        g.astrometric_params_solved > 3 
        
        AND g.bp_rp IS NOT NULL
        AND (1 = CONTAINS(POINT({0}, {1}), CIRCLE(g.ra, g.dec, {2}))) 
        AND g.phot_g_mean_mag < {3}
    """.format(ra, dec, search_radius, mag_limit)
    job = Gaia.launch_job_async(query, dump_to_file=False)
    gdr3_source = job.get_results()
    return gdr3_source

def parallax_zero_point_correction(source_table):
    gmag = np.array(source_table['phot_g_mean_mag'])
    nueffused = np.array(source_table['nu_eff_used_in_astrometry'])
    psc = np.array(source_table['pseudocolour'])
    ecl_lat = np.array(source_table['ecl_lat'])
    soltype = np.array(source_table['astrometric_params_solved'])
    
    zpvals = zpt.get_zpt(gmag, nueffused, psc, ecl_lat, soltype) 
    corrected_parallax = np.array(source_table['parallax']) - zpvals
    return corrected_parallax
    
obtained_stars = get_dr3_data(oc_ra, oc_dec, oc_search_radius, oc_mag_limit)
obtained_stars.write('star_table_{0}.csv'.format(oc_identifier), format='csv')
