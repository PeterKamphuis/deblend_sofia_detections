
import urllib.request
from urllib.parse import urlencode
from astropy import units as u
from astroquery import log
log.setLevel('WARN')


from astropy.coordinates import SkyCoord
from astropy.table import vstack


import numpy as np
import warnings
import os
import sys

def stopPrint(func, *args, **kwargs):
    with open(os.devnull,"w") as devNull:
        original = sys.stdout
        sys.stdout = devNull
        func(*args, **kwargs)
        sys.stdout = original 
stopPrint(lambda: __import__('astroquery.gaia'))
from astroquery.gaia import Gaia
    
class LW_GaiaQuery:
    def __init__(self,  credentials=['NONE','NONE']):
      
        self.coords = "<RA DEC>"
        self.size = 1.*u.arcmin
        self.verbose = False
        self.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
        self.equinox = 'J2000'
        self.credentials = credentials
       
    def query_region(self, coord,radius = None, verbose=False, maxrec = None, 
                equinox = None):
        
        if radius is None:
            radius = self.size
        if radius < 3.*u.arcmin:
            table = self.single_query(coord, size=radius, 
                verbose=verbose, maxrec=maxrec, equinox=equinox)
        else:
            n_chunks = int(np.ceil((radius.to(u.arcmin)/(3.*u.arcmin)).decompose().value))
            single_chunk_size = radius.to(u.arcmin)/ n_chunks
            if verbose:
                print(f"Dividing query region into {n_chunks}x{n_chunks} chunks of size {single_chunk_size}")
            for i in range(n_chunks):
                for j in range(n_chunks):
                    ra_offset = ((i - n_chunks/2 + 0.5) * single_chunk_size) \
                        / np.cos(coord.dec.to(u.rad).value)
                    dec_offset = (j - n_chunks/2 + 0.5) * single_chunk_size
                    ind_coord = SkyCoord(
                        ra=coord.ra + ra_offset.to(u.deg),
                        dec=coord.dec + dec_offset.to(u.deg),
                        frame='fk5')
                    table_chunk = self.single_query(ind_coord, size=single_chunk_size, 
                        verbose=verbose, maxrec=2000, equinox=equinox)
                    if i == 0 and j == 0:
                        table = table_chunk
                    else:
                        table = vstack([table, table_chunk])
            table.sort('phot_rp_mean_mag')
            if not maxrec is None:
                table = table[0:maxrec]
           

        return table
    def login(self, verbose=False):
        if self.credentials[0] != 'NONE' and self.credentials[1] != 'NONE':
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                Gaia.login(user=self.credentials[0], password=self.credentials[1],verbose=verbose)    
    def single_query(self, coord, size = 3*u.arcmin, 
                verbose=False, maxrec=2000, equinox=None):
        query  = self.construct_gaia_query(coord, size=size, maxrec=maxrec)
        Gaia.MAIN_GAIA_TABLE = self.MAIN_GAIA_TABLE
        Gaia.ROW_LIMIT = maxrec
        job = Gaia.launch_job(query,dump_to_file=False,verbose=verbose)
        table = job.get_results()
        return table
    
    def construct_gaia_query(self, coord, size = 1* u.arcmin, maxrec=2000):
        '''Construct a Gaia query for a given coordinate and size. Exclude galaxies and order by magnitude'''
        ramin = coord.ra.deg - size.to(u.deg).value/(2.*np.cos(coord.dec.to(u.rad).value))
        ramax = coord.ra.deg + size.to(u.deg).value/(2.*np.cos(coord.dec.to(u.rad).value))
        decmin = coord.dec.deg - size.to(u.deg).value/2.
        decmax = coord.dec.deg + size.to(u.deg).value/2.
        if maxrec > 2000:
            print(f"Warning: maxrec {maxrec} exceeds the recommended limit of 2000.")
            maxrec = 2000
        query = f"SELECT TOP {int(maxrec)} ra,dec,phot_rp_mean_mag,in_galaxy_candidates FROM gaiadr3.gaia_source"
        query += f" WHERE ra BETWEEN {ramin} AND {ramax} AND dec BETWEEN {decmin} AND {decmax}"
        query += f" AND phot_rp_mean_mag IS NOT NULL AND in_galaxy_candidates = 'False'"
        query += f" ORDER BY phot_rp_mean_mag "
        return query

class LW_NedQuery:
    # Lightweight NED query class for cone searches because astropquery hasn't updated the api's 
    # which leads to continuous failures when using astropquery.
    def __init__(self):
        
        self.coords = "<RA DEC>"
        self.radi = 1.*u.arcmin
        self.verbose = False
        self.equinox = 'J2000'
        self.z_constraint = "Unconstrained"

        self.payload = {
            "MAXREC": 0,
            "CSYS": "Equatorial",
            "EQUINOX": 'J2000',
            "RA": "<RA>",
            "DEC": "<DEC>",
            "RADIUS": self.radi.to(u.arcmin).value,
            "Z_CONSTRAINT": self.z_constraint
        }
        self.url =  "https://ned.ipac.caltech.edu/NED::API"
        '''Copyright 2019-2026 Caltech, Support from NASA is acknowledged. This information is taken from a collection created and curated by the NASA/IPAC Extragalactic Database (NED), operated by the California Institute of Technology (Caltech). If your research benefits from the use of NED, the following acknowledgement in your paper would be appreciated: "This research has made use of the NASA/IPAC Extragalactic Database (NED), which is funded by the National Aeronautics and Space Administration and operated by the California Institute of Technology." See: https://ned.ipac.caltech.edu/Documents/Overview/Acknowledgments'''
        self.timeout = 600

    def query_region(self, coord,radius = None, verbose=None, maxrec = None, 
            equinox = None, legacy =True): 
        self.url_addition = "ConeSearchByPosition"
        self.coords = coord
        if radius is not None:
            self.radi = radius
        if verbose is not None:
            self.verbose = verbose
        if maxrec is not None:
            self.payload["MAXREC"] = maxrec
        if equinox is not None:
            self.equinox = equinox
            self.payload["EQUINOX"] = self.equinox
        self.payload["RA"] = self.coords.ra.to_string(unit=u.hour, sep='hms')
        self.payload["DEC"] = self.coords.dec.to_string(unit=u.degree, sep='dms')
        self.payload["RADIUS"] = self.radi.to(u.arcmin).value
        request = self.make_url()
        temp = urllib.request.urlopen(request)
        table = self.xml_to_table(temp)
        if legacy:
            table = self.rename_columns(table)

        #with open(f"test_ned_Rad{self.radi.to(u.arcmin).value}.html", "w") as f:
        #    f.write(temp.read().decode('utf-8'))
        return table
    
    def query_object(self, name, verbose=None, maxrec=None, equinox=None, legacy=True):
        self.url_addition = f'OverviewOfObject'
        
        self.payload["TARGET"] = name
        if verbose is not None:
            self.verbose = verbose
        request = self.make_url()
        temp = urllib.request.urlopen(request)
        table = self.xml_to_table(temp)
        if legacy:
            table = self.rename_columns(table)
        return table
    
    def rename_columns(self, table):
        # Implement the logic to rename columns for legacy support
        column_names = table.colnames
        # Example renaming logic for legacy support
        ['row', 'prefname', 'equ_j2000_lon_s', 'equ_j2000_lat_s', 'ra', 'dec', 'z', 'velocity', 'zflag', 'ptype', 'emtype', 'sep', 'n_ref', 'n_notes', 'n_gphot', 'n_posd', 'n_zdf', 'n_ddf', 'n_dist', 'n_class', 'n_images', 'n_spectra']
        ['No.', 'Object Name', 'RA', 'DEC', 'Type', 'Velocity', 'Redshift', 'Redshift Flag', 'Magnitude and Filter', 'Separation', 'References', 'Notes', 'Photometry Points', 'Positions', 'Redshift Points', 'Diameter Points', 'Associations']
        rename_map = {
            "ra": "RA",
            "dec": "DEC",
            "row": "No.",
            "prefname": "Object Name",
            "velocity": "Velocity",
            "z" : "Redshift",
            "zflag": "Redshift Flag",
            "ptype": "Type",
            "sep": "Separation",
            "n_ref": "References",
            "n_notes": "Notes",
            "n_gphot": "Photometry Points",
            "n_posd": "Positions",
            "n_zdf": "Redshift Points",
            "n_ddf": "Diameter Points",
            "n_dist": "Associations"

        }
        for old_name, new_name in rename_map.items():
            if old_name in column_names:
                table.rename_column(old_name, new_name)

        return table

    def xml_to_table(self, xml_data):
        from astropy.io import votable
        from io import BytesIO
        votable_file = BytesIO(xml_data.read())
        table = votable.parse(votable_file).get_first_table().to_table()
        return table
    
    def make_url(self):
        payload_url = f'{self.url}/{self.url_addition}?'
        for key, value in self.payload.items():
            if key == 'MAXREC' and value == 0:
                continue
            payload_url += f'{urlencode({key: value})}&'
        payload_url = payload_url.rstrip('&')
        if self.verbose:
            print(f' Querying this url {payload_url}')
        return payload_url





class NullGate:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False
