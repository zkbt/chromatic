import h5py as h5

try:
    import cPickle as pickle
except:
    import _pickle as pickle
    
import numpy as np


def loadevent(filename, load=[], loadfilename=None): #from Eureka source code
    """
    Loads an event stored in .dat and .h5 files.

    Parameters
    ----------
    filename : String
               The string contains the name of the event file.
    load     : String tuple
               The elements of this tuple contain the parameters to read.
               We usually use the values: 'data', 'uncd', 'head', 'bdmskd',
               'brmskd' or 'mask'.

    Notes
    -----
    The input filename should not have the .dat nor the .h5 extentions.

    Returns
    -------
    This function return an Event instance.

    Examples
    --------
    See package example.

    Revisions
    ---------
    2010-07-10  patricio  Added documentation.     pcubillos@fulbrightmail.org
    """

    from astropy.io import fits
    handle = open(filename + '.dat', 'rb')
    event = pickle.load(handle, encoding='latin1')
    handle.close()

    if loadfilename == None:
        loadfilename = filename

    if load != []:
        handle = h5.File(loadfilename + '.h5', 'r')
        for param in load:
            exec('event.' + param + ' = handle["' + param + '"][:]')
      
    # calibration data:
            if event.havecalaor:
                exec('event.pre'  + param + ' = handle["pre'  + param + '"][:]')
                exec('event.post' + param + ' = handle["post' + param + '"][:]')

        handle.close()

    return event

def eureadka(filename,metaname):
    
    event = loadevent(filename,load=[])
    
    meta = loadevent(metaname,load=[])
       
    
    t = event.int_times['int_mid_BJD_TDB'] #UT time for the midpoint of exposure
    f = event.optspec[:, meta.xwindow[0]:meta.xwindow[1]]
    e = event.opterr[:, meta.xwindow[0]:meta.xwindow[1]]    
    w = event.subwave[meta.src_ypos, meta.xwindow[0]:meta.xwindow[1]]
    #this assumes the wavelength axis is the same for all exposures
    
    timelike = {}
    timelike['time'] = t
    
    wavelike = {}
    wavelike['wavelength'] = w
    
    fluxlike = {}
    fluxlike['flux'] = f
    fluxlike['error'] = e
    
    return timelike, wavelike, fluxlike
    