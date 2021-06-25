'''tool for reading a Eureka
output file and making a
Rainbow out of it'''

def eureadka(filename):
    
    eureka_data = np.load(filename, allow_pickle = True)
    
    eureka_data_times = 
    eureka_data_wave = 
    eureka_data_flux = 
    eureka_data_err = 
    
    timelike = {}
    timelike['time'] = eureka_data_times
    
    wavelike = {}
    wavelike['wavelength'] = eureka_data_wave
    
    fluxlike = {}
    fluxlike['flux'] = eureka_data_flux
    fluxlike['error'] = eureka_data_err
    
    return timelike, wavelike, fluxlike
    