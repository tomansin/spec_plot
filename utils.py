import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import specutils as sp
from astropy.time import Time


def read_fits(file_name):
    with fits.open(file_name) as hdul:
        header = hdul[0].header
        header['CUNIT1'] = 'Angstrom'
        
        if header['CTYPE1'] == 'LINEAR' or header['CTYPE1'] == 'wavelengthGTH':
            
            start, step, num, ref = (hdul[0].header['CRVAL1'],
                                     hdul[0].header['CDELT1'],
                                     hdul[0].header['NAXIS1'],
                                     hdul[0].header['CRPIX1'])
            
            wavelength = (np.arange(num, dtype=float) + 1 - ref) * step + start
            #wavelength = np.linspace(start=start, stop=start+num*step, num=num)
            
            # Para logitudes de onda en escala logaritmica
            if 'DC-FLAG' in header and header['DC-FLAG'] == 1:
                wavelength = 10.0**wavelength
            
            flux = (hdul[0].data) 
#            data = (hdul[0].data) * u.adu
            
#            spec_final = sp.Spectrum1D(flux=data, spectral_axis=wavelength*u.AA, meta=header)
        
        """ elif header['CTYPE1'] == 'MULTISPE':
            
            #Cantidad de ordenes
            nspec = header['NAXIS2']
            
            # Los parametros estan en los keywords WAT2...
            WAT2 = str()
            for linea in header['WAT2*']:
                v = header[linea]
                v = v + (" " * (68 - len(v)))
                WAT2 = WAT2 + v
            
            # Separamos las lineas
            specstr = [] 
            for i in range(nspec):
                sname = 'spec' + str(i+1)
                p1 = WAT2.find(sname)
                p2 = WAT2.find('"', p1)
                p3 = WAT2.find('"', p2+1)
                specstr.append(WAT2[p2 +1 :p3])
                
            # Separamos los parametros
            wparms = np.zeros((nspec, 9))
            for i in range(nspec):
                wparms[i,:] = np.asanyarray(specstr[i].split())
                
            # Armo el header
            header['NAXIS'] = 1
            header['DATE'] = Time(Time.now(), format='fits').value
            header['IRAF-TLM'] = Time(Time.now(), format='fits').value
            header['WCSDIM'] = 1
            
            header['WAT0_001'] = 'system=equispec'
            header['WAT1_001'] = 'wtype=linear label=wavelengthgth units=angstroms units_display=angstrom'
            header['CTYPE1']  = 'LINEAR'            
            header['DC-FLAG'] = 0            
            header['CRPIX1'] = 1.
        
            if 'NAXIS2' in header: del header['NAXIS2'] 
            if 'LTM2_2' in header: del header['LTM2_2'] 
            if 'WAT2_' in header:  del header['WAT2_*'] 
            if 'CTYPE2' in header: del header['CTYPE2'] 
            if 'CDELT2' in header: del header['CDELT2'] 
            if 'CD2_2' in header:  del header['CD2_2']  
            
            # Armar los espectros
            spec_final = []
            for i in range(nspec):                               
                
                star, step, num = wparms[i,3:6]
                apnum1 = str(wparms[i,0])+' '+str(wparms[i,1])+' '+str(wparms[i,7])+' '+str(wparms[i,8])
                
                # Agrego coefs                
                header['CDELT1']  = step
                header['CD1_1']  = step
                header['CRVAL1'] = star
                header['APNUM1'] = apnum1
                
                wavelength = np.linspace(star, star+num*step,int(num))
                # Para longitud de onda en escala logar
                if wparms[i,2] == 1:
                    wavelength = 10.0**wavelength
                    header['DC-FLAG'] = 1
                    
                # Si tiene más de una banda solo me quedo con la primera
                # Más adelante implementar levantar todas
                if 'NAXIS3' in header:
                    data = (hdul[0].data[0])[i] * u.adu
                else:
                    data = (hdul[0].data)[i] * u.adu
                    
                meta = header.copy(strip=True)
                spec = sp.Spectrum1D(flux=data, spectral_axis=wavelength*u.AA, meta=meta)
                spec_final.append(spec)
        
        else:
            raise ValueError('El keyword CTYPE1={} no se encuentra entre las opciones (Linear, wavelengthGTH y MULTISPE) que puede manejar esta función.'.format(header['CTYPE1'])) """
            
    return header, wavelength, flux
#    return spec_final


def spec_plot(wavelength, flux, header, wv_range = None, fl_range = None, fig_name = None):
        
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.plot(wavelength, flux)

    if wv_range is not None:
        shift_wv = (wv_range[1] - wv_range[0]) * .05
        ax.set_xlim(wv_range[0] - shift_wv, wv_range[1] + shift_wv)
    
    if fl_range is not None:
        shift_f = (fl_range[1] - fl_range[0]) * .05
        ax.set_ylim(fl_range[0] - shift_f, fl_range[1] + shift_f)
    
    ax.set_title(header['OBJECT'])
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Flux')
 
    if fig_name is not None:
        plt.savefig(fig_name)
    plt.show() 

# header, wv, flux = read_fits('./redsgs/fit/a12a.fit')
# spec_plot(wv, flux, header, wv_range=(7500, 8000), fl_range=(.6, .8), fig_name='cuidadoconlamesa')