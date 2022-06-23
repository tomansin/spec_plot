import numpy as np
import utils as ut

class spectrum:
    
    def __init__(self):        
        self.header = None
        self.wavelenght = None
        self.flux = None


    # def open(self, path):
    #     header, wavelenght, flux = ut.read_fits(path)
    #     self.header = header
    #     self.wavelenght = wavelenght
    #     self.flux = flux

    def plot(self, wv_range = None, fl_range = None, fig_name = None):
        ut.spec_plot(self.wavelenght, self.flux, self.header, wv_range, fl_range, fig_name)

def open_spec(path):
    output = spectrum()
    header, wavelenght, flux = ut.read_fits(path)
    output.header = header
    output.wavelenght = wavelenght
    output.flux = flux
    return output

# s = spectrum()
# s.open('./redsgs/fit/a12a.fit')
