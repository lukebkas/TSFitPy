import sys
import numpy as np
from math import log10

def calc_bandstrength_from_file(filename: str, doppler: float = 0.0):
    """
    Loads wave and value from a given file before running the regular calc_bandstrength
    """
    data = np.loadtxt(filename)
    wave = data[ : , 0]
    value = data[ : , 1]
    return calc_bandstrength(wave, value, doppler)

def calc_bandstrength(wave: np.ndarray, value: np.ndarray, doppler: float = 0.0):
    """
    Calculates the band strength of given spectra
    """
    try:
        dopl = (1+(doppler/300000.))
        # Mask wave and value with proper start and end points to
        # calculate fch, cbch, and crch, then use those to calculate ch
        tw, tv = apply_exact_mask(wave, value, 4285.*dopl, 4315*dopl)
        fch = np.trapz(tv, tw)
        tw, tv = apply_exact_mask(wave, value, 4240.*dopl, 4280*dopl)
        cbch = np.trapz(tv, tw)
        tw, tv = apply_exact_mask(wave, value, 4390.*dopl, 4460.*dopl)
        crch = np.trapz(tv, tw)
        ch = -2.5 * log10(fch/(0.5*cbch + 0.5*crch))

        return ch
    
    except ValueError:
        print("Given spectra does not fit within bands, skipping band strength calculation")
        return 0
    except ZeroDivisionError:
        print("cbch is equal to -1*crch, unable to calculate band strength (cannot divide by zero)")
        return 0

def apply_exact_mask(x: np.ndarray, y: np.ndarray, start: float, end: float):
    """
    Removes all (x,y) paris in x and y arrarys where x does not fall between start and end.
    Interoplates the exact value of the given start and end in between existing points.
    """
    # Creates a mask between values and extends by one in each direction
    # so that interpolation is possible
    mask = np.logical_and(x >= start, x <= end)
    if mask[0] != True:
        mask = np.roll(mask,-1)+mask
    if mask[-1] != True:
        mask = np.roll(mask,1)+mask
    newx = x[mask]
    newy = y[mask]
    # Interpolates the actual values of y where x = start and end, then
    # replace endpoints of arrys with those values
    ys = np.interp([start,end],newx[0:2],newy[0:2])
    newx[0] = start
    newx[-1] = end
    newy[0] = ys[0]
    newy[-1] = ys[1]
    return newx, newy

def main():
	print(calc_bandstrength_from_file("input_files/sample_spectrum/85-scaled.txt", 0))


main()
