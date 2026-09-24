import sys
import numpy as np
from math import log10

ch_constants = {"Martell": ((4212., 4242., 1.), (4297., 4317., 1.), (4330., 4375., 1)),
                "Harbeck_DDO_red": ((4285., 4315, 30.), (4240., 4280., 40.), (4479., 4555., 76.)),
                "DDO": ((4220.5, 4293.5, 1.), (4479., 4555., 1), None),
                "Harbeck_Weighted": ((4285., 4315., 30.), (4240., 4280., 40.), (4390., 4460., 70.)),
                "Harbeck_Modified": ((4280., 4300., 20.), (4240., 4280., 40.), (4390., 4460., 70.)),
                "Harbeck": ((4285., 4315., 1.), (4240., 4280., 1.), (4390., 4460., 1.)),
                "Lee": ((4270., 4320., 1.), (4230., 4260., 1.), (4390., 4420., 1)),
                "Smith": ((4285., 4315., 1.), (4180., 4250., 1.), (4380., 4430., 1.))
                }

def calc_bandstrength_from_file(filename: str, choice: str = "Harbeck", doppler: float = 0.0):
    """
    Loads wave and value from a given file before running the regular calc_bandstrength
    """
    data = np.loadtxt(filename)
    wave = data[ : , 0]
    value = data[ : , 1]
    return calc_bandstrength(wave, value, choice=choice, doppler=doppler)

def calc_bandstrength(wave: np.ndarray, value: np.ndarray, choice: str = "Harbeck", doppler: float = 0.0):
    """
    Calculates the band strength of given spectra.
    """
    try:
        constants = ch_constants[choice]
        dopl = (1+(doppler/300000.))
        # Mask wave and value with proper start and end points to
        # calculate fch, cbch, and crch, then use those to calculate ch
        tw, tv = apply_exact_mask(wave, value, constants[0][0]*dopl, constants[0][1]*dopl)
        fch = np.trapz(tv, tw)/constants[0][2]
        tw, tv = apply_exact_mask(wave, value, constants[1][0]*dopl, constants[1][1]*dopl)
        cbch = np.trapz(tv, tw)/constants[1][2]
        crch = 0
        if constants[2] != None:
            tw, tv = apply_exact_mask(wave, value, constants[2][0]*dopl, constants[2][1]*dopl)
            crch = np.trapz(tv, tw)/constants[2][2]
        denom = (0.5*cbch + 0.5*crch)
        if choice == "Martell":
            denom *= 2
        ch = -2.5 * log10(fch/denom)

        return ch
    
    except ValueError:
        print("Given spectra does not fit within bands, skipping band strength calculation")
        return 0
    except ZeroDivisionError:
        print("cbch is equal to -1*crch, unable to calculate band strength (cannot divide by zero)")
        return 0
    except KeyError:
        print("Choose from one of these for ch next time")
        print("Martell")
        print("Harbeck_DDO_red")
        print("DDO")
        print("Harbeck_Weighted")
        print("Harbeck Modified")
        print("Harbeck")
        print("Lee")
        print("Smith")
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
	print(calc_bandstrength_from_file("input_files/sample_spectrum/85-scaled.txt"))


main()