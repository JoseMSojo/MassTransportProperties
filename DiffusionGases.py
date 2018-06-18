# -*- coding: utf-8 -*-


def DiffGas_AB(A, B, T, P):
    """
    Function to compute difusivities in gas mixtures
    
    Correlation from Fuller, Schettler and Giddings 1966
    
    Database of coeficients from Fuller, Ensley and Giddings 1969
    
    E. N. Fuller, P. D. Schettler, and J. C. Giddings, “A New Method for Prediction of Binary
    Gas-Phase Diffusion Coefficients,” Industrial and Engineering Chemistry, 58, 1966, pp. 19-27.

    E. N. Fuller, K. Ensley, and J. C. Giddings, “Diffusion of Halogenated Hydrocarbons in Helium.
    The Effect of Structure on Collision Cross Sections,” Journal of Physical Chemistry, 73, 1969,
    pp. 3679-3685.
    
    Input
    -----
    A and B (string): Diffused gases 
            
    T (float): Temperature in [K]
    
    P (float): Pressure in [Pa]
    
    Output
    ------
    D_AB (float): Diffusion coefficient of A in B [m^2/s]
    
    
    ------
    Jose Manuel Sojo Gordillo
    Created on Tue May 8 2018
    
    """
    # Atomic weights [g/mol]
    M = {'H2':    2.01593,
         'D2':    4.028204,
         'He':    4.002598,
         'N2':   28.01403,
         'H2O':  18.01534,
         'Kr':   83.800,
         'Xe':  131.300,
         'Air':  28.963,
         'Ar':   39.948,
         'SF6': 146.05,
         'O2':   32.000,
         'CO2':  44.010,
         'CO':   28.020,
         'N2O':  44.010,
         'NH3':  17.030,
         'SO2':  64.0638,
         'Cl2':  70.90,
         'Br2':  159.808}
    
    # Atomic diffusion volume [cm^3]
    SIGMA_v = {'H2':    6.12,
               'D2':    6.84,
               'He':    2.67,
               'N2':   18.5,
               'H2O':  13.1,
               'Kr':   24.5,
               'Xe':   32.7,
               'Air':  19.7,
               'Ar':   16.2,
               'SF6':  71.3,
               'O2':   16.3,
               'CO2':  26.9,
               'CO':   18.0,
               'N2O':  35.9,
               'NH3':  20.7,
               'SO2':  41.8,
               'Cl2':  38.4,
               'Br2':  69.0}
    
    K = 0.0101325  # [m^2/s] if P in [Pa]
    #K = 1e-7        # [m^2/s] if P in [atm]
    
    try:
        return K*( T**1.75 * (1/M[A] + 1/M[B])**(1/2) ) / (P*(SIGMA_v[A]**(1/3) + SIGMA_v[B]**(1/3))**2)
    except KeyError:
        raise KeyError('Unsuported gas selected')

# Test the function
if (__name__ == '__main__'):
    
    D = DiffGas_AB('HO','Air', 60+273, 1.0*101325)
    