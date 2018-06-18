# -*- coding: utf-8 -*-


def MetalDiff(A, Q, T, W = True):
    """
    Function to calculate diffusivities of several isotopes in metallic alloys
    
    If the isotope requested is not in the database, Grham's aproximation law is
    used in order to estimate it from the H values (Ref 7: [Grahams1967])
    
    Available Alloys: SS304, SS316, Alpha-Fe, 75Pd-25Ag, EUROFER97, OPTIFER-IVb    
    Available Isotopes: H, D, T, He
    
    References
    ----------
    Ref 1 [Grant1987]: D. M. Grant, D. L. Cummings, and D. A. Blackburn, 
    “Hydrogen in 304 Steel - Diffusion, Permeation and Surface-Reaction” 
    Journal of Nuclear Materials, vol. 149, no. 2, pp. 180–191, 1987.
    
    Ref 2 [Lee2014]: S. K. Lee, S. H. Yun, H. G. Joo, and S. J. Noh, 
    “Deuterium transport and isotope effects in type 316L stainless steel at high
    temperatures for nuclear fusion and nuclear hydrogen technology applications”
    Current Applied Physics, vol. 14, no. 10, pp. 1385–1388, 2014.
    
    Ref 3 [Forcey1997]: K. S. Forcey, I. Iordanova, and M. Yaneva, 
    “The diffusivity and solubility of deuterium in a high chromium martensitic steel”
    Journal of Nuclear Materials, vol. 240,no. 2, pp. 118–123, 1997.
    
    Ref 4 [Serra1998]: E. Serra, G. Benamati, and O. Ogorodnikova, 
    “Hydrogen isotopes transport parameters in fusion reactor materials” 
    Journal of Nuclear Materials, vol. 255, no. 2-3, pp. 105–115, 1998.
    
    Ref 5 [Esteban2007]: G. A. Esteban, A. Pena, I. Urra, F. Legarda, and B. Riccardi, 
    “Hydrogen transport and trapping in EUROFER’97” 
    Journal of Nuclear Materials, vol. 367-370 A, no. SPEC.ISS., pp. 473–477, 2007.
    
    Ref 6 [Esteban2000]: G. A. Esteban, A. Perujo, K. Douglas, and L. A. Sedano, 
    “Tritium diffusive transport parameters and trapping effects in the reduced 
    activating martensitic steel OPTIFER-IVb” 
    Journal of Nuclear Materials, vol. 281, no. 1, pp. 34–41, 2000.
    
    Ref 7: [Grahams1967]: E. A. Mason and B. Kronstadt, 
    “Graham’s Laws of diffusion and effusion” 
    Journal of Chemical Education, vol. 44, no. 12, p. 740, 1967.
    
    Inputs
    ------
    A (string): Metal alloy
    
    Q (string): Diffused isotope
            
    T (float): Temperature in [K]
      
    Output
    ------
    D (float): Diffusion coefficient of isotope Q in alloy A [m^2/s]
    
    ------
    Jose Manuel Sojo Gordillo
    Created on Mon Jun 18 2018
    """
    from math import exp, sqrt
    R = 8.314  # [J/mol-K] Ideal gas constant
    
    # Atomic weights [g/mol]
    M = {'H':    1.00794,
         'D':    2.01410,
         'T':    3.01605,
         'He':   4.00260}
    
    # Pre-Exponential factors [m^2/s]
    D0 = {'SS304': {'H': 1.22e-6},   # Ref 1: [Grant1987]
                  
          'SS316': {'H': 1.24e-6,
                    'D': 1.38e-6},   # Ref 2: [Lee2014]
          
          'Alpha-Fe': {'H': 3.87e-8}, # Ref 3: [Forcey1997]
          
          '75Pd-25Ag': {'H': 3.07e-7,
                        'D': 1.87e-7}, # Ref 4: [Serra1998]
          
          'EUROFER97': {'H': 4.57e-7,
                        'D': 1.50e-7}, # Ref 5: [Esteban2007]
        
          'OPTIFER-IVb': {'H': 5.49e-8,
                          'D': 4.61e-8,
                          'T': 4.17e-8}} # Ref 6: [Esteban2000]
            
            
    # Activation energies [kJ/mol]
    Ea = {'SS304': {'H': 54.85},   # Ref 1: [Grant1987]

          'SS316': {'H': 55.10,
                    'D': 57.50},   # Ref 2: [Lee2014]
                    
          'Alpha-Fe': {'H': 4.50},  # Ref 3: [Forcey1997]
          
          '75Pd-25Ag': {'H': 25.90,
                        'D': 24.69}, # Ref 4: [Serra1998]
          
          'EUROFER97': {'H': 22.30,
                        'D': 14.50}, # Ref 5: [Esteban2007]
                    
          'OPTIFER-IVb': {'H': 10.60,
                          'D': 11.30,
                          'T': 12.00}} # Ref 6: [Esteban2000]
    
    # Temperature range [K]
    # Lower
    T_min = {'SS304': 645,   # Ref 1: [Grant1987]
                  
             'SS316': 623,    # Ref 2: [Lee2014]
          
             'Alpha-Fe': 573,    # Ref 3: [Forcey1997]
          
             '75Pd-25Ag': 323,  # Ref 4: [Serra1998]
          
             'EUROFER97': 373, # Ref 5: [Esteban2007]
                              
             'OPTIFER-IVb': 423} # Ref 6: [Esteban2000]
          
    # Upper
    T_max = {'SS304': 945,   # Ref 1: [Grant1987]
                  
             'SS316': 1123,    # Ref 2: [Lee2014]
          
            'Alpha-Fe': 873,   # Ref 3: [Forcey1997]
          
            '75Pd-25Ag': 773,  # Ref 4: [Serra1998]
          
            'EUROFER97': 723, # Ref 5: [Esteban2007]
                        
            'OPTIFER-IVb': 892} # Ref 6: [Esteban2000]
        
    try:
        if((T_min[A] > T)or(T_max[A] < T))and(W):
            print('Warning: Temperature used out of the correlation range [',T_min[A],'-',T_max[A],']K')
    except KeyError:
        raise KeyError('Unsuported Alloy selected')
    
    
    try:
        return D0[A][Q]*exp(-Ea[A][Q]*1000/(R*T)) 
    except:
        None
        
    try:   # If Q requested is not tabulated, compute D using Grahm's aproximation law
        D_H = D0[A]['H']*exp(-Ea[A]['H']*1000/(R*T))   
        return D_H*sqrt(M['H']/M[Q])                    # Ref 7: [Grahams1967]
        
    except KeyError:
        raise KeyError('Unsuported Isotope selected')
        
# Test the function
if (__name__ == '__main__'):
    
    D = MetalDiff('SS316','H', 700, W = False)      
