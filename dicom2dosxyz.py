###################################################################
# This code is part of PYCOM, The Python DICOM Processing Toolkit.
# PYCOM is distributed open source under GPL v3 license.
###################################################################
#
# Coordiante transform from DICOM to DOSXYZnrc
#    Input: gantry angle, couch angle and collimator angle from DICOM
#    Output: theta, phi, phicol for DOSXYZnrc
# 
#########################################################
# Developed by Lixin Zhan (lixinzhan AT gmail DOT com). #
#########################################################

import numpy as np

def dcm2dosxyz(AngleGantry, AngleCouch, AngleCollimator, Method='Zhan'):
    gamma = AngleGantry*np.pi/180.0
    col   = AngleCollimator*np.pi/180.0
    rho   = AngleCouch*np.pi/180
    # distort Couch and Gantry angles slightly to avoid the very special cases.
    if AngleCouch in (90.0,270.0) and AngleGantry in (90.0,270.0):
        rho   = rho   * 0.999999
        gamma = gamma * 0.999999
        
    if Method=='Zhan':  # Mapping couch rotation to collimator plane
        sgsr = np.sin(gamma)*np.sin(rho)
        sgcr = np.sin(gamma)*np.cos(rho)
        
        theta = np.arccos(-sgsr)
        phi = np.arctan2(-np.cos(gamma),sgcr)
        CouchAngle2CollPlane = np.arctan2(-np.sin(rho)*np.cos(gamma),np.cos(rho))
        phicol = (col-np.pi/2) + CouchAngle2CollPlane
        # coord. trans. for BEAMnrc generated phsp to DOSXYZnrc.
        phicol = np.pi - phicol
    elif Method=='Bush':  # Bush, Australas. Phys. Eng. Sci. (2010) 33:351
        rho = 2*np.pi-rho
        sgsr = np.sin(gamma)*np.sin(rho)
        sgcr = np.sin(gamma)*np.cos(rho)
        cgsr = np.cos(gamma)*np.sin(rho)
        
        theta = np.arctan2(np.sqrt(1.0-sgsr**2),sgsr)
        phi = np.arctan2(-np.cos(gamma),sgcr)
        col = col - np.pi/2
        phicol = np.arctan2( (-cgsr*np.cos(col)-np.cos(rho)*np.sin(col)),
            (cgsr*np.sin(col)-np.cos(rho)*np.cos(col)) )
        phicol = 2*np.pi-phicol
    elif Method=='Thebaut':  # Thebaut, Phys. Med. Biol. (2006) 51:N441
        rho = 2*np.pi-rho
        sgsr = np.sin(gamma)*np.sin(rho)
        sgcr = np.sin(gamma)*np.cos(rho)
        cgcr = np.cos(gamma)*np.cos(rho)
        
        theta = np.arccos(sgsr)
        phi = np.arctan2(-np.cos(gamma),sgcr)
        col = col - np.pi/2
        cos_phicol = np.cos(col)*cgcr*np.sin(phi) - \
                np.sin(col)*np.sin(rho)*np.sin(phi) - \
                np.cos(col)*np.sin(gamma)*np.cos(phi)
        if cos_phicol>1.0: # for possible binary express problem.
            cos_phicol = 1.0
        elif cos_phicol<-1.0:
            cos_phicol = -1.0
        phicol = np.arccos(cos_phicol) # return value within [0,pi]
        direct = np.cos(col)*np.sin(rho)*np.sin(phi) + \
                np.sin(col)*( cgcr*np.sin(phi)-np.sin(gamma)*np.cos(phi) )        
        if direct>0:
            phicol = -np.abs(phicol)
        else:
            phicol = np.abs(phicol)
                    
    else:
        raise ValueError, 'Incorrect c.s. transformation method!'
    
    ###########################################################
    # Range: theta - (0,180), phi - (0,360), phicol - (0,360) #
    ###########################################################

    return (theta*180/np.pi, np.mod(phi*180/np.pi,360), \
            np.mod(phicol*180/np.pi,360))

