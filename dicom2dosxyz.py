###################################################################
# This code is part of PYCOM, The Python DICOM Processing Toolkit.
# PYCOM is distributed open source under GPL v3 license.
###################################################################
#
# Coordinate transform from DICOM to DOSXYZnrc
#    Input: gantry angle, couch angle and collimator angle from DICOM
#    Output: theta, phi, phicol for DOSXYZnrc
#
#########################################################
# Developed by Lixin Zhan (lixinzhan AT gmail DOT com). #
#########################################################

import numpy as np

# Constants for angle conversions
DEG_TO_RAD = np.pi / 180.0
RAD_TO_DEG = 180.0 / np.pi
HALF_PI = np.pi / 2.0
TWO_PI = 2.0 * np.pi

def dcm2dosxyz(AngleGantry, AngleCouch, AngleCollimator, Method='Zhan'):
    # Convert angles to radians
    gamma = AngleGantry * DEG_TO_RAD
    col = AngleCollimator * DEG_TO_RAD
    rho = AngleCouch * DEG_TO_RAD

    # Distort Couch and Gantry angles slightly to avoid special cases
    if AngleCouch in (90.0, 270.0) and AngleGantry in (90.0, 270.0):
        rho *= 0.999999
        gamma *= 0.999999

    # Precompute trigonometric values
    sin_gamma = np.sin(gamma)
    cos_gamma = np.cos(gamma)
    sin_rho = np.sin(rho)
    cos_rho = np.cos(rho)
    sin_col = np.sin(col)
    cos_col = np.cos(col)

    if Method == 'Zhan':  # Mapping couch rotation to collimator plane
        sgsr = sin_gamma * sin_rho
        sgcr = sin_gamma * cos_rho

        theta = np.arccos(-sgsr)
        phi = np.arctan2(-cos_gamma, sgcr)
        couch_angle_2_coll_plane = np.arctan2(-sin_rho * cos_gamma, cos_rho)
        phicol = (col - HALF_PI) + couch_angle_2_coll_plane
        # Coordinate transformation for BEAMnrc generated phsp to DOSXYZnrc
        phicol = np.pi - phicol

    elif Method == 'Bush':  # Bush, Australas. Phys. Eng. Sci. (2010) 33:351
        rho = TWO_PI - rho
        sin_rho = np.sin(rho)  # Recalculate after rho change
        cos_rho = np.cos(rho)

        sgsr = sin_gamma * sin_rho
        sgcr = sin_gamma * cos_rho
        cgsr = cos_gamma * sin_rho

        theta = np.arctan2(np.sqrt(1.0 - sgsr**2), sgsr)
        phi = np.arctan2(-cos_gamma, sgcr)
        col_adj = col - HALF_PI
        phicol = np.arctan2(
            (-cgsr * np.cos(col_adj) - cos_rho * np.sin(col_adj)),
            (cgsr * np.sin(col_adj) - cos_rho * np.cos(col_adj))
        )
        phicol = TWO_PI - phicol

    elif Method == 'Thebaut':  # Thebaut, Phys. Med. Biol. (2006) 51:N441
        rho = TWO_PI - rho
        sin_rho = np.sin(rho)  # Recalculate after rho change
        cos_rho = np.cos(rho)

        sgsr = sin_gamma * sin_rho
        sgcr = sin_gamma * cos_rho
        cgcr = cos_gamma * cos_rho

        theta = np.arccos(sgsr)
        phi = np.arctan2(-cos_gamma, sgcr)
        col_adj = col - HALF_PI

        cos_phicol = (
            cos_col * cgcr * np.sin(phi) -
            sin_col * sin_rho * np.sin(phi) -
            cos_col * sin_gamma * np.cos(phi)
        )
        # Clamp to avoid numerical issues
        cos_phicol = np.clip(cos_phicol, -1.0, 1.0)
        phicol = np.arccos(cos_phicol)  # Returns value within [0, pi]

        direct = (
            cos_col * sin_rho * np.sin(phi) +
            sin_col * (cgcr * np.sin(phi) - sin_gamma * np.cos(phi))
        )
        if direct > 0:
            phicol = -np.abs(phicol)
        else:
            phicol = np.abs(phicol)

    else:
        raise ValueError('Incorrect c.s. transformation method!')

    ###########################################################
    # Range: theta - (0,180), phi - (0,360), phicol - (0,360) #
    ###########################################################

    return (
        theta * RAD_TO_DEG,
        np.mod(phi * RAD_TO_DEG, 360.0),
        np.mod(phicol * RAD_TO_DEG, 360.0)
    )

