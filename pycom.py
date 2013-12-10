###################################################################
# This code is part of PYCOM, The Python DICOM Processing Toolkit.
# PYCOM is distributed open source under GPL v3 license. 
###################################################################
#
# Python DICOM Plan interprator
#
# This script automatically generates BEAMnrc and DOSXYZnrc input files
# from a DICOM Plan (DICOM RP) file. The generated *.egsinp files and
# their necessary support files will be in a beamnrc folder for BEAMnrc
# and a dosxyznrc folder for DOSXYZnrc. 
# This script is generating input files based on DOSXYZnrc source 21,
# see Lobo and Popescu Phys. Med. Biol. 55:4431-4443, 2010 for detail.
# Hence, the input files generated are only compatible with 
# EGSnrc/BEAMnrc version 2011 and above.
#
# A typical clinical RT plan usually contains many fields. This script
# combines all fields into a single 'ARC' field. 'ARC' is in quation marks
# is due to the fact that it is actually not a real arc and can be of 
# non-coplanar. This technique simplifies the simulation procedure. 
# A single MC simulation job will be able to generate the 3d dose distribution
# of a plan, no matter how many treatment fields the plan has and no matter
# the fields are of standard geometry, IMRT, IMAT or their mixing.
# The only limitation is that they must be of the same energy.
#
# This script requires two templates for BEAMnrc and DOSXYZnrc egsinp files.
# 1. beamnrc_template.egsinp:
#    Use your regular beamnrc.egsinp file prepared for working with 
#    DOSXYZnrc source 21, 
#    (a) replace the regular jaw file position with SYNCJAWS_FILE, 
#    (b) replace the regular mlc file position with SYNCVMLC_FILE.
#    The parameters will be updated with the actual jaw_file and mlc_file
#    when running this script.
# 2. dosxyznrc_template.egsinp:
#    Use your regular dosxyznrc.egsinp file prepared for DOSXYZNRC source 21,
#    (a) replace the regular egsphant file position with EGSPHANT_FILE,
#    (b) replace the # of control points position with NUMBER_OF_CONTROL_POINTS,
#    (c) replace the list of gantry positions with GANTRY_POSITION_LIST,
#    (d) replace beam model position with BEAM_MODEL,
#    (e) replace the beamnrc input file position with BEAM_INPUT_FILE.
#    Both BEAM_MODEL and BEAM_INPUT_FILE are specific to source 21.
#    They will be updated with the actual values extracted from the 
#    DICOM plan files or corresponding to the beamnrc input file, when running
#    this script.
#
# NOTE: This script is only tested on plans of Eclipse 8.9 and Pinnacle3 9.0
#       with Varian high energy Linacs.
#       Some modification might be required for plans generated using 
#       other treatment planning systems and Linacs from other vendors.
#
# usage: python pycom.py [-h] -p DCMPLAN [-B BEAMDIR] [-D DOSXYZDIR]
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   -p DCMPLAN, -i DCMPLAN, --dcmplan DCMPLAN
#                         DICOM Plan input
#   -B BEAMDIR, --beamdir BEAMDIR
#                         Directory for BEAM model
#   -D DOSXYZDIR, --dosxyzdir DOSXYZDIR
#                         Directory for dosxyznrc
#
# Example:
# python pycom.py -p /path/to/dcm/plan.dcm 
#
#
#########################################################
# Developed by Lixin Zhan (lixinzhan AT gmail DOT com). #
#########################################################

import os
import argparse
import numpy as np
import dicom
from dicom2dosxyz import dcm2dosxyz

parser=argparse.ArgumentParser(description='Create MC Input from DICOM Plan')
parser.add_argument('-p', '-i', '--dcmplan', help='DICOM Plan input', required=True)
parser.add_argument('-B', '--beamdir', help='Directory for BEAM model')
parser.add_argument('-D', '--dosxyzdir', help='Directory for dosxyznrc')
args=parser.parse_args()

if args.beamdir:
    beamdir = args.beamdir
else:
    beamdir = os.getenv('EGS_HOME') + '/BEAM_ARC'
if args.dosxyzdir:
    dosxyzdir = args.dosxyzdir
else:
    dosxyzdir = os.getenv('EGS_HOME') + '/dosxyznrc'

dcmplanfile = args.dcmplan
dcmplan = dicom.read_file(dcmplanfile)

patientname = dcmplan.PatientsName
patientid = dcmplan.PatientID
planname = dcmplan.RTPlanLabel

# open the serial of files for output
outpath = os.path.dirname(os.path.abspath(dcmplanfile)) + '/'
outpathbeam = outpath+'beamnrc/'
outpathdose = outpath+'dosxyznrc/'
if not os.path.exists(outpathbeam):
    os.makedirs(outpathbeam)
if not os.path.exists(outpathdose):
    os.makedirs(outpathdose)

filemainpart = patientid + '_' + planname

egsinp_beam = outpathbeam + filemainpart + '.egsinp'
egsinp_dose = outpathdose + filemainpart + '.egsinp'

cp_gantry_file = outpath + patientid+'_gantry.txt'   # gantry
cpgantry = open(cp_gantry_file,'w')
cp_mlc_file = outpathbeam + filemainpart + '.mlc'    # mlc
cpmlc = open(cp_mlc_file,'w')
djaw_file = outpathbeam + filemainpart + '.djaw'     # jaws
djaw = open(djaw_file,'w')

dcm_gantry_file = outpath + patientid+'_dcmgantry.txt'
dcmgantry = open(dcm_gantry_file,'w')

# for verification of gantry angles
cpgantry_cmp = open(outpath + patientid+'_cmp.txt','w')
# in treatment log compatible format.
pdmlc_file = outpath + patientid+'_pdmlc.txt'
pdmlc = open(pdmlc_file,'w')

dsource = 45.0  # in cm, the source-iso distance

RefBeam = {}
MUtot = 0
for rbeam in dcmplan.FractionGroups[0].ReferencedBeams:
    if hasattr(rbeam,'BeamMeterset'):
        RefBeam[rbeam.ReferencedBeamNumber] = rbeam.BeamMeterset
        MUtot = MUtot + rbeam.BeamMeterset
print 'MU for each beam: ', RefBeam
print 'Total MU: ', MUtot

isIMRT = True
isARC = True
isALL6MV = True
ncptot = 0
ncpgantry = 0
for beam in dcmplan.Beams:
    if beam.TreatmentDeliveryType != 'TREATMENT': # tested for varian
        continue
    cp = beam.ControlPoints[0]
    NominalEnergy = cp.NominalBeamEnergy
    ncp = beam.NumberofControlPoints
    print '# of CP:', ncp, 'for beam #', beam.BeamNumber, '(', beam.BeamName,')', \
        ' with MU ', RefBeam[beam.BeamNumber], 'and MU Weight', RefBeam[beam.BeamNumber]/MUtot
    ncptot = ncptot + ncp
    if NominalEnergy-6.0 >= 1.0e-3:
        isALL6MV = False
    if cp.GantryRotationDirection == 'NONE':
        isIMRT = isIMRT * True
        isARC  = isARC  * False
    elif cp.GantryRotationDirection in ['CW', 'CC', 'CCW']:
        isIMRT = isIMRT * False
        isARC  = isARC  * True
    else:
        isIMRT = False
        isARC  = False

print 'Total CP:', ncptot
print    
print 'isIMRT = ', isIMRT, 'isARC = ', isARC, 'isALL6MV = ', isALL6MV

print >>djaw, 'SYNCJAW File for patient', patientid
print >>djaw, '%d' % (len(RefBeam)*2) # number of beams
print >>djaw, '2'                     # paired bars

# first two lines for dynvmlc input file.
if isIMRT:
    print >>cpmlc, 'IMRT MLC Control Points'
elif isARC:
    print >>cpmlc, 'RapidArc MLC Control Points'
print >>cpmlc, ncptot

# the head part for pdMLC file
print >>pdmlc, 'File Rev = G'
print >>pdmlc, 'Treatment = Dynamic Dose'
print >>pdmlc, 'Last Name = tst'
print >>pdmlc, 'First Name = tst'
print >>pdmlc, 'Patient ID =', patientid
print >>pdmlc, 'Number of Fields =', beam.NumberofControlPoints
print >>pdmlc, 'Number of Leaves = 120'
print >>pdmlc, 'Tolerance = 0.100000'

JawOpenMax = 0.0
PreBeamWeight = 0.0
for beam in dcmplan.Beams:
    if beam.TreatmentDeliveryType != 'TREATMENT': # tested for varian
        continue
        
    cp = beam.ControlPoints[0] # the first control pt for this beam
    iso = cp.IsocenterPosition
    BeamWeight = RefBeam[beam.BeamNumber]/MUtot

    mu_index1 = PreBeamWeight
    mu_index2 = PreBeamWeight + BeamWeight
    
    xjawmin = np.nan
    xjawmax = np.nan
    yjawmin = np.nan
    yjawmax = np.nan
    
    CurrentBeamLimitingDevicePositions = cp.BeamLimitingDevicePositions  # initialization for this beam
    for device in cp.BeamLimitingDevicePositions:
        if device.RTBeamLimitingDeviceType in ['X','ASYMX']: # or 'X' for symm fld
            xjawmin = device.LeafJawPositions[0]/10.0
            xjawmax = device.LeafJawPositions[1]/10.0
        elif device.RTBeamLimitingDeviceType in ['Y','ASYMY']: # or 'Y' for symm fld
            yjawmax = -device.LeafJawPositions[0]/10.0
            yjawmin = -device.LeafJawPositions[1]/10.0
    print >>djaw, mu_index1
    print >>djaw,'Y'
    print >>djaw,'%.4f, %.4f, %.4f, %.4f,' % \
            (yjawmax*0.28, yjawmax*0.358,yjawmin*0.28,yjawmin*0.358)
    print >>djaw,'X'
    print >>djaw,'%.4f, %.4f, %.4f, %.4f,' % \
            (xjawmax*0.367,xjawmax*0.445,xjawmin*0.367,xjawmin*0.445)
    print >>djaw, mu_index2
    print >>djaw,'Y'
    print >>djaw,'%.4f, %.4f, %.4f, %.4f,' % \
            (yjawmax*0.28, yjawmax*0.358,yjawmin*0.28,yjawmin*0.358)
    print >>djaw,'X'
    print >>djaw,'%.4f, %.4f, %.4f, %.4f,' % \
            (xjawmax*0.367,xjawmax*0.445,xjawmin*0.367,xjawmin*0.445)
    JawOpenMax = max(JawOpenMax,-xjawmin,xjawmax,-yjawmin,yjawmax)        
    
    # gamma: gantry angle, rho: couch angle, col: colimator angle
    # in DICOM RT Plan.
    col = cp.BeamLimitingDeviceAngle
    rho = cp.PatientSupportAngle
    if isIMRT:
        gamma = cp.GantryAngle
        out=dcm2dosxyz(gamma, rho, col, 'Zhan') # coord. transform
        print >>cpgantry, '%.6f, %.6f, %.6f, %.2f, %.2f, %.2f, %.6f, %.6f' % \
            (iso[0]/10.0, iso[1]/10.0, iso[2]/10.0, out[0], out[1], out[2], dsource, mu_index1)
        ncpgantry = ncpgantry + 1
        print >>dcmgantry,'%.6f, %.6f, %.6f, %.2f, %.2f, %.2f, %.6f, %.6f' % \
            (iso[0]/10.0, iso[1]/10.0, iso[2]/10.0, gamma, rho, col, dsource, mu_index1)
        print >>cpgantry, '%.6f, %.6f, %.6f, %.2f, %.2f, %.2f, %.6f, %.6f' % \
            (iso[0]/10.0, iso[1]/10.0, iso[2]/10.0, out[0], out[1], out[2], dsource, mu_index2)
        ncpgantry = ncpgantry + 1
        print >>dcmgantry,'%.6f, %.6f, %.6f, %.2f, %.2f, %.2f, %.6f, %.6f' % \
            (iso[0]/10.0, iso[1]/10.0, iso[2]/10.0, gamma, rho, col, dsource, mu_index2)
    
        output_bush=dcm2dosxyz(gamma,rho,col,'Bush')
        output_thebaut=dcm2dosxyz(gamma,rho,col,'Thebaut')
        print >>cpgantry_cmp, 'Zhan:    %.2f, %.2f, %.2f, %.6f' % \
            (out[0], out[1], out[2], mu_index1)
        print >>cpgantry_cmp, 'Thebaut: %.2f, %.2f, %.2f, %.6f' % \
            (output_thebaut[0], output_thebaut[1], output_thebaut[2], mu_index1)
        print >>cpgantry_cmp, 'Bush:    %.2f, %.2f, %.2f, %.6f' % \
            (output_bush[0], output_bush[1], output_bush[2], mu_index1)
        
    for cp in beam.ControlPoints:
        mu_index = cp.CumulativeMetersetWeight*BeamWeight/beam.FinalCumulativeMetersetWeight + PreBeamWeight
        if isARC: # is ARC --- other than ARC, gantry angle is only in cp[0]
            gamma = cp.GantryAngle
            out=dcm2dosxyz(gamma, rho, col, 'Zhan') # coord. transform
        
            # print the control points for gantry
            print >>cpgantry, '%.6f, %.6f, %.6f, %.2f, %.2f, %.2f, %.6f, %.6f' % \
                (iso[0]/10.0, iso[1]/10.0, iso[2]/10.0, out[0], out[1], out[2], dsource, mu_index)
            ncpgantry = ncpgantry + 1
            print >>dcmgantry,'%.6f, %.6f, %.6f, %.2f, %.2f, %.2f, %.6f, %.6f' % \
                (iso[0]/10.0, iso[1]/10.0, iso[2]/10.0, gamma, rho, col, dsource, mu_index)
        
            output_bush=dcm2dosxyz(gamma,rho,col,'Bush')
            output_thebaut=dcm2dosxyz(gamma,rho,col,'Thebaut')
            print >>cpgantry_cmp, 'Bush:    %.2f, %.2f, %.2f, %.6f' % \
                (out[0], out[1], out[2], mu_index1)
            print >>cpgantry_cmp, 'Thebaut: %.2f, %.2f, %.2f, %.6f' % \
                (output_thebaut[0], output_thebaut[1], output_thebaut[2], mu_index1)
            print >>cpgantry_cmp, 'Zhan:    %.2f, %.2f, %.2f, %.6f' % \
                (output_bush[0], output_bush[1], output_bush[2], mu_index1)

        # find mlcx
        try:
            # make sure there is BeamLimitingDevice defined for this control point.
            getattr(cp, 'BeamLimitingDevicePositions')
            CurrentBeamLimitingDevicePositions = cp.BeamLimitingDevicePositions
        except AttributeError: # using the pos from last control point
            pass
    
        # now the mlc part
        print >>cpmlc, '%.6f' % mu_index    
        
        print >>pdmlc, ''
        print >>pdmlc, 'Field = 1-%d' % (cp.ControlPointIndex)
        print >>pdmlc, 'Index = %.6f' % mu_index
        print >>pdmlc, 'Carriage Group = 1'
        print >>pdmlc, 'Operator ='
        print >>pdmlc, 'Collimator = 0.0'
        
        mlc_exist = False
        for device in CurrentBeamLimitingDevicePositions:
            if device.RTBeamLimitingDeviceType == 'MLCX': # valid for varian.
                mlc_exist = True
                leaves = device.LeafJawPositions
                for n in range(len(leaves)/2):
                    m=60-n-1
                    print >>cpmlc, '%.3f, %.3f, %d' % \
                        (0.51*leaves[m]/10.0, 0.51*leaves[len(leaves)/2+m]/10.0, 1)
                    
                for n in range(len(leaves)/2):
                    print >>pdmlc, 'Leaf %2dA = %.2f' % (n+1, leaves[n]/10.0)
                for n in range(len(leaves)/2):
                    m = len(leaves)/2+n
                    print >>pdmlc, 'Leaf %2dB = %.2f' % (n+1, leaves[m]/10.0)
        if not mlc_exist: # Case of 'MLCY' not checked. Possibly needed later.
            for n in range(60): 
                print >>cpmlc, '%.3f, %.3f, %d' % (-20, 20, 1)
            for n in range(60):
                print >>pdmlc, 'Leaf %2dA = %.2f' % (n+1, -20)
            for n in range(60):
                print >>pdmlc, 'Leaf %2dB = %.2f' % (n+1, 20)
        
        print >>pdmlc, 'Note = 0'
        print >>pdmlc, 'Shape = 0'
        print >>pdmlc, 'Magnification = 1.00'        
                
    PreBeamWeight = PreBeamWeight + BeamWeight
    
    print >>pdmlc, ''
    print >>pdmlc, 'CRC = 43B6'

print 'Jaw Max Opening: ', JawOpenMax
if JawOpenMax > 15:
    print 'WARNING: Jaw Opening is BIGGER than maximun allowance from DBS setting !!'
    
cpgantry.close()
dcmgantry.close()
cpmlc.close()
pdmlc.close()
djaw.close()

# generate the BEAMnrc input file
s = open('beamnrc_template.egsinp').read()
s = s.replace('SYNCJAWS_FILE', os.path.normpath(beamdir+'/'+filemainpart+'.djaw'))
s = s.replace('SYNCVMLC_FILE', os.path.normpath(beamdir+'/'+filemainpart+'.mlc'))
outfile = open(egsinp_beam,'w')
outfile.write(s)
outfile.close()

gpos = open(cp_gantry_file).read()[:-1]
# generate the dosxyznrc input file
s = open('dosxyznrc_template.egsinp').read()
s = s.replace('NUMBER_OF_CONTROL_POINTS', str(ncpgantry))
s = s.replace('EGSPHANT_FILE', os.path.normpath(dosxyzdir+'/'+filemainpart+'.egsphant'))
s = s.replace('GANTRY_POSITION_LIST', gpos)
s = s.replace('BEAM_MODEL',os.path.normpath(beamdir).split('/')[-1])
s = s.replace('BEAM_INPUT_FILE',filemainpart) 
outfile = open(egsinp_dose,'w')
outfile.write(s)
outfile.close()

