###################################################################
# This code is part of PYCOM, The Python DICOM Processing Toolkit.
# PYCOM is distributed open source under GPL v3 license. 
###################################################################
#
# Python DICOM RD generator
#
# This script generates DOCOM Dose file from the 3ddose file of DOSXYZnrc
#
# usage: python pycom_dose.py [-h] -p DCMPLAN -d DCMDOSE -m MCDOSE -c PPICDOSE
#                      [-B BEAMDIR] [-D DOSXYZDIR] [-O OUTPUTDIR]
#                      [-T DCMDOSETEMPLATEDIR]
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   -p DCMPLAN, --dcmplan DCMPLAN
#                         DICOM Plan input. If --dcmdosetemplatedir is not
#                         presented, path is required.
#   -d DCMDOSE, --dcmdose DCMDOSE
#                         DICOM Dose input. if --dcmdosetemplatedir is not
#                         presented, path is required
#   -m MCDOSE, --mcdose MCDOSE
#                         Monte Carlo 3ddose input, default to dosxyznrc folder,
#                         if no path specified.
#   -c PPICDOSE, --ppicdose PPICDOSE
#                         Parallel Plate Ion Chamber Dose, from MC simulation.
#   -B BEAMDIR, --beamdir BEAMDIR
#                         Directory for BEAM model
#   -D DOSXYZDIR, --dosxyzdir DOSXYZDIR
#                         Directory for dosxyznrc
#   -O OUTPUTDIR, --outputdir OUTPUTDIR
#                         Output directory for the converted DICOM RD file.
#   -T DCMDOSETEMPLATEDIR, -I DCMDOSETEMPLATEDIR, --dcmdosetemplatedir DCMDOSETEMPLATEDIR
#                         Directory that DICOM dose template exist. Once this is
#                         provided, only filename is required for both dcmplan
#                         and dcmdose; otherwise, paths for both dcmplan and
#                         dcmdose need to be specified.
#
# Example:
#
# python pycom_dose.py -p /path/to/dcm/plan/rp.dcm -d /path/to/dcm/dose/rd.dcm 
#                      -m /path/to/dosxyznrc/mc.3ddose -c 3.459e-15 
#                      -O /path/for/dcm/dose/output/
#
# Before using this script, you need to perform an absolute dose calibration
# first. The absolute dose calibration method is following 
# Popescu et al Phys. Med. Biol. 50:3375-3392, 2005.
#
# Suppose your calibration is using water phantom and results in 
# 1cGy/MU at depth d=5cm with calibration geometry SSD=100cm and FS=10x10cm2.
# Your MC simulation using the exact same calibration setup obtains
# scored dose in the parallel plate chamber of the Linac head is 3.459e-15,
# and dose in water phantom at d=5cm is 1.43915e-16. You can then find
# variables d_cal_ch, d_depcal, and D_cal_abs and change to their corresponding 
# values as below:
#
# d_cal_ch = 3.459e-15   # chamber dose from MC calib.
# d_depcal = 1.43915e-16 # dose @ d=5cm from MC calib.
# D_cal_abs = 0.01       # 1 cGy/MU = 0.01 Gy/MU
#
# Once you have finished absolute dose calibration, you can then convert 
# your DOSXYZnrc generated 3ddose file to DICOM RD format using this tool.
#
#########################################################
# Developed by Lixin Zhan (lixinzhan AT gmail DOT com). #
#########################################################

import os
import argparse
import struct
import StringIO as sio
import numpy as np
import dicom
from pycomuid import get_dicomuid

parser = argparse.ArgumentParser(description='Create DICOM RD from MC 3ddose file')
parser.add_argument('-p', '--dcmplan', help='DICOM Plan input. If --dcmdosetemplatedir is not presented, path is required.', required=True)
parser.add_argument('-d', '--dcmdose', help='DICOM Dose input. if --dcmdosetemplatedir is not presented, path is required', required=True)
parser.add_argument('-m', '--mcdose', help='Monte Carlo 3ddose input, default to dosxyznrc folder, if no path specified.', required=True)
parser.add_argument('-c', '--ppicdose', type=float, help='Parallel Plate Ion Chamber Dose, from MC simulation.', required=True)
parser.add_argument('-B', '--beamdir', help='Directory for BEAM model')
parser.add_argument('-D', '--dosxyzdir', help='Directory for dosxyznrc')
parser.add_argument('-O', '--outputdir', help='Output directory for the converted DICOM RD file.')
parser.add_argument('-T', '-I', '--dcmdosetemplatedir', help='Directory that DICOM dose template exist. Once this is provided, only filename is required for both dcmplan and dcmdose; otherwise, paths for both dcmplan and dcmdose need to be specified.')
args = parser.parse_args()

if args.beamdir:
    beamdir = args.beamdir
else:
    beamdir = os.getenv('EGS_HOME') + '/BEAM_ARC'
if args.dosxyzdir:
    dosxyzdir = args.dosxyzdir
else:
    dosxyzdir = os.getenv('EGS_HOME') + '/dosxyznrc'

d_ch = args.ppicdose

if os.path.normpath(args.mcdose)==os.path.basename(args.mcdose):
    mcdosefile  = dosxyzdir+'/'+args.mcdose
else:
    mcdosefile = os.path.normpath(args.mcdose)

if args.dcmdosetemplatedir:
    dcmplanfile = args.dcmdosetemplatedir+'/'+args.dcmplan
    dcmdosefile = args.dcmdosetemplatedir+'/'+args.dcmdose
else:
    dcmplanfile = args.dcmplan
    dcmdosefile = args.dcmdose

if args.outputdir:
    outfile = args.outputdir + '/mc_' + os.path.basename(dcmdosefile)
else:
    outfile = os.path.dirname(os.path.abspath(dcmdosefile)) + \
        '/mc_'+os.path.basename(dcmdosefile)

print
print '#################################'
print 'BEAMnrc Dir:       ' + os.path.normpath(beamdir)
print 'DOSXYZnrc Dir:     ' + os.path.normpath(dosxyzdir)
print 'DCM RP Input:      ' + os.path.normpath(dcmplanfile)
print 'DCM RD Template:   ' + os.path.normpath(dcmdosefile)
print 'MC 3ddose file:    ' + os.path.normpath(mcdosefile)
print 'New DCM RD output: ' + os.path.normpath(outfile)
print
print '#################################'
print

################################################################################
# This part is to be change by user for corresponding absolute dose calibration.
d_depcal = 1.43915e-16 # calibration @ d=5cm
d_cal_ch = 3.459e-15   # chamber dose
D_cal_abs = 0.01       # 1 cGy/MU = 0.01 Gy/MU
################################################################################

# Ion chamber absolute dose (in Gy) per MU. D_ch_MU are obtained from numbers above.
D_ch_MU = D_cal_abs * d_cal_ch/d_depcal

print 'MC Chamber Dose:    ', d_ch
print 'Calibration Factor: ', D_ch_MU
print
###### Info from DICOM RP file
dcmplan = dicom.read_file(dcmplanfile)
RefBeam = {}
MUtot = 0
#DoseTot = 0
for rbeam in dcmplan.FractionGroups[0].ReferencedBeams:
    if hasattr(rbeam,'BeamMeterset'):
        RefBeam[rbeam.ReferencedBeamNumber] = rbeam.BeamMeterset
        MUtot = MUtot + rbeam.BeamMeterset
        #DoseTot = DoseTot + rbeam.BeamDose
print 'MU for each beam: ', RefBeam
print 'Total MU: ', MUtot
#print 'Total Dose', DoseTot
nFrac = dcmplan.FractionGroups[0].NumberofFractionsPlanned
print 'nFrac = ', nFrac
print

# Info from the model DICOM RT Dose file
dcmdose = dicom.read_file(dcmdosefile)
xsizeold = dcmdose.Columns
ysizeold = dcmdose.Rows
zsizeold = dcmdose.NumberofFrames

DoseSumType = dcmdose.DoseSummationType
if DoseSumType=='FRACTION':
    MU = MUtot/nFrac
elif DoseSumType=='PLAN':
    MU = MUtot
elif DoseSumType=='BEAM':
    MU = MUtot
    print 'BEAM ONLY !!!'
else:
    print 'DoseSumType ',DoseSumType, ' not supported.'
print
print 'DoseSumType: ', DoseSumType, ' MUtot: ', MUtot, ' Fracions: ',nFrac, 'MU used for Generating DICOM RD: ', MU
print

dosescaling = dcmdose.DoseGridScaling

length = xsizeold*ysizeold*zsizeold
BitsAllocated = dcmdose.BitsAllocated
if BitsAllocated==32:
    format = 'I'*length
elif BitsAllocated==16:
    format = 'H'*length
else:
    print 'BitsAllocated: ', BitsAllocated, ' is not supported yet.'
pdata = struct.unpack(format, dcmdose.PixelData)
pdmax = max(pdata)

print '#################################'
print 'ImagePos:', dcmdose.ImagePositionPatient
print 'xsizeold:',xsizeold, 'ysizeold:',ysizeold, 'zsizeold:',zsizeold
print 'DoseGridScaling', dosescaling 
print 'MaxDose(INT):',pdmax, 'MaxDose(Gy):',pdmax*dosescaling, '<--'
print 'TotLength:',length, 'Bits:',length*dcmdose.BitsAllocated/8
print 'OffsetVector Length:', np.size(dcmdose.GridFrameOffsetVector)
print '-----------------------------'


# Info from the DOSXYZnrc 3ddose file
dose3d = open(mcdosefile)
xsize, ysize, zsize = np.array(dose3d.readline().split(),np.int)

x = np.array(dose3d.readline().split(),np.float)
y = np.array(dose3d.readline().split(),np.float)
z = np.array(dose3d.readline().split(),np.float)
d = np.array(dose3d.readline().split(),np.float) # dose
e = np.array(dose3d.readline().split(),np.float) # dose error

dx = x[1]-x[0]
dy = y[1]-y[0]
dz = z[1]-z[0]

dcmdose.SeriesDescription = 'DOSXYZnrc Doses'
dcmdose.ManufacturersModelName = 'DOSXYZnrc'
dcmdose.DeviceSerialNumber = '000000000'

dcmdose.Rows = ysize
dcmdose.Columns = xsize
dcmdose.NumberofFrames = zsize
dcmdose.ImagePositionPatient = [(x[0]+dx/2.0)*10.0, (y[0]+dy/2.0)*10.0,
    (z[0]+dz/2.0)*10.0]
z0 = float(dcmdose.ImagePositionPatient[2])
dcmdose.PixelSpacing = [(x[-1]-x[0])*10.0/xsize,(y[-1]-y[0])*10.0/ysize]
dcmdose.GridFrameOffsetVector = ((z[1:]-dz/2.0)*10.0-z0).tolist()

length = xsize * ysize * zsize
#format = 'I'*length
if BitsAllocated==32:
    format = 'I'*length
elif BitsAllocated==16:
    format = 'H'*length
else:
    print 'BitsAllocated: ', BitsAllocated, ' is not supported yet.'

# di = np.array(d*pdmax/d.max(),np.int)
# di = np.array(0.01*(d/d_cal5cm)*MU*nFrac/dosescaling,np.int)
di = np.array((d/d_ch)*nFrac*MU*D_ch_MU/dosescaling,np.int)
dcmdose.PixelData = struct.pack(format,*(di.tolist()))

print 'check dose volume sizes'
print 'Pixel Range:', x[0],x[-1],y[0],y[-1],z[0],z[-1]
print 'dx, dy, dz:', dx, dy, dz
print 'ImageOrigien:',dcmdose.ImagePositionPatient
print 'sizes:',xsize, ysize, zsize, np.size(x),np.size(y),np.size(z)
print 'Max D/Ne:',d.max(), 'Abs Dmax:',d.max()*MU*nFrac*D_ch_MU/d_ch, '<--'
print 'PixelSpacing:',dcmdose.PixelSpacing, 'MaxDose:',di.max()
print length, np.size(d), np.size(e)
print 'New Image Position:', dcmdose.ImagePositionPatient
print 'OffsetVector Length:', np.size(dcmdose.GridFrameOffsetVector)
print 

dcmdose.SOPInstanceUID = get_dicomuid()
dcmdose.StudyUID = get_dicomuid()
dcmdose.SeriesUID = get_dicomuid()
dcmdose.FrameUID = get_dicomuid()
dcmdose.SyncUID = get_dicomuid()
dcmdose.SrUID = get_dicomuid()
dcmdose.StudyInstanceUID = get_dicomuid()
dcmdose.SeriesInstanceUID = get_dicomuid()
        
dicom.write_file(outfile, dcmdose)
print 'New DCM RD file created: ',outfile

#print xsize, ysize, zsize, length, length*dcmdose.BitsAllocated/8, max(di)
