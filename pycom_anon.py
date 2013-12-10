###################################################################
# This code is part of PYCOM, The Python DICOM Processing Toolkit.
# PYCOM is distributed open source under GPL v3 license.
###################################################################
#
# This is a script for anonymizing DICOM CT files.
#
# usage: python pycom_anon.py [-h] -i DCMLIST [-n ANON_NAME] [-id ANON_ID]
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   -i DCMLIST, -l DCMLIST, --dcmlist DCMLIST
#                         List of DICOM files to be anonymized
#   -n ANON_NAME, --anon_name ANON_NAME
#                         Anonymous name used for replacing the original patient
#                         name. Default: anon
#   -id ANON_ID, --anon_id ANON_ID
#                         Anonymous ID used for replacing the original patient
#                         ID. Default: 00000000
# 
# Example:
# python pycom_anon.py -i dcmlist -n TEST -id 0001
#
#########################################################
# Developed by Lixin Zhan (lixinzhan AT gmail DOT com). #
#########################################################

import os
import argparse
import dicom
from pycomuid import get_dicomuid

parser = argparse.ArgumentParser(description='Anonymizing DICOM files')
parser.add_argument('-i', '-l', '--dcmlist', help='List of DICOM files to be anonymized', required=True)
parser.add_argument('-n', '--anon_name', help='Anonymous name used for replacing the original patient name. Default: anon')
parser.add_argument('-id', '--anon_id', help='Anonymous ID used for replacing the original patient ID. Default: 00000000')
args=parser.parse_args()

f = open(os.path.normpath(args.dcmlist))
dcm_list = f.readlines()

if args.anon_name:
    anon_name = args.anon_name
else:
    anon_name = 'anon'
if args.anon_id:
    anon_id = args.anon_id
else:
    anon_id = '00000000'

FOR_uid = get_dicomuid()

for dcmfile in dcm_list:    
    dcmfilename = dcmfile.strip()
    dcm = dicom.read_file(dcmfilename)

    dcm.FrameOfReferenceUID = FOR_uid        
    dcm.PatientsName = anon_name
    dcm.PatientID = anon_id
    dcm.PatientsSex = 'O'
    dcm.SOPInstanceUID = get_dicomuid()
    dcm.StudyUID = get_dicomuid()
    dcm.SeriesUID = get_dicomuid()
    dcm.FrameUID = get_dicomuid()
    dcm.SyncUID = get_dicomuid()
    dcm.SrUID = get_dicomuid()
    dcm.StudyInstanceUID = get_dicomuid()
    dcm.SeriesInstanceUID = get_dicomuid()
    
    outfile = 'anon_'+dcmfilename
    dicom.write_file(outfile, dcm)

