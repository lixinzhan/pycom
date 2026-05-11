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
#   --output-dir OUTPUT_DIR
#                         Directory for output files. Default: current directory
# 
# Example:
# python pycom_anon.py -i dcmlist -n TEST -id 0001
#
#########################################################
# Developed by Lixin Zhan (lixinzhan AT gmail DOT com). #
#########################################################

from pathlib import Path
import argparse
import pydicom
import pydicom.uid

parser = argparse.ArgumentParser(description='Anonymizing DICOM files')
parser.add_argument('-i', '-l', '--dcmlist', help='List of DICOM files to be anonymized', required=True)
parser.add_argument('-n', '--anon_name', default='anon', help='Anonymous name to be used. Default: anon')
parser.add_argument('-id', '--anon_id', default='00000000', help='Anonymous ID to be used. Default: 00000000')
parser.add_argument('--output-dir', default='.', help='Directory for output files. Default: current directory')
args=parser.parse_args()

dcmlist_path = Path(args.dcmlist).resolve()
with open(dcmlist_path, 'r') as f:
    dcm_list = [line.strip() for line in f.readlines()]

output_dir = Path(args.output_dir)
output_dir.mkdir(exist_ok=True)

anon_name = args.anon_name
anon_id = args.anon_id

FOR_uid = pydicom.uid.generate_uid()

for dcmfile in dcm_list:
    try:
        dcm = pydicom.read_file(dcmfile)
        dcm.FrameOfReferenceUID = FOR_uid
        dcm.PatientsName = anon_name
        dcm.PatientID = anon_id
        dcm.PatientsSex = 'O'
        dcm.SOPInstanceUID = pydicom.uid.generate_uid()
        dcm.FrameUID = pydicom.uid.generate_uid()
        dcm.SyncUID = pydicom.uid.generate_uid()
        dcm.SrUID = pydicom.uid.generate_uid()
        dcm.StudyInstanceUID = pydicom.uid.generate_uid()
        dcm.SeriesInstanceUID = pydicom.uid.generate_uid()

        outfile = output_dir / f"anon_{Path(dcmfile).name}"
        pydicom.write_file(outfile, dcm)
        print(f"Processed: {dcmfile} -> {outfile}")
    except Exception as e:
        print(f"Error processing {dcmfile}: {e}")


# Note: 
# 1. StudyUID and SeriesUID removed from the for loop. To be double checked later.
# 2. Using pycom.uid.generate_uid() now, which replaced the original pycomuid.py code.