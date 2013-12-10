###################################################################
# This code is part of PYCOM, The Python DICOM Processing Toolkit.
# PYCOM is distributed open source under GPL v3 license. 
###################################################################
#
# subroutine for generating DICOM UID, which is used for 
# DICOM file anonymizing.
#
#########################################################
# Developed by Lixin Zhan (lixinzhan AT gmail DOT com). #
#########################################################

import os
from datetime import datetime
from random import random

def get_dicomuid():
    uid_root = '9.9.9.9.9.9.9999.999.9.1.'

    uid = uid_root

    t = str(datetime.now())
    year = t[0:4]
    month = t[5:7]
    day = t[8:10]
    hr = t[11:13]
    min = t[14:16]
    sec = t[17:19]
    msec = repr(datetime.now().microsecond).zfill(6)
    rand = repr(int(random()*10**9)).zfill(9)

    uid = uid_root+year+month+day+hr+min+sec+msec+rand
    
    return uid
