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
import random

def get_dicomuid():
    uid_root = '9.9.9.9.9.9.9999.999.9.1.'

    now = datetime.now()
    time_str = now.strftime('%Y%m%d%H%M%S') + f"{now.microsecond:06d}"
    rand_str = f"{random.randint(0, 999999999):09d}"

    uid = uid_root + time_str + rand_str
    
    return uid
