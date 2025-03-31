import gzip
import pandas as pd
hdffile = pd.read_hdf('C:\MOPGA 2022-2023 Dodji\ARGO float Angola\dataNPP data Oregon state University\eppley.m.2021\eppley.2021001.hdf', keys='eppley.2021001')


import glob
import os
from zipfile import ZipFile

zipfiles = [os.path.basename(x) for x in glob.glob(r'C:\MOPGA 2022-2023 Dodji\ARGO float Angola\dataNPP data Oregon state University\*.tar')]

for zipfile in zipfiles:
    new_dirname = zipfile.replace('.bsw.', '_').replace('.tar', '')
    try:
        os.mkdir(f'output/{new_dirname}')
    except:
        print(f'output/{new_dirname} already exists')
    try:
        os.mkdir(f'output/{new_dirname}/subdir')
    except:
        print(f'output/{new_dirname}/subdir already exists')

    with ZipFile(f'C:\MOPGA 2022-2023 Dodji\ARGO float Angola\dataNPP data Oregon state University/{zipfile}') as zipObject:
        zipObject.extractall(path=f'output/{new_dirname}/subdir')