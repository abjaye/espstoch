#!/usr/bin/env python3
import os, sys
import glob, shutil
from datetime import timedelta, datetime
import argparse
import subprocess
import numpy as np
import xarray as xr 
import warnings
warnings.filterwarnings("ignore")

def parse_command_line(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--date",
                        help="Specify a start Date")
    parser.add_argument("--ens",default=0,
                        help="Specify the ensemble member")

    args = parser.parse_args()

    if args.date:
        try:
            date = datetime.strptime(args.date, '%Y-%m-%d')
        except ValueError as verr:
            raise ValueError("Incorrect data format, should be YYYY-MM-DD") from verr
    elif cdate:
        date = datetime.strptime(cdate, '%Y-%m-%d')
    else:
        date = datetime.today() - timedelta(days=1)

    return date.strftime("%Y-%m-%d"),args.ens

if __name__ == '__main__':
    date, ens = parse_command_line(sys.argv)

    var  = ["smois","tslb","sh2o","skintemp","snow","snowh"]
    var2 = ["SOIL_M","SOIL_T","SOIL_W","TG","SNEQV"]

    dated = str(datetime.strptime(date,'%Y-%m-%d').date() - timedelta(days=1))
    yyyy = dated.split("-")[0]; mm = dated.split("-")[1]; dd = dated.split("-")[2]
    date2 = yyyy + mm + dd

    noah_path = "/glade/campaign/ral/hap/zhezhang/offline_hrldas/LDASOUT_mesh/"
    noah_file = noah_path+date2+"06.LDASOUT_DOMAIN1"
    noah_in = xr.open_dataset(noah_file) 
    base_path = "/glade/derecho/scratch/espstoch/2012case/mpass2s_"+date+".0"+ens+"/"
    init_file = base_path+"x1.163842.init.nc"
    static = xr.open_dataset(init_file)

    init_spin = base_path+"x1.163842.init.nc_spinup"
    os.system("cp "+init_file+" "+init_spin)
    init_new = xr.open_dataset(init_spin)

    for s in range(3):
        data_array = np.zeros((163842,4),float)
        for l in range(4):
            data_array[:,l] = noah_in.variables[var2[s]][3,0,l,:]
            data_array[:,l] = np.where(data_array[:,l]>0.0,data_array[:,l],0.0)

        init_new.variables[var[s]][0,:,:] = data_array

    for s in range(4,5):
        data_array = np.zeros((163842),float)
        data_array[:] = noah_in.variables[var2[s]][3,:]
        data_array[:] = np.where(data_array[:]>0.0,data_array[:],0.0)

        init_new.variables[var[s]][0,:] = data_array

