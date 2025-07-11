#!/usr/bin/env python3
import os, sys
import glob, shutil
from datetime import timedelta, datetime
import argparse
import subprocess
import xarray as xr
import numpy as np
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
            raise ValueError("Incorrect data format, should be YYYY-MM-DD or YYYY-MM") from verr
    elif cdate:
        date = datetime.strptime(cdate, '%Y-%m-%d')
    else:
        date = datetime.today() - timedelta(days=1)

    return date.strftime("%Y-%m-%d"),args.ens

if __name__ == '__main__':
    date, ens = parse_command_line(sys.argv)

    workflow = "mpass2s_120km"
    if workflow == "mpass2s_120km":
        cells = "40962"
    elif workflow == "mpass2s":
        cells = "163842"

    month_abbr = ["","jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"]
    yyyy = date.split("-")[0]; mm = date.split("-")[1]; dd = date.split("-")[2]
    outdate = dd + month_abbr[int(mm)] + yyyy

    base_path = "/glade/derecho/scratch/espstoch/2012case/"+workflow+"_"+date+".0"+ens+"/"
    files = base_path + "diag_sfc*.nc"
    files2 = base_path + "diag_p1p2*.nc"
    grid_path = base_path + "x1."+cells+".grid.nc"
    uxds = xr.open_mfdataset(files,combine="nested",concat_dim="Time")
    uxds2 = xr.open_mfdataset(files2,combine="nested",concat_dim="Time")

    # Make p1 calculations

    prec_c = uxds["prec_acc_c"].resample(Time="D",closed="right").sum()[1::]
    prec_nc = uxds["prec_acc_nc"].resample(Time="D",closed="right").sum()[1::]
    pr_sfc = (prec_c + prec_nc)/86400
    pr_sfc.attrs ={"units": "kg m^-2 s^-1", "long_name": "Total (convective and large-scale) precipitation rate (liq + ice)"}
    pr_sfc = pr_sfc.to_dataset(name="pr_sfc")

    rlut = uxds2["olrtoa"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"olrtoa": "rlut"})
    tas_2m = uxds["t2m"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"t2m": "tas_2m"})
    ts = uxds["skintemp"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"skintemp": "ts"})

    zg_200 = uxds2["height_200hPa"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"height_200hPa": "zg_200"})
    zg_500 = uxds2["height_500hPa"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"height_500hPa": "zg_500"})

    ua_200 = uxds2["uzonal_200hPa"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"uzonal_200hPa": "ua_200"})
    ua_850 = uxds2["uzonal_850hPa"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"uzonal_850hPa": "ua_850"})
    va_200 = uxds2["umeridional_200hPa"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"umeridional_200hPa": "ua_200"})
    va_850 = uxds2["umeridional_850hPa"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"umeridional_850hPa": "ua_850"})

    datasetsp1 = [rlut,pr_sfc,tas_2m,ts,zg_200,zg_500,ua_200,ua_850,va_200,va_850]
    varsp1 = ["rlut","pr_sfc","tas_2m","ts","zg_200","zg_500","ua_200","ua_850","va_200","va_850"]

    for name, ds in zip(varsp1, datasetsp1):
        print(name)
        path1 = "/glade/derecho/scratch/espstoch/2012case/"+workflow+"_native/p1/"+name+"/"+yyyy+"/"+mm+"/"
        if not os.path.exists(path1):
            os.system("mkdir -p "+path1)
        path2 = "/glade/derecho/scratch/espstoch/2012case/"+workflow+"_latlon/p1/"+name+"/"+yyyy+"/"+mm+"/"
        if not os.path.exists(path2):
            os.system("mkdir -p "+path2)
        filename = path1+f"{name}_{workflow}_{outdate}_00z_unst_d46_m0{ens}.nc"
        filename_regrid = path2+f"{name}_{workflow}_{outdate}_00z_d01_d46_m0{ens}.nc"
        ds.to_netcdf(filename,unlimited_dims="Time")
        os.system("cdo -P 1 -f nc5 remapcon,r360x181 -setgrid,mpas:"+grid_path+" "+filename+" "+filename_regrid)

    # Make p2 calculations

    tc850 = uxds2["temperature_850hPa"].resample(Time="D",closed="right").mean()[1::] - 273.15
    rh850 = uxds2["temperature_850hPa"].resample(Time="D",closed="right").mean()[1::]
    pressure = 850.
    e_s = 6.112 * np.exp((17.67 * tc850) / (tc850 + 243.5))
    e = rh850 / 100 * e_s
    q850 = (0.622 * e) / (pressure - (0.378 * e))
    q850.attrs ={"units": "kg/kg", "long_name": "Specific Humidity at 850 mbar pressure surface"}
    huss_850 = q850.to_dataset(name="huss_850")

    q2m = uxds["q2"].resample(Time="D",closed="right").mean()[1::]
    sfcp = uxds["surface_pressure"].resample(Time="D",closed="right").mean()[1::]/100
    tdc = q2m*sfcp/(.622 + q2m)
    tdc = tdc.where(tdc>0.001, 0.001)
    tdps = ((243.5*np.log(tdc) - 440.8)/(19.48 - np.log(tdc))) + 273.15
    tdps.attrs ={"units": "K", "long_name": "Dewpoint temperature at 2 meters"}
    tdps = tdps.to_dataset(name="tdps")

    swdnb = uxds["swdnb"].resample(Time="D",closed="right").mean()[1::]
    swupb = uxds["swupb"].resample(Time="D",closed="right").mean()[1::]
    lwdnb = uxds["lwdnb"].resample(Time="D",closed="right").mean()[1::]
    lwupb = uxds["lwupb"].resample(Time="D",closed="right").mean()[1::]
    rad_sfc = lwupb-lwdnb+swupb-swdnb
    rad_sfc.attrs ={"units": "W/m2", "long_name": "Net surface radiation (lwupb-lwdnb+swupb-swdnb)"}
    rad_sfc = rad_sfc.to_dataset(name="rad_sfc")

    smois = uxds["smois"].resample(Time="D",closed="right").mean()[1::]
    mrso = smois.sum(dim="nSoilLevels")#.rename({"smois": "mrso"})
    mrso.attrs ={"units": "kg/m2", "long_name": "Vertically integrated soil moisture"}
    mrso = mrso.to_dataset(name="mrso")

    snow = uxds2["snow"].resample(Time="D",closed="right").mean()[1::]
    snowh = uxds2["snowh"].resample(Time="D",closed="right").mean()[1::]
    snd = snow/snowh
    snd.attrs ={"units": "kg/m3", "long_name": "Snow Density (snow/snowh)"}
    snd = snd.to_dataset(name="snd")

    snc = uxds2["snowc"].resample(Time="D",closed="right").mean()[1::]*100
    snc.attrs ={"units": "%", "long_name": "Snow Cover"}
    snc = snc.to_dataset().rename({"snowc": "snc"})
    sic = uxds2["xice"].resample(Time="D",closed="right").mean()[1::]*100
    sic.attrs ={"units": "%", "long_name": "Sea Ice Concentration"}
    sic = sic.to_dataset().rename({"xice": "sic"})

    tasmax_2m = uxds2["tasmax"].resample(Time="D",closed="right").max()[1::].to_dataset().rename({"tasmax": "tasmax_2m"})
    tasmin_2m = uxds2["tasmin"].resample(Time="D",closed="right").min()[1::].to_dataset().rename({"tasmin": "tasmin_2m"})
    hfss_sfc = uxds["hfx"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"hfx": "hfss_sfc"})
    hfls_sfc = uxds["lh"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"lh": "hfls_sfc"})
    wap_500 = uxds2["w_500hPa"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"w_500hPa": "wap_500"})
    ua_100 = uxds2["uzonal_100hPa"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"uzonal_100hPa": "ua_100"})
    va_100 = uxds2["umeridional_100hPa"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"umeridional_100hPa": "va_100"})
    uas = uxds["u10"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"u10": "uas"})
    vas = uxds["v10"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"v10": "vas"})
    psl = uxds["mslp"].resample(Time="D",closed="right").mean()[1::].to_dataset().rename({"mslp": "psl"})
    swe = uxds["snow_acc_nc"].resample(Time="D",closed="right").sum()[1::].to_dataset().rename({"snow_acc_nc": "swe"})
    cape = uxds["cape"].resample(Time="D",closed="right").mean()[1::].to_dataset()

    datasetsp2 = [huss_850,tasmax_2m,tasmin_2m,hfss_sfc,hfls_sfc,wap_500,ua_100,va_100,uas,vas,tdps,psl,swe,rad_sfc,snd,snc,mrso,sic,cape]
    varsp2 = ["huss_850","tasmax_2m","tasmin_2m","hfss_sfc","hfls_sfc","wap_500","ua_100","va_100","uas","vas","tdps","psl","swe","rad_sfc","snd","snc","mrso","sic","cape"]


    for name, ds in zip(varsp2, datasetsp2):
        print(name)
        path1 = "/glade/derecho/scratch/espstoch/2012case/"+workflow+"_native/p2/"+name+"/"+yyyy+"/"+mm+"/"
        if not os.path.exists(path1):
            os.system("mkdir -p "+path1)
        path2 = "/glade/derecho/scratch/espstoch/2012case/"+workflow+"_latlon/p2/"+name+"/"+yyyy+"/"+mm+"/"
        if not os.path.exists(path2):
            os.system("mkdir -p "+path2)
        filename = path1+f"{name}_{workflow}_{outdate}_00z_unst_d46_m0{ens}.nc"
        filename_regrid = path2+f"{name}_{workflow}_{outdate}_00z_d01_d46_m0{ens}.nc"
        ds.to_netcdf(filename,unlimited_dims="Time")
        os.system("cdo -P 1 -f nc5 remapcon,r360x181 -setgrid,mpas:"+grid_path+" "+filename+" "+filename_regrid)
