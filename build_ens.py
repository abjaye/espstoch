#!/usr/bin/env python3
import os, sys
import glob, shutil
from datetime import timedelta, datetime
import argparse
import subprocess

def parse_command_line(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--date",
                        help="Specify a start Date")
    parser.add_argument("--ensemble-start",default=0,
                        help="Specify the first ensemble member")
    parser.add_argument("--ensemble-end",default=9,
                        help="Specify the last ensemble member")

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

    return date.strftime("%Y-%m-%d"),int(args.ensemble_start),int(args.ensemble_end)

if __name__ == '__main__':
    date, ensemble_start, ensemble_end = parse_command_line(sys.argv)

    workflow = "mpass2s"
    baseroot = "/glade/derecho/scratch/espstoch/2012case"
    caseroot = os.path.join(baseroot,"{}_{}".format(workflow,date)+".00")

    basemonth = int(date[5:7])
    baseyear = int(date[0:4])
    baseday = int(date[8:9])

    base_date_str = date  # starting date
    num_days = 46  # number of days

    startday = datetime.strptime(base_date_str, '%Y-%m-%d').date()
    endday = startday + timedelta(days=num_days)
    startday_str = startday.strftime('%Y-%m-%d')
    endday_str = endday.strftime('%Y-%m-%d')

    startval = "00"
    nint = len(startval)

    for i in range(ensemble_start, ensemble_end+1):
        member_string = '{{0:0{0:d}d}}'.format(nint).format(i)
        caseroot = caseroot[:-nint] + member_string
        print(caseroot)

        os.system("rm -rf "+caseroot)

        if not os.path.isdir(baseroot):
            os.mkdir(baseroot)
        if not os.path.isdir(caseroot):
            os.mkdir(caseroot)

        subprocess.run(["python", "/glade/work/espstoch/era5_to_int/era5_to_int_abby.py", \
                        "--path", "/glade/work/espstoch/cds_test","--po",caseroot,"-i","-e", \
                        str(i),date+"_00"])

        #now need to link the SST update, static and grid files to each ens directory

        sst_update = "/glade/work/espstoch/mpas_grid_input/x1.163842.sfc_update.nc"
        os.system("ln -sf "+sst_update+" "+caseroot+"/.")
        grid_file = "/glade/work/espstoch/mpas_grid_input/x1.163842.grid.nc"
        os.system("ln -sf "+grid_file+" "+caseroot+"/.")
        static_file = "/glade/work/espstoch/mpas_grid_input/x1.163842.static.nc"
        os.system("ln -sf "+static_file+" "+caseroot+"/.")

        # copy model run dir to ens dirs

        model_dir = "/glade/work/espstoch/mpas_workflow/mpas_base/"
        os.system("cp "+model_dir+"* "+caseroot+"/.")

        # change dates in namelists

        with open(caseroot+"/namelist.init_atmosphere_temp") as fin, open(caseroot+"/namelist.init_atmosphere","w") as fout:
            input_lines = fin.readlines()
            for line in input_lines:
                if "    config_start_time = '" in line:
                    fout.write("    config_start_time = '{}_00:00:00'\n".format(startday_str))
                else:
                    fout.write(line)

        with open(caseroot+"/namelist.atmosphere_temp") as fin, open(caseroot+"/namelist.atmosphere","w") as fout:
            input_lines = fin.readlines()
            for line in input_lines:
                if "    config_start_time = '" in line:
                    fout.write("    config_start_time = '{}_00:00:00'\n".format(startday_str))
                else:
                    fout.write(line)


