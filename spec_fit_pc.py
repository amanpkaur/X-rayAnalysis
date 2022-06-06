# coding: utf8
# To find all the cleaned event files and then do the spectral fits for all the grouped spectra.
# this script saves two files : 1) log file with all the fitting results+tstart for all the files and 2) the final flux and error
import os
import xspec
from xspec import *
import glob
from astropy.io import fits
from astropy.table import Table
import sys
import os.path
import pandas as pd
import subprocess
import fileinput
import numpy as np
path="/Users/aman/Desktop/3C273/2017/final"
res=[]
nores=[]
os.chdir(path)
name="3C273"
for dirpath, dirnames, filenames in os.walk(path):
    for filename in [f for f in filenames if f.endswith(("xpcw2po_cl.evt","xpcw3po_cl.evt","xpcw4po_cl.evt"))]:
        os.chdir(dirpath)
        data=fits.getdata(filename,2)
        df=pd.DataFrame(data=data)
        
        os.system("cp /Users/aman/Desktop/scripts/tclex_pc.xcm .")
        for line2 in fileinput.FileInput("tclex_pc.xcm", inplace=1):
            line2=line2.replace("fit_result.dat","pc_fitres"+str(filename[2:13])+".dat")
            print(line2) #required to replace
        os.chdir(path)
        dir=open("log"+str((filename[2:13]))+".txt","w")
        
        os.chdir(dirpath)
        times=subprocess.check_output("fkeyprint "+str(filename)+"[1] TST\n",shell=True)
        mjd=subprocess.check_output("fkeyprint "+str(filename)+"[1] MJD-OBS\n",shell=True)
        #     print(len(df),file=dir)
        # print(times,file=dir)
        
        cmd="xspec <<EOF\n@tclex_pc.xcm\nexit\nEOF"
        os.system(cmd)
        out=subprocess.check_output(cmd,shell=True).decode()
        lines=out.split("\n")
        for line in lines:
            print(line,file=dir)
        dir.close()
        
        if os.stat("pc_fitres"+str(filename[2:13])+".dat").st_size==0:
            nores.append((str(filename[2:13])))
            continue
        a=pd.read_csv("pc_fitres"+str(filename[2:13])+".dat",sep=" ",header=None)
        flux=10**(a[a.columns[4]]) # convert from log
        lerr=flux-10**(a[a.columns[5]])
        uerr=10**(a[a.columns[6]])-flux
        spec_ind=a[a.columns[1]]
        lerr_spec_ind=(spec_ind-a[a.columns[2]])/len(flux)
        uerr_spec_ind=(a[a.columns[3]]-spec_ind)/len(flux)
        dof=sum(a[a.columns[9]])
        chi2=sum(a[a.columns[8]])
        spec_ind_mean=np.mean(spec_ind)
        lerr_spec_ind_mean=np.mean(np.sqrt(sum(lerr_spec_ind**2)))
        uerr_spec_ind_mean=np.mean(np.sqrt(sum(uerr_spec_ind**2)))
        logflux_mean=(np.log10(np.mean(flux)))
        
        loglerr_mean=((np.sqrt(sum(lerr**2)))/len(flux))/np.mean(flux) # error propagation del(log f) = (del f)/f
        loguerr_mean=((np.sqrt(sum(uerr**2)))/len(flux))/np.mean(flux)
        res.append((int(filename.strip()[2:13]) ,(float(str(mjd).split("n")[5][10:31]),float(str(times).split("n")[5][10:31]),float(str(times).split("n")[6][10:31]),spec_ind_mean,lerr_spec_ind_mean,uerr_spec_ind_mean,logflux_mean,loglerr_mean,loguerr_mean,chi2,dof,chi2/dof)))

os.chdir(path)

with open(f'{name}_results_pc.csv',"w") as f:
    f.write("#Model employed tbabs*cflux*powerlaw\n#NH fixed =0.022cm-2\n#Emin=0.3 keV\n#Emax=10.0 keV\n \n#################################\n#MJD : Modified Julian Date\n#TSTART : Observation start time (sec)\n#TSTOP :Observation stop time (sec)\n#Spec_Ind : Spectral Index\n#Spec_Ind_lowErr : Lower bound error on Spectral Index\n#Spec_Ind_upErr : Upper bound error on Spectral Index\n#log10Flux : Log10(unabsorbed flux (erg/cm2/s))\n#lowErr_log10Flux : Lower bound error on log10(flux)\n#upErr_log10Flux : Upper bound error on log10(flux)\n#chi2: chi-squared\n#dof: degrees of freedom\n#red-chi2: reduced chi2 per degrees of freedom\n")
    f.write("MJD,TSTART,TSTOP,Spec_Ind,Spec_Ind_lowErr,Spec_Ind_upErr,log10Flux,lowErr_log10Flux,upErr_log10Flux,chi2,dof,redchi2\n")
    for i in range(len(res)):
        f.write(str(res[i])+'\n')

print(nores)
