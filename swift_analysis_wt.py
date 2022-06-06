# run heainit import os
import os.path
import fileinput
import re
import subprocess
from astropy import coordinates as coord
import sys
from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord
from astropy.time import Time
import numpy as np
import glob
from astropy.io import fits
from astropy.table import Table
import pandas as pd
heasarc = Heasarc()
counts_xsel=[]
file_name=[]
mjd=[]
fluence=[]
counts_xspec=[]
os.system("export HEADASNOQUERY= \nexport HEADASPROMPT=/dev/null")
#### Change these values as needed #####
name="3C273"
mission = 'swiftmastr'
coords=coord.SkyCoord.from_name(name,frame="icrs")
ra=coords.ra.deg
dec=coords.dec.deg
pix_Scale=2.36
box_l=40*pix_Scale # pixels
box_w=20*pix_Scale
back_in=100*pix_Scale
back_out=200*pix_Scale
path="/Users/test"

##### THE DATA DOWNLOAD SCRIPT IS IN THE "swift_analysis.py" script ################

for dirpath, dirnames, filenames in os.walk(path):
    for filename0 in [f for f in filenames if f.endswith("xwtw2po_cl.evt.gz")]:
        os.chdir(path)
        dir=filename0[2:13]
        filename=filename0[0:27]
        swfilename=filename0[0:13]
        outpath=f'{path}/final/{dir}'
        if not os.path.exists(outpath):
            cmd=f'xrtpipeline indir={path}/{dir} outdir={path}/final/{dir} srcra={ra} srcdec={dec} steminputs={dir} createexpomap=yes clobber=yes'
            os.system(cmd)
            os.chdir(outpath)
        else:
            os.chdir(outpath)
            
 
            data=fits.getdata(filename,2) # 2 means we read GTIs
            df=pd.DataFrame(data=data)
        # roll angle and time info are in s.mkf.gz files under auxil directory
        # match the time given in GTIs in the cl_evt file for WT mode to the time (nearest) in the auxil file and read off the roll angle at that time
            hdulist=fits.open(f'{dir}s.mkf')
            roll=hdulist['FILTER'].data['ROLL']
            time=hdulist['FILTER'].data['TIME']
            count=0
            mjd_i=[]
            counts_xsel_i=[]
            fluence_i=[]
            counts_xspec_i=[]
            for l in range(len(df)):
                df.iloc[l:l+1].to_csv('gti'+str(l+1)+'.txt',index=None,header=None,mode="w")
                gti=np.loadtxt('gti'+str(l+1)+'.txt',delimiter=",")
                t=int(gti[0])
                i=np.argwhere(time.astype(int) == t)[0] # index where time corresponds to gti start..
                ang=roll[i]-90
                ang=ang[0]
                file=open("src"+str(l+1)+".reg","w")
                file.write(f'global color=magenta dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+"\n")
                file.write("fk5"+"\n")
                file.write(f'box({ra},{dec},{box_l}",{box_w}",{ang})') # for some reason the correct orientation is rollangle-90 in the displayed image..
                file.close()
        
                # extract the image and spectrum using the extracted info. Please select only 0.3 to 10 kev for our work
                cmd="xselect <<EOF\nxel\nread event\n./\n"+str(filename)+"\nyes\nfilter pha_cutoff 30 1000\nfilter time file gti"+str(l+1)+".txt\nextract image\nsave image image"+str(l+1)+".img\nyes\nexit\nno\nEOF\n"
                os.system(cmd)
                # remove the first 150 seconds of data in WT mode and extract light curve and spectrum and check for pile up
            
                pileup_check="xselect <<EOF\nxel\nread event\n./\n"+str(filename)+"\nyes\ncpd /xw\nfilter region src"+str(l+1)+".reg\nfilter time file gti"+str(l+1)+".txt\nset binsize 10\nextract curve\nfilter time scc\nquit\n150,r\nx\nextract spectrum\nexit\nno\nEOF\n"
                os.system(pileup_check)
                    #save spectrum spec"+str(l+1)+"\nyes\nexit\nno\nEOF\n"
                with open("xselect.log","r") as f:
                    for line in f:
                        if line.startswith(" Spectrum         has"):
                                # print(line)
                            rate=(line.strip().split()[5])
                            counts_xsel_i.append(float(rate))
                            pix_Scale=2.36
                        
                            if (300 >= float(rate) >= 100): # based on Romani's work these pile up values are provided as follows and replace these in source region files:
                                inner=1*pix_Scale
                                for line in fileinput.FileInput("src"+str(l+1)+".reg", inplace=1):
                                    line=line.replace(f'box({ra},{dec},{box_l}",{box_w}",{ang})',f'box({ra},{dec},{box_l}",{box_w}",{ang},{inner},{box_w}",{box_l}",{box_w}",{ang})')
                                    print(line)
                            elif (400 >= float(rate) > 300.):
                                inner=2*pix_Scale
                                for line in fileinput.FileInput("src"+str(l+1)+".reg", inplace=1):
                                    line=line.replace(f'box({ra},{dec},{box_l}",{box_w}",{ang})',f'box({ra},{dec},{box_l}",{box_w}",{ang},{inner},{box_w}",{box_l}",{box_w}",{ang})')
                                    print(line)
                            elif (float(rate) > 400.):
                                inner=4*pix_Scale
                                for line in fileinput.FileInput("src"+str(l+1)+".reg", inplace=1):
                                    line=line.replace(f'box({ra},{dec},{box_l}",{box_w}",{ang})',f'box({ra},{dec},{box_l}",{box_w}",{ang},{inner},{box_w}",{box_l}",{box_w}",{ang})')
                                    print(line)
                            elif (float(rate)<100.):
                                continue
                    f.close()
                bfile=open("back"+str(l+1)+".reg","w") # create the background file
                bfile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+"\n")
                bfile.write("fk5"+"\n")
                bfile.write(f'box({ra},{dec},{back_in}",{box_w}",{back_out}",{box_w}",{ang})')
                bfile.close()
                # now we have modified the source region file according to pile up correction and we are ready to extract the actual light curve and spectrum and save these curves and spectra
                cmdf="xselect <<EOF\nxel\nread event\n./\n"+str(filename)+"\nyes\ncpd /xw\nfilter region src"+str(l+1)+".reg\nfilter time file gti"+str(l+1)+".txt\nset binsize 10\nextract curve\nfilter time scc\nquit\n150,r\nx\nextract curve\nsave curve lc"+str(l+1)+"\nyes\nextract spectrum\nsave spectrum spec"+str(l+1)+"\nyes\n"
                cmdf+="clear region\nfilter region back"+str(l+1)+".reg\nextract curve\nsave curve back_lc"+str(l+1)+"\nyes\nextract spectrum\nsave spectrum back"+str(l+1)+"\nyes\nexit\nno\nEOF\n"
                print(cmdf)
                os.system(cmdf)
                # note on BACKSCAL = ratio of area of extracted region/total detector area
                os.system("fparkey 800 spec"+str(l+1)+".pha BACKSCAL") # no. of pixels where the source is 40*20
                os.system("fparkey 2000 back"+str(l+1)+".pha BACKSCAL")# no. of total included pixels 20*50*2
                # create arf file to be used for spectral fitting later
                os.system("xrtmkarf outfile=WT"+str(l+1)+".arf  expofile="+str(filename0[0:21])+"ex.img phafile=spec"+str(l+1)+".pha srcx=-1 srcy=-1 psfflag=yes clobber=yes > arf"+str(l+1)+".log")
                for line in open("arf"+str(l+1)+".log","r"):
                    if "ON AVERAGE" in line:
                        a=[float(x) for x in re.findall(re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\*[0-9]+)?'), line)]
                        print(a)
                        if a[0]>0.0:
                            fluence_i.append(a[0])
                            mjd_in=subprocess.check_output("fkeyprint image"+str(l+1)+".img MJD-OBS\n",shell=True)
                            mjd_i.append(float(str(mjd_in).split("n")[5][10:31]))
                
                    if "/data/swift/xrt/cpf/rmf" in line:
                        link=line[34:125]
                        print(link)
                        os.system(f'wget {link}')# copy the rmf to the current folder which was used to create the arf
                        os.system("mv sw*.rmf WT.rmf")
                        # group all the required files to create a final spectrum which contains info about arf, rmf, background as well as any binning to be applied
                        times=subprocess.check_output("fkeyprint spec"+str(l+1)+".pha[0] TOTCTS\n",shell=True)
                        if [float(s) for s in times.split() if s.isdigit()][1]>0.:
                            count+=1
                            grp_cmd="grppha <<EOF\nspec"+str(l+1)+".pha\n!fspec"+str(count)+".pha\nchkey BACKFILE back"+str(l+1)+".pha\nchkey RESPFILE WT.rmf\nchkey ANCRFILE WT"+str(l+1)+".arf\ngroup min 20\nexit\nEOF\n"
                            os.system(grp_cmd)
                
                            xspec_cts="xspec <<EOF\ndata fspec"+str(l+1)+".pha\nexit\nEOF"
                            out=subprocess.check_output(xspec_cts,shell=True).decode()
                            for line in out.split("\n"):
                                if "cts/s" in line:
                                    counts_xspec_i.append(float(line.split(" ")[7]))
            
                        else:
                            grp_cmd="grppha <<EOF\nspec"+str(l+1)+".pha\n!noct_finalspec"+str(l+1)+".pha\nchkey BACKFILE back"+str(l+1)+".pha\nchkey RESPFILE WT.rmf\nchkey ANCRFILE WT"+str(l+1)+".arf\ngroup min 20\nexit\nEOF\n"
                            os.system(grp_cmd)
        
        

        
                    if len(mjd_i)>0.0:
                        mjd.append((mjd_i[-1]+mjd_i[0])/2.)
                        file_name.append(str(filename[2:13]))
                        counts_xsel.append(np.mean(counts_xsel_i))
                        counts_xspec.append(np.mean(counts_xspec_i))
                        fluence.append(np.mean(fluence_i))
                    else:
                        continue
                os.chdir(path)

                with open("allcounts_wt_test.txt","w") as f:
                    f.write("ID,MJD,fluence,xspec_counts\n")
                    for i in range(len(fluence_i)):
                        f.write(str(mjd_i[i])+","+str(fluence_i[i])+","+str(counts_xspec_i[i])+"\n")
                f.close()

                with open("pileup_counts_wt_test.txt","w") as f:
                    f.write("ID,MJD_center,fluence,xselect_cts,xspec_counts,count_rate\n")
                    for i in range(len(fluence)):
                        f.write(str(file_name[i])+","+str(mjd[i])+","+str(fluence[i])+","+str(counts_xsel[i])+","+str(counts_xspec[i])+","+str(counts_xspec[i]*100./fluence[i])+"\n")
                f.close()

