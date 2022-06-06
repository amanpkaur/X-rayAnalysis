#!/usr/bin/env python3
#!/usr/bin/python3
import os
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
heasarc = Heasarc()
os.system("export HEADASNOQUERY= \nexport HEADASPROMPT=/dev/null")

counts=[]
count_rate=[]
bkg_counts=[]
bkg_count_rate=[]
fluence=[]
countrate_xspec=[]
file_name=[]
mjd=[]
unc=[]
res=[]
counts_xspec=[]
### Run this part when you want to query heasrc for

#### Change these values#####
name="Mrk421"
mission = 'swiftmastr'
path="/Users/aman/Desktop/Mrk421/2021"
ap=45 #aperture
back_in=170
back_out=200
#############################

coords=coord.SkyCoord.from_name(name,frame="icrs")
ra=coords.ra.deg
dec=coords.dec.deg
cols = heasarc.query_mission_cols(mission=mission) # check what cols are available for a mission
table=heasarc.query_region(coords,mission=mission,radius='8 arcmin')
 # add fields='OBSID,START_TIME,NAME,RA,DEC,POINTING_MODE,EXPOSURE')
#table=table_tot[np.where((table_tot["START_TIME"] > 57177.0)&(table_tot["START_TIME"] < 57997.5))]
#table_tot[np.where((table_tot["START_TIME"])>57000.)]
#convert MJD to year_mm format#

os.chdir(path)
dates=table['START_TIME']
obsid=table['OBSID']
#type(dates.data) # lists the astropy column as an array
date=Time(dates.data,format='mjd')
date.format='fits'
yrmm=[x[0:7] for x in date.value] # retrieve the yyyy-mm part from the string
date=[x.replace("-","_") for x in yrmm] # the wget command needs "_" not "-"
"""
filename="download_{}.sh".format(name) #name of the file
f=open(filename,"w")
for i in range(len(table)):
    f.write(("wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks"
           " https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/{}//{}/xrt/".format(date[i],obsid[i]))+"\n")
    f.write(("wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks"
           " https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/{}//{}/auxil/".format(date[i],obsid[i]))+"\n")
f.close()

os.system("chmod a+x {}".format(filename)) # change the permissions
os.system("./{}".format(filename)) # run the script
"""

## running the pipeline and analysis now

for dirpath, dirnames, filenames in os.walk(path):
    for filename0 in [f for f in filenames if f.endswith(("xpcw2po_cl.evt.gz","xpcw3po_cl.evt.gz","xpcw4po_cl.evt.gz"))]:
        os.chdir(path)
        dir=filename0[2:13]
        filename=filename0[0:27]
        print(filename)
        cmd=f'xrtpipeline indir={path}/{dir} outdir={path}/final/{dir} srcra={ra} srcdec={dec} steminputs={dir} createexpomap=yes clobber=yes'
        #print(cmd)
        os.system(cmd)
        outpath=f'{path}/final'
        os.chdir(f'{outpath}/{dir}')
        file=open("pc_src.reg","w")
        file.write(f'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n fk5 \n')
        file.write(f'circle({ra},{dec},20")')
        file.close()
        filename=filename0[0:27]
        mjd_i=subprocess.check_output(f'fkeyprint {filename}[1] MJD-OBS \n',shell=True)
        mjd.append(str(mjd_i).split("n")[5][10:31])
        cmd="xselect <<EOF\nxel\nread event\n./\n"+str(filename)+"\nyes\nfilter pha_cutoff 30 1000\nextract image\nsave image image.img\nyes\nexit\nno\nEOF\n"
        os.system(cmd)

        """
        %%only for centroid finding%%
        center="xrtcentroid infile=image.img outfile=pos.txt outdir=./ interactive=no boxra='18 44 41.6' boxdec='-03 05 35.97' boxradius=4 calcpos=yes clobber=yes"
        os.system(center)
        xpos=[]
        ypos=[]
        with open("pos.txt","r") as pos:
            for line in pos:
                if line.startswith("X"):
                    xpos.append(float((line.strip().split()[2])))
                elif line.startswith("Y"):
                    ypos.append(float((line.strip().split()[2])))
        for line in fileinput.FileInput("pc_src.reg", inplace=1):
            line=line.replace('circle(500,500,45)','circle('+str(xpos[0])+','+str(ypos[0])+',45)')
            print(line)
        """
        
        pileup_check="xselect <<EOF\nxel\nread event\n./\n"+str(filename)+"\nyes\ncpd /xw\nfilter region pc_src.reg\nextract spectrum\nexit\nno\nEOF\n"
        os.system(pileup_check)
        with open("xselect.log","r") as f:
            for line in f:
                if line.startswith(" Spectrum         has"):
                    counts=(line.strip().split()[5])
                    res.append(float(counts))
                    pix_Scale=2.36 # arcsec 2.36 arcsec =1 pixel
                    if (1.0>= float(counts) > 0.5):
                        inner=2*pix_Scale
                        for line in fileinput.FileInput("pc_src.reg", inplace=1):
                            line=line.replace(f'circle({ra},{dec},{ap}")',f'annulus({ra},{dec},{inner}",{ap}")')
                            print(line)
                    elif (3.0>= float(counts) > 1.0):
                        inner=4*pix_Scale
                        for line in fileinput.FileInput("pc_src.reg", inplace=1):
                            line=line.replace(f'circle({ra},{dec},{ap}")',f'annulus({ra},{dec},{inner}",{ap}")')
                            print(line)
                    elif (4.0>= float(counts) > 3.0):
                        inner=5*pix_Scale
                        for line in fileinput.FileInput("pc_src.reg", inplace=1):
                            line=line.replace(f'circle({ra},{dec},{ap}")',f'annulus({ra},{dec},{inner}",{ap}")')
                            print(line)
                    elif (6.0>= float(counts) > 4.0):
                        inner=6*pix_Scale
                        for line in fileinput.FileInput("pc_src.reg", inplace=1):
                            line=line.replace(f'circle({ra},{dec},{ap}")',f'annulus({ra},{dec},{inner}",{ap}")')
                            print(line)
                    elif (8.0>= float(counts) > 6.0): #4.5 to 9 - 7 pixel
                        inner=7*pix_Scale
                        for line in fileinput.FileInput("pc_src.reg", inplace=1):
                            line=line.replace(f'circle({ra},{dec},{ap}")',f'annulus({ra},{dec},{inner}",{ap}")')
                            print(line)

                    elif (10.>= float(counts) > 8.0):# 10 and 12 - 11 pix
                        inner=9*pix_Scale
                        for line in fileinput.FileInput("pc_src.reg", inplace=1):
                            line=line.replace(f'circle({ra},{dec},{ap}")',f'annulus({ra},{dec},{inner}",{ap}")')
                            print(line)
                    elif (float(counts) > 10.): #
                        inner=12*pix_Scale
                        for line in fileinput.FileInput("pc_src.reg", inplace=1):
                            line=line.replace(f'circle({ra},{dec},{ap}")',f'annulus({ra},{dec},{inner}",{ap}")')
                            print(line)
                
                    elif (0.5>= float(counts) > 0.0):
                            continue
        f.close()
        #no time filtering needed for pc mode
        bfile=open("pc_back.reg","w")
        bfile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+"\n")
        bfile.write("image"+"\n")
        bfile.write(f'annulus({ra} {dec}),{back_in}",{back_out}")')
        bfile.close()
        cmdf="xselect <<EOF\nxel\nread event\n./\n"+str(filename)+"\nyes\nfilter region pc_src.reg\nextract spectrum \nsave spectrum spec\nyes\nextract curve\nset binsize 10\nextract curve\nsave curve pc_lc\nyes\n"
        cmdf+="clear region\nfilter region pc_back.reg\nextract curve\nset binsize 10\nextract curve\nsave curve pc_lcback\nyes\nextract spectrum\nsave spectrum back\nyes\nexit\nno\nEOF\n"
        os.system(cmdf)
        os.system(f'xrtmkarf outfile=PC.arf expofile={filename[0:20]}_ex.img phafile=spec.pha srcx=-1 srcy=-1 psfflag=yes clobber=yes > arf.log')
        
        for line in open("arf.log","r"):
            if "ON AVERAGE" in line:
                a=[float(x) for x in re.findall(re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\*[0-9]+)?'), line)]
                fluence.append(a[0])
            
            if "/data/swift/xrt/cpf/rmf" in line:
                print("wget "+str(line[34:126])+"")
                os.system("wget "+str(line[34:126])+"")
                os.system("mv sw*.rmf PC.rmf") # same rmf for this
        grp_cmd="grppha <<EOF\nspec.pha\n!fspec.pha\nchkey BACKFILE back.pha\nchkey RESPFILE PC.rmf\nchkey ANCRFILE PC.arf\ngroup min 20\nexit\nEOF\n"
        os.system(grp_cmd)
        grp_cmd_nogroup="grppha <<EOF\nspec.pha\n!fspec_cstat.pha\nchkey BACKFILE back.pha\nchkey RESPFILE PC.rmf\nchkey ANCRFILE PC.arf\ngroup min 1\nexit\nEOF\n"
        os.system(grp_cmd_nogroup)
        file_name.append(filename)
#write counts from xspec
        new=os.getcwd()
        xspec_cts="xspec <<EOF\ndata fspec.pha\nexit\nEOF"
        out=subprocess.check_output(xspec_cts,shell=True).decode()
        for line in out.split("\n"):
            if "cts/s" in line:
                counts_xspec.append(float(line.split(" ")[7]))

os.chdir(path)
with open("pileup_counts_new.txt","w") as f:
    f.write("ID,MJD,piled_up_cts,fluence,xspec_counts\n")
    for i in range(len(res)):
        f.write(f'{file_name[i]},{mjd[i]},{res[i]},{fluence[i]},{counts_xspec[i]} \n')

    
        
        #f.write(str(file_name[i])+","+str(mjd[i])+","+str(res[i])+","+str(fluence[i])+","+str(counts_xspec[i])+"\n")

        
