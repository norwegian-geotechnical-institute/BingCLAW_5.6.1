#!/usr/bin/python
#
# Original python2 code by Sylfest Glimsdal
# Converted to python 3 by Steven Gibbons (2023)
#
 
import sys
import getopt
import netCDF4
import os
import numpy
import re
import math

from sct3 import io 
from scipy.io import netcdf

dd=io.io()

def help():
    print ("""
    This script reads fort.q files from GeoClaw.

    The content behind each column (VoellmyClaw): 
    1. col is the slidethickness, 
    2. col is the u component, 
    3. col is the v component, 
    4. col is the topo

    Three possible choices:

    1. read a single fort.q file and saves the wanted data in
    column <col> into a netcdf-file:

        claw2file.py -f <fort.q file> -c <column> -s <out netcdf filename>

    
    2. read column <col> of a set of fort.q files, and store it in a netcdf
    propagation file
    
        claw2file.py -d <dir to fort.q files> -c <column> -p <out netcdf filename>

    3. same as 2. but store the result as sink-source data (velocity of 
    surface/slidethickness) for GloBouss and not as netcdf. Computes also a 
    "time" file for gglo:

        claw2file.py -d <dir to fort.q files> -c <column> -g <out sink-source>

    Can change lengthscale by using [-x scale]
    
    4. reads the vales of the veloicity components (column numbers -u <col> and
    -v <col>) of a set of fort.q files, and store the total velocity in a 
    netcdf propagation file.
    
        claw2file.py -d <dir to fort.q files> -u <column> -v <column> -p <out netcdf>

    If the inputdata is from BingClaw use flag [-t] (Will for results form 
    VoellmyClaw divide hu and hv with slidethickness to get out the 
    velocities u and v. In the BingClaw the output at col. 2 and 3 is pure u and v). 

    Removing dry land (only tested for BoussClaw):

        Add -l <col number of eta+h, normally col 1>

    

    """)


def write_netcdf(filename,x,y,z):
    out=netcdf.netcdf_file(filename, 'w') #, format='NETCDF3_CLASSIC')
    out.createDimension("lon", len(x[:,0]))
    out.createDimension("lat", len(y[0,:]))
    lonv = out.createVariable("lon",'d',("lon",))
    lonv.unit='east'
    #lonv.point_spacing="even"
    latv = out.createVariable("lat",'d',("lat",))
    latv.unit='north'
    #latv.point_spacing="ueven"
    field = out.createVariable('Field','d',("lat","lon",))
    lonv[:]=x[:,0]
    latv[:]=y[0,:]
    field[:]=z.transpose()
    out.sync()
    out.close()

def write_netcdf_propf(inputdir, outfilename, model,remove_land, remove_land_col):
    #z is a 4D array x,y,z,t-levels
    field_unit="XXX"
    field="Field"
    field_longn="fn"
    xname='X'
    yname='Y'

    #count number of fort.t* files = nt, find start- and stop-time, dt
    line=open(os.path.join(inputdir,"claw.data"),'r').readlines()
    line=' '.join(line)
    print (line)
    p=re.compile(r'(-?\d*?\.\d+|-?\d+)\s*=:\s*num_output_times')
    nt=int(p.findall(line)[0])+1
    p=re.compile(r'(-?\d*?\.\d+|-?\d+)\s*=:\s*t0')
    tstart=float(p.findall(line)[0])
    print ("tstart",tstart)
    p=re.compile(r'(-?\d*?\.\d+|-?\d+)\s*=:\s*tfinal')
    tend=float(p.findall(line)[0])
    dt=(float(tend)-tstart)/(nt-1)


    #read first data-file:
    #[x,y,z]=d.read(list[0])
    basename=os.path.join(inputdir,"fort.q%04d")
    print (basename)
    [x,y,z]=read_fortq(col,model,basename %(0),remove_land, remove_land_col)

        
    ## Open the file for writing
    fname=outfilename
#    file = NetCDF.NetCDFFile(fname, 'w')
    file=netcdf.netcdf_file(fname, 'w')
    file.createDimension(xname, len(x[:,0]))
    file.createDimension(yname, len(y[0,:]))
    file.createDimension('TIME', nt)
    #Longitue
    var = file.createVariable(xname, 'd',(xname,))
    setattr(var, 'units', 'degrees_east')
    setattr(var, 'points_spacing', 'even')
    # Put some data in the variable
    var[:] = x[:,0]    

    #Latitude
    var = file.createVariable(yname, 'd', (yname,))
    setattr(var, 'units', 'degrees_north')
    setattr(var, 'points_spacing', 'even')
    # Put some data in the variable
    var[:] = y[0,:]

    #Time
    time = file.createVariable('TIME', 'd', ('TIME',))
    setattr(time, 'units', 'seconds')
    setattr(time, 'points_spacing', 'even')
    #Field
    var = file.createVariable(field, 'f', ('TIME',yname,xname))
    setattr(var, 'units', field_unit)
    setattr(var, 'long_name', field_longn)
    # Put some data in the variable

    k=0
    T=tstart
    maxval=z.copy()*0.0
    maxval[numpy.isnan(maxval)]=0

    for n in range(nt):
        [x,y,z]=read_fortq(col,model,basename %(n),remove_land, remove_land_col)
        time[k]=T
        if n==0:
            z[z==-150]=0
        for j in range(var.shape[2]):
            for i in range(var.shape[1]):
                var[k,i,j]=z[j,i]
                if not numpy.isnan(z[j,i]):
                    if maxval[j,i]<z[j,i]:
                        maxval[j,i]=z[j,i]

        k+=1
        T+=dt
        print ("T",T)
    if remove_land:
        maxval[maxval<0.00001]=numpy.nan
    #write out maxfield to asci-esri:
    #dd.write(propfile+"_max.asc",x,y,maxval, type="esri")
    write_netcdf(propfile+"_max.nc",x,y,maxval)
    ## Close the netCDF file
    
    file.close()
    print ("NetCDF-file",fname,"has been created")
    
def write_sinksource(inputdir, outfilename, model,remove_land, remove_land_col):

    #count number of fort.t* files = nt, find start- and stop-time, dt
    line=open(os.path.join(inputdir,"claw.data"),'r').readlines()
    line=' '.join(line)
    print (line)
    p=re.compile(r'(-?\d*?\.\d+|-?\d+)\s*=:\s*num_output_times')
    nt=int(p.findall(line)[0])+1
    p=re.compile(r'(-?\d*?\.\d+|-?\d+)\s*=:\s*t0')
    tstart=float(p.findall(line)[0])
    print ("tstart",tstart)
    p=re.compile(r'(-?\d*?\.\d+|-?\d+)\s*=:\s*tfinal')
    tend=float(p.findall(line)[0])
    dt=(float(tend)-tstart)/(nt-1)


    #read first data-file:
    basename=os.path.join(inputdir,"fort.q%04d")
    [x,y,zp]=read_fortq(col,model,basename %(0),remove_land, remove_land_col)
    zp[zp==-150]=0
        
    ## Open the file for writing
    out=outfilename+"F"
    file=open(out,'w')
    file.close()
    outtime=open(outfilename+".tim",'w')
    k=0
    T=tstart
    maxval=zp.copy()*0.0
    for n in range(1,nt):
        print ("n",n)
        T+=dt
        [x,y,z]=read_fortq(col,model,basename %(n),remove_land, remove_land_col)
        dd.write("tmp",x*lengthscale,y*lengthscale,lengthscale*(z-zp)/dt,type="gphov")
        os.system("cat tmp  >> "+out)
        outtime.write("%f\n"%(T))
        zp=z.copy()
    outtime.close()

    #adjust heading in gphov (include the number of timelevels):
    lines=open(out,'r').readlines()
    firstline=lines[0].strip()
    firstline=firstline[:-1]+str(nt)
    file=open(out,'w')
    file.write(firstline+"\n")
    for line in lines[1:]:
        file.write(line)
    file.close()

    print ("Files written with sink-source distribution "+str(out)+" and "+str(outfilename)+".tim.")


    
def read_topo_type3(file):
    from numpy import mod,zeros
    infile1=file
    
    fid= open(infile1,'r')

    mx=int(fid.readline().split()[0])
    my=int(fid.readline().split()[0])
    xlow=float(fid.readline().split()[0])
    ylow=float(fid.readline().split()[0])
    dx=float(fid.readline().split()[0])
    nodata=fid.readline()
    #nodata=-9999
    #print "herr",grid,amr,mx,my,xlow,ylow,dx,dy
    dy=dx
    val=zeros((mx,my))
    X=zeros((mx,my))
    Y=zeros((mx,my))
    
    for j in range(my):
        line=fid.readline().split()
        for i in range(mx):
            Y[i,j] = ylow+(my-j)*dy
            #Y[i,j] = ylow+j*dy
            X[i,j] = xlow+i*dx
            val[i,j]=float(line[int(i)])
    fid.close()

    return X,Y,val

def read_fortq(col,model,file,remove_land,remove_land_col): #model=voell or bing
    from numpy import mod,zeros
    
    infile=file
    #dir=dir
    infile1=infile #os.path.join(dir,infile)

    print ("read fort.q",infile1)
    fid= open(infile1,'r')
    
    grid=fid.readline()
    amr=fid.readline()
    mx=int(fid.readline().split()[0])
    my=int(fid.readline().split()[0])
    xlow=float(fid.readline().split()[0])
    ylow=float(fid.readline().split()[0])
    dx=float(fid.readline().split()[0])
    dy=float(fid.readline().split()[0])
    nodata=fid.readline()
    #nodata=-9999
    #print "herr",grid,amr,mx,my,xlow,ylow,dx,dy
    val=zeros((mx,my))
    X=zeros((mx,my))
    Y=zeros((mx,my))
    
    for j in range(my):
        for i in range(mx):
            #Y[i,j] = ylow+(my-j-.5)*dy
            Y[i,j] = ylow+(j+.5)*dy
            X[i,j] = xlow+(i+.5)*dx
            line=fid.readline().split()

            if vcol:
                #calculate velocity
                u=v=0.
                h=float(line[0])
                if abs(h)>0.000000000000000001:
                    if model=="voell":
                        u=float(line[int(col)])/float(line[0])
                        v=float(line[int(vcol)])/float(line[0])
                    else: #bing
                        u=float(line[int(col)])
                        v=float(line[int(vcol)])
                val[i,j]=math.sqrt(u*u+v*v)
            else:
                val[i,j]=float(line[int(col)])
            if remove_land:
                #print "line[0]",line[0]
                if float(line[remove_land_col])< 0.0001:
                    val[i,j]=numpy.nan
        fid.readline()
    fid.close()

    return X,Y,val





#---------------------------------------------------------------------------
#start of program:
#---------------------------------------------------------------------------


try: 
    arg=sys.argv[1];
except:
    help();
    sys.exit();
    
#default values:
fortq=False
col=1
vcol=False
outfile="test.nc"
propfile=False
sinksource=False
datadir="."
dtopo=False
lengthscale=1.0
model="voell"
remove_land=False
remove_land_col=0   #column of eta+h

#read input
options,args=getopt.getopt\
              (sys.argv[1:],'hf:c:s:p:d:b:u:v:g:x:tl:',\
               ['help','fortq=','col=', 'save=', 'prop=', 'datadir=','dtopo=','u=','v=','sinksource=','lengthscale=','bing','remove_land='])
for option,value in options:
    if option in ('-h','--help'):
        help()
        sys.exit(0)
    elif option in ('-f','--fortq'):
        fortq=str(value)
    elif option in ('-c','--col'):
        col=str(int(value)-1)
    elif option in ('-u'):
        col=str(int(value)-1)
    elif option in ('-v'):
        vcol=str(int(value)-1)
    elif option in ('-s','--save'):
        outfile=str(value)
    elif option in ('-p','--prop'):
        propfile=str(value)
    elif option in ('-d','--datadir'):
        datadir=str(value)
    elif option in ('-b','--dtopo'):
        dtopo=str(value)
    elif option in ('-g','--sinksource'):
        sinksource=str(value)
    elif option in ('-x','--lengthscale'):
        lengthscale=float(value)
    elif option in ('-t','--bing'):
        model="bing"
    elif option in ('-l','--remove_land'):
        remove_land=True
        remove_land_col=int(value)-1

print ("dtopo",dtopo,remove_land)




if fortq:
    #
    
    x,y,z=read_fortq(col,model,fortq,remove_land,remove_land_col)
    #read in x,y vectors and z array
    write_netcdf(outfile,x,y,z)

        
if propfile:
    print ("Write propfiles from fort.q files")
    #x,y,z=read_topo_type3(dtopo)
    #write_netcdf("heihei.nc",x,y,z)
    write_netcdf_propf(datadir,propfile,model,remove_land,remove_land_col)


if sinksource:
    print ("Write sinksource file for GloBouss")
    write_sinksource(datadir,sinksource, model,remove_land,remove_land_col)
