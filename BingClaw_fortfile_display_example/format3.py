#! /usr/bin/python
#
# Original python2 code by Sylfest Glimsdal
# Converted to python 3 by Steven Gibbons (2023)
#


import sys
import getopt
import os.path
import numpy
import netCDF4
from sct3 import io
from sct3.Grid2D import Grid2D as G

#netcdf:
#from Scientific.IO import NetCDF
#from Numeric import *
# from scipy.io import netcdf
# import netcdf
from numpy import *

def help():
    print ("""
    Reads, writes, and manipulates  3D files:
    -----------------------------------------
    
    mtv, xyz, gphov, geoclaw (esri), ESRI-ascii (also extracting fields from
    gphov), and most grid files.
    
    Can also read,extract fields at given times and write Most propagation
    files (netcdf. By giving the extention '.nc' the data will be written
    into netcdf format.

    Useage:
    -------
    
    format.py -f|--infile filename -s|--save output [-h|--help] 

     
    Options:
    
    [-r|--resolution extrp_factor - for output]
    [-t|--type output type: xyz|mtv|gphov|most|esri|gcl or by extension of -s
          option] 
    [-k|--skredp 'stem4depth stem4eta'] only for gphov (first field is depth)
    [-[x|y|z] or --[x|y|z]sc scalingfactor]
    [-[X|Y|Z] value or --[x|y|z]add value]
    [-c|--plotcmds 'cmds' (only for mtv)]
    [-d|--domain ] xstart/xend/ystart/yend
          scrinking the domain or expanding it giving zero as field value
          in the new parts of the domain
    [-g|--gauges file] two (posx,posy)  or three column (name,posx,posy) 
          file with positions of points to be plotted. Numbering of points: 
          tow column file - points are numbered subsequently starting with 
          1,2,3,4 ... For three column file, the first column must be the 
          name of the point. (Only for the mtv format)
          
    [-i|--propfield] name. Name of field to be saved in (or read from)
          the propagation files in NetCDF format for ComMIT/MOST.

             EG:

                 HA (eta - centimeters)
                 UA (u - centimeters/second)
                 VA (v - centimeters/second)
    [-m|--ncxname] name. Name of x-coordinate of netcdf files to be
          read or written
    [-n|--ncyname] name. Name of y-coordinate of netcdf files to be
          read or written
    [-q|--sift] flag for ComMIT - sift files (time range but not time varying)
             
    [-v|--variables] if the inputfile is a NetCDF-file (-f option), this option
          writes out the variables of the NetCDF-file
    [-a|--arrivaltimes] computes the arrivaltimes based on a netCDF file
          (propagation file). Arrivaltime is defined as a surface elevation
          above the threshold given by the -b (--treshold_arrival)  option.
          Default value is 0.1.
    [-F|--sortcol 0|1] defines the sorting of xyz output files. 
          0 sorted along x, 1 sorted along y.

    Scaleing and adding values -x/-X ... option:
    --------------------------------------------
    If a coord. is both scaled and added to a value the coord. is first
    scaled then the value is added by default. If you want to do it in
    oposite order use this option:
                [-u|--add_first]

   
    
    Adding a value on land (only for -k/--skredp option):
    -----------------------------------------------------

    [-l|--landvalue value] adds a value for eta where there is no water ...

    Producing propagation files for MOST (NetCDF):
    ------------------------------------

    [-f <file>] files (xyz or gphov format) are listed in <file> in
                correct order
    [-T|--timeinstr_propf 'tstart tstop step'] time for the input files
    [-i|--field <type>] type of input (eta: HA, velocities: u: UA, v: VA)
    [-p|--strech <value>] the grid is stretched by chaning the coordinates for
          the first and last five coordinates in both x and y direction with a
          factor = <value>.
    
    Scaling (domain and value) is possible through the  -x,-y,-z (scaling)
    -X,-Y,-Z (adding values) and -d (changing the domain)

    The name of the NetCDF-file: the stem is given by the -s option, added
    the type of input (HA,UA,VA) and given an extention '.nc'. E.g. 'aatHA.nc'
    
    NetCDF I/O:
    -----------

    Extracting a data set at a given time from a netcdf file::
        -e|--extr_time 'time'

    For writing out the content of a NetCDF-file the -s option gives the
    stem of the outputfiles, while an integer is added (timelevel) as well
    as the extension (-t type of output e.g. 'xyz', test001.xyz)

    NetCDF:
    -------

    Give the extention '.nc' to the outputfile. By default the data is
    supposed to be given in lon/lat. By using the flag -C the data will
    be saved as x/y. (e.g. utm)

    """)


def read_netcdf(xname,yname,fname,prntvar,pftype,sortcol):
    #####################################################
    #
    #
    #   5/5-2011: bgynt aa legge inn mulighet for aa velge
    #             navn paa parametre som skal hentes ut.
    #             f.eks. for ComMIT sift filer ...
    #
    #             ikke ferdig!!!! (xname, yname)...
    #
    ######################################################
    tmvar=1 #time is a variable?
    file = netCDF4.Dataset(fname, 'r')
    string=""
    for f in file.variables.keys():
        string+=f+" "
    print ("This is the variables of the NetCDF-file", fname,":",string)
    if prntvar:
        sys.exit()
    #find correct pftype:
    keys=file.variables.keys()
    if pftype not in keys:
        if "HA" in keys:
            pftype="HA"
        elif "z" in keys:
            pftype="z"
        elif "Z" in keys:
            pftype="Z"
        elif "Field" in keys:
            pftype="Field"
    var=file.variables[pftype]
    #check variable name for time
    if 'time' in file.variables:
        tname = "time"
    elif 'TIME' in file.variables:
        tname = "TIME"
    else:
        tname = False
        #check for upper/lowercase:
    if not xname or not yname:
        lat='LAT'
        lon='LON'
        if 'lat' in file.variables:
            lat='lat'
            lon='lon'
        elif 'LAT' in file.variables:
            lat='LAT'
            lon='LON'
        elif 'x' in file.variables:
            lat='y'
            lon='x'
        elif 'X' in file.variables:
            lat='Y'
            lon='X'
    else:
        lon=xname
        lat=yname
    nlat=file.dimensions[lat].size
    nlon=file.dimensions[lon].size
    print ("nlat,nlon",nlat,nlon,lat,lon)
    if not tname or commit_sift:
        print ("TIME is not a variable")
        tmvar=0
    if tmvar:
        print ("TIME: ",file.variables[tname][0]," to ",file.variables[tname][-1],", dt=",file.variables[tname][1]-file.variables[tname][0])

        ntim=file.dimensions[tname]

    x=numpy.zeros([nlon,nlat],float)
    print ("hh")
    y=numpy.zeros([nlon,nlat],float)
    z=numpy.zeros([nlon,nlat],float)
    data=file.variables[lon][:]
    for i in range(nlat):
        x[:,i]=data
    data=file.variables[lat][:]
    for i in range(nlon):
        y[i,:]=data
    data=var[:]

#    from sct.Grid2D import Grid2D as G


    k=0
    if tmvar:
        if extrtm:
            #write out data at given time:
            i=0
            for t in file.variables[tname][:]:
                if t >=extrtmval:
                    #if t-1 is closer to extrtmval:
                    diff1=abs(t-extrtmval)
                    tprev=file.variables[tname][i-1]
                    diff2=abs(tprev-extrtmval)
                    if diff1<diff2:
                        z=data[i,:,:]
                    else:
                        z=data[i-1,:,:]
                        t=tprev
                    z=numpy.transpose(z)
                    break
                i+=1
            print ("Data at time",t,"is now ready to be written to file...")
            if coarse:
                g=G(x,y,z)
                g.extrap(coarse)
                x=g.x
                y=g.y
                z=g.z
            g=G(x,y,z)
            #z=g.replace_value(g.z,-500,500,0)
            d.write(save,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt)
            # d.write(save,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt,sortcol=sortcol)
        #tsunami arrival time
        elif arrival:
            #mask:
            mask=data[k,:,:]*0.0+1.0
            #arrivaltime:
            arrivalarr=data[k,:,:]*float("NaN")
            #t=file.variables[tname].getValue()[k]
            
            for t in file.variables[tname][:]:
                #find indices where the value is above the threshold
                #mask out old areas:
                indices=numpy.greater(abs(data[k,:,:])*mask,threshold_arrival)
                #update mask:
                mask=numpy.where(indices,0.0,mask)
                #update arrivaltimes:
                arrivalarr=numpy.where(indices,t,arrivalarr)

                k+=1

            z=numpy.transpose(arrivalarr)
            z[numpy.isnan(z)]=t
            d.write(save,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt)
            # d.write(save,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt,sortcol=sortcol)

        else:
            #write out all fields
            print ("else",file.variables[tname][:])
            for t in file.variables[tname][:]:
                print ("TIME", t)
                for i in arange(nlat):
                    for j in arange(nlon):
                        #print ("heiii",i,j,k,len(z))
                        #print ("z",z)
                        #print ("data",data)
                        z[j,i]=data[k,i,j]
                k+=1
                if not extrtm:
                    print ("type",typ)
                    out=save+"%(k)04d" % vars() +typ
                    if coarse:
                        g=G(x,y,z)
                        g.extrap(coarse)
                        x=g.x
                        y=g.y
                        z=g.z
                        d.write(out,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt)
                        # d.write(out,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt,sortcol=sortcol)
                    else:
                        d.write(out,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt)
                        # d.write(out,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt,sortcol=sortcol)
    elif commit_sift:
        z=data[:,:]
        z=numpy.transpose(z)
        d.write(save,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt)
        # d.write(save,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt,sortcol=sortcol)
        
    else:
        print ("no time-varying field...")
        for i in arange(nlat):
            for j in arange(nlon):
                z[j,i]=data[i,j]
#        from sct.Grid2D import Grid2D as G
        g=G(x,y,z)
        z=g.replace_value(z,-99999,99999,float("NaN"))
        if coarse:
            g=G(x,y,z)
            g.extrap(coarse)
            x=g.x
            y=g.y
            z=g.z
 
        d.write(save,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt)
        # d.write(save,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt,sortcol=sortcol)




def mtv_gauges(no,x,y):
    
    entry="\n@point x1=%(x)f y1=%(y)f label=\"%(no)s\"" %vars()
    return entry
    
def write_netcdf(list,start,stop,step,field,xname,yname):
    #print (list,start,stop,step,field)
    
    field_unit="centimeters/seconds"
    if field=="HA":
        field_unit="centimeters"
        field_longn="Wave Amplitude"
    elif field=="UA":
        field_longn="Velocity Component along Longitude"
    else:
        field_longn="Velocity Component along Latitude"

    #from sct.Grid2D import Grid2D as G
    #read first data-file:
    [x,y,z]=d.read(list[0])
    if coarse:
        g=G(x,y,z)
        x,y,z=g.extrap(coarse)
    if x0: #extend domain with zeros or make a zoom (-d option)
#        from sct.Grid2D import Grid2D as G
        g=G(x,y,z)
        x,y,z=g.zoom(x0,y0,x1,y1)
        print ("zooming propfiles")
        
    ## Open the file for writing
    fname=save+field.lower()+".nc"
    file = netcdf.netcdf_file(fname, 'w')

    #Dimensions
    if not xname:
        xname='LON'
    if not yname:
        yname='LAT'
    file.createDimension(xname, len(x[:,0]))
    file.createDimension(yname, len(y[0,:]))
    file.createDimension('TIME', len(list))
    #Longitue
    var = file.createVariable(xname, 'd',(xname,))
    setattr(var, 'units', 'degrees_east')
    setattr(var, 'points_spacing', 'even')
    # Put some data in the variable
    var[:] = (x[:,0]*xsc)+xadd
    var[0:4]=var[0:4]-stretch
    var[-5:]=var[-5:]+stretch
    

    #Latitude
    var = file.createVariable(yname, 'd', (yname,))
    setattr(var, 'units', 'degrees_north')
    setattr(var, 'points_spacing', 'even')
    # Put some data in the variable
    var[:] = (y[0,:]*ysc)+yadd
    var[0:4]=var[0:4]-stretch
    var[-5:]=var[-5:]+stretch

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
    for f in list:
        [x,y,z]=d.read(list[k])
        if coarse:
            g=G(x,y,z)
            x,y,z=g.extrap(coarse)
        if x0: #extend domain with zeros or make a zoom (-d option)
            g=G(x,y,z)
            x,y,z=g.zoom(x0,y0,x1,y1)
        time[k]=T
        for i in range(var.shape[1]):
           for j in range(var.shape[2]):
               if z[j,i]<-50000000:
                   var[k,i,j]=0
               else:
                   var[k,i,j]=z[j,i]*zsc+zadd
        #var[k,:,:]=z[:,:]
        k+=1
        T+=tstep

    ## Close the netCDF file
    
    file.close()
    print ("NetCDF-file",fname,"has been created")


      
#---------------------------------------------------------------------------
#start of program:
#---------------------------------------------------------------------------


try: 
    arg=sys.argv[1];
except:
    help();
    sys.exit();
    
#default values:
depth=False
eta=False
coarse=False
xsc=1
ysc=1
zsc=1
xadd=0
yadd=0
zadd=0
typ=None
txt=""
skredp=0
landvalue=0
gauges=False
#netcdf prop. files for most:
propf  = [] #propagation files (xyz,mtv,gphov...)
pftype = "HA"
xname = False
yname = False
tname = False
#flag for commit sift files:
commit_sift=False
#time for extracting data from a netcdf-file:
extrtm=False
#numbering
step   = 1
start  = 1
stop   = 51
#corresponding time
tstart = 0 #0.0   403.030261993  6.01537704468'
tstep  = 6.01537704468 #12.031 #31.1688308716
tstop  = 403.030261993 #601.6 #14399.9998627
#write out the variable of a netcdf-file
printvar=False
readnc=False
#adjust domain:
x0=False
#stretching first and last five coord. in both spatial direction by
#adding/subtracting this value:
stretch=0
#scaling/adding or adding/scaling:
add_first=False
#arrival times:
arrival=False
#arrival defined as first values above threshold:
threshold_arrival=0.1
degr=True
sortcol=0  #For output xyz: if 0 sort first along x, if 1 first along y
#read input
options,args=getopt.getopt\
              (sys.argv[1:],'hf:s:k:t:r:x:y:z:c:X:Y:Z:l:g:i:vd:T:e:p:uCm:n:ab:qF:',\
               ['help','infile=','save=','skredp:','type=','resolution=','xsc=','ysc=','zsc=','plotcmds=','xadd=','yadd=','zadd=','landvalue=','gauges=','propfield=','printvars','domain=','timeinstr_propf=','extr_time=','stretch=','add_first','cartesian','ncxname=','ncyname=','arrivaltime','threshold_arrival=','sift=','sortcol='])
for option,value in options:
    if option in ('-h','--help'):
        help()
        sys.exit(0)
    elif option in ('-f','--infile'):
        infilename=str(value)
        print ("infilename", infilename)
        #file = netcdf.netcdf_file(infilename, 'r')
        #is it a nc-file?
        #file = netcdf.netcdf_file(infilename, 'r')
        try:
            print ("readnc...",readnc)
            print ("Is ",infilename,"a NetCDF-file?")
            #file = netcdf.netcdf_file(infilename, 'r')
            file = netCDF4.Dataset(infilename, 'r')
            readnc=True
        except:
            pass
        #test if input-file is a list of files -> makeing propagation file
        #for ComMIT/MOST
        propf=[] #list of filenames for propagation file
        if not readnc:
            try:
                for f in open(infilename).readlines():
                    f=f.strip().split()[0]
                    if os.path.isfile(f):
                        propf.append(f.strip().split()[0])
            except:
                prof=[]
                pass

    elif option in ('-s','--save'):
        save=str(value)
    elif option in ('-t','--type'):
        typ=str(value)
        print ("typ",typ)
    elif option in ('-k','--skredp'):
        skredp=1
        depth=str(value.strip().split()[0])
        print ("depth",depth)
        eta=str(value.strip().split()[1])
        print ("eta",eta)
    elif option in ('-r','--resolution'):
        coarse=int(value)
    elif option in ('-X','--xadd'):
        xadd=float(value)
    elif option in ('-Y','--yadd'):
        yadd=float(value)
    elif option in ('-Z','--zadd'):
        zadd=float(value)
    elif option in ('-x','--xsc'):
        xsc=float(value)
    elif option in ('-y','--ysc'):
        ysc=float(value)
    elif option in ('-z','--zsc'):
        zsc=float(value)
    elif option in ('-u','--add_first'):
        add_first=True
    elif option in ('-c','--plotcmds'):
        txt=str(value)
    elif option in ('-l','--landvalue'):
        landvalue=float(value)
    elif option in ('-g','--gauges'):
        gauges=str(value)
    elif option in ('-p','--stretch'):
        stretch=float(value)
    elif option in ('-i','--propfield'):
        pftype=str(value)
    elif option in ('-v','--printvars'):
        printvar=True
    elif option in ('-C','--cartesian'):
        degr=False
    elif option in ('-q','--sift'):
        commit_sift=True
    elif option in ('-T','--timeinstr_propf'):
        tm=value.strip().split()
        tstart=float(tm[0])
        tstop=float(tm[1])
        tstep=float(tm[2])
        print ("TIME",tstart,tstop,tstep)
    elif option in ('-e','--extr_time'):
        extrtm=True
        extrtmval=float(value) #time for data to be extracted (from netcdf-file)

    elif option in ('-d','--domain'):
        domain=value.strip().split("/")
        x0=float(domain[0])
        x1=float(domain[1])
        y0=float(domain[2])
        y1=float(domain[3])
        print ("new domain:",x0,x1,y0,y1)
    
    elif option in ('-m','--ncxname'):
        xname=str(value)
    
    elif option in ('-n','--ncyname'):
        yname=str(value)
    elif option in ('-a','--arrivaltime'):
        arrival=True
    elif option in ('-b','--threshold_arrival'):
        threshold_arrival=float(value)
    elif option in ('-F','--sortcol'):
        sortcol=int(value)
        print ("sortcol",sortcol)
        
XA="x*xsc+xadd"
YA="y*ysc+yadd"
ZA="z*zsc+zadd"

if add_first:
    XA="(x+xadd)*xsc"
    YA="(y+yadd)*ysc"
    ZA="(z+zadd)*zsc"
    
#XA="(x+xadd)*xsc"
d=io.io()
if len(propf)==0 and not readnc:
    [x,y,z]=d.read(infilename)

d=io.io()



if readnc:
    print ("Reading the netcdf-file ...")
    read_netcdf(xname,yname,infilename,printvar,pftype,sortcol)
    sys.exit()



#from sct.Grid2D import Grid2D as G

#if the outputfile is on the mtv-format, gauges can be plotted:
if gauges:
    i=0
    for line in open(gauges).readlines():
        line=line.strip().split()
        i+=1
        print (line)
        if len(line)==2:
            name=i
            txt+=mtv_gauges(name,float(line[0]),float(line[1]))
        else:
            name=line[0]
            txt+=mtv_gauges(name,float(line[1]),float(line[2]))
        

if len(propf)>0:
    print ("propagation files for most")
    write_netcdf(propf,tstart,tstop,tstep,pftype,xname,yname)
    
elif skredp:
    for i in range(len(x)):
        print (i)
        g=G(x[i],y[i],z[i])
        if coarse:
            g.extrap(coarse)
        if i==0:
            dp=g
            d.write(depth+"."+typ,
                    g.x*xsc+xadd,
                    g.y*ysc+yadd,
                    g.z*zsc+zadd,
                    txt=txt)
        else:
            if landvalue:
                g.replace_value(dp.z, let=0, value=landvalue)
            x=g.x
            y=g.y
            z=g.z
            d.write(eta+str(i)+"."+typ,eval(XA),eval(YA),eval(ZA),type=typ,txt=txt)
else:
    print ("x[0]",x[0][0])
    g=G(eval(XA),eval(YA),eval(ZA))
    x=g.x
    y=g.y
    z=g.z

    if x0: #extend domain with zeros or make a zoom (-d option)
        x,y,z=g.zoom(x0,y0,x1,y1)
    if coarse:
        g.update(x,y,z)
        g.extrap(coarse)
        x=g.x
        y=g.y
        z=g.z
    print ("sortcol end",sortcol)
    d.write(save,x,y,z,type=typ,txt=txt,degr=degr)
    # d.write(save,x,y,z,type=typ,txt=txt,degr=degr,sortcol=sortcol)

