#!/usr/bin/python
"""
This module reads 1D (two column) and 2D data into numpy arrays, and writes
data from numpy arrays into specified format. Currently supported formats

   1D:  xy 
   2D:  xyz, mtv, gphov, most (ascii grid), grd (surfer), esri (ascii)
   
Example:

>>> from sct3.io import io 
>>> file="Verify/test.xyz"
>>> d=io()
>>> [x,y,z]=d.read(file)
Your input is a XYZ-file
>>> print z[:,9]
[  4.07247000e+03   4.03103000e+03   3.67675000e+03   2.74900000e+03
   3.34292000e+03   4.25419000e+03   2.38240000e+03   1.99097000e+03
   1.04734000e+03  -1.11137000e+02  -1.84614000e+03   3.69400000e+00
   6.23330000e+01   2.60550000e+01  -5.07870000e+02  -1.15211000e+02
  -4.59420000e+01   5.98900000e+01]
>>> file="Verify/test.gphov"
>>> [x,y,z]=d.read(file)
Your input is a GPHOV-file
>>> print z[:,9] 
[  4.07247000e+03   4.03103000e+03   3.67675000e+03   2.74900000e+03
   3.34292000e+03   4.25419000e+03   2.38240000e+03   1.99097000e+03
   1.04734000e+03  -1.11137000e+02  -1.84614000e+03   3.69400000e+00
   6.23330000e+01   2.60550000e+01  -5.07870000e+02  -1.15211000e+02
  -4.59420000e+01   5.98900000e+01]
>>> d.write("tmpgphov",x,y,z,type="gphov",NC=-11)
2D data written into gphov, filename: tmpgphov
>>> import os
>>> os.system("rm -f tmpgphov")
0
>>> file="Verify/test.1D"
>>> [x,y]=d.read(file)
Your input is a 1D-file
>>> print x,y
[ 0.   3.4  5.3  4.3  4.   5. ] [ 1.2  2.   3.2  3.4  3.   4. ]

"""

import os, sys, re, scipy, scipy.io, numpy, time, string, netCDF4

XYZ=0; MTV=1; GPHOV=2; ONEDIM=3; GRD=4; MOST=5; COMCDEP=5; COMCSLIDE=6; NC=-12;
ESRI=7; GCL=8;

class io:
    def __init__(self,quiet=0):
        """init"""
        #gphov regex:
        i=r"\s+(\d+)"                 #integer
        f=r"\s+(-?\d*?\.\d+|-?\d+)"   #float or integer
        F=r"\s*?(-?\d{5}\.\d{6}|-?\d*?\.\d{6}|-?\d*?\.\d+|-?\d+)"
        self.gpexp=r"\s*(-?\d*)"+i+i+F+F+F+F+i
        #integer, float and scientific notation:
        self.real=r"[+/-]?(\d+|\d+(\.\d*)?|\d*\.\d+)([eE][+/-]?\d+)?"
        self.real2=r"("+self.real+r"\s+"+self.real+")"
        self.quiet=quiet

    def _test(self):
        import doctest
        doctest.testmod()

    ########################################################################
    #                                                                      #
    #                               READ                                   #
    #                                                                      #
    ########################################################################

    def pscr(self,txt):
        if not self.quiet:
            print (txt)
 

    def read(self,file,nds=False,unsrt=False):
        """reads  1D data (two column file) or uniform 2D data
        into numpy arrays.
        1D: strictly two columns with numbers, no comments etc.
        2D: available formats:
             xyz   - sorted first in x then in y
             gphov -
             mtv   -
             most gridfiles -
             esrii (ascii)
        
        If the file is a gphov type, you may override the number of data set
        in the heading of the gphov by spec. "nds". Only the nds numbers of
        data sets will be extracted from the data file.

        """
        self.nds=nds
        #do the file exist?
        if not os.path.isfile(file):
            print ("The file '"+str(file)+"', do not exist! Bailing out ...")
            sys.exit(0)
        #read file into one line:
        self.first=open(file,'r').readline()
        self.start=open(file,'r').readlines()[:20]
        self.infile=file

        if str(file)[-3:]=="grd":
            type=GRD
        else:
            type=self.check_type()
        #number of dataset in a file (gphov):
        self.no=0
        #ts=time.time()
        if type==ONEDIM:
            self.pscr(file+" is a 1D-file")
            #her kan man lage en bedre rutine som ikke tar med linjer
            #med bare ett tall, eller linjer med tegn etc. som ikke
            #er tall
            #ts=time.time()
            self.inarray=numpy.loadtxt(self.infile)
            #self.input=open(self.infile,'r').read()
            self._parse1D()
            #ts=time.time()-ts
            #self.pscr "File parsed in %4.2f seconds." % ts

            
            return [self.x,self.y]
        else:
            if type==GPHOV:
                print (file+" is a GPHOV-file")
                self.input=open(self.infile,'r').read()
                self._parse_gphov()
            elif type==MTV:
                print (file+" is a MTV-file")
                self.input=open(self.infile,'r').read()
                self._parse_mtv()
            elif type==XYZ:
                print (file+" is a XYZ-file")
                self.inarray=numpy.loadtxt(self.infile)
                self._parse_xyz(unsrt=unsrt)
            elif type==GRD:
                print (file+" is a grd-file (surfer)")
                self.input=open(self.infile,'r').read()
                self._parse_grd()
            elif type==MOST:
                print (file+" is a most grid-file (ascii)")
                self.input=open(self.infile,'r').read()
                self._parse_most()
            elif type==ESRI:
                print (file+" is a ESRI ASCII file")
                self._parse_esri(file)
            elif type==GCL:
                print (file+" is a GeoClaw type ESRI ASCII file")
                self._parse_gcl_esri(file)
            else:
                print ("The format of your data-file is not supported...")
                print ("...bailing out!")
                sys.exit()
            #ts=time.time()-ts
            #print "File read in %4.2f seconds." % ts
            #print "read: ",self.x0,self.y0,self.dx,self.dy,self.nx,self.ny
            return [self.x, self.y, self.z]

    def _list2array(self,list):
        arr=numpy.zeros([len(list)],float)
        #print "len(string)", len(string)
        count=0
        for item in list:

            try:
                arr[count]=float(item)
                #print "item",float(item)
                count+=1
            except:
                #remove an entry
                arr=arr[:-1]
        #check number of entries:
        if self.no>0 and len(arr)/self.no!=self.nx*self.ny:
            print ( len(arr) )
            print ( self.nx*self.ny )
            print ( "in sct3.io._list2array: error! Incorrect length of array,"+str(len(arr)),"vs.",str(self.nx*self.ny)+"." )
        return arr
        
    def _string2array(self, string):
        string=string.split()
        arr=numpy.zeros([len(string)],float)
        remove=arr.copy()
        #if a joined gphov, mark entries to remove (headings ...)
        for i in range(self.no):
            print ( (i+1)*self.nx*self.ny+i*8 )
            print ( (i+1)*self.nx*self.ny+(i+1)*8 )
            remove[(i+1)*self.nx*self.ny+i*8:(i+1)*self.nx*self.ny+(i+1)*8]=1
        #print "len(string)", len(string)
        count=0
        cc=0
        for item in string:
            if remove[cc]:
                arr=arr[:-1]
            else:
                arr[count]=float(item)
                count+=1
            
##             try:
##                 arr[count]=float(item)
##                 #print "item",float(item)
##                 count+=1
##             except:
##                 #remove an entry
##                 arr=arr[:-1]
            cc+=1

        #check number of entries:
        if self.no>0 and len(arr)/self.no!=self.nx*self.ny:
            print ( len(arr),self.no )
            print ( self.nx*self.ny )
            print ( "in sct3.io._string2array: error! Incorrect length of array,"+str(len(arr)),"vs.",str(self.nx*self.ny)+"." )
        return arr


    def _data_to_arrays(self):
        #make arrays for x,y, and z.
        if self.no > 1:
            self.x=numpy.zeros([self.no,self.nx,self.ny],float)
            self.y=numpy.zeros([self.no,self.nx,self.ny],float)
            self.z=numpy.zeros([self.no,self.nx,self.ny],float)
    
            for n in range(self.no):
                y=self.y0
                for j in range(self.ny):
                    x=self.x0
                    for i in range(self.nx):
                        self.y[n,i,j]=y
                        self.x[n,i,j]=x
                        self.z[n,i,j]=self.data[n*self.nx*self.ny+j*self.nx+i]
                        x+=self.dx
                    y+=self.dy
        else:
            self.x=numpy.zeros([self.nx,self.ny],float)
            self.y=numpy.zeros([self.nx,self.ny],float)
            self.z=numpy.zeros([self.nx,self.ny],float)
            y=self.y0
            for j in range(self.ny):
                x=self.x0
                for i in range(self.nx):
                    self.y[i,j]=y
                    self.x[i,j]=x
                    self.z[i,j]=self.data[j*self.nx+i]
                    x+=self.dx
                y+=self.dy
            
    def _parse_gphov(self):

        flb=self.input.find("\n")     #first linebreak
        firstline=self.input[:flb]    #first line
        data=self.input[flb:]         #content of file without first line
        param=re.search(self.gpexp,firstline)
        if param:
            self.NC=int(param.group(1))
            self.nx=int(param.group(2))
            self.ny=int(param.group(3))
            self.dx=float(param.group(4))
            self.dy=float(param.group(5))
            self.x0=float(param.group(6))
            self.y0=float(param.group(7))
            if self.nds:
                self.no=self.nds
            else:
                self.no=int(param.group(8))
        print ( self.NC,self.nx,self.ny,self.dx,self.dy,self.x0,self.y0,self.no )
        #replace '-' with ' -' if there is no space between numbers :
        self.data=self._string2array(data.replace("-"," -"))
        self._data_to_arrays()

    def _parse_mtv(self):
        #first data is a number on a new line:
        #print "her",self.input
        s=re.search(r"\n\s*"+self.real2,self.input)
        #print self.input
        #print "s.group(1)",s.group(1)
        endmtvinstr=self.input.find(s.group(1))
        #print "endmtv...",endmtvinstr
        mtvinstr=self.input[:endmtvinstr]     #first lines: mtv instructions
        #print mtvinstr
        data=self.input[endmtvinstr-1:]         #followed by the data   
        self.nx=int(re.search(r"nx\s*=\s*(\d+)",mtvinstr).group(1))
        self.ny=int(re.search(r"ny\s*=\s*(\d+)",mtvinstr).group(1))
        self.x0=float(re.search(r"xmin\s*=\s*(-?\d*\.\d*|-?\d+)"\
                                ,mtvinstr).group(1))
        self.y0=float(re.search(r"ymin\s*=\s*(-?\d*\.\d*|-?\d+)"\
                                ,mtvinstr).group(1))
        xm=re.search(r"xmax\s*=\s*(-?\d*\.\d*|-?\d+)",mtvinstr).group(1)
        ym=re.search(r"ymax\s*=\s*(-?\d*\.\d*|-?\d+)",mtvinstr).group(1)
        self.dx=(float(xm)-float(self.x0))/(int(self.nx-1))
        self.dy=(float(ym)-float(self.y0))/(int(self.ny-1))
        #print "parse mtv", self.x0,self.y0,self.dx,self.dy,self.nx,self.ny
        #make array of data:
        self.data=self._string2array(data)
        self._data_to_arrays()

    def _parse_grd(self):
        input=self.input.strip().split()
        self.nx=int(input[1])
        self.ny=int(input[2])
        self.x0=float(input[3])
        self.y0=float(input[5])
        xm=input[4]
        ym=input[6]
        self.dx=(float(xm)-float(self.x0))/(int(self.nx-1))

        self.dy=(float(ym)-float(self.y0))/(int(self.ny-1))
        self.data=self._list2array(input[7:])
        self._data_to_arrays()
        
    def _parse1D(self):
        self.x=self.inarray[:,0]
        self.y=self.inarray[:,1]
        #s=re.findall(r""+self.number+r"\s+"+self.number,self.input)
        #self.x=numpy.zeros([len(s)])
        #self.y=numpy.zeros([len(s)])
        #for item in s:
        #    print item
        
    def _parse_most(self):
        input=self.input.strip().split()
        self.nx=int(input[0])
        self.ny=int(input[1])
        self.x=numpy.zeros([self.nx,self.ny],float)
        self.y=numpy.zeros([self.nx,self.ny],float)
        self.z=numpy.zeros([self.nx,self.ny],float)

        for i in range(self.nx):
            start=2
            self.x[i,:]+=float(input[i+2])
            start+=self.nx
        for j in range(self.ny):
            r=float(input[j+start])
            self.y[:,self.ny-1-j]+=float(r)

        count=start+self.ny

        for j in range(self.ny):
            for i in range(self.nx):
                self.z[i,self.ny-1-j]=float(input[count])
                count+=1

       
    def _parse_esri(self,file):
        self.input=open(self.infile,'r').read()
        #find first datapoint
        s=re.search(r"\n\s*(-?\d*\.\d*|-?\d+)",self.input)
        enddata=self.input.find(s.group(1))
        head=self.input[:enddata]

        #extract data from headingilnes
        nx=int(re.search(r"ncols\s*(\d+)",head).group(1))
        ny=int(re.search(r"nrows\s*(\d+)",head).group(1))
        xll=float(re.search(r"xllcorner\s*(-?\d*\.\d*|-?\d+)"\
                            ,head).group(1))
        yll=float(re.search(r"yllcorner\s*(-?\d*\.\d*|-?\d+)"\
                            ,head).group(1))
        dx=float(re.search(r"cellsize\s*(-?\d*\.\d*|-?\d+)"\
                           ,head).group(1))
        dy=dx

        #split data
        data=self.input[enddata:].strip().split()

        self.nx=nx
        self.ny=ny
        self.x=numpy.zeros([self.nx,self.ny],float)
        self.y=numpy.zeros([self.nx,self.ny],float)
        self.z=numpy.zeros([self.nx,self.ny],float)


        #fill x,y,z arrays
        count=0
        for j in range(self.ny-1,-1,-1):
            for i in range(self.nx):
                self.x[i,j]=xll+dx*i
                self.y[i,j]=yll+dy*j
                self.z[i,j]=data[count]
                count+=1
                
    def _parse_gcl_esri(self,file):
        """Esri-asci ala Geoclaw"""
        self.input=open(self.infile,'r').read()
        #find first datapoint
         
        fid=open(self.infile,'r')
        line=fid.readline().split()
        mx = int(line[0])
        line=fid.readline().split()
        my = int(line[0])
        line=fid.readline().split()
        xll = float(line[0])
        line=fid.readline().split()
        yll = float(line[0])
        line=fid.readline().split()
        dx = float(line[0])
        line=fid.readline().split()
        dy=dx
        nodata = float(line[0])
        #nodata=-9999

        print ( "geoclaw",mx,my,xll,yll,dx,nodata )
        

        ## #split data
        data=self.input.strip().split()[6:]
        #print "data",data

        self.nx=mx
        self.ny=my
        self.x=numpy.zeros([self.nx,self.ny],float)
        self.y=numpy.zeros([self.nx,self.ny],float)
        self.z=numpy.zeros([self.nx,self.ny],float)


        ## #fill x,y,z arrays
        count=0
        for j in range(self.ny-1,-1,-1):
            for i in range(self.nx):
                self.x[i,j]=xll+dx*i
                self.y[i,j]=yll+dy*j
                self.z[i,j]=data[count]
                count+=1
         
        
        
    def _parse_xyz(self,unsrt=False):
        #find structure, suppose sorted in x then y:
        #sort if not sorted!!

        
        #find ny (i.e. number of equal x entry in one block)
        self.x0=self.inarray[0,0]
        for i in range(len(self.inarray)):
            if self.inarray[i,0]!=self.x0:
                break
        self.ny=i
        #find nx, limits of domain, dx,dy:
        self.nx=len(self.inarray)/self.ny
        self.y0=self.inarray[0,1]
        xm=self.inarray[len(self.inarray)-1,0]
        ym=self.inarray[len(self.inarray)-1,1]
        self.dx=(xm-self.x0)/(self.nx-1)
        self.dy=(ym-self.y0)/(self.ny-1)
        #construct arrays:
        self.x=numpy.zeros([self.nx,self.ny],float)
        self.y=numpy.zeros([self.nx,self.ny],float)
        self.z=numpy.zeros([self.nx,self.ny],float)
        #exit if incorrecte nx,ny,len(data)
        if self.nx*self.ny!=len(self.inarray) and not unsrt:
            print ( "Your xyz-file is not uniform or sorted correctly!" )
            print ( "... bailing out!" )
            sys.exit(0)
        #read data
        for i in range(self.nx):
            for j in range(self.ny):
                self.x[i,j]=float(self.inarray[i*self.ny+j,0])
                self.y[i,j]=float(self.inarray[i*self.ny+j,1])
                self.z[i,j]=float(self.inarray[i*self.ny+j,2])

                
    def check_type(self):
        #start is the first 20 lines of file:
        xyzcount=0
        mtvcount=0
        mostcount=0
        gclcount=0
        ONEDIMcount=0
        TYPE=None
        for line in self.start:
            if len(line.strip())>0:
                if line.strip()[0]==r"%":          #test mtv
                    mtvcount+=1
                if len(line.split())==3:           #test xyz
                    xyzcount+=1
                if len(line.split())==2:
                    try:
                        float(line.split()[0])
                        float(line.split()[1])
                        ONEDIMcount+=1
                    except:
                        pass
                if len(line.split())==1:
                    mostcount+=1
                    gclcount+=1
                gphov=re.search(self.gpexp,line)   #test gphov
                try:
                    nx=int(re.search(r"ncols\s*(\d+)",line).group(1))
                    print ( "ESRI",nx )
                    TYPE=ESRI
                    break
                except:
                    pass
                if gphov:
                    TYPE=GPHOV
                    break
                elif xyzcount>4:
                    TYPE=XYZ
                    break
                elif mtvcount>2:
                    TYPE=MTV
                    break
                elif ONEDIMcount>=2:
                    TYPE=ONEDIM
                    break
                elif mostcount>15:
                    TYPE=MOST
                    break
                elif gclcount==6:
                    TYPE=GCL
        print ( "check type", TYPE )
        return TYPE


    #######################################################################
    #                                                                     #
    #                             WRITE                                   #
    #                                                                     #
    #######################################################################


                   
    def write(self,filename,x,y,z=None,type=None, NC=NC,txt="",degr=True,nds=False):
        """writing data in arrays to file (filename).
        
        1D: x,y is specified (z=None) and is 1D arrays,
            written as two-column file
            
        2D: x,y,z is specified and is 2D arrays:

            type="mtv"|"gphov"|"xyz"|"most"|"esri"|"gcl". If no type is given, the type 
            is taken as the extension to the filename.
            For format gphov NC is an index for lengthscales
            (see below). Default: NC=-12:

            NC: depth grid
            -11   m     m
            -12   m     km
            -21   km    m
            -22   km    km

        Plotting commands for mtv may be given as the string 'txt'.
        
        """
        self.nds=nds
        #1D:
        
        print ("io.py : z = ", z)
        if z.all()==None:
            self._write1D(filename,x,y)
            self.pscr("1D data written to a two-column file: "+filename)
        #2D
        else:
            self.x=x
            self.y=y
            self.z=z

            #find parameters:
            self.x0=self.x[0][0]
            self.y0=self.y[0][0]
            self.dx=self.x[1][0]-self.x[0][0]
            self.dy=self.y[0][1]-self.y[0][0]
            self.nx=len(self.z[:,0])
            self.ny=len(self.z[0,:])
            
            self.fn=filename
            #NC is first entry in gphov:
            self.NC=NC  
            #if type is not given - use the extention of the filename:
            if type==None:
                type=filename.split('.')[-1]
                print ("TYPE",type)
            if type=="mtv" or type=="MTV":
                self._write_mtv(txt)
                self.pscr("2D data written into mtv, filename: "+filename)
            elif type=="gphov":
                self._write_gphov()
                self.pscr("2D data written into gphov, filename: "+filename)
            elif type=="xyz" or type=="XYZ":
                self._write_xyz()
                self.pscr("2D data written into xyz, filename: "+filename)
            elif type=="esri" or type=="ESRI":
                self._write_esri()
                self.pscr("2D data written into ESRI ASCII format, file:"+filename)
            elif type=="most" or type=="MOST":
                self._write_most()
                self.pscr("2D data written into grid-most format, file:"+filename)
            elif type=="dep" or type=="COMCDEP":
                self._write_comcot_dep()
                self.pscr("2D data written into COMCOT depth format, file:"+filename)
            elif type=="csl" or type=="COMCSLIDE":
              self._write_comcot_slide()
              self.pscr("2D data written into COMCOT slide format, file:"+filename)
            elif type=="gcl" or type=="GCL":
              self._write_gcl_esri()
              self.pscr("2D data written into Geoclaw/ESRI format, file:"+filename)
            elif type=="nc":
                x="lon"
                y="lat"
                if not degr:
                    x="x"
                    y="y"
                self._write_netcdf(x,y)
                self.pscr("2D data written into netCDF format, file:"+filename)
                
            else:
                self.pscr("The type of output (\""+type+"\") is not valid")
                sys.exit()

    def _write1D(self,filename,x,y):
        self.x=x
        self.y=y
        file=open(filename,'w')
        for i in range(len(self.x)):
            file.write("%15.8f %15.8f\n"% (self.x[i],self.y[i]))
        file.close()
        
    def _write_mtv(self,txt):
        file=open(self.fn,'w')
        file.write(self._mtv_cmds()+"\n\n")
        file.write("%"+txt+"\n")
        count=0
        for j in range(self.ny):
            for i in range(self.nx):
                count+=1
                file.write("%g " % self.z[i,j])
                if count % 50 ==0  and i!=0:
                    file.write("\n")
        file.write("\n")
        file.close()

    def _write_netcdf(self,x="lon",y="lat"):
       out=netCDF4.Dataset(self.fn, 'w', format='NETCDF3_CLASSIC')
       out.createDimension(x, len(self.z[:,0]))
       out.createDimension(y, len(self.z[0,:]))
       lonv = out.createVariable(x,'d',(x,))
       lonv.unit='east'
       #lonv.point_spacing="even"
       latv = out.createVariable(y,'d',(y,))
       latv.unit='north'
       #latv.point_spacing="ueven"
       field = out.createVariable('Field','d',(y,x,))
       lonv[:]=self.x[:,0]
       latv[:]=self.y[0,:]
       field[:]=self.z.transpose()
       out.sync()
       out.close()



        
    def _write_gphov(self):
        file=open(self.fn,'w')
        file.write(self._gphov_index()+"\n")
        count=0
        for j in range(self.ny):
            for i in range(self.nx):
                file.write("%12.6f"% self.z[i,j])
                #file.write("%9.3f"% self.z[i,j])
                count+=1
                if (count==10 or i==self.nx-1):
                    file.write("\n")
                    count=0
        file.close()

    def _write_comcot_dep(self):
        file=open(self.fn,'w')
        print ( "Comcot - depth file", self.nx, self.ny )
        count=0
        for j in range(self.ny):
            for i in range(self.nx):
                file.write("%9.3f"% self.z[i,j])
                count+=1
                if (count==10 or i==self.nx-1):
                    file.write("\n")
                    count=0
        file.close()

    def _write_comcot_slide(self):
        file=open(self.fn,'w')
        print ( "Comcot - slide file", self.nx, self.ny )
        for j in range(self.ny):
            for i in range(self.nx):
                file.write("%10.5f"% self.z[i,j])
                if (i==self.nx-1):
                    file.write("\n")
        file.close()

    def _write_most(self):
        file=open(self.fn,'w')
        file.write(" %12d %12d\n" %(self.nx,self.ny))
        for i in range(self.nx):
            #file.write(" %11.4f\n" %(self.x[i,0]))
            file.write(" %11.8f\n" %(self.x[i,0]))
        for j in range(self.ny-1,-1,-1):
#       for j in range(self.ny):
            #file.write(" %11.4f\n" %(self.y[0,j]))
            file.write(" %11.8f\n" %(self.y[0,j]))
        for j in range(self.ny-1,-1,-1):
            for i in range(self.nx):
                file.write(" %9.3f" %(self.z[i,j]))
            file.write("\n")
        file.close()

    def _write_esri(self):
        file=open(self.fn,'w')
        #write heading
        nx=self.nx
        ny=self.ny
        x0=self.x0
        y0=self.y0
        dx=self.dx
        file.write("ncols %(nx)d\n" %vars())
        file.write("nrows %(ny)d\n" %vars())
        file.write("xllcorner %(x0).10f\n" %vars())
        file.write("yllcorner %(y0).10f\n" %vars())
        file.write("cellsize %(dx).10f\n" %vars())
        file.write("NODATA_value -9999\n")
        

        for j in range(self.ny-1,-1,-1):
            for i in range(self.nx):
                file.write(" %9.3f" %(self.z[i,j]))
            file.write("\n")
        file.close()

    def _write_gcl_esri(self):
        file=open(self.fn,'w')
        #write heading
        nx=self.nx
        ny=self.ny
        x0=self.x0
        y0=self.y0
        dx=self.dx
        nodata=-9999
        file.write("%(nx)10d\n" %vars())
        file.write("%(ny)10d\n" %vars())
        file.write("%(x0)20.10f\n" %vars())
        file.write("%(y0)20.10f\n" %vars())
        file.write("%(dx)20.10f\n" %vars())
        file.write("%(nodata)20.10f\n" %vars())
        

        for j in range(self.ny-1,-1,-1):
            for i in range(self.nx):
                file.write(" %20.10f" %(self.z[i,j]))
            file.write("\n")
        file.close()

    def _write_xyz(self):
        file=open(self.fn,'w')
        for i in range(self.nx):
            for j in range(self.ny):
                file.write("%f %f %f\n"% (self.x[i,j],self.y[i,j],float(self.z[i,j])))
        file.close()
        

    def _mtv_cmds(self):
        """constructing the mtv instructions (first lines in a mtv-file)"""
        cmds="$DATA=CONTOUR \n% contstyle=2 \n% nsteps=30\n"
        cmds+=r"%"+' toplabel="" subtitle="" xlabel="" ylabel=""\n'
        cmds+=r"%"+' nx=%d xmin=%s xmax=%s\n' \
               %(self.nx,self.x0,self.x0+(self.nx-1)*self.dx)
        cmds+=r"%"+' ny=%d ymin=%s ymax=%s\n' \
               %(self.ny,self.y0,self.y0+(self.ny-1)*self.dy)
        return cmds

    def _gphov_index(self):
        """initializing grid specifications (first line) for gphov"""
        no=1
        if self.nds:
            no=self.nds
        index="%d %d %d %12.6f%12.6f%12.6f%12.6f %d"\
               % (self.NC,self.nx,self.ny,self.dx,self.dy,self.x0,self.y0,no)
        return index


    def writeDP(self,grid,filename):
        """writing the one-dimensional array grid, into file 'filename' in
        Diffpacks grid format"""
        nnodes=len(grid)
        file=open(filename,'w')
        dp="""
Finite element mesh (GridFE):

  Number of space dim. =   1  embedded in physical space with dimension 1
  Number of elements   =  %d
  Number of nodes      =  %d

  All elements are of the same type : true
  Max number of nodes in an element: 2
  Only one material                : true
  Lattice data                     ? 0



  2 boundary indicators: 
    x1=%f
    x2=%f


  Nodal coordinates and nodal boundary indicators,
  the columns contain:
   - node number
   - coordinates
   - no of boundary indicators that are set (ON)
   - the boundary indicators that are set (ON) if any.
#
""" %(len(grid)-1,len(grid),grid[0],grid[-1])

        for i in range(1,len(grid)+1,1):
            dp+="   "+str(i)+"  ( %e)" % grid[i-1]
            if i==1:
                dp+="  [1] 1\n"
            elif i==nnodes:
                dp+="  [1] 1\n"
            else:
                dp+="  [0]\n"
                
        dp+="""    
  Element types and connectivity
  the columns contain:
   - element number
   - element type
   - material number
   - the global node numbers of the nodes in the element.
#
"""
        for i in range(1,len(grid)+1,1):
            dp+=" "+str(i)+"  \tElmB2n1D \t1 \t"+str(i)+"\t"+str(i+1)+"\n"    
        file.write(dp)
        print ( "\nGrid is converted into Diffpack-format found in file \""
              +filename+"\"!" )
        print ( "("+str(nnodes)+" nodes totally)\n" )
        
        file.close()

    #######################################################################
    #                                                                     #
    #                             LOGFILE                                 #
    #                                                                     #
    #######################################################################



    def log(self):
        """creates a logfile wherever this function is called, appends
        to the logfile if already present. """
        #create at log-file if not present filename.log add to svn
        #write to logfile (filename, path, time, commands, revision number etc.)
        #check logfile into svn
        dir=os.getcwd()
        tm=time.asctime()
        script=sys.argv[0]
        scriptrev=script
        log=script+".log"
        opt=""
        prev=""
        for o in range(1,len(sys.argv)):
            opt+=o+" "
        if not os.path.isfile(log):
            file=open(log,'w')
            os.system("svn add "+log)
            os.system("svn ci -m update "+log)

        else:
            prev=open(log,'r').read()
            file=open(log,'w')

 
        #revision number:
        #if not in svn, add and check in

        os.system("svn add "+script)
        os.system("svn ci -m update "+script)
        rev=" - "
        try:
            pop=os.popen("svn log "+script,'r').read()
            pat=r"r(\d*)"
            rev=re.findall(pat,pop)[0]  #first occurence (last revision)
            scriptrev+="_rev"+rev
            print ( rev )
        except:
            pass
        import socket
        host = socket.gethostname()
        print ( host )
        file.write("="*80)
        file.write("\nTim: "+tm+"\n")
        file.write("Scr: "+script+"\n")
        file.write("Hos: "+host+"\n")
        file.write("Rev: "+rev+"\n")
        file.write("Opt: "+opt+"\n")
        file.write("Dir: "+dir+"\n")
        file.write("\n\n")
        file.write(prev)
        
        
        file.close()

        os.system("svn ci -m update "+log)
        os.system("cp "+script+" "+scriptrev) #make a copy of current rev
        
        return scriptrev
    

if __name__ == '__main__':

    #1D:
    file="etamax.dat"
    file="test.gphov"
    d=io()
    [x,y,z]=d.read(file)
    print ( x )
    #x=arr[0]
    #y=arr[1]
    #d.write("tmp1D",x,y)
    #2D:
    #file="Verify/test.mtv"
    #[x,y,z]=d.read(file)
    #print "d.z",d.z
    #d.write("tmp2D.xyz",x,y,z)
    #d.write("tmp2D",x,y,z,type="gphov")
    
    #running doc-test:
    print ( "running test:" )
    d._test()

    
