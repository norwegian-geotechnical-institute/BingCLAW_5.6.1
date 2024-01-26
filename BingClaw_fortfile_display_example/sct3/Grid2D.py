#!/usr/bin/python 
"""
Module for manipulating a 2D grid. The module reads arrays for x,y, and z,
and manipulates the grid according to the parameters and functions called.
"""

import os, sys, math,numpy
import scipy.interpolate as intp
import scipy
from sct3.io import io

# from sct3.Grid1D import Grid1D as G



class Grid2D:
    def __init__(self,x,y,z):
        #read the arrays:
        self.x=x
        self.y=y
        self.z=z
        self._update_grid_parameters()

        
    def _update_grid_parameters(self):
        #find parameter x0,y0,xm,ym,nx,ny,dx,dy:
        self.x0=self.x[0][0]
        self.y0=self.y[0][0]
        self.xm=self.x[-1][0]
        self.ym=self.y[0][-1]
        self.dx=self.x[1][0]-self.x[0][0]
        self.dy=self.y[0][1]-self.y[0][0]
        self.nx=len(self.z[:,0])
        self.ny=len(self.z[0,:])


    def update(self,x,y,z):
        """exchange the arrays of the instance with new one"""
        self.x=x
        self.y=y
        self.z=z

    
    def extrap(self,factor=2):
        """making the grid coarser with an integer factor. factor=2 means
        that every second gridpoint in x and y direction is used as the new
        grid. factor=3 -> every third point is used etc."""


        self.x=self.x[::factor,::factor]
        self.y=self.y[::factor,::factor]
        self.z=self.z[::factor,::factor]
        self._update_grid_parameters()
        
        return self.x,self.y,self.z

    
    #def interp(self,xnew,ynew,x=self.x,y=self.y,z=self.z):
    def interp_bspl(self,xnew,ynew,x=None,y=None,z=None):
        """defining a new grid based on the arrays xnew and ynew. Field
        values are interpolated using B-spline. The function is too slow
        for large grid, but is ok for grids less than 50 x 50 points. This
        function may be used as standalone, you must then specify the
        grid to interpolate the new grid from. If x=y=z=None, we use
        x=self.x, y=self.y, and z=self.z.
        """
    
        try:
            if x==None or y==None or z==None:
                x=self.x
                y=self.y
                z=self.z
            
        except:
            pass
#        print "len(x)",len(x)
            
        self._update_grid_parameters()
        tck=intp.bisplrep(x,y,z,s=0)
        znew=intp.bisplev(xnew[:,0],ynew[0,:],tck)
        self.update(xnew,ynew,znew)
        self._update_grid_parameters()
        
        return self.x,self.y,self.z

    
    def interp_lin(self,xnew,ynew,outsidevalue=0):
        """interpolate values of self.z onto grid xnew,ynew.
        Field values are interpolated using linear interpolation.
        For points outside the old grid are given the
        'outsidevalue' """

        #check if new grid is outside or idenstical to the old:

        #return only the overlap between the grids!!!!!!!!!!!!
        ######################################################
        znew=xnew.copy()
        znew.fill(0.0)

        for i in range(len(xnew)):
            for j in range(len(ynew[0,:])):
                x=xnew[i,j]
                y=ynew[i,j]
                ni=0
                if self._inside(x,y):
                    while self.x[ni,0]<x:
                        ni+=1
                    iv=ni-1
                    znew[i,j]=self._interp_col(iv,x,y)
                else:
                    znew[i,j]=outsidevalue
                



        self.update(xnew,ynew,znew)
        self._update_grid_parameters()
        
        return self.x,self.y,self.z
    
    def _inside(self,x,y):

        if self.x[0,0]<x<self.x[len(self.x[:,0]),0] \
               and self.y[0,0]<y<self.x[0,len(self.y[0,:])]:
            return 1
        else:
            return 0
            

    def _interp_col(self,i,x,y):
        #x,y is a point inside the y-columns i and i+1:
        
        X   = [self.x[i,0],self.x[i+1,0]]
        Y   = self.y[0,:]
        #interpolate y-direction, x-index i:
        f   = intp.interp1d(Y,self.z[i,:])
        vl  = f(y)
        #interpolate y-direction, x-index i+1:
        f   = intp.interp1d(Y,self.z[i+1,:])
        vr  = f(y)
        #interpolate x-direction:
        f   = intp.interp1d(X,[vl[0],vr[0]])

        return f(x)

            


    def _interp_cell(self,i,j,x,y):
        #x,y is a point inside the cell i,i+1,j,j+1
        X   = [self.x[i,0],self.x[i+1,0]]
        Y   = [self.y[0,j],self.y[0,j+1]]
        #interpolate x-direction, bottom:
        f   = intp.interp1d(X,[self.z[i,j],self.z[i+1,j]])
        vj  = f(x)
        #interpolate x-direction, top:
        f   = intp.interp1d(X,[self.z[i,j+1],self.z[i+1,j+1]])
        vjp = f(x)
        #interpolate y-direction:
        f   = intp.interp1d(Y,[vj[0],vjp[0]])

        return f(y)
            

    def smooth(self, steps=1):
        """smoothing data on a 2D uniform grid using simple 5-pt scheme.
        Updates the data values in array self.z"""
        self._update_grid_parameters()
        d=self.z.copy()
        s=d.copy()
        nx=self.nx
        ny=self.ny
        print ("Starting smoothing",steps,"steps!")
        for no in range(steps):
            d=s
            #BOUNDARIES:
            #i=0
            s[0,1:ny-2]=d[0,1:ny-2]*0.5+(d[0,2:ny-1]+d[0,0:ny-3]\
                                         +d[1,1:ny-2])*(1/3.0)*0.5
            #i=nx-1
            s[nx-1,1:ny-2]=d[nx-1,1:ny-2]*0.5+(d[nx-1,2:ny-1]+d[nx-2,1:ny-2]\
                                               +d[nx-1,0:ny-3])*(1/3.0)*0.5
            #j=0
            s[1:nx-2,0]=d[1:nx-2,0]*0.5+(d[2:nx-1,0]+d[0:nx-3,0]\
                                         +d[1:nx-2,1])*(1/3.0)*0.5
            #j=ny-1
            s[1:nx-2,ny-1]=d[1:nx-2,ny-1]*0.5+(d[2:nx-1,ny-1]+d[1:nx-2,ny-2]\
                                               +d[0:nx-3,ny-1])*(1/3.0)*0.5
            #INTERIOR
            s[1:nx-2,1:ny-2]=d[1:nx-2,1:ny-2]*0.5+\
                              (d[0:nx-3,1:ny-2]+d[2:nx-1,1:ny-2]\
                               +d[1:nx-2,0:ny-3]+d[1:nx-2,2:ny-1])*0.125
        self.z=s
        return self.x,self.y,self.z




    def scale(self,xsc=1,ysc=1,zsc=1):
        """scaleing x,y,z by multiplying the arrays"""
        self.x*=xsc
        self.y*=ysc
        self.z*=zsc
        self._update_grid_parameters()
        
        return self.x,self.y,self.z
        

    def gradient(self,dir=None, abs=True, degrees=False):
        """find the gradient of the field defined by the arrays x,y,and z.
        Default is the absolute gradient while:
        dir=x: abs. value of gradient in x direction (or real values abs=False)
        dir=y: abs. value of gradient in y direction (or real values abs=False)
        The function returns the grid as three arrays (same size):

        x,y,grad(z).

        If degrees=True the z-values returned are given in degrees."""
        

        ydir=1
        xdir=1
        if dir=="x":
            ydir=0
        elif dir=="y":
            xdir=0
        elif dir==None:
            pass
        else:
            print ("Wrong direction:",dir,"Inside Grid2d.gradient()")
            sys.exit(0)
        grad=self.z.copy()
        grad.fill(0.0)
        for i in range(1,self.nx-1):
            for j in range(1,self.ny-1):
                gry=(self.z[i,j+1]-self.z[i,j-1])/(2*self.dy)
                grx=(self.z[i+1,j]-self.z[i-1,j])/(2*self.dx)
                if not abs:
                    grad[i,j]=xdir*grx+ydir*gry
                else:
                    grad[i,j]=math.sqrt(xdir*grx*grx+ydir*gry*gry)
                if degrees:
                    grad[i,j]=math.atan(grad[i,j])*180/math.pi
        self.update(self.x,self.y,grad)
        self._update_grid_parameters()
        return self.x,self.y,self.z

    def max(self):
        """maximum z-value and the position in the grid is returned"""
        print ("not implemented")
        
    def min(self):
        """minimum z-value and the position in the grid is returned"""
        print ("not implemented")
        

    def steepest_descends(self,xs=0,ys=0,dr=1,totlen=100):
        """Returns a pair of a x array and y array of the path defined by
        the steepest descends with start at (xs,ys)"""
        print ("steepest descends, path length", totlen)
        
        Gx=Grid2D(self.x,self.y,self.z)
        Gy=Grid2D(self.x,self.y,self.z)
        tmpx,tmpy,tmpz=Gx.gradient(dir="x",abs=False)
        tmpx,tmpy,tmpz=Gy.gradient(dir="y",abs=False)
        dx=tmpx[1,0]-tmpx[0,0]
        dy=tmpy[0,1]-tmpy[0,0]

        length=0
        x=xs
        y=ys
        xl=[xs]
        yl=[ys]
        while length<totlen:
            #interpolate values for gradients
            for i in range(len(self.x[:,0])):
                if self.x[0,0]+i*dx>=x:
                    i-=1
                    break

            for j in range(len(self.y[0,:])):
                if self.y[0,0]+j*dy>=y:
                    j-=1
                    break
            #find direction
            grx=-Gx._interp_cell(i,j,x,y)[0]
            gry=-Gy._interp_cell(i,j,x,y)[0]
            add=0

            if grx>0 and gry<0:
                add=math.pi
                
            alpha=math.atan(gry/grx)+add
            
            #find new point
            x+=math.cos(alpha)*dr
            y+=math.sin(alpha)*dr
            
            xl.append(x)
            yl.append(y)

            length+=dr
        xval=numpy.array(xl,float)
        yval=numpy.array(yl,float)

        return xval,yval

            
    def crossection(self,xs=0,ys=0,xe=10,ye=10,dr=1,start=0,arg="r",\
                    pts=None, spline=False, steepest_desc=False,\
                    pathlength=100):
        """interpolate z-values along a given path defined by the start
        point (xs,ys) and end point (xe,ye) or the points defined in x,y pairs
        in the file 'pts'. The path defined by these points may be smoothed
        using spline.
        
        The resolution of the crossection is given by 'dr'. 

        \'arg\' is the spatial argument given in the output.  Possible choices
        of \'arg\' is

        \'x\' (x-axis),
        '\y'\ (y-axis),
        or
        \'r\' (distance along crossection - default)

        If arg='r', then 'start' is the
        argument of the first point of the crossection, the second point
        have argument start+dr. The function returns one 1D array for
        x values and one for the interpolated field values.

        If 'follow_slope' is True, the interpolation of the path-values starts
        at point (xs,ys) (the values for xe/ye make no sense) and for each
        segment dr the direction of the path is chosen to follow the direction
        of the slope until the length of the path is reaching the 'pathlength'.
        
        """
        
        if steepest_desc:
            print ("crossection - steepest descends")
            G=Grid2D(self.x,self.y,self.z)
            [x,y]=G.steepest_descends(xs,ys,dr,pathlength)
            length=pathlength
            xy=numpy.zeros([len(x),2],float)
            xy[:,0]=x[:]
            xy[:,1]=y[:]
            d=io()
            d.write(".track",xy[:-1,0],xy[:-1,1])
            print ("The position of the track is plotted in file \'.track\'.")
            
            
        #find elevation ratio and dx,dy
        if xe!=xs:
            a=(float(ye)-ys)/(xe-xs)
            dx=dr*((1+a**2)**(-0.5))
            dy=abs(a)*dr*((1+a**2)**(-0.5))
        else:
            dx=0
            dy=dr
        if xs>xe:
            dx=-dx
        if ys>ye:
            dy=-dy
            
        length=0 #length of path
        if pts:  #points on file for defining the path
            print ("points on file for defining the path")
            d=io()
            [x,y]=d.read(pts)
            if spline:
                print ("no spline implemented")
##                 xspl=scipy.mgrid[x[0]:x[-1]:10j]
##                 print "xspl",xspl
##                 print "x",x
##                 g=G(x,y)
##                 [x,y]=g.interp_spl(xspl)
##                 print "x",x
                
            for i in range(1,len(x)):
                length+=math.sqrt((x[i]-x[i-1])**2+(y[i]-y[i-1])**2)
            xy=numpy.zeros([int(1+length/dr),2],float) #x,y points along path.
            xy[0,0]=x[0]
            xy[0,1]=y[0]
            #dR=0
            #rlength=dr #running length along path
            #dRold=dR
            drl=dr
            xsign=1
            ysign=1
            ii=1  #index for xy
            for i in range(0,len(x)-1):
                if (x[i]==x[i+1]):
                    if y[i+1]-y[i]>=0:
                        sign=1
                        ysign=1
                    else:
                        sign=-1
                        ysign=-1
                    a=1000000000*sign
                else:
                    a=(y[i+1]-y[i])/(x[i+1]-x[i]) #slope angle
                    if y[i+1]-y[i]<0:
                        ysign=-1
                    if x[i+1]-x[i]<0:
                        xsign=-1
                    
                seglen=math.sqrt((\
                    x[i+1]-x[i])**2+(y[i+1]-y[i])**2) #length betw.
                                                      #tow neighb. nodes
                
                drest=seglen
                if drest>=dr:
                    while drest>=dr:
                        #x- and y-value of new point on path
                        if drl==dr:
                            xy[ii,0]=xsign*math.sqrt(drl*drl/(a*a+1))+xy[ii-1,0]
                            xy[ii,1]=ysign*abs(a*(xy[ii,0]-xy[ii-1,0]))+xy[ii-1,1]
                        else:
                            xy[ii,0]=xsign*math.sqrt(drl*drl/(a*a+1))+x[i]
                            xy[ii,1]=ysign*abs(a*(xy[ii,0]-x[i]))+y[i]
                        drl=dr
                        drest-=drl
                        ii+=1
                while drest<dr:
                    drest+=dr
                    drl=dr-drest
            d.write(".track",xy[:-1,0],xy[:-1,1])
            print ("The position of the track is plotted in file \'.track\'.")
                    
        elif not pts and not steepest_desc:
            length=math.sqrt((xe-xs)**2+(ye-ys)**2)

        
        
        #initialize an array for holding the values along the crossection:
        if steepest_desc:
            cross=xy.copy()
            #cross.fill(0.0)
                
        else:
            cross=numpy.zeros(\
                [int(1+length/dr),2],float)

        count=0
        if pts or steepest_desc:
            x=xy[count,0]
            y=xy[count,1]
        else:
            x=xs
            y=ys
        r=start
        for c in range(len(cross)):
            if arg:
                if arg=="x":
                    r=x
                elif arg=="y":
                    r=y
                elif arg=="r":
                    pass
                else:
                    print ("Your input to Grid2D.crossection, \'arg\':\
                    %s is not valid" %arg)
            #is point outside, value=0?
            if x<self.x[0,0]:
                cross[c,0]=r
                cross[c,1]=0
            elif x>self.x[self.nx-1,0]:
                cross[c,0]=r
                cross[c,1]=0
            elif y<self.y[0,0]:
                cross[c,0]=r
                cross[c,1]=0
            elif y>self.y[0,self.ny-1]:
                cross[c,0]=r
                cross[c,1]=0
            #inside, interpolate from neighbours
            else:
                i=0
                while self.x[i,0]<=x:
                    i+=1
                iv=i-1
                j=0
                while self.y[0,j]<=y:
                    j+=1
                jv=j-1
             
                cross[c,0]=r
                #interpolate values from the cell of 4 neighb. points:
                cross[c,1]=self._interp_cell(iv,jv,x,y)
            count+=1
            r+=dr
            if pts or steepest_desc:
                if count>=len(xy):
                    x=xy[len(xy)-1,0]
                    y=xy[len(xy)-1,1]
                else:
                    x=xy[count,0]
                    y=xy[count,1]
            else:
                x+=dx
                y+=dy
 
            

        return cross[:,0],cross[:,1]

        

    def zoom(self,x0,y0,xm,ym,value=0):
        """changing the physical domain of a grid by given the coordinates
        of lower left corner (x0,y0) and upper right corner (xm,ym) of the
        new domain. If the new domain is larger than the old, the array
        will be filled with <value> in new points"""

        #check that x0,y0,xm,ym is inside the original grid:
#        print "zoom", x0,self.x[0,0], xm,self.x[-1,0] , y0,self.y[0,0],ym,self.x[0,-1]
#        if x0>=self.x[0,0] and xm<=self.x[-1,0] and y0>=self.y[0,0] and ym<=self.x[0,-1]:
#            print "zoom"
        #find indexes for lower left corner, i0, j0:
        i0=0
        for i in range(self.nx):
            if self.x[i,0]>=x0:
                i0=i
                break
        im=self.nx-1
        for i in range(i0,self.nx):
            if self.x[i,0]>=xm:
                im=i
                break

        #find indexes for upper right corner, im, jm:
        j0=0
        for j in range(self.ny):
            if self.y[0,j]>=y0:
                j0=j
                break
        jm=self.ny-1
        for j in range(j0,self.ny):
            if self.y[0,j]>=ym:
                jm=j
                break
        if i0==im or j0==jm:
            print ("""The definition of the smaller area is outside
        data, lower left=(%f,%f), upper right=(%f,%f).""" %(x0,y0,xm,ym) )
            print ( "Terminating ..." )
            sys.exit()

        I0=i0
        Im=im
        J0=j0
        Jm=jm
        #extending data area with zeros?
        dx=self.x[1,0]-self.x[0,0]
        dy=self.y[0,1]-self.y[0,0]
        #adjust x0,xm,y0,ym to match a integer numbers of gridcells:
        nx1=int((self.x[0,0]-x0)/dx)
        nx2=int((self.x[-1,0]-xm)/dx)
        ny1=int((self.y[0,0]-y0)/dy)
        ny2=int((self.y[0,-1]-ym)/dy)
        nx=len(self.x[:,0])
        ny=len(self.y[0,:])
#        print "nx1,nx2....",nx1,nx2,ny1,ny2,nx,ny
        NX=nx1-nx2+nx
        NY=ny1-ny2+ny

#        print "her",nx,ny,dx,dy,NX,NY
        #if domain is larger than grid (self.x and self.y), the arrays
        #x,y,z includes zones with fieldvalue zero
        
        x=numpy.zeros([NX,NY],float)
        y=numpy.zeros([NX,NY],float)
        z=numpy.zeros([NX,NY],float)+value
        #indices i,j is for orig. data, I,J is for new 
        if nx1>=0:
            i0=0
            I0=nx1
        elif nx1<0:
            i0=-nx1
            I0=0
        if nx2>=0:
            im=nx-nx2
            Im=nx1+nx-nx2
        elif nx2<0:
            im=nx
            Im=nx+nx1
        if ny1>=0:
            j0=0
            J0=ny1
        elif ny1<0:
            j0=-ny1
            J0=0
        if ny2>=0:
            jm=ny-ny2
            Jm=ny1+ny-ny2
        elif ny2<0:
            jm=ny
            Jm=ny+ny1
#        print "indexes",i0,im,I0,Im,j0,jm,J0,Jm
        x[I0:Im,J0:Jm]=self.x[i0:im,j0:jm]
        y[I0:Im,J0:Jm]=self.y[i0:im,j0:jm]
        z[I0:Im,J0:Jm]=self.z[i0:im,j0:jm]
        #print z[I0:Im,J0:Jm]
        #fill up grid outside original data area (x and y).

        #x: #if I0>0:
        for i in range(I0-1,-1,-1):
            x[i,:]=x[I0,J0]+dx*(i-I0)
            #print "1",i,x[I0,J0],I0,dx,x[I0,J0]-dx*(I0-i)
        for i in range(I0,Im):
            x[i,:]=x[i,J0]
            #print "2",i
        for i in range(Im,NX):
            x[i,:]=x[Im-1,Jm-1]+dx*(i-Im)
            #print "3",i
        #y:
        for j in range(J0-1,-1,-1):
            y[:,j]=y[I0,J0]+dy*(j-J0)
            #print "4",j
        for j in range(J0,Jm):
            y[:,j]=y[I0,j]
            #print "5",j
        for j in range(Jm,NY):
            #print "6",j
##             print Jm,NY,j,Im,Jm,dy,j,Jm
##             print y
            y[:,j]=y[Im-1,Jm-1]+dy*(j-Jm)



        
        #print "x",x,self.x
        #print "y",y,self.y
        #print "z",z,self.z
        #makeing slices:
        #self.x=self.x[i0:im,j0:jm]
        #self.y=self.y[i0:im,j0:jm]
        #self.z=self.z[i0:im,j0:jm]
        #self.update(self.x,self.y,self.z)
        self.update(x,y,z)
        self._update_grid_parameters()
        return self.x,self.y,self.z

    def _inside(self,x,y):
        """checks if point (x,y) is inside grid"""
        if self.x0<=x<=self.xm and self.y0<=y<=self.ym:
            return True
        else:
            return False
    
    def value_pt(self,x,y):
        """returns the z-value of point (x,y)"""
        #check id x,y is outside
        if self._inside(x,y):
            #interpolated value at x,y
            ni=0
            while self.x[ni,0]<=x:
                ni+=1
            iv=ni-1
            #interpolate values from the y-columns i and i+1:
            value=self._interp_col(iv,x,y)
            return value[0]
        else:
            print ("Point to be evaluated ("+str(x)+","+str(y)+\
                  ") is outside grid, returning zero." )
            return None
        
    def replace_value(self,arr,let=None,lat=None,value=-666):
        """modifies values of self.z if value of same point in
        'arr' is less than 'let' or larger than 'lat'. The value at this
        point in 'arr' is then replaced by 'value'. self.z and 'arr' are the
        values of a field and must be on identical grids.  A typical use of
        this functionality is to replace surface elevation over land by
        a value outside the plotting range. 'arr' is then the depth matrix."""
        
        #may add functionality so that arr2 is interpolated into the
        #same grid as arr1 if they differ...

        if let!=None:
            code="arr[i,j]<%s" % let
        elif lat!=None:
            code="arr[i,j]>%s" % lat

        for i in range(self.nx):
            for j in range(self.ny):
                if eval(code):
                    self.z[i,j]=value
        
        return self.z

            
if __name__ == '__main__':

    
    from sct3.io import io

    import scipy,time
    file="Verify/test.xyz"
    file="mjolnir/ampH_w_pb.mtv"
    #file="/home/sylfest/SEAtsu/gebco1_3D/gphov/gbc1m_new.dat"
    #file="/home/sylfest/transf/xing/depth.mtv"
    pd=io()
    [x,y,z]=pd.read(file)
    g=Grid2D(x,y,z)
    xt,yt=scipy.mgrid[-800:800:70j,-800:800:70j]
    #print xt
    #print yt
    #z=xt.copy()
    #z*=0.0
    #z+=1
    #z[1,1]=5
    #z[0,2]=2
    #z[1,2]=0
    #print z
    #x=xt
    #y=yt
    #g.update(x,y,z)
    #x,y,z=g.zoom(-800,-800,800,800)
    print ( "xt",xt[69,69] )
    #wd=W(x,y,z)
    #wd.write("Azoom.mtv")
    #x,y,z=g.extrap(3)
    #print "len(x)",len(x)
    #wd=W(x,y,z)
    #wd.write("Bextrap.mtv")
    #w=W(x,y,z)
    #w.write("TESTIN-before.mtv")
    #x,y,z=g.scale(zsc=0.001)
    #wd=W(x,y,z)
    #wd.write("Cscale.mtv")
    #g.crossection(-300,-300,1000,1,dr=10, start=-666)
    #x,y,grad=g.gradient()
    #wd=W(x,y,grad)
    #wd.write("Dgradient.mtv")
    #g.update(x,y,grad)
    #t=time.time()

    x,y,z=g.interp_lin(xt,yt)
    #t=time.time()-t
    #print "Time for interp_lin in seconds:%4.2f" % t
    #g.update(X,Y,Z)
    #t=time.time()
    #x,y,z=g.interp_bspl(xt,yt,x,y,z)
    #t=time.time()-t
    #print "Time for interp_lin in seconds:%4.2f" % t
    #x,y,z=g.smooth(3)
    #print z

    pd.write("new.mtv",x,y,z)


#    pd.write("SFtest.mtv",x,y,z)
#    pd.write("SVtest.gphov",x,y,z)
    #x,y,z=g.zoom(-1000,-1000,1000,1000)
    #x,y,z=g.extrap(8)

    #wd=W(x,y,z)
    #wd.write("TESTING.mtv")
    

    #file="eta_1h.xyz"
    #file="Verify/test.xyz"
    #d=io()
    #[x,y,z]=d.read(file)
    #g=Grid2D(x,y,z)
    #file="depth.dat.xyz"
    #[x,y,z]=d.read(file)
    #g2=Grid2D(x,y,z)
    #g2.zoom(g.x0,g.y0,g.xm,g.ym)
    #x,y,z=g2.interp_lin(g.x,g.y)
    #z=g.replace_value(g2.z,let=0,value=-666)
    #d.write("tmp.mtv",x,y,z)
    
    
    
    
