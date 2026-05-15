"""
Module for dtopo file

by Jihwan Kim

"""
import os

import numpy as np
import shutil

def dtopo_tt3(indir,outfile):
   """
   Generated tt3 file from fort.q output
   """
   from numpy import mod,zeros

   fid_data = indir+'/claw.data'

   fid=open(fid_data,'r')
   for i in range(15):
      line=fid.readline()

   line=fid.readline().split()
   t0  =float(line[0])

   for i in range(2):
      line=fid.readline()

   line=fid.readline().split()
   mt  =int(line[0])

   line=fid.readline().split()
   tf  =float(line[0])

   fid.close()

   dt =(tf-t0)/mt

   dir = indir

   infile="fort.qxxxx"

   infile1=dir+"/fort.q0000"
   fid= open(infile1,'r')
   fout=open(outfile,'w')

   for i in range(2):
      line=fid.readline()

   line=fid.readline().split()
   mx = int(line[0])
   line=fid.readline().split()
   my = int(line[0])
   line=fid.readline().split()
   xlow = float(line[0])
   line=fid.readline().split()
   ylow = float(line[0])
   line=fid.readline().split()
   dx = float(line[0])
   line=fid.readline().split()
   dy = float(line[0])

   line=fid.readline()

   fout.write("%10d \n" % mx)
   fout.write("%10d \n" % my)
   fout.write("%10d \n" % mt)
   fout.write("%20.10e \n" % xlow)
   fout.write("%20.10e \n" % ylow)
   fout.write("%20.10e \n" % (t0+dt))
   fout.write("%20.10e \n" % dx)
   fout.write("%20.10e \n" % dy)
   fout.write("%20.10e \n" % dt)

   h_orig=zeros(mx*my)
   k=0
   for j in range(my):
      for i in range(mx):
         line=fid.readline().split()
         h_orig[k]=float(line[0])
         k+=1
      line=fid.readline()

   fid.close()

   infile1=infile
   dz=zeros(mx*my)

   for k in range(1,mt+1):
      for ipos in range(10,6,-1):
         idigit= int(np.mod(k,10))
         infile1=infile1[0:ipos-1]+str(idigit)+infile1[ipos:10]
         k=k/10

      fid1= open(dir+'/'+infile1,'r')
      print(fid1)

      for i in range(9):
         line=fid1.readline().split()

      k=0
      for j in range(my):
         for i in range(mx):
            line1=fid1.readline().split()
            hnew=float(line1[0])
            dz[i+(my-j-1)*mx]=hnew-h_orig[k]
            k+=1
         line=fid1.readline()

      for j in range(my):
         for i in range(mx):
            k=i+j*mx

            fout.write("%20.10e " % dz[k])

         fout.write("\n")

      fid1.close()

   fout.close()

def topo_w_slide_tt3(outfile):
   """
   Generated tt3 file from fort.q0000 output
   """
   from numpy import mod,zeros

   infile="fort.q0000"
   dir="_output/"
   infile1=dir+infile

   fid= open(infile1,'r')
   fout=open(outfile,'w')

   grid=fid.readline()
   amr=fid.readline()
   mx=int(fid.readline().split()[0])
   my=int(fid.readline().split()[0])
   xlow=float(fid.readline().split()[0])
   ylow=float(fid.readline().split()[0])
   dx=float(fid.readline().split()[0])
   dy=fid.readline()
   nodata=fid.readline()
   nodata=-9999

   fout.write("%10d \n" % mx)
   fout.write("%10d \n" % my)
   fout.write("%20.10e \n" % xlow)
   fout.write("%20.10e \n" % ylow)
   fout.write("%20.10e \n" % dx)
   fout.write("%20.10e \n" % nodata)

   h=zeros((mx,my))
   eta=zeros((mx,my))

   for j in range(my):
      for i in range(mx):
         line=fid.readline().split()
         h[i,j]=float(line[0])
         eta[i,j]=float(line[6])
      fid.readline()
   fid.close()

   for j in range(my):
      for i in range(mx):
         fout.write("%20.10e " % eta[i,my-j-1])

      fout.write("\n")

   fout.close()

if __name__=='__main__':
    indir = '../BingCLAW_sim_GEO'
    outdir = "./makedtopo_output"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outfile1 = os.path.join(outdir,"dtopo_test0.tt3" )
    dtopo_tt3(indir,outfile1)
    #outfile2 = os.path.join(outdir,"bathymetry_with_slide.tt3")
    #topo_w_slide_tt3(outfile2)
