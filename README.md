# BingCLAW_5.6.1
BingClaw software for simulating landslide dynamics (requires ClawPACK version 5.6.1)

BingCLAW is licensed under the Creative Commons BY-NC-ND license.

https://creativecommons.org/licenses/by-nc-nd/4.0/

Release v1.0.1 from 2023/09/18 is permanently available on Zenodo with DOI 10.5281/zenodo.8354763  

https://dx.doi.org/10.5281/zenodo.8354763  

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8354763.svg)](https://doi.org/10.5281/zenodo.8354763)

This repository contains the source code for BingClaw version 5.6.1 for simulating landslide dynamics, with an example case study.
We describe here how to set up and compile BingClaw and run a simulation of the Storegga landslide.

A comprehensive description of the simulation of the Storegga landslide is provided in the following publication:  
https://doi.org/10.1029/2018JC014893  
*(Jihwan Kim, Finn Løvholt, Dieter Issler, and Carl Fredrik Forsberg,  
Landslide Material Control on Tsunami Genesis—The Storegga Slide and Tsunami (8,100 Years BP),  
JGR Oceans, vol 124(6), June 2019, Pages 3607-3627)*   

See also https://www.ngi.no/eng/Services/Technical-expertise/Tsunamis/Model-for-simulating-dynamics-of-cohesive-landslides

This directory contains four files: BingClaw5.6.1_singlesource.tar setrun.py PaleoNorthAtlantic.tt3.gz storegga_ini.tt3

Here we outline the steps required to perform a simulation of the Storegga slide using BingClaw on a Linux system.

We will need clawpack (version 5.6.1), gfortran, python3, LAPACK and BLAS.

If you do not have LAPACK and BLAS, the code can be obtained from https://github.com/Reference-LAPACK/lapack Follow the instructions for compilation and you should end up with libraries liblapack.a and librefblas.a. You will need to modify a Makefile to point to these files.
You will also need HDF5. If you have superuser rights, then LAPACK, BLAS, and HDF5 are all installed using  

```
sudo apt-get install libblas-dev liblapack-dev  
sudo apt-get install libhdf5-serial-dev
```
If you get a message of the form *E: Unable to locate package libhdf5-dev* then you will most likely have to update apt-get, e.g.  
```
sudo apt-get update -y
sudo apt-get install -y libhdf5-dev
```
Should gfortran be lacking from your system, this is installed using
```
sudo apt install gfortran
```

If you do not have python3 and do not have superuser rights to your Linux, you can type

```
wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
```

followed by

```
bash Anaconda3-2020.11-Linux-x86_64.sh
```

and follow the instructions.

We then need to install and activate clawpack5.6.1.

Choose a directory in which to place it: e.g. ${HOME}/clawpack_src

Obtain the tar file *clawpack-v5.6.1.tar.gz* using  

```
wget https://github.com/clawpack/clawpack/files/3781620/clawpack-v5.6.1.tar.gz
```

then do  

```
tar xvf clawpack-v5.6.1.tar.gz
```

then  

```
cd clawpack-v5.6.1  
pwd
```

Then, add the following lines to ${HOME}/.bashrc file (or equivalent)

```
export CLAW=$HOME/clawpack_src/clawpack-v5.6.1
```

(setting this variable to what is shown on typing "pwd")  

```
export PYTHONPATH=$CLAW:$PYTHONPATH
```

We then need to compile BingClaw. If you then do

```
tar xvf BingClaw5.6.1_singlesource.tar  
cd BingClaw5.6.1_singlesource
```

you need to alter the files src/Makefile* such that the variables LAPACK and BLAS point to the correct places:

e.g.

```
LAPACK=/usr/lib/x86_64-linux-gnu/lapack/liblapack.a  
BLAS=/usr/lib/x86_64-linux-gnu/blas/libblas.a
```

and with the correct path for the HDF5 library, e.g.:

```
-I/usr/include/hdf5/serial
```

You should then be able to compile BingClaw5.6.1 using the following sequence:  

```
sh compile_src.sh
```

This generates the executable *bin/BingClaw5.6.1*

Finally, create a directory in which you want to run BingClaw Save the file setrun.py in it and make sure that the lines with "PaleoNorthAtlantic.tt3" and "storegga_ini.tt3" point to the correct location. (You will have to (g)unzip the file PaleoNorthAtlantic.tt3.gz before running the program.)    

Execute the setrun.py file using

```
python setrun.py
```

We then execute the BingClaw program you have just compiled, and this should be generate the required output.

There are many different ways of processing the output from the BingClaw simulations.
In the directory *BingClaw_fortfile_display_example* we present a simple solution for displaying the output from the Storegga simulation.

# Specification of a spatially varying initial yield strength    

It is possible to specify a file providing a spatial distribution of the initial yield strength.  
You can now set an integer flag "i_use_var_tau_y" to 1 to impose a spatially varying initial yield strength using absolute values, and to 2 to impose a spatially varying initial yield strength using "relative values" (i.e. a map of values 0 < val < 1 such that tauy = tauy_r + val*(tauy_i - tauy_r) ... ) The file containing the map of initial yield strength is specified using the variable name "fname_tau_y". Examples of setrun.py files and distribution files are found in the directories "absolute_initial_tauy_example" and "relative_initial_tauy_example".
Note that the the global values of tauy_i and tauy_r still apply, the spatially varying map just alters the total shear deformation at the start of the simulation. The value of tauy must always lie between tauy_i and tauy_r.




# Note on compiling 2023/03/09

The file Makefile.nonmodules requires the following line to allow it to compile in situations described on this page:
https://forum.mmm.ucar.edu/threads/how-to-fix-error-rank-mismatch-between-actual-argument-at-1-and-actual-argument-at-2.11998/#post-27521 

```
OPTIM=    -O3 -cpp -fallow-argument-mismatch
```


