# BingCLAW_5.6.1
BingClaw software for simulating landslide dynamics (requires ClawPACK version 5.6.1)

BingCLAW is licensed under the Creative Commons BY-NC-ND license.

https://creativecommons.org/licenses/by-nc-nd/4.0/

This repositoey contains the source code for BingClaw version 5.6.1 for simulating landslide dynamics, with an example case study.
We describe here how to set up and compile BingClaw and run a simulation of the Storegga landslide.

See https://www.ngi.no/eng/Services/Technical-expertise/Tsunamis/Model-for-simulating-dynamics-of-cohesive-landslides

This directory contains four files: BingClaw5.6.1_singlesource.tar setrun.py PaleoNorthAtlantic.tt3.gz storegga_ini.tt3

Here we outline the steps required to perform a simulation of the Storegga slide using BingClaw on a Linux system.

We will need clawpack (version 5.6.1), python3, LAPACK and BLAS.

If you do not have LAPACK and BLAS, the code can be obtained from https://github.com/Reference-LAPACK/lapack Follow the instructions for compilation and you should end up with libraries liblapack.a and librefblas.a. You will need to modify a Makefile to point to these files.

If you do not have python3 and do not have superuser rights to your Linux, you can type

wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh

followed by

bash Anaconda3-2020.11-Linux-x86_64.sh

and follow the instructions.

We then need to install and activate clawpack5.6.1.

Choose a directory in which to place it: e.g. ${HOME}/clawpack_src

Run the command: pip install --src=$HOME/clawpack_src --user -e git+https://github.com/clawpack/clawpack.git@v5.6.1#egg=clawpack-v5.6.1

Then, add the following lines to ${HOME}/.bashrc file (or equivalent)

export CLAW=$HOME/clawpack_src/clawpack-v5-6-1 export PYTHONPATH=$CLAW:$PYTHONPATH

We then need to compile BingClaw. If you then do

tar xvf BingClaw5.6.1_singlesource.tar cd BingClaw5.6.1_singlesource

you need to alter the Makefile such that the variables LAPACK and BLAS point to the correct places:

e.g.

LAPACK=/usr/lib/x86_64-linux-gnu/lapack/liblapack.a BLAS=/usr/lib/x86_64-linux-gnu/blas/libblas.a

You should then be able to compile BingClaw5.6.1 using the following sequence:

make -f Makefile.modules make -f Makefile.nonmodules make BingClaw5.6.1

This generates the executable BingClaw5.6.1

Finally, create a directory in which you want to run BingClaw Save the file setrun.py in it and make sure that the lines with "PaleoNorthAtlantic.tt3" and "storegga_ini.tt3" point to the correct location. (You will have to (g)unzip the file PaleoNorthAtlantic.tt3.gz before running the program.)

Execute the setrun.py file using

python setrun.py

We then execute the BingClaw program you have just compiled, and this should be generate the required output.

Further scripts on visualizing and analyzing the output will be provided at a later date.
