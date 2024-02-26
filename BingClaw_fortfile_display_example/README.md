**BingCLAW fort file display example**


The BingClaw program writes out files with names fort.*

There are many different approaches to processing and visualizing these files; a familiarity with ClawPACK is useful here.

This directory contains one ad-hoc solution to converting and visualizing an output from the Storegga simulation using GMT software https://www.generic-mapping-tools.org/ 

It consists of a bash script *one_fortfile_2_png.sh* and three python scripts:

```
claw2file3.py
format3.py
sct3/io.py
sct3/Grid2D.py
```

The execution takes the form:  

```
USAGE: 
./one_fortfile_2_png.sh   ./fort.q0179 
```

but you will need to adjust paths for GMT and install all of the python dependencies.  

If successful, you should end up with a figure resembling this one:  

![png file generated of BingClaw simulation output](example_q0179.png)  

In the case of running GMT6, where the program gmt is in the user's path (e.g. under /usr/bin), then the script *gmt6_one_fortfile_2_png.sh* should work.
Running 

```
bash make_all_png_files.sh
```

should make one png file for each output file from BingCLAW (../fort.q*) and an animated gif file can be made by running

```
python3 make_animated_gif_file.py
```
