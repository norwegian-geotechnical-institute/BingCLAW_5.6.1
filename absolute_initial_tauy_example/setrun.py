

""" 
Module to set up run time parameters for Clawpack -- AMRClaw code.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    
"""

#Parameters

rho_a = 1000.0        # Density of ambient fluid
rho_s = 1860.0        # Density of slide
n_param = 0.5        # Bing rheology parameter
gamma_r = 100.0        # Reference strain rate
c_mass = 0.1          # Added Mass
hydrodrag = True      # Use hydrodrag force?
cF_hyd = 0.01
cP_hyd = 1.0
remolding = True     # remolding? if not, only tauy_i is used
tauy_i =   12000.0   # initial yield strength (Pa)
tauy_r =   3000.0    # residual yield strength (Pa)
remold_coeff = 5.e-4  # remolding parameter

#Initial Condition
qinit_style = 3       # topotype style initial conditions
i_use_var_tau_y = 1                    # Initial yield strength file with absolute tau values.
fname_tau_y = 'absolute_tau_init.tt3'  # File for the variable yield stress
#
# The file specified above must contain, as a function of location, initial values of tau
# from which the initial shear deformation is calculated.
# These values of tauy must lie between the (global) initial yield stress tauy_i and
# the (global) residual yield stress tauy_r.
#
# i_use_var_tau_y = 0    - do not set a spatially varying iniital yield stress
# i_use_var_tau_y = 1    - set a spatially varying iniital yield stress with absolute values
# i_use_var_tau_y = 2    - set a spatially varying iniital yield stress with relative values
#                          (i.e. with 0 < val < 1 such that tauy = tauy_r + val*(tauy_i - tauy_r)


# Variables for variables tau


import os
import numpy as np

#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """ 
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData 
    
    """

    from clawpack.clawutil import data


    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------
    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('rho_a'       ,  rho_a       ,'Density of ambient fluid')
    probdata.add_param('rho_s'       ,  rho_s       ,'Density of slide')
    probdata.add_param('n_param'     ,  n_param     ,'Bing rheology parameter')
    probdata.add_param('gamma_r'     ,  gamma_r     ,'Reference strain rate')
    probdata.add_param('c_mass'      ,  c_mass      ,'Added mass coeff.')
    probdata.add_param('hydrodrag'   ,  hydrodrag   ,'Hydrodrag?')
    probdata.add_param('cF_hyd'      ,  cF_hyd      ,'Hydrodrag coeff.')
    probdata.add_param('cP_hyd'      ,  cP_hyd      ,'Hydrodrag coeff.')
    probdata.add_param('remolding'   ,  remolding   ,'Remolding?')
    probdata.add_param('tauy_i'      ,  tauy_i      ,'Initial yield strength')
    probdata.add_param('tauy_r'      ,  tauy_r      ,'Residual yield strength')
    probdata.add_param('remold_coeff',  remold_coeff,'Remolding param. gamma')
    probdata.add_param('qinit_style' ,  qinit_style ,'Data style of the IC')
    probdata.add_param('i_use_var_tau_y',i_use_var_tau_y ,'Use variablt tau_y (0,1,2)?')
    probdata.add_param('fname_tau_y' ,  fname_tau_y      ,'File name for variable tau_y')


    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = -5.0
    clawdata.upper[0] = 7.0
    clawdata.lower[1] = 62.0
    clawdata.upper[1] = 68.0

    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] =          313
    clawdata.num_cells[1] =          334

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 6

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 7

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2


    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False              # True to restart from prior results
    clawdata.restart_file = 'fort.chk00028'  # File to use for restart data


    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output ntimes frames at equally spaced times up to tfinal:
        # Can specify num_output_times = 0 for no output
        clawdata.num_output_times = 1
        clawdata.tfinal = 60
        # clawdata.num_output_times = 180
        # clawdata.tfinal = 60.*60*6.
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times = np.linspace(20,60,81)*60

    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 3
        clawdata.output_t0 = True  # output at initial (or restart) time?


    clawdata.output_format == 'ascii'      # 'ascii', 'binary', 'netcdf'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'all'  # could be list
    clawdata.output_aux_onlyonce = False    # output aux arrays only at t0


    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==True:  variable time steps used based on cfl_desired,
    # if dt_variable==Falseixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True
    # Initial time step for variable dt.  
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 0.015

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used 
    clawdata.cfl_desired = 0.4
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = .45

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 50000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 1

    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'

    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 0

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 6

    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = [3,3,3,3,3,3]

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms

    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 1


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'   # at xlower
    clawdata.bc_upper[0] = 'extrap'   # at xupper

    clawdata.bc_lower[1] = 'extrap'   # at ylower
    clawdata.bc_upper[1] = 'extrap'   # at yupper


    # ---------------
    # Gauges:
    # ---------------
    #rundata.gaugedata.gauges = []

    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
      # Do not checkpoint at all
      pass

    elif clawdata.checkpt_style == 1:
      # Checkpoint only at tfinal.
      pass

    elif clawdata.checkpt_style == 2:
      # Specify a list of checkpoint times.  
      clawdata.checkpt_times = [0.1,0.15]

    elif clawdata.checkpt_style == 3:
      # Checkpoint every checkpt_interval timesteps (on Level 1)
      # and at the final time.
      clawdata.checkpt_interval = 5

    # ---------------
    # AMR parameters:   (written to amr.data)
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 1

    # List of refinement ratios at each level (length at least amr_level_max-1)
    amrdata.refinement_ratios_x = [1]
    amrdata.refinement_ratios_y = [1]
    amrdata.refinement_ratios_t = [1]


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length num_aux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).
    amrdata.aux_type = ['center', 'capacity', 'yleft', 'center', 'center','center','center']


    # Flag for refinement based on Richardson error estimater:
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.0  # Richardson tolerance

    # Flag for refinement using routine flag2refine:
    amrdata.flag2refine = True      # use this?
    amrdata.flag2refine_tol = 0.5  # tolerance used in this routine
    # Note: in geoclaw the refinement tolerance is set as wave_tolerance below 
    # and flag2refine_tol is unused!

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.7

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0


    # ---------------
    # Regions:
    # ---------------
    regions = rundata.regiondata.regions
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    regions.append([1, 3, 0., 1e9, -5000, 5000, -6000, 6000])

    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting

    return rundata

    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    """
    try:
        geo_data = rundata.geo_data
    except:
        print ("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system =  2
    geo_data.earth_radius = 6367500.0

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.
    geo_data.dry_tolerance = 0.001
    geo_data.friction_forcing = False
    geo_data.manning_coefficient = 0.0
    #geo_data.manning_coefficient = [0.025,0.04]
    #geo_data.manning_break = [0.]
    geo_data.friction_depth = 0.0

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.02
    refinement_data.deep_depth = 100.0
    refinement_data.max_level_deep = 4

    # == settopo.data values ==
    rundata.topo_data.topofiles = []
    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]
    fname = 'PaleoNorthAtlantic.tt3'
    # fname = "/home/stg/NGI/P/2021/05/20210583/Background-Others/Bathymetry/PaleoNorthAtlantic.tt3"
    topofiles.append([-3, 1, 1, 0., 1.e10, fname])

    # == setdtopo.data values ==
    dtopofiles = rundata.dtopo_data.dtopofiles
    # for moving topography, append lines of the form :  
    #   [topotype, minlevel,maxlevel,fname]

    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 1
    # 1 => height, 4 => eta
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]
    fname = 'storegga_ini.tt3'
    # fname = "/home/stg/NGI/P/2021/05/20210583/Background-Others/Bathymetry/storegga_ini.tt3"
    # rundata.qinit_data.qinitfiles.append([1, 2, fname])
    rundata.qinit_data.qinitfiles.append([0, 0, fname])

    # == fixedgrids.data values ==
    rundata.fixed_grid_data.fixedgrids = []
    fixedgrids = rundata.fixed_grid_data.fixedgrids
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,    #  ioutarrivaltimes,ioutsurfacemax]
    #fixedgrids.append([0.,3600.*6,301,-4.,6.,62,67,10*30+1,5*30+1,0,0])

    # == fgmax.data values ==
    fgmax_files = rundata.fgmax_data.fgmax_files
   # for fixed grids append to this list names of any fgmax input files
    #fgmax_files.append('fgmax_grid.txt')

    return rundata
    # end of function setgeo
    # ----------------------

if __name__ == '__main__':
    """
    if os.path.exists('fgmax_grid.txt'):
        print ("File fgmax_grid.txt exists, not regenerating")
    else:
        try:    
            fname = 'make_fgmax_grid.py'
            execfile(fname)
            print ("Created fixed grid data by running ",fname)
        except:
            print ("Did not find fixed grid specification ", fname)
    """

    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
