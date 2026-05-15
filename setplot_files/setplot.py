""" 
Set up the plot figures, axes, and items to be done for each frame.
This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters. 
""" 

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from clawpack.visclaw import colormaps
    from clawpack.geoclaw import fgmax_tools, geoplot
    import numpy,pylab

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from pyclaw.plotters import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)
    def topo_q(current_data):
        q=current_data.q
        b=current_data.aux
        var = q[0,:,:]+b[0,:,:]
        return var

    def slide_only(current_data):
        q=current_data.q
        var = numpy.ma.masked_where(q[0,:,:]>0.,q[0,:,:])
        return var

    def topo_wo_slide(current_data):
        h=current_data.q[0,:,:]
        b=current_data.aux[0,:,:]
        var = numpy.ma.masked_where(h>0.,b)
        return var
        
    def topo(current_data):
        b=current_data.aux
        var = b[0,:,:]
        return var

    def surface_or_depth(current_data):
        from numpy import ma, where
        q = current_data.q
        h = q[0,:,:]
        eta = q[6,:,:]
        topo = eta - h
        surface = ma.masked_where(h<=drytol, eta)
        depth = ma.masked_where(h<=drytol, h)
        surface_or_depth =  where(topo<0, surface, depth)
        return surface_or_depth

    def hillshade(array, azimuth, angle_altitude):  

        from numpy import gradient  
        from numpy import pi  
        from numpy import arctan  
        from numpy import arctan2  
        from numpy import sin  
        from numpy import cos  
        from numpy import sqrt  
        #from numpy import zeros  
        #from numpy import uint8  
          
        x, y = gradient(array)#,dx,dy)  
        slope = pi/2. - arctan(sqrt(x*x + y*y))  
        aspect = arctan2(-x, y)  
        azimuthrad = azimuth*pi / 180.  
        altituderad = angle_altitude*pi / 180.  
        shaded = sin(altituderad) * sin(slope) + cos(altituderad) * cos(slope) * cos(azimuthrad - aspect)  

        return 255*(shaded + 1)/2  

    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)
    #plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    #plotfigure.kwargs = {'figsize': (18.5,10.5)}
    plotaxes.scaled = True
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    #plotaxes.afteraxes = addgauges
    
    drytol = 1e-4

    def fixup(current_data):
        t = current_data.t
        t = t / 60.  # minutes
        x = current_data.x
        y = current_data.y
        h = current_data.q[0,:,:]
        b = current_data.aux[0,:,:]
        tmp = numpy.where(h==0.,b,numpy.nan)
        plt.title('Landslide at %4.2f minutes' % t, fontsize=20)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        #colors = plt.cm.get_cmap('RdPu', 10)
        #plt.pcolor(x,y,h,cmap=colors)
        #hs_array = hillshade(tmp,20., 1.)
        #plt.pcolor(x,y,hs_array,cmap='Greys')

    plotaxes.afteraxes = fixup

    # Add contour lines of topography:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = topo
    from numpy import arange, linspace
    plotitem.contour_levels = linspace(0,100,21)
    plotitem.patchedges_show = 0
    plotitem.amr_contour_colors = ['.5']  # color on each level
    plotitem.kwargs = {'linestyles':'solid', 'linewidth': 1} #, 'alpha': 0.5}
    plotitem.amr_contour_show = [1,0,0]  # show contours only on coarsest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    # Landslide
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    colors = plt.cm.get_cmap('RdPu', 10)
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colors
    plotitem.pcolor_cmin =  0.
    plotitem.pcolor_cmax =  5.
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0
    plotitem.show = True
    
    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = topo_wo_slide
    plotitem.pcolor_cmap = 'gist_earth' #geoplot.land_colors
    plotitem.pcolor_cmin = 40.0
    plotitem.pcolor_cmax = 140.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0
    plotitem.show = True

    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Yield Strength', figno=2)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Max velocity'
    #plotfigure.kwargs = {'figsize': (18.5,10.5)}
    plotaxes.scaled = True
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'

    drytol = 1e-4
    
    # Landslide
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    def ys(current_data):
        a=current_data.q[5,:]
        var=600. + (10000.-600.)*numpy.exp(-a*5e-3)
        return var

    colors = plt.cm.get_cmap('RdPu', 10)
    plotitem.plot_var = ys
    plotitem.pcolor_cmap = 'Reds'
    plotitem.pcolor_cmin =  600.
    plotitem.pcolor_cmax =  10000.
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    def fixup(current_data):
        t = current_data.t
        t = t / 60.  # minutes
        pylab.title('Maximum speed at %4.2f minutes' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)

    plotaxes.afteraxes = fixup

    # Add contour lines of topography:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = topo
    from numpy import arange, linspace
    plotitem.contour_levels = linspace(0,100,21)
    plotitem.patchedges_show = 0
    plotitem.amr_contour_colors = ['.5']  # color on each level
    plotitem.kwargs = {'linestyles':'solid', 'linewidth': 1} #, 'alpha': 0.5}
    plotitem.amr_contour_show = [1,0,0]  # show contours only on coarsest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = topo_wo_slide
    plotitem.pcolor_cmap = 'gist_earth' #geoplot.land_colors
    plotitem.pcolor_cmin = 40.0
    plotitem.pcolor_cmax = 140.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0
    plotitem.show = True

    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='max velocity', figno=1)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Max velocity'
    #plotfigure.kwargs = {'figsize': (18.5,10.5)}
    plotaxes.scaled = True
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'

    drytol = 1e-4
    
    # Landslide
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    def upmax(current_data):
        b=current_data.aux
        var=b[7,:,:]
        return var

    colors = plt.cm.get_cmap('RdPu', 10)
    plotitem.plot_var = upmax
    plotitem.pcolor_cmap = 'Reds'
    plotitem.pcolor_cmin =  0.
    plotitem.pcolor_cmax =  10.
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0

    def fixup(current_data):
        t = current_data.t
        t = t / 60.  # minutes
        pylab.title('Maximum speed at %4.2f minutes' % t, fontsize=20)
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)

    plotaxes.afteraxes = fixup

    # Add contour lines of topography:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = topo
    from numpy import arange, linspace
    plotitem.contour_levels = linspace(0,100,21)
    plotitem.patchedges_show = 0
    plotitem.amr_contour_colors = ['.5']  # color on each level
    plotitem.kwargs = {'linestyles':'solid', 'linewidth': 1} #, 'alpha': 0.5}
    plotitem.amr_contour_show = [1,0,0]  # show contours only on coarsest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = topo_wo_slide
    plotitem.pcolor_cmap = 'gist_earth' #geoplot.land_colors
    plotitem.pcolor_cmin = 40.0
    plotitem.pcolor_cmax = 140.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0
    plotitem.show = True

    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='slide and topo', figno=3)
    plotfigure.show = False

    def x_slice(current_data):
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        q = current_data.q
        m = len(y)
        m = int(m/2)
        eta = q[5,:,m]
        return x,eta

    def y_slice(current_data):
        y = current_data.y[0,:]
        q = current_data.q
        eta = q[5,1,:]
        return y,eta

    def B_slice(current_data):
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        q = current_data.q
        m = len(y)
        m = int(m/2)
        h = q[0,:,m]
        eta = q[5,:,m]
        B = eta - h
        return x,B
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('water surface')
    plotaxes.title = 'line'
    plotaxes.xlimits = 'auto'#[20000,30000]
    plotaxes.ylimits = 'auto'#[-1500,-1450]
    
    # Topography
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = B_slice
    plotitem.color = 'k'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':2}

    # Slide Surface
    #plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    #plotitem.outdir = '../_output'
    #plotitem.map_2d_to_1d = x_slice
    #plotitem.color = 'b'
    #plotitem.plotstyle = '-'
    #plotitem.kwargs = {'linewidth':2}

    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='x slice', figno=5)
    #plotfigure.kwargs = {'figsize': (8,4)}
    plotfigure.show = False

    def x_slice(current_data):
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        q = current_data.q
        m = len(y)
        m = int(m/2)
        h = q[0,:,m]
        return x,h
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('landslide')
    plotaxes.title = 'landslide'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [0,5]
    
    # Slide Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'b'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':1}

    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='y slice', figno=6)
    #plotfigure.kwargs = {'figsize': (8,4)}
    plotfigure.show = False

    def y_slice(current_data):
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        q = current_data.q
        m = len(x)
        m = int(m/2)
        h = q[0,m,:]
        return y,h
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('landslide')
    plotaxes.title = 'landslide'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [0,5]
    
    # Slide Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = y_slice
    plotitem.color = 'b'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':1}

    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='cross_xdiff', figno=10)
    #plotfigure.kwargs = {'figsize': (8,4)}
    plotfigure.show = False

    def x_slice(current_data):
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        q = current_data.q
        m = len(y)
        m = int(m/2)
        h = q[0,:,m]
        return x,h
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('x slice')
    plotaxes.title = 'landslide'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [0,25]
    
    # Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = '../_output_dx20/'
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'b'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':1}

    # Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.outdir = '../_output_dx10/'
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'r'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':1}

    def add_legend(current_data):
        from pylab import plot, legend
        legend(['dx=20','dx=10'],'upper right')

    plotaxes.afteraxes = add_legend

    #-----------------------------------------
    # Figure for Cross section plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='combined', figno=7)
    plotfigure.kwargs = {'figsize': (8,4)}
    plotfigure.show = False

    def x_slice(current_data):
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        q = current_data.q
        m = len(y)
        m = int(m/2)
        eta = q[0,:,m]
        return x,eta
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('landslide')
    plotaxes.axescmd = 'axes([.1,.6,.8,.3])'
    plotaxes.title = 'landslide'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    
    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = x_slice
    plotitem.color = 'b'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':1}

    def x_slice_speed(current_data):
        x = current_data.x[:,0]
        y = current_data.y[0,:]
        q = current_data.q
        m = len(y)
        m = int(m/2)
        u = q[2,:,m]/q[0,:,m]
        return x,u
    
    # Figure for speed:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'axes([.1,.1,.8,.3])'
    plotaxes.title = 'velocity'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    
    # Water Surface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = x_slice_speed
    plotitem.color = 'b'
    plotitem.plotstyle = '-'
    plotitem.kwargs = {'linewidth':1}

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=300, \
                    type='each_gauge')
                    
    plotfigure.show = True
    plotfigure.kwargs = {'figsize': (12,2)}

    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-0.01, .2]
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    #plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    #plotitem.outdir='../_output'
    #plotitem.plot_var = 3
    #plotitem.plotstyle = 'b-'
    #plotitem.kwargs = {'linewidth':1}
    
    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'
    plotitem.kwargs = {'linewidth':1}

    # Plot topo as green curve:
    #plotitem = plotaxes.new_plotitem(plot_type='1d_plot')

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    #plotitem.plot_var = gaugetopo
    #plotitem.plotstyle = 'g-'
    
    def add_zeroline(current_data):
        from pylab import plot, legend
        t = current_data.t
        plot(t, 0*t, 'k')


    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = []         # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.format = 'ascii'                # Format of output
    # plotdata.format = 'netcdf'             

    return plotdata

    
