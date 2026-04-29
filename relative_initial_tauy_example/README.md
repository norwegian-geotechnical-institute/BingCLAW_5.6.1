A directory containing a sample setrun.py file to implement the specification of a spatially varying initial yield stress using values (val) between 0 and 1 in order to specify an initial tauy between tauy_i (val)=1 and tauy_r (val=0).  
This is to say that at the start of the simulation, the value tauy is given by
```
tauy = tauy_r + val*(tauy_i - tauy_r)
```

Please look at the file **setrun.py** in this directory.
It is *almost* identical to the **setrun.py** in the folder above but there are two crucial differences.  

Firstly, we need to set the following lines near the top:  
```
i_use_var_tau_y = 2                    # Initial yield strength file with absolute tau values.
fname_tau_y     = 'relative_tau_init.tt3'  # File for the variable yield stress
```

(this is lines 28 and 29 in the **setrun.py** file in this folder).  

but then, you also need to add the following two lines further down:  
```
    probdata.add_param('i_use_var_tau_y',i_use_var_tau_y ,'Use variablt tau_y (0,1,2)?')
    probdata.add_param('fname_tau_y' ,  fname_tau_y      ,'File name for variable tau_y')
```

(this is lines 90 and 91 in the **setrun.py** file in this folder).  
