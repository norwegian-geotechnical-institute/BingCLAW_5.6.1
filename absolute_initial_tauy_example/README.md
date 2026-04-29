A directory containing a sample setrun.py file to implement the specification of a spatially varying initial yield stress using absolute values (between tauy_i and tauy_r)  

Please look at the file **setrun.py** in this directory.
It is *almost* identical to the **setrun.py** in the folder above but there are two crucial differences.  

Firstly, we need to set the following lines near the top:  
```
i_use_var_tau_y = 1                    # Initial yield strength file with absolute tau values.
fname_tau_y     = 'absolute_tau_init.tt3'  # File for the variable yield stress
```

(this is lines 28 and 29 in the **setrun.py** file in this folder).  

but then, you also need to add the following two lines further down:  
```
    probdata.add_param('i_use_var_tau_y',i_use_var_tau_y ,'Use variablt tau_y (0,1,2)?')
    probdata.add_param('fname_tau_y' ,  fname_tau_y      ,'File name for variable tau_y')
```

(this is lines 88 and 89 in the **setrun.py** file in this folder).  
