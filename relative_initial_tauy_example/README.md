A directory containing a sample setrun.py file to implement the specification of a spatially varying initial yield stress using values (val) between 0 and 1 in order to specify an initial tauy between tauy_i (val)=1 and tauy_r (val=0).  
This is to say that at the start of the simulation, the value tauy is given by
```
tauy = tauy_r + val*(tauy_i - tauy_r)
```
