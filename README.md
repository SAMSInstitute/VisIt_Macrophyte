## Analyze IBAMR data via VisIt/Python

### Parsing IBAMR data into average velocities
- You must be running Python 2.7 to parse data from IBAMR.
- You must have VisIt installed
- You must edit the Python file with the correct system path to the visit 
  lib/site-packages folder

collect_avgs_py27.py is for parsing a single IBAMR simulation
collect_avgs_loop.py is for looping over several IBAMR simulations, parsing 
    each in turn

All averages and dp/dx are stored in the folder MacrophyteAvgs.

### Visualizing the averages
This can be done in Python 3 using plot_avgs.py. Python 2 may work too, but I 
    haven't tried it.

This is a bit hacky at the moment, but you should import plot_avgs into an 
    IPython session and then run its functions from there. The good news is 
    that each individual function is decently documented, so it shouldn't be 
    too difficult to find your way from there. Loading of avg data is done 
    automatically based on the paramters you pass.
