# To run the bioreactor test-case

see folder VirtualEngineering/bioreactor/bubble_column
To run on eagle, you have just submit the ofoamjob script.
It takes about 6 hours to run to get to steady-state.

You can use the pv_extract_analyze_script.py script (also invoked through pvjob script) to 
analyse oxygen concentration and gas hold-up. This script uses paraview python, which is part of 
paraview module on eagle. This script writes a file called volume_avg.dat which stores the 
time history (first column) of oxygen concentration in mol/m3 (second column). 
The fifth column stores gas hold-up.
