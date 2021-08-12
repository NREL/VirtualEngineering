# Compilation

1) compile openfoam-dev using instructions in VirtualEngineering/README.md

2) set openfoam's bashrc - source <path-to-openfoam-folder>/etc/bashrc

3) Go to EHFoam folder and do wmake

4) Try out any of the cases in tests, for e.g. 
RushtonReact is a production case for lignocellulose digestion. 
* The biomass inputs are in constant/globalVars. 
* constant/EHProperties has transport coefficients and some solver knobs/parameters. 
* constant/MRFproperties has the rotational speed

5) The coupled runs creates a file called integrated_quantities.dat which 
has integrated values of solids volume fractions and averaged dissolved sugar concentrations.
This file can be postprocessed to obtain conversion data
