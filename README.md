# l2rn
Random Ecological Reaction Network Creator

This program generates a random ecological network in Matlab from [com_assambly](https://github.com/pmaldona/com_assembly)
repository (LVG model ecological community matrix and growth rate vector). 
This is then used as input for the repository program [crn](https://github.com/pmaldona/crn) in R
in order to make structural calculations of these reaction networks. 


The function ``lvm2mat.m``, requires the functions present in the repository
[com_assambly](https://github.com/pmaldona/com_assembly/src). This function generates a ``.mat``
file which contains the communities matrix and growth rate vectors that are used by the ``lvm2rn`` function. 

The ``lvm2rn.R`` function requires the installation of the [crn](https://github.com/pmaldona/crn) package.
This function takes as input a .mat file generated by lvm2mat.m and returns a reaction network
object that can be used to calculate the different structural properties from the repository 
[crn](https://github.com/pmaldona/com_assembly/src).  
