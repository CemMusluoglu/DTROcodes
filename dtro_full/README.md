This folder contains the code to obtain the figures shown in [1].  
The following scripts contain implement independent Monte-Carlo runs, and can therefore take hours to finish computing.  
`dtro_full_mc_K.m` : Script to run the DTRO algorithm for various number of nodes *K*.  
`dtro_full_mc_Q.m` : Script to run the DTRO algorithm for various number of projection dimensions *Q*.  
`dtro_full_mc_topo.m` : Script to run the DTRO algorithm for various network topologies.  
`dtro_full_mc.m` : Script to run the DTRO algorithm for one specific network.  

Each script outputs a '_.mat' file, which can then be used to plot the results shown in [1].  

`plot_func.m` : Script in the *plots* folder which plots the figures in [1]. Assumes the '_.mat' data files are in the folder above the *plots* folder.  

**EXTERNAL DEPENDENCIES**:  

* The *Random trees* package https://www.mathworks.com/matlabcentral/fileexchange/2516-random-trees .  

Used to create randomly generated trees.  

* The Graph Signal Processing Toolbox https://epfl-lts2.github.io/gspbox-html/ [2].  

Used for creating random graphs using the Erd&#337-R&#233nyi model.  

* The 'legendflex' package https://www.mathworks.com/matlabcentral/fileexchange/31092-legendflex-m-a-more-flexible-customizable-legend .  

Used for the legends in the plots.


## References ##

[1] C. A. Musluoglu and A. Bertrand, “Distributed adaptive trace ratio optimization in wireless sensor networks”.

[2] Perraudin Nathanaël, Johan Paratte, David Shuman, Lionel Martin, Vassilis Kalofolias, Pierre Vandergheynst and David K. Hammond}, GSPBOX: A toolbox for signal processing on graphs. Arxiv e-print, 08-2014.
