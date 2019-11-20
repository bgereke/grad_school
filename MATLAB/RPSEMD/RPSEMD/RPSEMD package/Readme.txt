
MAIN Matlab CODEs.

main_Rpsemd.m
	The example of running the Regenerated Phase-shifted Sinusoids assisted Empirical Mode Decomposition (RPSEMD) code. 
	The way of computing the indexs of Orthogonality of Successive IMFs is also given. 

Rpsemd.m 
	The main function of the algorithm, which is the outer loop of RPSEMD to decompose signal Y into a series of pure IMFs.



SUBFUNCTIONS:

Clustering.m
	To execute the hierarchical clustering.

Cluster_IMF.m
	To cluster the modes in an IMF.

Design_Asin_Fsin.m
	To design the amplitude and frequency of the auxiliary sinusoid through the clustered IMs.

Phaseshift.m
	One stage of decomposition to produce a pure IMF through phase shift the sinusoid.



----------------------------------------------------------------------------------------------------------------------------------------

CITED AND REFERRED Matlab CODEs 
We thank to G. Rilling, P. Flandrin and their groups because a part of the free EMD toolbox is provided them and downloaded from http://perso.ens-lyon.fr/patrick.flandrin/emd.html. 

Compu_instA_instF.m
	To compute instantaneous ampitudes instA and frequencies instF of an IMF.

Detect_extrema.m
	To detect the indexes of extrema, including minima and maxima.

EMD File:
	EMD Codes downloaded from http://perso.ens-lyon.fr/patrick.flandrin/emd.html.



-----------------------------------------------------------------------------------------------------------------------------------------

DATA

ECG:  It is an ECG from the MIT-BIH Normal Sinus Rhythm Database, which is available at http://www.physionet.org/cgi-bin/atm/ATM.

simulation_generation.m
	To generate the simulation data in the manuscript.


----------------------------------------------------------------

Written by Wang Chenxing, Jan 13, 2016, w.chenxing@gmail.com
Copyright Preserved.