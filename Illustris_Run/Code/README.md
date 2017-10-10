# All Masses are in 10^10 solar masses
# All distances are in kpc

Remove all data and start over:

	rm ../data/*.csv
	rm ../data/phaseSpace/*.csv
	rm ../data/nnEnv/*.data
	rm ../data/twEnv/*.csv

-----------------------------------------------------------------------------------------------------------------------------

Obtain files for ../data/phaseSpace from Illustris simulations 1,2,3:

	python ./Py_libs/Illustris_Cluster.py

Format of ../data/phaseSpace files:

	x,y,z,vx,vy,vz,mg,mdm,mstellar,mbh

-----------------------------------------------------------------------------------------------------------------------------

Obtain files for ../data/nnEnv (nearest neighbor environment) from ../data/phaseSpace:

	./C_libs/env.out ../data/phaseSpace/Illustris1.csv ../data/nnEnv/env1.data ../data/nnEnv/specs1.data
	./C_libs/env.out ../data/phaseSpace/Illustris2.csv ../data/nnEnv/env2.data ../data/nnEnv/specs2.data
	./C_libs/env.out ../data/phaseSpace/Illustris3.csv ../data/nnEnv/env3.data ../data/nnEnv/specs3.data

Format of ../data/nnEnv/env*.data (Needed for MCMCfits.py to make graphics):
	
	d^2,mg,mdm

	(squared distance to the 3rd nearest neighbor), (mass gas of center galaxy), (mass dark matter of center galaxy)

_Format of ../data/nnEnv/specs*.data complete data to compare numeric and mass densities:_
	
	mg,mdm,mst,mbh,dist proyected,dist tridimensional

	masses include all 3 nearest neighbors.
	distances are not squared.

----------------------------------------------------------------------------------------------------------------------------

Another way to create needed files in ../data/ :

	make

----------------------------------------------------------------------------------------------------------------------------

_Obtain graphics of Schechter fits:_
	
	cd Py_libs
	python MCMCfits.py
	
	



