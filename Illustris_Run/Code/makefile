all: ../data/nnenv/env1.data ../data/nnenv/env2.data ../data/nnenv/env3.data

../data/nnenv/env1.data: ../data/phaseSpace/Illustris1.txt ./C_libs/env.out
	./C_libs/env.out ../data/nnenv/env1.data env1.data

../data/nnenv/env2.data: ../data/phaseSpace/Illustris2.txt ./C_libs/env.out
	./C_libs/env.out ../data/nnenv/env2.data env2.data

../data/nnenv/env3.data: ../data/phaseSpace/Illustris3.txt ./C_libs/env.out
	./C_libs/env.out ../data/nnenv/env3.data env3.data

./C_libs/env.out: ./C_libs environment.c
	gcc -Wall environment.c -o env.out -lm

./data/phaseSpace/Illustris1.txt: ./Py_libs/Illustris_Cluster.py
	python Illustris_Cluster.py

./data/phaseSpace/Illustris2.txt: ./Py_libs/Illustris_Cluster.py
	python Illustris_Cluster.py

./data/phaseSpace/Illustris3.txt: ./Py_libs/Illustris_Cluster.py
	python Illustris_Cluster.py