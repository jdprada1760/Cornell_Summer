all: ../data/nnEnv/env1.data ../data/nnEnv/env2.data ../data/nnEnv/env3.data ../data/twEnv/Tweb1.csv ../data/twEnv/Tweb2.csv ../data/twEnv/Tweb3.csv

../data/twEnv/Tweb1.csv ../data/twEnv/Tweb2.csv ../data/twEnv/Tweb3.csv: ../data/phaseSpace/Illustris1.csv ../data/phaseSpace/Illustris2.csv ../data/phaseSpace/Illustris3.csv ./Py_libs/read_tweb.py
	python ./Py_libs/read_tweb.py

../data/twEnv/Tweb1.csv: ../data/phaseSpace/Illustris1.csv ./Py_libs/read_tweb.py

../data/nnEnv/env1.data: ../data/phaseSpace/Illustris1.csv ./C_libs/env.out
	./C_libs/env.out ../data/phaseSpace/Illustris1.csv ../data/nnEnv/env1.data

../data/nnEnv/env2.data: ../data/phaseSpace/Illustris2.csv ./C_libs/env.out
	./C_libs/env.out ../data/phaseSpace/Illustris2.csv ../data/nnEnv/env2.data

../data/nnEnv/env3.data: ../data/phaseSpace/Illustris3.csv ./C_libs/env.out
	./C_libs/env.out ../data/phaseSpace/Illustris3.csv ../data/nnEnv/env3.data

./C_libs/env.out:
	gcc -Wall ./C_libs/environment.c -o ./C_libs/env.out -lm -fopenmp

../data/phaseSpace/Illustris1.csv ../data/phaseSpace/Illustris2.csv ../data/phaseSpace/Illustris3.csv:
	python ./Py_libs/Illustris_Cluster.py

