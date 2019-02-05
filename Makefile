all:
	g++ -o3 -march=native -std=c++17 -Wall src/infoli.cpp src/infoli_sim_impl.cpp -o infoli_sim_impl -lm -ffast-math -fopenmp -mavx2

cluster:
	icpc -O3 -march=native -std=c++17 -Wall src/infoli.cpp src/infoli_sim_impl.cpp -o infoli_sim_impl -lm -qopenmp -ffast-math
