all:
	g++ -o3 -march=native -std=c++17 -Wall src/infoli.cpp src/infoli_sim_impl.cpp -o infoli_sim_impl -lm -ffast-math -fopenmp

cluster:
	icpc -O3 -inline-max-size=350 -qopt-report=5 -march=knl -std=c++17 -Wall src/infoli.cpp src/infoli_sim_impl.cpp -o infoli_sim_impl -lm -qopenmp -ffast-math
