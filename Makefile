all:
	g++ -o3 -march=native -std=c++14 -Wall infoli.cpp infoli_sim_impl.cpp -o infoli_sim_impl -lm -ffast-math -fopenmp -mavx2 
