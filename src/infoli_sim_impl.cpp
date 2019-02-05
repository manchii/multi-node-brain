/*
Note: This license has also been called the "New BSD License" or "Modified BSD License". See also the 2-clause BSD License.

Copyright 2017 Erasmus Brain Project

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#******************************************************************************
#* Vivado_ZedBoard/hls/infoli_tb.cpp
#*
#* Written by: George Smaragdos.
#* Modified by : Carlos Salazar-García, 2017
#* This code is the testbech when you compile using gcc
#*
#*
#******************************************************************************/
#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <memory>
#include <array>
#include <vector>
#include <omp.h>
#include "infoli.h"

void initialize_cellstates(std::array<cellState,MAX_TIME_MUX> &IniArray){
  for(cellState &cell:IniArray){
    cell=InitState();
  }
}
void generate_conn_matrix(mod_prec *Connectivity_Matrix){
  for(int i=0; i<MAX_TIME_MUX;i++){
    Connectivity_Matrix[i]=CONDUCTANCE;
  }
}
void stimulus_iApp(auto &iAppArray, int simStep){
  mod_prec iApp = (simStep>20000-1 && simStep<20500-1)? 6.0f : 0.0f;
  for(auto &iA: iAppArray){
      iA = iApp;
  }
}
void appendToFileStim(std::ofstream &outFile, int simStep, mod_prec Stim){
  outFile <<simStep<<" "<<simStep*0.05<<" "<<Stim<<" ";
}
void appendToFileAxon(std::ofstream &outFile, const returnState &cellOut){
  const int max_print = std::min(MAX_TIME_MUX,4);
  for(int j=0;j<max_print;j++){
    outFile <<j<<": "<<cellOut.axonOut[j]<<" ";
  }
}

int main(int argc, char *argv[]){

  std::array<mod_prec,MAX_TIME_MUX> iAppArray;
  std::array<cellState,MAX_TIME_MUX> IniArray;
  std::array<mod_prec,MAX_TIME_MUX> Connectivity_Matrix;
  returnState cellOut;

  std::cout << "Inferior Olive Model ("<<TIME_MUX_FACTOR<<" cell network)"<< std::endl;
  std::ofstream outFile;
  outFile.open("InferiorOlive_Output.txt");
  outFile << "#simSteps Time(ms) Input(Iapp) Output(V_axon)" << std::endl;

  initialize_cellstates(IniArray);
  generate_conn_matrix(Connectivity_Matrix.data());

  const int simTime = SIMTIME; // in miliseconds
  const int simSteps = std::ceil(simTime/DELTA);
  const int N_Size = NUM_NEIGH_CELLS;
  const int Mux_Factor = MAX_TIME_MUX;

  for(int i=0;i<simSteps/10;i++){
    stimulus_iApp(iAppArray,i);
    appendToFileStim(outFile,i,iAppArray[0]);
    auto start = std::chrono::system_clock::now();
    ComputeNetwork( IniArray.data(), iAppArray.data(), N_Size, Mux_Factor,Connectivity_Matrix.data(), cellOut.axonOut);
    //ComputeNetwork_iC( IniArray.data(), iAppArray.data(), N_Size, Mux_Factor,Connectivity_Matrix.data(), cellOut.axonOut);
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time elapsed: " << elapsed.count() << std::endl;
    appendToFileAxon(outFile,cellOut);
    outFile<<std::endl;
  }
  std::cout << "End simulation" << std::endl;
  outFile.close();

  return 0;
}

cellState InitState(){
    cellState initState;
    //Initial dendritic parameters
    initState.dend.V_dend = -60;
    initState.dend.Calcium_r = 0.0112788;// High-threshold calcium
    initState.dend.Potassium_s = 0.0049291;// Calcium-dependent potassium
    initState.dend.Hcurrent_q = 0.0337836;// H current
    initState.dend.Ca2Plus = 3.7152;// Calcium concentration
    initState.dend.I_CaH   = 0.5;// High-threshold calcium current
    //Initial somatic parameters
    initState.soma.g_CaL = 0.68; //default arbitrary value but it should be randomized per cell
    initState.soma.V_soma = -60;
    initState.soma.Sodium_m = 1.0127807;// Sodium (artificial)
    initState.soma.Sodium_h = 0.3596066;
    initState.soma.Potassium_n = 0.2369847;// Potassium (delayed rectifier)
    initState.soma.Potassium_p = 0.2369847;
    initState.soma.Potassium_x_s = 0.1;// Potassium (voltage-dependent)
    initState.soma.Calcium_k = 0.7423159;// Low-threshold calcium
    initState.soma.Calcium_l = 0.0321349;
    // Initial axonal parameters
    initState.axon.V_axon = -60;
    //sisaza: Sodium_m_a doesn't have a state, therefore this assignment doesn'thave any effect
    initState.axon.Sodium_m_a = 0.003596066;// Sodium (thalamocortical)
    initState.axon.Sodium_h_a = 0.9;
    initState.axon.Potassium_x_a = 0.2369847;// Potassium (transient)

    return (initState);
}

/*
#******************************************************************************
#* Vivado_ZedBoard/hls/infoli_tb.cpp
#*
#* Written by: George Smaragdos.
#* Modified by : Carlos Salazar-García, 2017
#* This code is the testbech when you compile using gcc
#*
#*
#******************************************************************************/
