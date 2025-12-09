#pragma once
#include <bits/stdc++.h>
#include "../../system/Header.hpp"
#include "../Lattice/Lattice.hpp"

class Configuration
{
    Clock& ck;
    IOControl& io;
    RandomNumGen& rn;
    Parameter& para;
    Observable& obs;
    Histogram& his;

public:
    //--- Configuration
    std::vector<std::vector<double>> Site;
    Hypercubic Latt;   
    // Triangle Latt;         
    int Dim, Nspin, L, NNb;
    long Vol;
    double Beta;

    //--- Basic parameter for algorithms
    std::vector<int> Mem, Que;
    double P_metro, P_sw;

    //--- For spin flip and k-space
    std::vector<double> tempSpin;
    std::vector<double> k_vec, k_cos, k_sin;

    //--- For long-range interactions
    int Nlist;
    std::vector<std::vector<int>> range_list;
    std::vector<double> dis_list, strength_list;
    std::vector<long> clock_order;
    std::vector<double> clock_rho;

    //--- Basic Observables
    std::vector<int> x_max, x_min, x_now;
    long NCluster;                  // Number of clusters
    double C1, C2, S2, S4;          // Cluster size defined by particle number

public:
std::string infoConfig()
{
    return 
        "====================\n"
        "This program simulates the long-range O(N) spin model, with 1/r^{d+σ} interaction strength.\n"
        "It can simulate the O(N) model in any spatial dimension d and with any spin dimension N.\n"
        "====================\n";
};


    Configuration(Clock& _ck, IOControl& _io, RandomNumGen& rn, Parameter& _para, Observable& _obs, Histogram& _his);
    void initialConf();
    void initialRang();
    void initialAlgo();
    void initialMeas();
    void initialObsr();

    void updateCnf();
    bool measureOrNot();
    void measure();
    void writeCnf();

    void printConfig(int _index);
    void squarePrint();
    // void corrFunPrint();

    void checkCnf();

    // Basic Operation
    void randSpin();
    void randSpin_();
    void flipSpin(long _site, std::vector<double> _axis);
    
    // Update Algorithm
    void Metropolis();
    void Metropolis_LR();
    void Metropolis_LR_Clock_Thinning();

    void Wolff();
    void Wolff_LR();
    void Wolff_LR_Clock_Thinning(); 

    void SwendsenWang();
    void SwendsenWang_LR();
    void SwendsenWang_LR_Clock_Thinning(); 


};
