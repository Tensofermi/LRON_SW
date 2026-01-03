#pragma once
#include <bits/stdc++.h>
#include "../../system/Header.hpp"
#include "../Lattice/Lattice.hpp"
#include "../BondTable/BondTable.hpp"
#include "../BondTable/BondTable_v2.hpp"
#include "../Wrapping/Cluster_Wrapping.hpp"


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
    // std::vector<std::vector<double>> Site;
    std::vector<int> Site;
    Hypercubic Latt;   
    // Triangle Latt;  
    BondTable_v2 Bond;       
    int Dim, Nspin, L, NNb;
    long Vol;
    double Beta;

    //--- Basic parameter for algorithms
    std::vector<int> Mem, Que;
    bool bondOrNot = false;    // check if 
    // ClusterFind cluster;
    ClusterFind_Wrapping cluster;
    std::vector<int> cluster_flip_tag;
    double P_metro, P_sw;

    //--- For spin flip and k-space
    double tempSpin;
    std::vector<double> k_vec, k_cos, k_sin;

    //--- For long-range interactions
    int Nlist;
    std::vector<std::vector<int>> range_list;
    std::vector<double> dis_list, strength_list;
    std::vector<long> clock_order;
    std::vector<double> clock_rho;
    std::vector<int> shift_, k_;

    //--- Basic Observables
    std::vector<int> x_max, x_min, x_now;
    std::vector<int> Tx, Ty, listWrap_x, listWrap_y;
    long NCluster;                  // Number of clusters
    double C1, C2, S2, S4;          // Cluster size defined by particle number
    double S1;                      // Shortest path
    double Rx, Ry, R0, R1, Rd, R2;
    void BFS_Bond();

public:
std::string infoConfig()
{
    return 
        "====================\n"
        "This program simulates the long-range O(1) spin model, with 1/r^{2+σ} interaction strength.\n"
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
    void SwendsenWang_LR_Clock_Bond(); 

};
