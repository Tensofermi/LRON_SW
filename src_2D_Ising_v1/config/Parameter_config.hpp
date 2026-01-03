#pragma once

struct Parameter
{
    // Simulation Parameters
    int Seed;
    unsigned long N_Measure;
    unsigned long N_Each;
    unsigned long N_Therm;
    unsigned long N_Total;
    unsigned long NBlock;
    unsigned long MaxNBin;
    unsigned long NperBin;

    // Model Parameters
    int Nspin, D, L;
    double sigma, beta;

    // Observable label
    int i_M, i_absM, i_M2, i_M4, i_E, i_E2;
    int i_NCluster, i_S2, i_S4, i_SM4, i_C1, i_C2;
    int i_Qm, i_Qm_fk, i_Cv, i_C12;
    int i_Mk2, i_corr_L;
    int i_M2_E, i_Mk2_E, i_K, i_Kk;
    int i_R0, i_R1, i_Rd, i_R2, i_Rp;
    int i_S1;

    int i_R_E, i_R_K, i_Rx;
    

    // Distribution label

};