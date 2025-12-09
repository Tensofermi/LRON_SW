#include "Configuration.hpp"

Configuration::Configuration(Clock& _ck, IOControl& _io, RandomNumGen& _rn, Parameter& _para, Observable& _obs, Histogram& _his) 
: ck(_ck), io(_io), rn(_rn), para(_para), obs(_obs), his(_his)
{
    initialConf();  // for site spin 
    initialRang();  // for long-range interactions
    initialAlgo();  // for algorithms
    initialMeas();  // for k-space and cluster measurement
    initialObsr();  // for basic observables
}

void Configuration::initialConf()
{
    //--- Initialize Basic Parameters
    Beta = para.beta;
    Dim = para.D;
    L = para.L;
    Nspin = para.Nspin;

    //--- Initialize Lattice
    Latt.set(Dim, L);
    // Latt.set({L,L});
    Vol = Latt.getVol();
    NNb = Latt.getNNb();
    Site.resize(Vol, std::vector<double>(Nspin));
    
    //--- Initialize Temporary Spin
    tempSpin.resize(Nspin);

    //--- Initialize Spin Configuration
    for (unsigned int i = 0; i < Vol; i++)
    {
        randSpin();
        Site[i] = tempSpin;
    }

}

void Configuration::initialRang()
{
    //--- Long-Range Interaction
    double norm = 0;   // check the normalization
    std::vector<int> dx, central_point;
    dx.resize(Dim);
    central_point.resize(Dim);

    // Set Central Point
    for (unsigned int i = 0; i < Dim; i++)
        central_point[i] = ceil(L / 2.0) - 1;
    
    // Set Interaction Lists
    for (unsigned int i = 0; i < Vol; i++)
    {
        double dis_square = 0;
        for (unsigned int j = 0; j < Dim; j++)
        {
            dx[j] = Latt.getComponent(i, j);
            dx[j] = dx[j] - central_point[j];

            dis_square += dx[j] * dx[j];
        }

        // for triangle lattice
        if(Dim==2 && NNb==6) dis_square = (0.5 * dx[0] + dx[1]) * (0.5 * dx[0] + dx[1]) + dx[0] * dx[0] * 3.0 / 4.0;
        
        if (dis_square == 0) continue;

        // storage vector from site to central point
        range_list.push_back(dx);   
        
        // storage distance from site to central point
        double r = sqrt(dis_square);
        dis_list.push_back(r);      

        // storage the power-law interacting strength as 1/r^{d+\sigma}
        double J = 1 / pow(r, (double)Dim + (double)para.sigma);
        strength_list.push_back(J);
        
        norm += strength_list.back();
    }
    Nlist = range_list.size();

    // // Normalize the Strength List
    double factor = NNb / norm;
    for (unsigned int i = 0; i < Nlist; i++)
        strength_list[i] *= factor;

    // // Check the Normalization
    double sum = 0;
    for (unsigned int i = 0; i < Nlist; i++)
        sum += strength_list[i];

    // io.exportInfo(io.OuputInfo, "NNb = " + toStr(NNb) + ", factor = " + toStr(factor) + ", Sum = " + toStr(sum) +"\n");

    // Dynamic Thinning Method
    clock_order.resize(Nlist);
    clock_rho.resize(Nlist);

    // initialize clock index
    std::iota(clock_order.begin(), clock_order.end(), 0);
    
    // sort the clock_order[] according to the order of dis_list[]
    std::stable_sort(clock_order.begin(), clock_order.end(),
                     [this](size_t i1, size_t i2)
                     { return dis_list[i1] < dis_list[i2]; });

    // generate upper bound probability 
    for (unsigned int i = 0; i < Nlist; i++)
    {
        clock_rho[i] = 1 - exp(- 2 * Beta * strength_list[clock_order[i]]);
        // io.exportInfo(io.OuputInfo, toStr(i) + ": " + toStr(clock_rho[i]) + "\n");
    }

}

void Configuration::initialAlgo()
{
    //--- Initialize Queue and Memory
    Que.resize(Vol);
    Mem.resize(Vol);
    for (int i = 0; i < Vol; i++)
    {
        Que[i] = 0;
        Mem[i] = 0;
    }

    //--- Initialize Probabilities
    P_metro = 0.0;
    P_sw = 0.0;

}

void Configuration::initialMeas()
{
    //--- k-space Measurement
    k_vec.resize(Dim);
    for (int i = 0; i < Dim; i++) 
        k_vec[i] = 0;
    
    k_vec[0] = 2 * M_PI / L;    // kx = 2 Pi / L

    k_cos.resize(Vol);
    k_sin.resize(Vol);

    std::vector<double> coor;
    coor.resize(Dim);
    for (int i = 0; i < Vol; i++)
    {
        for (int j = 0; j < Dim; j++)
            coor[j] = Latt.getComponent(i, j);

        k_cos[i] = cos(k_vec * coor);
        k_sin[i] = sin(k_vec * coor);
    }
    coor.clear();

    //--- Geometric Measurement
    x_max.resize(Dim);
    x_min.resize(Dim);
    x_now.resize(Dim);

    for (int i = 0; i < Dim; i++)
    {
        x_max[i] = 0;
        x_min[i] = 0;
        x_now[i] = 0;
    }

}

void Configuration::initialObsr()
{
    //--- Initialize Cluster Observables
    NCluster = 0;
    C1 = 0;
    C2 = 0;
    S2 = 0;
    S4 = 0;

    // para.Corr_Fun.resize(L);
}

#include "../Measurement_config.hpp"

bool Configuration::measureOrNot()
{
    return true;
}

void Configuration::writeCnf()
{
    
}

void Configuration::checkCnf()
{

}
