#include "../Configuration.hpp"

/* ------------------------------------------------------------------------------------*/
/* ------------------------------ Short-Range Interaction -----------------------------*/
/* ------------------------------------------------------------------------------------*/

void Configuration::Metropolis()
{
    long i_Site, j_Site;
    double dE;
    std::vector<double> dV;

    for (unsigned int i = 0; i < Vol; i++)
    {
        i_Site = rn.getRandomNum(Vol);
        randSpin();
        
        dV = tempSpin - Site[i_Site];

        //--- Calculate Delta Energy
        dE = 0.0;
        for (unsigned int j = 0; j < NNb; j++)
        {
            j_Site = Latt.getNNSite(i_Site, j);
            dE = dE - dV * Site[j_Site];
        }

        //--- Attempt to Update
        P_metro = exp(- Beta * dE);
        if(rn.getRandomDouble() < P_metro)
            Site[i_Site] = tempSpin;
    }

}


/* ------------------------------------------------------------------------------------*/
/* ------------------------------ Long-Range Interaction ------------------------------*/
/* ------------------------------------------------------------------------------------*/

void Configuration::Metropolis_LR()
{
    long i_Site, j_Site;
    std::vector<int> shift_, j_;
    std::vector<double> dV;
    double dE;

    shift_.resize(Dim);
    j_.resize(Dim);

    for (unsigned int i = 0; i < Vol; i++)
    {
        i_Site = rn.getRandomNum(Vol);
        randSpin();
        
        dV = tempSpin - Site[i_Site];

        //--- Calculate Delta Energy
        dE = 0.0;
        for (unsigned int j = 0; j < Nlist; j++)  // cause O(N) computational cost
        {
            // obtain the interaction j_Site
            for (int k = 0; k < Dim; k++)
            {
                shift_[k] = range_list[j][k];
                j_[k] = Latt.getComponent(i_Site, k) + shift_[k];
                j_[k] = pbc_mod(j_[k], L);
            }
            j_Site = Latt.getSite(j_);
            
            dE += - strength_list[j] * dV * Site[j_Site];
        }

        //--- Attempt to Update
        P_metro = exp(- Beta * dE);
        if(rn.getRandomDouble() < P_metro)
            Site[i_Site] = tempSpin;
    }

}


/* ------------------------------------------------------------------------------------*/
/* ----------- Long-Range Interaction with Clock and Dynamic Thinning Method ----------*/
/* ------------------------------------------------------------------------------------*/

void Configuration::Metropolis_LR_Clock_Thinning()
{
    long i_Site, j_Site;
    std::vector<int> shift_, j_;
    std::vector<double> dV;
    double dE;

    shift_.resize(Dim);
    j_.resize(Dim);

    for (unsigned int i = 0; i < Vol; i++)
    {
        i_Site = rn.getRandomNum(Vol);
        randSpin();
        
        dV = tempSpin - Site[i_Site];

        // ----- Clock Method ----- //
        int n_box = Nlist;
        int clock_index = 0;
        while (true)
        {
            // 1. calculate number of skip steps
            double P_mb = clock_rho[clock_index];   // here, P_mb means 1 - P_metro_filter
            if (P_mb > 0)
            {
                double n_skip = log(1.0 - rn.getRandomDouble()) / log(1.0 - P_mb);
                if (n_skip < n_box)
                    clock_index += floor(1 + n_skip);
                else
                {
                    Site[i_Site] = tempSpin;    // update
                    break;
                }
            }
            else
            {
                Site[i_Site] = tempSpin;    // update
                break;
            }

            // 2. check whether to reject.
            if (clock_index <= n_box)
            {
                double hat = clock_rho[clock_index - 1];

                if (rn.getRandomDouble() < hat / P_mb)
                {
                    long j_index = clock_order[clock_index - 1];
                    // obtain the interacting j_Site
                    for (int k = 0; k < Dim; k++)
                    {
                        shift_[k] = range_list[j_index][k];
                        j_[k] = Latt.getComponent(i_Site, k) + shift_[k];
                        j_[k] = pbc_mod(j_[k], L);
                    }
                    j_Site = Latt.getSite(j_);
                    dE = - strength_list[j_index] * dV * Site[j_Site];
                    P_sw = 1.0 - exp(- Beta * dE);
                    P_sw = P_sw / hat;
                    if(rn.getRandomDouble() < P_sw)
                        break;  // reject

                }
            }

            // 3. judge whether exit
            if (clock_index >= n_box)
            {
                Site[i_Site] = tempSpin;    // update
                break;
            }

        }
        // ----- Clock Method ----- //
    }
    
}