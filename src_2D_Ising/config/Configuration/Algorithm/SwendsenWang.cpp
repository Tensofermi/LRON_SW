#include "../Configuration.hpp"

void Configuration::SwendsenWang_LR_Clock_Bond()
{
    long i_Site, j_Site, k_Site;
    // std::vector<int> shift_, k_;
    double project;
    int rear, front;

    bondOrNot = true;  // flag for measurement

    //--- Clear Observables
    Bond.reset();      // reset the Bond configuration
    cluster.reset();

    rn.getRandomNum(2);

    //--- Generate clusters
    for (long i_Bond = 0; i_Bond < Nlist; i_Bond++)
    {
        int k_index = clock_order[i_Bond];
        double P_mb = clock_rho[i_Bond];

        // consider bond vectors pointing in the upper plane (dy >= 0) and with negative axis (dy = 0 and dx <= 0) removed
        if (range_list[k_index][1] < 0 || ((range_list[k_index][1] == 0 || range_list[k_index][1] == int((L + 1) / 2)) && range_list[k_index][0] < 0))
        {
            continue;
        }

        int n_box = Vol;
        int clock_index = 0;
        while (true)
        {
            // 1. calculate number of skip steps
            if (P_mb > 0)
            {   
                double n_skip = log(1.0 - rn.getRandomDouble()) / log(1.0 - P_mb);
                if (n_skip < n_box)
                    clock_index += floor(1 + n_skip);
                else
                    break;
            }
            else
                break;

            // 2. check whether to make a bond
            if (clock_index <= n_box)
            {
                j_Site = clock_index - 1;
                
                // avoiding double counting
                if(range_list[k_index][1] == 0 && range_list[k_index][0] == int((L+1)/2) && Latt.getComponent(j_Site, 0) + 1 > int((L + 1)/2))
                {
                    continue;
                }
                if(range_list[k_index][0] == 0 && range_list[k_index][1] == int((L+1)/2) && Latt.getComponent(j_Site, 1) + 1 > int((L + 1)/2))
                {
                    continue;
                }
                if(range_list[k_index][0] == int((L+1)/2) && range_list[k_index][1] == int((L+1)/2) && Latt.getComponent(j_Site, 0) + 1 > int((L + 1)/2))
                {
                    continue;
                } 

                // obtain k_Site
                shift_[0] = range_list[k_index][0];
                k_[0] = Latt.getComponent(j_Site, 0) + shift_[0];
                k_[0] = pbc_mod(k_[0], L);

                shift_[1] = range_list[k_index][1];
                k_[1] = Latt.getComponent(j_Site, 1) + shift_[1];
                k_[1] = pbc_mod(k_[1], L);

                k_Site = Latt.getSite(k_);

                project = Site[j_Site] * Site[k_Site];
                P_sw = 1.0 - exp(-2 * Beta * strength_list[k_index] * project);
                P_sw = P_sw / P_mb;
                
                if (rn.getRandomDouble() < P_sw)
                {
                    // unite clusters
                    cluster.unite(j_Site, k_Site);
                    // add double Bonds
                    Bond.addDouble(j_Site, k_Site, shift_[0], shift_[1]);
                }

            }

            // 3. judge whether exit
            if (clock_index >= n_box)
                break;
        }

    }

    // -1 means not examined yet; 0 means not flip; 1 means flip
    std::fill(cluster_flip_tag.begin(), cluster_flip_tag.end(), -1); 
    //--- Update spins
    for (long i_Site = 0; i_Site < Vol; i_Site++)
    {
        // find the root of the cluster
        int root = cluster.find(i_Site);

        // if cluster root is not examined yet, assign 0 or 1 to the root
        if (cluster_flip_tag[root] == -1)
            cluster_flip_tag[root] = rn.getRandomDouble() < 0.5 ? 0 : 1;

        // flip the site based on the root
        if (cluster_flip_tag[root] == 1)
        {
            // flipSpin(i_Site, tempSpin);
            Site[i_Site] = - Site[i_Site];
        }

    }

}