#include "../Configuration.hpp"

/* ------------------------------------------------------------------------------------*/
/* ------------------------------ Short-Range Interaction -----------------------------*/
/* ------------------------------------------------------------------------------------*/

void Configuration::SwendsenWang()
{
    long i_Site, j_Site, k_Site;
    double project;
    int rear, front;

    //--- Clear Observables
    NCluster = 0;
    C1 = 0;
    C2 = 0;
    S2 = 0;
    S4 = 0;

    //--- Clear Queue and Memory
    for (unsigned int i = 0; i < Vol; i++)
    {
        Que[i] = 0;
        Mem[i] = 0;
    }

    //--- Randomly choose an axis to project
    randSpin();

    //--- BFS Search
    for (long i_Site = 0; i_Site < Vol; i_Site++)
    {
        //--- Check Whether Visited
        if(Mem[i_Site] == 1) continue;
        Mem[i_Site] = 1;

        NCluster++;

        rear = 0;
        front = 1;
        Que[front] = i_Site;     // start from 1 !

        while(front > rear)
        {
            rear++;
            j_Site = Que[rear];
            for (unsigned int j = 0; j < NNb; j++)
            {
                // obtain the interacting k_Site
                k_Site = Latt.getNNSite(j_Site, j);
                project = (tempSpin * Site[j_Site]) * (tempSpin * Site[k_Site]);
                P_sw = 1.0 - exp(-2 * Beta * project);
                if(Mem[k_Site] == 0 && rn.getRandomDouble() < P_sw)
                {
                    front++;
                    Mem[k_Site] = 1;
                    Que[front] = k_Site;
                }
            }
        }

        //--- Cluster Measurement
        S2 = S2 + pow(front, 2);
        S4 = S4 + pow(front, 4);

        if(front > C1)
        {
            C2 = C1;
            C1 = front;
        }
        else if(front > C2)
        {
            C2 = front;
        }

        //--- Flip Cluter Spins
        if(rn.getRandomDouble() < 0.5)
        {
            for (unsigned j = 1; j <= front; j++)    // start from 1 !
                flipSpin(Que[j], tempSpin);
        }

    }

}


/* ------------------------------------------------------------------------------------*/
/* ------------------------------ Long-Range Interaction ------------------------------*/
/* ------------------------------------------------------------------------------------*/

void Configuration::SwendsenWang_LR()
{
    long i_Site, j_Site, k_Site;
    std::vector<int> shift_, k_;
    double project;
    int rear, front;

    //--- Clear Observables
    NCluster = 0;
    C1 = 0;
    C2 = 0;
    S2 = 0;
    S4 = 0;

    shift_.resize(Dim);
    k_.resize(Dim);

    //--- Clear Queue and Memory
    for (unsigned int i = 0; i < Vol; i++)
    {
        Que[i] = 0;
        Mem[i] = 0;
    }

    //--- Randomly choose an axis to project
    randSpin();

    //--- BFS Search
    for (long i_Site = 0; i_Site < Vol; i_Site++)
    {
        //--- Check Whether Visited
        if(Mem[i_Site] == 1) continue;
        Mem[i_Site] = 1;

        NCluster++;

        rear = 0;
        front = 1;
        Que[front] = i_Site;     // start from 1 !

        while(front > rear)
        {
            rear++;
            j_Site = Que[rear];
            for (unsigned int j = 0; j < Nlist; j++)
            {
                // obtain the interacting k_Site
                for (int k = 0; k < Dim; k++)
                {
                    shift_[k] = range_list[j][k];
                    k_[k] = Latt.getComponent(j_Site, k) + shift_[k];
                    k_[k] = pbc_mod(k_[k], L);
                }
                k_Site = Latt.getSite(k_);

                project = (tempSpin * Site[j_Site]) * (tempSpin * Site[k_Site]);
                P_sw = 1.0 - exp(-2 * Beta * strength_list[j] * project);
                if(Mem[k_Site] == 0 && rn.getRandomDouble() < P_sw)
                {
                    front++;
                    Mem[k_Site] = 1;
                    Que[front] = k_Site;
                }
            }
        }

        //--- Cluster Measurement
        S2 = S2 + pow(front, 2);
        S4 = S4 + pow(front, 4);

        if(front > C1)
        {
            C2 = C1;
            C1 = front;
        }
        else if(front > C2)
        {
            C2 = front;
        }

        //--- Flip Cluter Spins
        if(rn.getRandomDouble() < 0.5)
        {
            for (unsigned j = 1; j <= front; j++)    // start from 1 !
                flipSpin(Que[j], tempSpin);
        }

    }

}


/* ------------------------------------------------------------------------------------*/
/* ----------- Long-Range Interaction with Clock and Dynamic Thinning Method ----------*/
/* ------------------------------------------------------------------------------------*/

void Configuration::SwendsenWang_LR_Clock_Thinning()
{
    long i_Site, j_Site, k_Site;
    std::vector<int> shift_, k_;
    double project;
    int rear, front;

    //--- Clear Observables
    NCluster = 0;
    C1 = 0;
    C2 = 0;
    S2 = 0;
    S4 = 0;
    shift_.resize(Dim);
    k_.resize(Dim);

    //--- Clear Queue and Memory
    for (unsigned int i = 0; i < Vol; i++)
    {
        Que[i] = 0;
        Mem[i] = 0;
    }

    //--- Randomly choose an axis to project
    randSpin();

    //--- BFS Search
    for (long i_Site = 0; i_Site < Vol; i_Site++)
    {
        //--- Check Whether Visited
        if(Mem[i_Site] == 1) continue;
        Mem[i_Site] = 1;

        NCluster++;

        rear = 0;
        front = 1;
        Que[front] = i_Site;     // start from 1 !

        while(front > rear)
        {
            rear++;
            j_Site = Que[rear];

            // ----- Clock Method ----- //
            int n_box = Nlist;
            int clock_index = 0;
            while (true)
            {
                // 1. calculate number of skip steps
                double P_mb = clock_rho[clock_index];
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

                // 2. check whether to make a bond between j_Site and the chosen site.
                if (clock_index <= n_box)
                {
                    double hat = clock_rho[clock_index - 1];

                    if (rn.getRandomDouble() < hat / P_mb)
                    {
                        long k_index = clock_order[clock_index - 1];
                        // obtain the interacting k_Site
                        for (int k = 0; k < Dim; k++)
                        {
                            shift_[k] = range_list[k_index][k];
                            k_[k] = Latt.getComponent(j_Site, k) + shift_[k];
                            k_[k] = pbc_mod(k_[k], L);
                        }
                        k_Site = Latt.getSite(k_);
                        
                        project = (tempSpin * Site[j_Site]) * (tempSpin * Site[k_Site]);
                        P_sw = 1.0 - exp(-2 * Beta * strength_list[k_index] * project);
                        P_sw = P_sw / hat;
                        if(Mem[k_Site] == 0 && rn.getRandomDouble() < P_sw)
                        {
                            front++; 
                            Mem[k_Site] = 1;
                            Que[front] = k_Site; 
                        }
                    }

                }

                // 3. judge whether exit
                if (clock_index >= n_box)
                    break;

            }
            // ----- Clock Method ----- //

        }

        //--- Cluster Measurement
        S2 = S2 + pow(front, 2);
        S4 = S4 + pow(front, 4);

        if(front > C1)
        {
            C2 = C1;
            C1 = front;
        }
        else if(front > C2)
        {
            C2 = front;
        }

        //--- Flip Cluter Spins
        if(rn.getRandomDouble() < 0.5)
        {
            for (unsigned j = 1; j <= front; j++)    // start from 1 !
                flipSpin(Que[j], tempSpin);
        }

    }


}