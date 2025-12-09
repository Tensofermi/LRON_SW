#include "../Configuration.hpp"

/* ------------------------------------------------------------------------------------*/
/* ------------------------------ Short-Range Interaction -----------------------------*/
/* ------------------------------------------------------------------------------------*/

void Configuration::Wolff()
{
    long i_Site, j_Site;
    double project;
    int rear, front;

    //--- Clear Queue and Memory
    for (unsigned int i = 0; i < Vol; i++)
    {
        Que[i] = 0;
        Mem[i] = 0;
    }

    //--- Randomly choose an axis to project
    randSpin();
    
    //--- BFS Search
    i_Site = rn.getRandomNum(Vol);
    Mem[i_Site] = 1;

    rear = 0; 
    front = 1;
    Que[front] = i_Site;    // start from 1 !

    while(front > rear)
    {
        rear++;
        i_Site = Que[rear];
        for (unsigned int i = 0; i < NNb; i++)
        {
            j_Site = Latt.getNNSite(i_Site, i);
            project = (tempSpin * Site[i_Site]) * (tempSpin * Site[j_Site]);
            P_sw = 1.0 - exp(-2 * Beta * project);
            if(Mem[j_Site] == 0 && rn.getRandomDouble() < P_sw)
            {
                front++; 
                Mem[j_Site] = 1;
                Que[front] = j_Site; 
            }
        }

    }

    //--- Flip Cluter Spins
    for (unsigned i = 1; i <= front; i++)    // start from 1 !
        flipSpin(Que[i], tempSpin);


}


/* ------------------------------------------------------------------------------------*/
/* ------------------------------ Long-Range Interaction ------------------------------*/
/* ------------------------------------------------------------------------------------*/

void Configuration::Wolff_LR()
{
    long i_Site, j_Site;
    std::vector<int> shift_, j_;
    double project;
    int rear, front;


    shift_.resize(Dim);
    j_.resize(Dim);

    //--- Clear Queue and Memory
    for (unsigned int i = 0; i < Vol; i++)
    {
        Que[i] = 0;
        Mem[i] = 0;
    }

    //--- Randomly choose an axis to project
    randSpin();
    
    //--- BFS Search
    i_Site = rn.getRandomNum(Vol);
    Mem[i_Site] = 1;

    rear = 0; 
    front = 1;
    Que[front] = i_Site;    // start from 1 !

    while(front > rear)
    {
        rear++;
        i_Site = Que[rear];
        for (unsigned int i = 0; i < Nlist; i++)
        {
            // obtain the interaction j_Site
            for (int k = 0; k < Dim; k++)
            {
                shift_[k] = range_list[i][k];
                j_[k] = Latt.getComponent(i_Site, k) + shift_[k];
                j_[k] = pbc_mod(j_[k], L);
            }
            j_Site = Latt.getSite(j_);

            project = (tempSpin * Site[i_Site]) * (tempSpin * Site[j_Site]);
            P_sw = 1.0 - exp(-2 * Beta * strength_list[i] * project);
            if(Mem[j_Site] == 0 && rn.getRandomDouble() < P_sw)
            {
                front++; 
                Mem[j_Site] = 1;
                Que[front] = j_Site; 
            }
        }

    }

    //--- Flip Cluter Spins
    for (unsigned i = 1; i <= front; i++)    // start from 1 !
        flipSpin(Que[i], tempSpin);

    
}


/* ------------------------------------------------------------------------------------*/
/* ----------- Long-Range Interaction with Clock and Dynamic Thinning Method ----------*/
/* ------------------------------------------------------------------------------------*/

void Configuration::Wolff_LR_Clock_Thinning()
{
    long i_Site, j_Site;
    std::vector<int> shift_, j_;
    double project;
    int rear, front;


    shift_.resize(Dim);
    j_.resize(Dim);

    //--- Clear Queue and Memory
    for (unsigned int i = 0; i < Vol; i++)
    {
        Que[i] = 0;
        Mem[i] = 0;
    }

    //--- Randomly choose an axis to project
    randSpin();
    
    //--- BFS Search
    i_Site = rn.getRandomNum(Vol);
    Mem[i_Site] = 1;

    rear = 0; 
    front = 1;
    Que[front] = i_Site;    // start from 1 !

    while(front > rear)
    {
        rear++;
        i_Site = Que[rear];

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

            // 2. check whether to make a bond between i_Site and the chosen site.
            if (clock_index <= n_box)
            {
                double hat = clock_rho[clock_index - 1];

                if (rn.getRandomDouble() < hat / P_mb)
                {
                    long j_index = clock_order[clock_index - 1];
                    // obtain the interaction j_Site
                    for (int k = 0; k < Dim; k++)
                    {
                        shift_[k] = range_list[j_index][k];
                        j_[k] = Latt.getComponent(i_Site, k) + shift_[k];
                        j_[k] = pbc_mod(j_[k], L);
                    }
                    j_Site = Latt.getSite(j_);
                    
                    project = (tempSpin * Site[i_Site]) * (tempSpin * Site[j_Site]);
                    P_sw = 1.0 - exp(-2 * Beta * strength_list[j_index] * project);
                    P_sw = P_sw / hat;
                    if(Mem[j_Site] == 0 && rn.getRandomDouble() < P_sw)
                    {
                        front++; 
                        Mem[j_Site] = 1;
                        Que[front] = j_Site; 
                    }
                }

            }
            
            // 3. judge whether exit
            if (clock_index >= n_box)
                break;

        }
        // ----- Clock Method ----- //

    }

    //--- Flip Cluter Spins
    for (unsigned i = 1; i <= front; i++)    // start from 1 !
        flipSpin(Que[i], tempSpin);

    

}