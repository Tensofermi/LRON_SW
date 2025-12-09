#include "Configuration.hpp"

void Configuration::updateCnf()
{
//--- Short-Range Algorithm // you need to keep input sigma > 10000
    // Metropolis();
    // Wolff();
    // SwendsenWang();

//--- Long-Range Algorithm
    // Metropolis_LR();
    // Wolff_LR();
    // SwendsenWang_LR();

    // Metropolis_LR_Clock_Thinning();
    // Wolff_LR_Clock_Thinning();
    SwendsenWang_LR_Clock_Thinning();
}
