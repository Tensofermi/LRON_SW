#pragma once
#include <bits/stdc++.h>

// This version is only applicable to 2D case.
class BondTable
{
public:
    unsigned long Vol;                 // Volume for lattice
    unsigned long H_index;             // hash index (event index)
    std::vector<unsigned long> V2H;    // Input site, tell you its hash index
    std::vector<unsigned long> H2V;    // Input hash index, tell you the ending site for bond
    std::vector<unsigned long> Lnk;    // Input hash index, tell you the previous hash index for this site
    std::vector<int> Dx, Dy;           // Displacement in the x and y directions for bonds

    void ini(unsigned long _Vol)
    {
        Vol = _Vol;
        unsigned long size = 10 * Vol;

        V2H.resize(size);  
        H2V.resize(size);   
        Lnk.resize(size);  

        Dx.resize(size);   
        Dy.resize(size);   
    }

    void reset()
    {
        H_index = 0;

        std::fill(V2H.begin(), V2H.end(), 0); 
        std::fill(H2V.begin(), H2V.end(), 0); 
        std::fill(Lnk.begin(), Lnk.end(), 0); 

        std::fill(Dx.begin(), Dx.end(), 0); 
        std::fill(Dy.begin(), Dy.end(), 0); 
    }

    void addBond(long _i, long _j, long _x, long _y)
    {
        /* This function adds a bond from site _i to site _j.
           The relative coordinate shift is represented by _x and _y
           for the x and y directions, respectively.
           
           Bond representation:
           _i ========> _j
           start       end
        */

        H_index++;                // Update the hash index (event index)
        Lnk[H_index] = V2H[_i];   // Save the previous hash index for starting site _i
        V2H[_i] = H_index;        // Update the hash index for starting site _i
        H2V[H_index] = _j;        // Record the ending site _j for the current hash index (event index)

        // Store the displacement for x and y directions
        Dx[H_index] = _x;
        Dy[H_index] = _y;
    }

    void addDouble(long _i, long _j, long _x, long _y)
    {
        addBond(_i, _j, _x, _y);        // add Bond: _i ========> _j
        addBond(_j, _i, - _x, - _y);    // add Bond: _j ========> _i 
    }

    // void listBond();                // output all the bonds
};

// void BondTable::listBond()
// {
//     std::string str = "";
//     str += "/---\n";
//     str += l_jf("# hash_index", 20) + l_jf("# site_index", 20) + l_jf("# link_hash", 20) + l_jf("# shift_x", 20) + l_jf("# shift_y", 20) + "\n";
//     for (unsigned long i_hash = 0; i_hash <= H_index; i_hash++)
//     {
//         str += l_jf(toStr(i_hash), 20) + l_jf(toStr(H2V[i_hash]), 20) + l_jf(toStr(Lnk[i_hash]), 20) + l_jf(toStr(Dx[i_hash]), 20) + l_jf(toStr(Dy[i_hash]), 20) + "\n";
//     }


//     str += "/---\n";
//     str += l_jf("# site_index", 20) + l_jf("# hash_index", 20) + "\n";
//     for (unsigned long i_site = 0; i_site < Vol; i_site++)
//     {
//         str += l_jf(toStr(i_site), 20) + l_jf(toStr(V2H[i_site]), 20) + "\n";
//     }


//     wrt("listBond.txt", str);
//     // hisWrt("listBond.txt", str);

// }