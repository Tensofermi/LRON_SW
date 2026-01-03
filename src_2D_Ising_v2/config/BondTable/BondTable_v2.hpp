#pragma once
#include <bits/stdc++.h>

// This version is only applicable to 2D case.
class BondTable_v2
{
public:
    unsigned long Vol;         // Volume for lattice
    unsigned long H_index;     // Hash index (event index)

    unsigned long* V2H;        // Input site, tell you its hash index
    unsigned long* H2V;        // Input hash index, tell you the ending site for bond
    unsigned long* Lnk;        // Input hash index, tell you the previous hash index for this site
    int* Dx;                   // Displacement in the x direction for bonds
    int* Dy;                   // Displacement in the y direction for bonds

    unsigned long capacity;    // Capacity to store bonds

    BondTable_v2() : V2H(nullptr), H2V(nullptr), Lnk(nullptr), Dx(nullptr), Dy(nullptr), capacity(0), H_index(0) {}

    // Destructor to free allocated memory
    ~BondTable_v2() {
        delete[] V2H;
        delete[] H2V;
        delete[] Lnk;
        delete[] Dx;
        delete[] Dy;
    }

    // Initialize the BondTable
    void ini(unsigned long _Vol)
    {
        Vol = _Vol;
        capacity = Vol;   // Start with initial capacity equal to Vol

        // Allocate memory for all arrays
        V2H = new unsigned long[Vol]();  // Initialize to 0
        H2V = new unsigned long[Vol]();
        Lnk = new unsigned long[Vol]();
        Dx = new int[Vol]();
        Dy = new int[Vol]();

        H_index = 0;
    }

    // Reset the BondTable
    void reset()
    {
        H_index = 0;

        // Reset all the arrays to zero
        std::memset(V2H, 0, capacity * sizeof(unsigned long));
        std::memset(H2V, 0, capacity * sizeof(unsigned long));
        std::memset(Lnk, 0, capacity * sizeof(unsigned long));
        std::memset(Dx, 0, capacity * sizeof(int));
        std::memset(Dy, 0, capacity * sizeof(int));
    }

    // Dynamically grow the capacity of the arrays when needed
    void grow()
    {
        unsigned long new_capacity = capacity * 2;

        // Allocate new arrays with double the capacity
        unsigned long* new_V2H = new unsigned long[new_capacity]();
        unsigned long* new_H2V = new unsigned long[new_capacity]();
        unsigned long* new_Lnk = new unsigned long[new_capacity]();
        int* new_Dx = new int[new_capacity]();
        int* new_Dy = new int[new_capacity]();

        // Copy existing data to new arrays
        std::memcpy(new_V2H, V2H, capacity * sizeof(unsigned long));
        std::memcpy(new_H2V, H2V, capacity * sizeof(unsigned long));
        std::memcpy(new_Lnk, Lnk, capacity * sizeof(unsigned long));
        std::memcpy(new_Dx, Dx, capacity * sizeof(int));
        std::memcpy(new_Dy, Dy, capacity * sizeof(int));

        // Free old arrays
        delete[] V2H;
        delete[] H2V;
        delete[] Lnk;
        delete[] Dx;
        delete[] Dy;

        // Point to new arrays
        V2H = new_V2H;
        H2V = new_H2V;
        Lnk = new_Lnk;
        Dx = new_Dx;
        Dy = new_Dy;

        // Update capacity
        capacity = new_capacity;

        // std::cout << "capacity update to: " << capacity << std::endl;

    }

    // Add a bond from site _i to site _j
    void addBond(long _i, long _j, long _x, long _y)
    {
        if (H_index >= capacity - 1) {
            grow();  // Grow the arrays when capacity is exceeded
        }

        /* This function adds a bond from site _i to site _j.
           The relative coordinate shift is represented by _x and _y
           for the x and y directions, respectively.
           
           Bond representation:
           _i ========> _j
           start       end
        */

        H_index++;                   // Increment the hash index
        Lnk[H_index] = V2H[_i];      // Store the previous bond's hash index for site _i
        V2H[_i] = H_index;           // Update the current bond's hash index for site _i
        H2V[H_index] = _j;           // Store the ending site _j for the current bond

        // Store displacement in x and y directions
        Dx[H_index] = _x;
        Dy[H_index] = _y;
    }

    // Add a double bond between two sites
    void addDouble(long _i, long _j, long _x, long _y)
    {
        addBond(_i, _j, _x, _y);       // Add bond: _i ========> _j
        addBond(_j, _i, -_x, -_y);     // Add bond: _j ========> _i
    }

    // Output the bonds (implementation placeholder)
    // void listBond();
};

// void BondTable_v2::listBond()
// {
//     std::string str = "";
//     str += "/---\n";
//     str += "capacity = " + toStr(capacity) + "\n";
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