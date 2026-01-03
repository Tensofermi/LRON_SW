#pragma once
#include <bits/stdc++.h>

class ClusterFind_Wrapping {

private:
    int size, L;   // set size
    std::vector<int> Dx, Dy;           // Displacement in the x and y directions for bonds
    std::vector<int> root; // Tracks the root of each element.
public:
    int wrp;
  // 0 -- no wrapping
  // 1 -- 1D wrapping (only x)
  // 2 -- 1D wrapping (only y)
  // 3 -- 1D spiraling wrapping
  // 4 -- 2D wrapping (x and y)


public:
    void ini(unsigned long _L)
    {
        L = _L;
        size = L * L;
        root.resize(size, -1);
        Dx.resize(size, 0);
        Dy.resize(size, 0);
        reset();
    }

    // Reset the cluster
    void reset(){
        std::fill(root.begin(), root.end(), -1);
        std::fill(Dx.begin(), Dx.end(), 0);
        std::fill(Dy.begin(), Dy.end(), 0);
        wrp = 0;
    }

    // Finds and returns the root of the set that element 'x' belongs to.
    // Implements path compression to make subsequent queries faster.
    int find(int x) {
        if (root[x] < 0) {
            return x;   // negtive value means root node
        }

        while (root[root[x]] >= 0) {
            Dx[x] += Dx[root[x]];
            Dy[x] += Dy[root[x]];
            x = root[x] = root[root[x]];
            if (root[x] < 0)
            {
                return x;
            }
        }

        return root[x]; 
    }
        
    std::vector<int> find_root_displacement(int x) {
        // Find the displacement of node x relative to its root in a 2D disjoint-set structure.
        // root[i] < 0 indicates that i is a root node.
        // Dx[i], Dy[i] represent the displacement from node i to its parent root[i].
        
        if (root[x] < 0) {
            return {0, 0}; // x is root
        }
        
        int dx = 0, dy = 0;
        while (root[root[x]] >= 0) {
            Dx[x] += Dx[root[x]];
            Dy[x] += Dy[root[x]];
            dx += Dx[x]; 
            dy += Dy[x];
            x = root[x] = root[root[x]];
            if (root[x] < 0) { // v is root
            return {dx, dy};
            }
        }
        return {dx + Dx[x], dy + Dy[x]}; // ds[v] is root
    }

    // Connects the sets containing elements 'x' and 'y' together.
    // Uses union by rank to choose the new root, which keeps the tree flat.
    void unite(int x, int y, int shift_x, int shift_y) {
        // shift: displacement from x to y
        // X----->Y
        //  shift

        int rootX = find(x);
        int rootY = find(y);
        
        if (rootX != rootY) { 
            
            std::vector<int> dispX = find_root_displacement(x);
            std::vector<int> dispY = find_root_displacement(y);

            if (root[rootX] > root[rootY]) // due to negtive value for root node, rootY is large
            {
                // X==>Y
                root[rootY] += root[rootX]; // rootX merge into rootY
                root[rootX] = rootY;

                if (Dx[rootX] != 0 || Dy[rootX] != 0)
                {
                    Dx[rootY] = Dx[rootX];
                    Dy[rootY] = Dy[rootX];
                }

                Dx[rootX] = dispY[0] - dispX[0] + shift_x;
                Dy[rootX] = dispY[1] - dispX[1] + shift_y;
            } 
            else  // otherwise, rootX is large or equal
            {
                // Y==>X
                root[rootX] += root[rootY]; // rootY merge into rootX
                root[rootY] = rootX;

                if (Dx[rootY] != 0 || Dy[rootY] != 0)
                {
                    Dx[rootX] = Dx[rootY];
                    Dy[rootX] = Dy[rootY];
                }

                Dx[rootY] = dispX[0] - dispY[0] - shift_x;
                Dy[rootY] = dispX[1] - dispY[1] - shift_y;
            }
        }
        else
        {
            // Already in the same cluster, check wrapping
            wrapping(x, y, shift_x, shift_y);
        }

    }

    // Returns true if elements 'x' and 'y' are part of the same set, false otherwise.
    bool clusterOrNot(int x, int y) {
        return find(x) == find(y);
    }

    // Return root value
    int getval(int x) { return root[x]; }

    int wrapping(int x, int y, int shift_x, int shift_y)
    {
        if (wrp == 4) return 4;
        // Checking wrapping cluster for the insertion of bond x-y and update wrp
        std::vector<int> dispX = find_root_displacement(x);
        std::vector<int> dispY = find_root_displacement(y);
        int total_dx = dispX[0] - dispY[0] - shift_x;
        int total_dy = dispX[1] - dispY[1] - shift_y;
        if(total_dx != 0 || total_dy != 0)
        {
            if (total_dx == 0)
            {
                wrp = 2;    // y wrapping
            }
            else if (total_dy == 0)
            {
                wrp = 1;    // x wrapping
            }
            else
            {
                wrp = 3;    // spiraling wrapping
            }

            int rootX = find(x);

            if (Dx[rootX] * total_dy - Dy[rootX] * total_dx != 0)
            {
                wrp = 4;    // 2D wrapping
            }
            
            Dx[rootX] = total_dx;
            Dy[rootX] = total_dy;
        }

        return wrp;

    }
};