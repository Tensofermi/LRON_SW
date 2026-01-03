#pragma once
#include <bits/stdc++.h>

class Triangle
{
    //                        [x]
    //             (2)        (1)
    //              o ______ o ______ o             
    //             /  \     /  \     /                  
    //            /    \   /    \   /      
    //           /      \ /      \ /        
    //  [-y](5) o ______⭐______ o (0) [y]     
    //        /  \      / \      / 
    //       /    \    /   \    /       
    //      /      \  /     \  /       
    //     o ______ o ______ o   
    //             (4)       (3)
    //            [-x]

protected:
    const int Dim = 2;
    const int NNb = 6;

    std::vector<int> L; // 不再初始化默认值
    long Vol = 9;
    int* x_ = nullptr;

public:
    void set(std::initializer_list<int> initList);
    
    long getVol();
    long getNBond();
    int getNNb();
    int getDim();
    int getOpsDir(int _Dir);
    int getComponent(long _Site, int _Dir);
    long getSite(int *_Component);
    long getSite(std::vector<int> &_Component);
    int getDir(long _Site, long _NNSite);
    long getNNSite(long _Site, int _Dir);
    long getNNBond(long _Site, int _Dir);
    int getParity(long _Site);
    long getDirSite(long _Site, int _Dir, int _Length);
};
