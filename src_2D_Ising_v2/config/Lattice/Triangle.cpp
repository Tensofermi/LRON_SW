#include "Triangle.hpp"

// Initialization
void Triangle::set(std::initializer_list<int> initList)
{
    L.assign(initList); // 使用assign将初始化列表赋给成员变量L
    Vol = 1;
    for (int i = 0; i < Dim; ++i)
        Vol *= L[i];

    x_ = new int[Dim];
}

long Triangle::getVol()
{
    return Vol;
}

long Triangle::getNBond()
{
    return 3 * Vol;
}

int Triangle::getNNb()
{
    return NNb;
}

int Triangle::getDim()
{
    return Dim;
}

int Triangle::getOpsDir(int _Dir)
{
    return NNb - _Dir - 1;
}

//-----------------------------------------------------------
//-----------------------------------------------------------

int Triangle::getComponent(long _Site, int _Dir)      // Site --> (y, x)
{
    for (int i = 1; i > _Dir; i--) _Site = _Site / L[i];
    return _Site % L[_Dir];
}

long Triangle::getSite(int *_Component)               // (y, x) --> Site
{
    long Vc;
    Vc = _Component[0];
    for (int i = 1; i < Dim; i++)
    {
        Vc = _Component[i] + Vc * L[Dim-i];
    }
    return Vc;
}

long Triangle::getSite(std::vector<int> &_Component)  // (y, x) --> Site
{
    long Vc;
    Vc = _Component[0];
    for (int i = 1; i < Dim; i++)
    {
        Vc = _Component[i] + Vc * L[Dim-i];
    }
    return Vc;
}

//-----------------------------------------------------------
//-----------------------------------------------------------

long Triangle::getNNSite(long _Site, int _Dir)
{
    for (int i = 0; i < Dim; i++) x_[i] = getComponent(_Site, i);

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


    if(_Dir == 0) // y
    {
        x_[0]++;    // y++
        if(x_[0]>=L[0]) x_[0] -= L[0];

    }
    else if(_Dir == 1) // x
    {
        x_[1]++;    // x++
        if(x_[1]>=L[1]) x_[1] -= L[1];

    }
    else if(_Dir == 2) // effective x
    {
        x_[0]--;    // y--
        if(x_[0]<0)     x_[0] += L[0];
        x_[1]++;    // x++
        if(x_[1]>=L[1]) x_[1] -= L[1];
    
    }
    else if(_Dir == 3) // effective x
    {
        x_[0]++;    // y++
        if(x_[0]>=L[0]) x_[0] -= L[0];
        x_[1]--;    // x--
        if(x_[1]<0)     x_[1] += L[1];
    
    }
    else if(_Dir == 4) // x
    {
        x_[1]--;    // x--
        if(x_[1]<0)     x_[1] += L[1];
    
    }
    else if(_Dir == 5) // y
    {
        x_[0]--;    // y--
        if(x_[0]<0)     x_[0] += L[0];

    }

    return getSite(x_);
}

long Triangle::getNNBond(long _Site, int _Dir)
{
    if(_Dir < 3) {return _Site + _Dir * Vol;}
    else {return getNNSite(_Site, _Dir) + getOpsDir(_Dir) * Vol;}
}

int Triangle::getDir(long _Site, long _NNSite)
{
    for (int Dir = 0; Dir < NNb; Dir++)
    {
        if(getNNSite(_Site, Dir)==_NNSite) return Dir;
    }
    return -1;
}

int Triangle::getParity(long _Site)
{
	long x = 0;
	for (int i = 0; i < Dim; i++)
	{
		x += getComponent(_Site, i);
	}
	if(x % 2) return -1;
	return 1;
}

long Triangle::getDirSite(long _Site, int _Dir, int _Length)
{
    for (int i = 0; i < Dim; i++) x_[i] = getComponent(_Site, i);
    
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


    if(_Dir == 0) // y
    {
        x_[0] += _Length;    // y++
        if(x_[0]>=L[0]) x_[0] -= L[0];

    }
    else if(_Dir == 1) // x
    {
        x_[1] += _Length;    // x++
        if(x_[1]>=L[1]) x_[1] -= L[1];

    }
    else if(_Dir == 2) // effective x
    {
        x_[0] -= _Length;    // y--
        if(x_[0]<0)     x_[0] += L[0];
        x_[1] += _Length;    // x++
        if(x_[1]>=L[1]) x_[1] -= L[1];
    
    }
    else if(_Dir == 3) // effective x
    {
        x_[0] += _Length;    // y++
        if(x_[0]>=L[0]) x_[0] -= L[0];
        x_[1] -= _Length;    // x--
        if(x_[1]<0)     x_[1] += L[1];
    
    }
    else if(_Dir == 4) // x
    {
        x_[1] -= _Length;    // x--
        if(x_[1]<0)     x_[1] += L[1];
    
    }
    else if(_Dir == 5) // y
    {
        x_[0] -= _Length;    // y--
        if(x_[0]<0)     x_[0] += L[0];

    }
    return getSite(x_);

}
