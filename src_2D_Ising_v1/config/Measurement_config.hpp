# pragma once

// ------------------------------------------------------------------------
// This function passes the config data into the observable for collection.
// ------------------------------------------------------------------------
void Configuration::measure()
{
//--- Spin Observables
		tempSpin = 0;
	
	obs.Ob[para.i_E] = 0;
	for (unsigned int i = 0; i < Vol; i++)
	{
		tempSpin = tempSpin + Site[i];
		for (unsigned int j = 0; j < NNb / 2; j++)	// avoid double counting
		{
			obs.Ob[para.i_E] += strength_list[clock_order[j]] * (Site[i] * Site[Latt.getNNSite(i,j)]);
		}
	}
	obs.Ob[para.i_M] = abs(tempSpin);
	obs.Ob[para.i_absM] = abs(tempSpin);
	obs.Ob[para.i_M2] = obs.Ob[para.i_M] * obs.Ob[para.i_M];
	obs.Ob[para.i_M4] = obs.Ob[para.i_M2] * obs.Ob[para.i_M2];
	obs.Ob[para.i_E] = -obs.Ob[para.i_E];	// add minus
	obs.Ob[para.i_E2] = obs.Ob[para.i_E] * obs.Ob[para.i_E];

//--- k-space Observables
	// std::vector<double> mk_c, mk_s;
	// mk_c.resize(Nspin);
	// mk_s.resize(Nspin);
    double mk_c, mk_s;
	
	// for (int i = 0; i < Nspin; i++)
	// {
	// 	mk_c[i] = 0;
	// 	mk_s[i] = 0;
	// }

    mk_c = 0;
    mk_s = 0;
	
	for (int i = 0; i < Vol; i++)
	{
		mk_c = mk_c + k_cos[Latt.getComponent(i,1)] * Site[i];
		mk_s = mk_s + k_sin[Latt.getComponent(i,1)] * Site[i];
	}
	
	// double Mk2_temp = 0;
	// // for (int i = 0; i < Nspin; i++)
	// 	// Mk2_temp += mk_c[i] * mk_c[i] + mk_s[i] * mk_s[i];

    // Mk2_temp += mk_c * mk_c + mk_s * mk_s;

	obs.Ob[para.i_Mk2] = mk_c * mk_c + mk_s * mk_s;

	obs.Ob[para.i_M2_E] = obs.Ob[para.i_M2] * obs.Ob[para.i_E];
	obs.Ob[para.i_Mk2_E] = obs.Ob[para.i_Mk2] * obs.Ob[para.i_E];

//--- Cluster Observables
	// if (bondOrNot)
	// {
	// 	NCluster = 0; C1 = 0; C2 = 0;
	// 	S2 = 0; S4 = 0;
	// 	for (unsigned int i = 0; i < Vol; i++)
	// 	{
	// 		double label = cluster.getval(i);

	// 		if(label < 0)	// find root node
	// 		{
	// 			double cluster_size = -label;  // cluster size
				
	// 			NCluster++;  // number of clusters

	// 			// largest and second-largest cluster sizes
	// 			if (C1 < cluster_size)
	// 			{
	// 				C2 = C1;
	// 				C1 = cluster_size;
	// 			}
	// 			else if (C2 < cluster_size)
	// 			{
	// 				C2 = cluster_size;
	// 			}

	// 			S2 += pow(cluster_size, 2);
	// 			S4 += pow(cluster_size, 4);

	// 		}
	// 	}
	// }

	// obs.Ob[para.i_NCluster] = NCluster;
	// obs.Ob[para.i_S2] = S2;
	// obs.Ob[para.i_S4] = S4;
	// obs.Ob[para.i_SM4] = 3 * S2 * S2 - 2 * S4;
	// obs.Ob[para.i_C1] = C1;
	// obs.Ob[para.i_C2] = C2;

//--- Wrapping Probability
	if (bondOrNot)
	{
		BFS_Bond();
		// if(C1 != obs.Ob[para.i_C1]) io.exportInfo(io.TestiInfo, "C1_bond = " + toStr(C1) + "C1_site = " + toStr(obs.Ob[para.i_C1]) + "\n");
		// if(C2 != obs.Ob[para.i_C2]) io.exportInfo(io.TestiInfo, "C2_bond = " + toStr(C2) + "C2_site = " + toStr(obs.Ob[para.i_C2]) + "\n");
		// if(S2 != obs.Ob[para.i_S2]) io.exportInfo(io.TestiInfo, "S2_bond = " + toStr(S2) + "S2_site = " + toStr(obs.Ob[para.i_S2]) + "\n");
		// if(S4 != obs.Ob[para.i_S4]) io.exportInfo(io.TestiInfo, "S4_bond = " + toStr(S4) + "S4_site = " + toStr(obs.Ob[para.i_S4]) + "\n");
		// if(NCluster != obs.Ob[para.i_NCluster]) io.exportInfo(io.TestiInfo, "NCluster_bond = " + toStr(NCluster) + "NCluster_site = " + toStr(obs.Ob[para.i_NCluster]) + "\n");

        obs.Ob[para.i_NCluster] = NCluster;
        obs.Ob[para.i_S2] = S2;
        obs.Ob[para.i_S4] = S4;
        obs.Ob[para.i_SM4] = 3 * S2 * S2 - 2 * S4;
        obs.Ob[para.i_C1] = C1;
        obs.Ob[para.i_C2] = C2;

        obs.Ob[para.i_R0] = R0;
        obs.Ob[para.i_R1] = R1;
        obs.Ob[para.i_Rd] = Rd;
        obs.Ob[para.i_R2] = R2;

        obs.Ob[para.i_Rx] = Rx;

        obs.Ob[para.i_S1] = S1;

        obs.Ob[para.i_R_E] = (Rx) * obs.Ob[para.i_E];

        // Bond.listBond();
	}

}

void Configuration::BFS_Bond()
{

	long hash_index, j_Site, k_Site;
    int wrap_x, wrap_y, Nwrap, wrap_flag;
    int project;
    int rear, front;
    int depth_cluster;


    //--- Clear Observables
    NCluster = 0;
    C1 = 0;
    C2 = 0;
    S2 = 0;
    S4 = 0;

    S1 = 0;

    Nwrap = 0;

    //--- Clear Queue and Memory
    for (unsigned int i = 0; i < Vol; i++)
    {
        Que[i] = 0;
        Mem[i] = 0;
        Tx[i] = 0;
        Ty[i] = 0;
    }

    //--- BFS Search
    for (long i_Site = 0; i_Site < Vol; i_Site++)
    {
        //--- Check Whether Visited
        if(Mem[i_Site] > 0) continue;

        NCluster++;
        Mem[i_Site] = NCluster;

        rear = 0;
        front = 1;
        Que[front] = i_Site;
        depth_cluster = 0;
        
        int current_level_count = 1;  // 当前层的节点数
        int next_level_count = 0;     // 下一层的节点数
        
        while(front > rear)
        {
            depth_cluster++;
            for(int i = 0; i < current_level_count; ++i)
            {
                rear++;
                j_Site = Que[rear];
                hash_index = Bond.V2H[j_Site];
        
                while (hash_index != 0)
                {
                    k_Site = Bond.H2V[hash_index];
                    x_now[0] = Tx[j_Site] + Bond.Dx[hash_index];
                    x_now[1] = Ty[j_Site] + Bond.Dy[hash_index];
                    hash_index = Bond.Lnk[hash_index];
        
                    if (Mem[k_Site] == 0)
                    {
                        front++;
                        Mem[k_Site] = NCluster;
                        Que[front] = k_Site;
                        Tx[k_Site] = x_now[0];    Ty[k_Site] = x_now[1];
                        next_level_count++;
                    }
                    else if (x_now[0] != Tx[k_Site] || x_now[1] != Ty[k_Site])
                    {
                        wrap_x = abs(x_now[0] - Tx[k_Site]) > 0 ? L : 0;
                        wrap_y = abs(x_now[1] - Ty[k_Site]) > 0 ? L : 0;
        
                        if(Nwrap == 0)
                        {
                            Nwrap = 1;
                            listWrap_x[1] = wrap_x;
                            listWrap_y[1] = wrap_y;
                        }
                        else
                        {
                            wrap_flag = 1;
                            for (int i = 1; i <= Nwrap; i++)
                            {
                                if(wrap_x == listWrap_x[i] && wrap_y == listWrap_y[i]) wrap_flag = 0;
                            }
                            if(wrap_flag)
                            {
                                Nwrap++;
                                listWrap_x[Nwrap] = wrap_x;
                                listWrap_y[Nwrap] = wrap_y;
                            }
                        }
                    }
                }
            }
            current_level_count = next_level_count;
            next_level_count = 0;
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

        if(S1 < depth_cluster) S1 = depth_cluster;

    }


    //--- Wrapping Probabilities
    R0 = 0; R1 = 0; R2 = 0; Rd = 0;
    Rx = 0; Ry = 0; 
    if (Nwrap == 0) R0 = 1;
    else if (Nwrap == 1)
    {
        R1 = 1;
        if (listWrap_x[1]!=0 && listWrap_y[1]==0) Rx = 1;
        else if (listWrap_x[1]==0 && listWrap_y[1]!=0) Ry = 1;
        else if (listWrap_x[1]!=0 && listWrap_y[1]!=0) {Rd = 1; Rx = 1; Ry = 1;};
    }
    else if (Nwrap > 1) 
    {
        R2 = 1;
        Rx = 1;
        Ry = 1;
    }

}