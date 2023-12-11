//
//  Copyright (C) 2018 Satoshi Shimizu
//
//
//  This file is part of MECQ shown in the following paper:
//  - Satoshi Shimizu, Kazuami Yamaguchi, Sumio Masuda,
//    ``A Maximum Edge-Weight Clique Extraction Algorithm Based on Branch-and-Bound,''
//    https://arxiv.org/abs/1810.10258, 2018.
//
//  Note that Phased Local Search used by MECQ is originally proposed in the following paper:
//  - Wayne Pullan, ``Approximating the maximum vertex/edge weighted clique using local search,''
//    Journal of Heuristics 14 (2) (2008) 117–134.
//
//  MECQ is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  MECQ is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with MECQ.  If not, see <http://www.gnu.org/licenses/>.
//
#include"weighted_graph.h"
#include"phased_local_search_adjmatrix.h"
#include<cstdlib>
#include<cstdio>
#include<memory.h>
#include<algorithm>
#include<functional>
#include <cassert>
#include <limits>

namespace PHASED_LOCAL_SEARCH
{
    phased_local_search_adjmatrix::phased_local_search_adjmatrix(weighted_graph *graph) : phased_local_search(graph)
    {
        edge_weight=graph->edge_weight;
        nonadjacency_list=graph->nonadjacency_list;

        num_of_nonadj_in_clique=new int[n];
        for(int i=0;i<n; i++)
        {
            num_of_nonadj_in_clique[i]=0;
        }

        c0list=new int[n+1];
        c1list=new int[n+1];
        c0index=new int[n];
        c1index=new int[n];
        for(int i=0; i<n;i++)
        {
            c0list[i]=i;
            c0index[i]=i;
            c1index[i]=n;
        }
        c0list[n]=n; //dummy
        c1list[n]=n; //dummy
    }

    //
    // add vertex to clique and update num_of_nonadj_in_clique
    //
    void phased_local_search_adjmatrix::add_vertex_to_clique(add_info *info)
    {
        int v=info->v_add;
        remove_from_c0(v);
        clique[clique_size++]=v;
        clique_weight += info->weight_add;
        num_of_nonadj_in_clique[v]=n; //dummy number

        //update num_of_nonadj_in_clique for nonadj vertices of v
        int *nadj=nonadjacency_list[v];
        int nd=n-degree[v];
        for(int i=0; i<nd; i++)
        {
            int nv=nadj[i];
            num_of_nonadj_in_clique[nv]++;
            if(num_of_nonadj_in_clique[nv]==1)
            {
                remove_from_c0(nv);
                add_to_c1(nv);
            }
            else if(num_of_nonadj_in_clique[nv]==2)
            {
                remove_from_c1(nv);
            }
        }
        num_of_nonadj_in_clique[v]=n; //dummy number
    }

    //
    // remove clique[index] from clique
    // and update num_of_nonadj_in_clique
    //
    void phased_local_search_adjmatrix::drop_vertex_from_clique(drop_info *info)
    {
        int index=info->index_drop;
        int v_out=clique[index];
        clique_weight+=info->weight_drop;
        clique_size--;
        clique[index]=clique[clique_size];

        // update num_of_nonadj_in_clique for nonadj vertices of v_out
        {
            int *nadj=nonadjacency_list[v_out];
            int nd=n-degree[v_out];
            for(int i=0; i<nd; i++)
            {
                int nv=nadj[i];
                num_of_nonadj_in_clique[nv]--;
                if(num_of_nonadj_in_clique[nv]==0)
                {
                    remove_from_c1(nv);
                    add_to_c0(nv);
                }
                else if(num_of_nonadj_in_clique[nv]==1)
                {
                    add_to_c1(nv);
                }

            }
        }
        //add v_out to c0
        add_to_c0(v_out);
        num_of_nonadj_in_clique[v_out]=0;
    }

    add_info phased_local_search_adjmatrix::select_from_C0_random()
    {
        int *set = new int[c0size];
        int size=0;
        int max_increase = std::numeric_limits<int>::min();

        for(int i=0; i<c0size; i++)
        {
            int c0i=c0list[i];
            int increase = vertex_weight[c0i] + calc_ewsum(c0i);

            if(increase > max_increase)
            {
                size = 1;
                set[0] = c0i;
                max_increase = increase;
            }
            else if (increase == max_increase)
            {
                set[size++] = c0i;
            }
            else
            {
                //Nothing to do.
            }
        }

        int r = (int)( rand() * ((size-1) + 1.0) / (1.0 + RAND_MAX) );
        add_info info;
        info.v_add=set[r];
        info.weight_add=max_increase;
        delete[] set;
        return info;
    }

    swap_info phased_local_search_adjmatrix::select_from_C1_random()
    {
        int *set = new int[c1size-c1andUsize];
        int *increase = new int[c1size-c1andUsize];
        int *decrease = new int[c1size-c1andUsize];
        int *index = new int[c1size-c1andUsize];
        int size=0;
        int max_diff = std::numeric_limits<int>::min();

        for(int i=0; i<c1size-c1andUsize; i++)
        {
            int c1i = c1list[i];
            int* ewi = edge_weight[c1i];
            increase[i] = vertex_weight[c1i] + calc_ewsum(c1i);
            for(int j=0; j<clique_size; j++)
            {
                int cj=clique[j];
                if(ewi[cj] == NOT_ADJACENT)
                {
                    index[i]=j;
                    decrease[i] = -vertex_weight[cj] - calc_ewsum(cj);
                    break;
                }
            }

            if(increase[i] + decrease[i] > max_diff)
            {
                size = 1;
                set[0] = c1i;
                max_diff = increase[i] + decrease[i];
            }
            else if (increase[i] + decrease[i] == max_diff)
            {
                set[size++] = c1i;
            }
            else
            {
                //Nothing to do.
            }
        }
        int r = (int)( rand() * ((size-1) + 1.0) / (1.0 + RAND_MAX) );

        swap_info info;
        info.ainfo.v_add = set[r];
        info.ainfo.weight_add = increase[c1index[set[r]]];
        info.dinfo.index_drop = index[c1index[set[r]]];
        info.dinfo.weight_drop = decrease[c1index[set[r]]];

        delete[] set;
        delete[] increase;
        delete[] decrease;
        delete[] index;
        return info;
    }

    add_info phased_local_search_adjmatrix::select_from_C0_degree()
    {
        int *set = new int[c0size];
        int size=0;
        int max_degree = std::numeric_limits<int>::min();

        for(int i=0; i<c0size; i++)
        {
            int c0i = c0list[i];
            int di = degree[c0i];

            if(di > max_degree)
            {
                size = 1;
                set[0] = c0i;
                max_degree = di;
            }
            else if (di == max_degree)
            {
                set[size++] = c0i;
            }
            else
            {
                //Nothing to do.
            }
        }

        int v=set[0];
        int increase = vertex_weight[v] + calc_ewsum(v);
        // find a vertex with max-increase
        for(int i=1; i<size; i++)
        {
            int v2 = set[i];
            int increase2 = vertex_weight[v2] + calc_ewsum(v2);
            if( increase < increase2 )
            {
                v=v2;
                increase=increase2;
            }
        }

        add_info info;
        info.v_add = v;
        info.weight_add = increase;

        delete[] set;
        return info;
    }

    swap_info phased_local_search_adjmatrix::select_from_C1_degree()
    {
        int *set = new int[c1size-c1andUsize];
        int size=0;
        int max_degree = std::numeric_limits<int>::min();

        for(int i=0; i<c1size-c1andUsize; i++)
        {
            int c1i = c1list[i];
            int di = degree[c1i];

            if(di > max_degree)
            {
                size = 1;
                set[0] = c1i;
                max_degree = di;
            }
            else if (di == max_degree)
            {
                set[size++] = c1i;
            }
            else
            {
                //Nothing to do.
            }
        }

        int max_diff = std::numeric_limits<int>::min();
        int v=-1;
        int w_add=0;
        int w_drop=0;
        int index = -1;

        for(int i=0; i<size; i++)
        {
            int seti = set[i];
            int* ewi = edge_weight[seti];
            int increase = vertex_weight[seti] + calc_ewsum(seti);
            for(int j=0; j<clique_size; j++)
            {
                int cj=clique[j];
                if(ewi[cj] == NOT_ADJACENT)
                {
                    int decrease = -vertex_weight[cj] - calc_ewsum(cj);
                    if(increase + decrease > max_diff)
                    {
                        v=seti;
                        index=j;
                        max_diff=increase+decrease;
                        w_add = increase;
                        w_drop = decrease;
                    }
                    break;
                }
            }
        }

        swap_info info;
        info.ainfo.v_add = v;
        info.ainfo.weight_add = w_add;
        info.dinfo.index_drop = index;
        info.dinfo.weight_drop = w_drop;

        delete[] set;
        return info;
    }

    add_info phased_local_search_adjmatrix::select_from_C0_penalty()
    {
        int *set = new int[c0size];
        int size=0;
        int min_penalty = std::numeric_limits<int>::max();

        for(int i=0; i<c0size; i++)
        {
            int c0i = c0list[i];
            int pi = penalty[c0i];

            if(pi < min_penalty)
            {
                size = 1;
                set[0] = c0i;
                min_penalty = pi;
            }
            else if (pi == min_penalty)
            {
                set[size++] = c0i;
            }
            else
            {
                //Nothing to do.
            }
        }

        int v=set[0];
        int increase = vertex_weight[v] + calc_ewsum(v);
        // find a vertex with max-increase
        for(int i=1; i<size; i++)
        {
            int v2 = set[i];
            int increase2 = vertex_weight[v2] + calc_ewsum(v2);
            if( increase < increase2 )
            {
                v=v2;
                increase=increase2;
            }
        }

        add_info info;
        info.v_add = v;
        info.weight_add=increase;

        delete[] set;
        return info;
    }

    swap_info phased_local_search_adjmatrix::select_from_C1_penalty()
    {
        int *set = new int[c1size-c1andUsize];
        int size=0;
        int min_penalty = std::numeric_limits<int>::max();

        for(int i=0; i<c1size-c1andUsize; i++)
        {
            int c1i = c1list[i];
            int pi = penalty[c1i];

            if(pi < min_penalty)
            {
                size = 1;
                set[0] = c1i;
                min_penalty = pi;
            }
            else if (pi == min_penalty)
            {
                set[size++] = c1i;
            }
            else
            {
                //Nothing to do.
            }
        }

        int max_diff = std::numeric_limits<int>::min();
        int v=-1;
        int index=-1;
        int w_add=0;
        int w_drop=0;

        for(int i=0; i<size; i++)
        {
            int seti = set[i];
            int* ewi = edge_weight[seti];
            int increase = vertex_weight[seti] + calc_ewsum(seti);
            for(int j=0; j<clique_size; j++)
            {
                int cj=clique[j];
                if(ewi[cj] == NOT_ADJACENT)
                {
                    int decrease = -vertex_weight[cj] - calc_ewsum(cj);
                    if(increase + decrease > max_diff)
                    {
                        v=seti;
                        index=j;
                        max_diff=increase + decrease;
                        w_add = increase;
                        w_drop = decrease;
                    }
                    break;
                }
            }
        }

        swap_info info;
        info.ainfo.v_add = v;
        info.ainfo.weight_add = w_add;
        info.dinfo.index_drop = index;
        info.dinfo.weight_drop = w_drop;

        delete[] set;
        return info;
    }

    void phased_local_search_adjmatrix::initialize()
    {
        int v = (int)( rand() * ((n-1) + 1.0) / (1.0 + RAND_MAX) );
        clique[0]=v;
        clique_size=1;
        clique_weight=vertex_weight[v];
        int *ewv=edge_weight[v];
        num_of_nonadj_in_clique[v]=n;
        c0index[v]=n;
        c1index[v]=n;

        c0size=0;
        c1size=0;
        c1andUsize=0;

        for(int i=0; i<n; i++)
        {
            if(i==v) continue;
            if(ewv[i] != NOT_ADJACENT)
            {
                num_of_nonadj_in_clique[i]=0;
                c0index[i]=c0size;
                c0list[c0size++]=i;
                c1index[i]=n;
            }
            else
            {
                num_of_nonadj_in_clique[i]=1;
                c1index[i]=c1size;
                c1list[c1size++]=i;
                c0index[i]=n;
            }
            U[i]=false;
        }
    }


    void phased_local_search_adjmatrix::reinitialize()
    {
        //initialize U
        c1andUsize=0;
        for(int i=0; i<n; i++)
        {
            U[i]=false;
        }

        // select a new vertex to add to clique
        int v_new = (int)( rand() * ((n-1) + 1.0) / (1.0 + RAND_MAX) );

        // increment v_new while v_new is already in clique
        while(num_of_nonadj_in_clique[v_new]>=n)
        {
            v_new++;
            if(v_new>=n)
            {
                v_new=0;
            }
        }

        //delete vertices nonadjacent to v_new from clique
        int *ewv_new=edge_weight[v_new];
        for(int i=0;i<clique_size;i++)
        {
            int v_in_clique=clique[i];
            if(ewv_new[v_in_clique] == NOT_ADJACENT)
            {
                drop_info info={i,-vertex_weight[v_in_clique]-calc_ewsum(v_in_clique)};
                drop_vertex_from_clique(&info);
                i--;
            }
        }
        //add v_new to clique
        add_info info;
        info.v_add = v_new;
        info.weight_add = vertex_weight[v_new] + calc_ewsum(v_new);
        add_vertex_to_clique(&info);
    }

    inline void phased_local_search_adjmatrix::add_to_c0(int v)
    {
        c0list[c0size]=v;
        c0index[v]=c0size;
        c0size++;
    }

    inline void phased_local_search_adjmatrix::remove_from_c0(int v)
    {
        int v_c0last=c0list[c0size-1];
        int c0index_v=c0index[v];
        c0list[c0index_v]=v_c0last;
        c0index[v_c0last]=c0index_v;
        c0index[v]=n;
        c0size--;
    }

    inline void phased_local_search_adjmatrix::add_to_c1(int v)
    {
        c1list[c1size]=v;
        c1index[v]=c1size;
        c1size++;

        if(U[v])
        {
            c1andUsize++;
        }
        else
        {
            int v2=c1list[c1size-c1andUsize-1];
            c1index[v]=c1size-c1andUsize-1;
            c1list[c1size-c1andUsize-1]=v;
            c1index[v2]=c1size-1;
            c1list[c1size-1]=v2;
        }
    }

    inline void phased_local_search_adjmatrix::remove_from_c1(int v)
    {
        if(U[v])
        {
            int v2=c1list[c1size-1];
            int c1index_v=c1index[v];
            c1list[c1index_v]=v2;
            c1index[v2]=c1index_v;
            c1index[v]=n;
            c1size--;
            c1andUsize--;
        }
        else
        {
            int v2=c1list[c1size-1];
            int v3=c1list[c1size-c1andUsize-1];
            int c1index_v=c1index[v];
            c1list[c1index_v]=v3;
            c1list[c1size-c1andUsize-1]=v2;
            c1index[v2]=c1size-c1andUsize-1;
            c1index[v3]=c1index_v;

            c1index[v]=n;
            c1size--;
        }
    }


    int phased_local_search_adjmatrix::calc_ewsum(int v)
    {
        int ewsum=0;
        int *ewv=edge_weight[v];
        for(int i=0; i<clique_size; i++)
        {
            int ewvci=ewv[clique[i]];
            if(ewvci != NOT_ADJACENT)
            {
                ewsum += ewvci;
            }
        }
        return ewsum;
    }


    bool phased_local_search_adjmatrix::check_status()
    {
        if(!graph->is_clique(clique,clique_size,clique_weight)) return false;
        for(int i=0; i<clique_size; i++)
        {
            int v=clique[i];
            if(c0index[v] != n) return false;
            if(c1index[v] != n) return false;
            if(num_of_nonadj_in_clique[v] != n) return false;
        }
        for(int i=0; i<c0size; i++)
        {
            int v=c0list[i];
            if(c0index[v] != i) return false;
            if(c1index[v] != n) return false;
            if(num_of_nonadj_in_clique[v] != 0) return false;
        }
        for(int i=0; i<c1size; i++)
        {
            int v=c1list[i];
            if(c0index[v] != n) return false;
            if(c1index[v] != i) return false;
            if(num_of_nonadj_in_clique[v] != 1) return false;
        }
        for(int i=0; i<n; i++)
        {
            int v=i;
            if((num_of_nonadj_in_clique[v] == 1) && (c1index[v] == n)) return false;
            if((num_of_nonadj_in_clique[v] == 0) && (c0index[v] == n)) return false;
        }
        return true;
    }

    phased_local_search_adjmatrix::~phased_local_search_adjmatrix()
    {
        delete[] num_of_nonadj_in_clique;
        delete[] c0list;
        delete[] c1list;
        delete[] c0index;
        delete[] c1index;
    }
}
