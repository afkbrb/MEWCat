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
#include"phased_local_search.h"
#include<cstdlib>
#include<cstdio>
#include<memory.h>
#include<algorithm>
#include<functional>
#include <cassert>
#include <limits>

namespace PHASED_LOCAL_SEARCH
{
    phased_local_search::phased_local_search(weighted_graph *graph)
    {
        this->graph=graph;
        n=graph->n;
        m=graph->m;
        degree=graph->degree;
        vertex_weight=graph->vertex_weight;

        start_time=clock();

        c0size=n;
        c1size=0;
        c1andUsize=0;
        U=new bool[n+1];
        for(int i=0; i<n; i++)
        {
            U[i]=false;
        }
        U[n]=true; //dummy

        penalty=new int[n];
        for(int i=0; i<n; i++)
        {
            penalty[i]=0;
        }
        penalty_deley=2; //set 2 initially
        penalty_cycle_count=0;
        num_of_penalised_vertices=0;

        clique=new int[n];
        clique_size=0;
        clique_weight=0;

        best_clique=new int[n];
        best_clique_size=0;
        best_clique_weight=0;
    }

    void phased_local_search::search(int iterations)
    {
        //Randomly select a vertex v and add it to clique
        {
            int v = (int)( rand() * ((n-1) + 1.0) / (1.0 + RAND_MAX) );
            add_info info={v,vertex_weight[v]};
            add_vertex_to_clique(&info);
        }

        for(int i=0; i<iterations; i++)
        {
            // phase random
            phase(50,
                    std::bind( &phased_local_search::select_from_C0_random, this),
                    std::bind( &phased_local_search::select_from_C1_random, this),
                    std::bind( &phased_local_search::reinitialize, this)
                 );

            // phase penalty
            phase(50,
                    std::bind( &phased_local_search::select_from_C0_penalty, this),
                    std::bind( &phased_local_search::select_from_C1_penalty, this),
                    std::bind( &phased_local_search::initialize, this)
                 );

            // phase degree
            phase(100,
                    std::bind( &phased_local_search::select_from_C0_degree, this),
                    std::bind( &phased_local_search::select_from_C1_degree, this),
                    std::bind( &phased_local_search::reinitialize, this)
                 );
        }
        clock_t time2=clock();
        elapsed_time_sec=((double)(time2-start_time)/CLOCKS_PER_SEC);
    }

    void phased_local_search::phase(int iterations,
            std::function<add_info()> select_from_C0,
            std::function<swap_info()> select_from_C1,
            std::function<void()> perturb)
    {
        while(iterations-->0)
        {
            while( c0size>0 || c1size-c1andUsize>0 )
            {
                // maximalize
                if(c0size>0)
                {
                    // add one vertex
                    {
                        add_info info=select_from_C0();
                        add_vertex_to_clique(&info);
                    }

                    init_U();

                    // add vertices until c0 is empty
                    while(c0size>0)
                    {
                        add_info info=select_from_C0();
                        add_vertex_to_clique(&info);
                    }
                }
                //update best clique
                if(clique_weight > best_clique_weight)
                {
                    best_clique_size=clique_size;
                    best_clique_weight=clique_weight;
                    memcpy(best_clique,clique,sizeof(int)*clique_size);
                    //clock_t time2=clock();
                    //double elapsed_time_sec=((double)(time2-start_time)/CLOCKS_PER_SEC);
                    //printf("%.2f,%d\n",elapsed_time_sec,best_clique_weight);
                }
                //swap
                if(c1size-c1andUsize>0)
                {
                    swap_info info=select_from_C1();
                    U[clique[info.dinfo.index_drop]]=true;
                    drop_vertex_from_clique(&info.dinfo);
                    add_vertex_to_clique(&info.ainfo);
                }
            }
            update_penalties();
            perturb();
        }
    }

    void phased_local_search::init_U()
    {
        c1andUsize=0;
        for(int i=0; i<n; i++)
        {
            U[i]=false;
        }
    }

    void phased_local_search::update_penalties()
    {
        // increment penalties for all vertices in clique
        for(int i=0; i<clique_size; i++)
        {
            if(penalty[clique[i]]++ == 0)
            {
                num_of_penalised_vertices++;
            }
        }

        // decrement all non-zero penalties for every penalty_deley
        if(++penalty_cycle_count > penalty_deley)
        {
            penalty_cycle_count=0;
            for(int i=0; i<n; i++)
            {
                int pi=penalty[i];
                if(pi)
                {
                    if(pi==1)
                    {
                        num_of_penalised_vertices--;
                    }
                    penalty[i]--;
                }
            }
        }

        // update penalty_deley
        int n75=n*0.75;
        if( num_of_penalised_vertices < n75 )
        {
            penalty_deley++;
        }
        else if( num_of_penalised_vertices > n75 && penalty_deley>1 )
        {
            penalty_deley--;
        }
    }

    bool phased_local_search::check_status()
    {
        return true;
    }

    phased_local_search::~phased_local_search()
    {
        delete[] clique;
        delete[] best_clique;
        delete[] U;
        delete[] penalty;
    }
}
