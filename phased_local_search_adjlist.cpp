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
#include"phased_local_search_adjlist.h"
#include<cstdlib>
#include<cstdio>
#include<memory.h>
#include<algorithm>
#include<functional>
#include <cassert>
#include <limits>

namespace PHASED_LOCAL_SEARCH
{
    phased_local_search_adjlist::phased_local_search_adjlist(weighted_graph *graph) : phased_local_search(graph)
    {
        adjacency_list=graph->adjacency_list;
        edge_weight_list=graph->edge_weight_list;

        num_of_adj_in_clique=new int[n];
        for(int i=0;i<n; i++)
        {
            num_of_adj_in_clique[i]=0;
        }

        sk_list_size = new int[n+1];
        sk_list_prev = new int[n+1];
        sk_list_next = new int[n+1];
        sk_list_top = new int[n+1];
        sk_list_last = new int[n+1];
        sk_and_U_size = new int[n+1];

        {
            for(int i=1; i<n+1; i++)
            {
                sk_list_top[i] = n; //dummy
                sk_list_last[i] = n; //dummy
                sk_list_size[i] = 0; //dummy
                sk_and_U_size[i] = 0;
            }
            sk_list_top[0] = 0;
            sk_list_last[0] = n-1;
            sk_list_size[0] = n;
            sk_and_U_size[0] = 0;

            for(int i=0; i<n;i++)
            {
                sk_list_prev[i] = i-1;
                sk_list_next[i] = i+1;
            }
            sk_list_prev[0] = n; //dummy
            sk_list_next[n-1] = n; //dummy
        }

        ewsum=new int[n];
        for(int i=0;i<n; i++)
        {
            ewsum[i]=0;
        }

        adjArray=new int[n];
        for(int i=0;i<n; i++)
        {
            adjArray[i]=-1;
        }
    }

    void phased_local_search_adjlist::init_U()
    {
        phased_local_search::init_U();
        for(int i=0; i<n; i++)
        {
            sk_and_U_size[i]=0;
        }

    }

    //
    // add vertex to clique and update num_of_nonadj_in_clique
    //
    void phased_local_search_adjlist::add_vertex_to_clique(add_info *info)
    {
        int v=info->v_add;
        remove_from_sk(v, clique_size);
        clique[clique_size++]=v;
        clique_weight += info->weight_add;
        num_of_adj_in_clique[v]=2*n; //dummy number

        //update num_of_adj_in_clique for adj vertices of v
        int *adj=adjacency_list[v];
        int *ew=edge_weight_list[v];
        int d=degree[v];
        for(int i=0; i<d; i++)
        {
            int u=adj[i];
            ewsum[u]+=ew[i];
            if(num_of_adj_in_clique[u] > n) // if u is in clique
            {
                continue;
            }

            remove_from_sk(u, num_of_adj_in_clique[u]);
            num_of_adj_in_clique[u]++;
            add_to_sk(u, num_of_adj_in_clique[u]);

        }
        num_of_adj_in_clique[v]=2*n; //dummy number

        c0size = sk_list_size[clique_size];
        c1size = sk_list_size[clique_size-1];
        c1andUsize = sk_and_U_size[clique_size-1];
    }

    //
    // remove clique[index] from clique
    // and update num_of_nonadj_in_clique
    //
    void phased_local_search_adjlist::drop_vertex_from_clique(drop_info *info)
    {
        int index=info->index_drop;
        int v_out=clique[index];
        clique_weight+=info->weight_drop;
        clique_size--;
        clique[index]=clique[clique_size];

        num_of_adj_in_clique[v_out]=clique_size;
        add_to_sk(v_out, clique_size);

        // update num_of_adj_in_clique for adj vertices of v_out
        {
            int *adj=adjacency_list[v_out];
            int *ew=edge_weight_list[v_out];
            int d=degree[v_out];
            for(int i=0; i<d; i++)
            {
                int u=adj[i];
                ewsum[u]-=ew[i];
                if(num_of_adj_in_clique[u] > n) // if u is in clique
                {
                    continue;
                }
                remove_from_sk(u, num_of_adj_in_clique[u]);
                num_of_adj_in_clique[u]--;
                add_to_sk(u, num_of_adj_in_clique[u]);
            }
        }

        c0size = sk_list_size[clique_size];
        c1size = sk_list_size[clique_size-1];
        c1andUsize = sk_and_U_size[clique_size-1];
    }

    add_info phased_local_search_adjlist::select_from_C0_random()
    {
        int *set = new int[sk_list_size[clique_size]];
        int size=0;
        int max_increase = std::numeric_limits<int>::min();

        {
            int v = sk_list_top[clique_size];
            while( v != n )
            {
                int increase = vertex_weight[v] + ewsum[v];
                if(increase > max_increase)
                {
                    size = 1;
                    set[0] = v;
                    max_increase = increase;
                }
                else if (increase == max_increase)
                {
                    set[size++] = v;
                }
                else
                {
                    //Nothing to do.
                }
                v = sk_list_next[v];
            }
        }

        int r = (int)( rand() * ((size-1) + 1.0) / (1.0 + RAND_MAX) );
        add_info info;
        info.v_add=set[r];
        info.weight_add=max_increase;
        delete[] set;
        return info;
    }

    swap_info phased_local_search_adjlist::select_from_C1_random()
    {
        int *set = new int[sk_list_size[clique_size-1] - sk_and_U_size[clique_size-1]];
        int *index= new int[sk_list_size[clique_size-1] - sk_and_U_size[clique_size-1]];
        int size=0;
        int max_diff = std::numeric_limits<int>::min();

        {
            int v = sk_list_top[clique_size-1];
            while( U[v] == false )
            {
                set_adjArray(v);
                int increase = vertex_weight[v] + ewsum[v];
                int drop_index = -1;
                for(int i=0; i<clique_size; i++)
                {
                    int ci=clique[i];
                    if(adjArray[ci] != v) //if nonadjacent
                    {
                        drop_index=i;
                        increase -= vertex_weight[ci] + ewsum[ci];
                        break;
                    }
                }

                if(increase > max_diff)
                {
                    size = 1;
                    set[0] = v;
                    index[0] = drop_index;
                    max_diff = increase;
                }
                else if (increase == max_diff)
                {
                    set[size] = v;
                    index[size] = drop_index;
                    size++;
                }
                else
                {
                    //Nothing to do.
                }
                v = sk_list_next[v];
            }
        }
        int r = (int)( rand() * ((size-1) + 1.0) / (1.0 + RAND_MAX) );
        int v_add = set[r];
        int index_drop = index[r];
        int v_drop = clique[index_drop];

        swap_info info;
        info.ainfo.v_add = v_add;
        info.ainfo.weight_add = vertex_weight[v_add] + ewsum[v_add];
        info.dinfo.index_drop = index_drop;
        info.dinfo.weight_drop = -vertex_weight[v_drop] - ewsum[v_drop];

        delete[] set;
        delete[] index;
        return info;
    }

    add_info phased_local_search_adjlist::select_from_C0_degree()
    {
        int *set = new int[sk_list_size[clique_size]];
        int size=0;
        int max_degree = std::numeric_limits<int>::min();

        {
            int v = sk_list_top[clique_size];
            while( v != n )
            {
                int d = degree[v];
                if(d > max_degree)
                {
                    size = 1;
                    set[0] = v;
                    max_degree = d;
                }
                else if (d == max_degree)
                {
                    set[size++] = v;
                }
                else
                {
                    //Nothing to do.
                }
                v = sk_list_next[v];
            }
        }

        int v=set[0];
        int increase = vertex_weight[v] + ewsum[v];
        // find a vertex with max-increase
        for(int i=1; i<size; i++)
        {
            int v2 = set[i];
            int increase2 = vertex_weight[v2] + ewsum[v2];
            if( increase < increase2 )
            {
                v=v2;
                increase=increase2;
            }
        }

        delete[] set;

        add_info info;
        info.v_add=v;
        info.weight_add=increase;

        return info;
    }

    swap_info phased_local_search_adjlist::select_from_C1_degree()
    {
        int *set = new int[sk_list_size[clique_size-1] - sk_and_U_size[clique_size-1]];
        int size=0;
        int max_degree = std::numeric_limits<int>::min();

        {
            int v = sk_list_top[clique_size-1];
            while( U[v] == false )
            {
                int d = degree[v];

                if(d > max_degree)
                {
                    size = 1;
                    set[0] = v;
                    max_degree = d;
                }
                else if (d == max_degree)
                {
                    set[size++] = v;
                }
                else
                {
                    //Nothing to do.
                }
                v = sk_list_next[v];
            }
        }

        int max_diff = std::numeric_limits<int>::min();
        int v=-1;
        int index = -1;

        for(int i=0; i<size; i++)
        {
            int seti = set[i];
            set_adjArray(seti);
            int increase = vertex_weight[seti] + ewsum[seti];
            for(int j=0; j<clique_size; j++)
            {
                int cj=clique[j];
                if(adjArray[cj] != seti) //if nonadjacent
                {
                    increase -= vertex_weight[cj] + ewsum[cj];
                    if(increase > max_diff)
                    {
                        v=seti;
                        index=j;
                        max_diff=increase;
                    }
                    break;
                }
            }
        }

        delete[] set;

        swap_info info;
        info.ainfo.v_add = v;
        info.ainfo.weight_add = vertex_weight[v] + ewsum[v];
        info.dinfo.index_drop = index;
        int v_drop = clique[index];
        info.dinfo.weight_drop = -vertex_weight[v_drop] - ewsum[v_drop];

        return info;
    }

    add_info phased_local_search_adjlist::select_from_C0_penalty()
    {
        int *set = new int[sk_list_size[clique_size]];
        int size=0;
        int min_penalty = std::numeric_limits<int>::max();

        {
            int v = sk_list_top[clique_size];
            while( v != n )
            {
                int p = penalty[v];

                if(p < min_penalty)
                {
                    size = 1;
                    set[0] = v;
                    min_penalty = p;
                }
                else if (p == min_penalty)
                {
                    set[size++] = v;
                }
                else
                {
                    //Nothing to do.
                }
                v = sk_list_next[v];
            }
        }

        int v=set[0];
        int increase = vertex_weight[v] + ewsum[v];
        // find a vertex with max-increase
        for(int i=1; i<size; i++)
        {
            int v2 = set[i];
            int increase2 = vertex_weight[v2] + ewsum[v2];
            if( increase < increase2 )
            {
                v=v2;
                increase=increase2;
            }
        }
        delete[] set;

        add_info info;
        info.v_add=v;
        info.weight_add=increase;

        return info;
    }

    swap_info phased_local_search_adjlist::select_from_C1_penalty()
    {
        int *set = new int[sk_list_size[clique_size-1] - sk_and_U_size[clique_size-1]];
        int size=0;
        int min_penalty = std::numeric_limits<int>::max();

        {
            int v = sk_list_top[clique_size-1];
            while( U[v] == false )
            {
                int p = penalty[v];

                if(p < min_penalty)
                {
                    size = 1;
                    set[0] = v;
                    min_penalty = p;
                }
                else if (p == min_penalty)
                {
                    set[size++] = v;
                }
                else
                {
                    //Nothing to do.
                }
                v = sk_list_next[v];
            }
        }

        int max_diff = std::numeric_limits<int>::min();
        int v=-1;
        int index=-1;

        for(int i=0; i<size; i++)
        {
            int seti = set[i];
            set_adjArray(seti);
            int increase = vertex_weight[seti] + ewsum[seti];
            for(int j=0; j<clique_size; j++)
            {
                int cj=clique[j];
                if(adjArray[cj] != seti) //if nonadjacent
                {
                    increase -= vertex_weight[cj] + ewsum[cj];
                    if(increase > max_diff)
                    {
                        v=seti;
                        index=j;
                        max_diff=increase;
                    }
                    break;
                }
            }
        }

        delete[] set;

        swap_info info;
        info.ainfo.v_add = v;
        info.ainfo.weight_add = vertex_weight[v] + ewsum[v];
        info.dinfo.index_drop = index;
        int v_drop = clique[index];
        info.dinfo.weight_drop = -vertex_weight[v_drop] - ewsum[v_drop];

        return info;
    }

    void phased_local_search_adjlist::initialize()
    {
        clique_size=0;
        clique_weight=0;

        for(int i=0; i<n; i++)
        {
            num_of_adj_in_clique[i]=0;
            ewsum[i]=0;
            U[i]=false;
        }

        {
            for(int i=1; i<n+1; i++)
            {
                sk_list_top[i] = n; //dummy
                sk_list_last[i] = n; //dummy
                sk_list_size[i] = 0; //dummy
                sk_and_U_size[i] = 0;
            }
            sk_list_top[0] = 0;
            sk_list_last[0] = n-1;
            sk_list_size[0] = n;
            sk_and_U_size[0] = 0;

            for(int i=0; i<n;i++)
            {
                sk_list_prev[i] = i-1;
                sk_list_next[i] = i+1;
            }
            sk_list_prev[0] = n; //dummy
            sk_list_next[n-1] = n; //dummy
        }

        // add a vertex to clique randomly
        int v = (int)( rand() * ((n-1) + 1.0) / (1.0 + RAND_MAX) );
        add_info info={v,vertex_weight[v]};
        add_vertex_to_clique(&info);

        c0size = sk_list_size[clique_size];
        c1size = sk_list_size[clique_size-1];
        c1andUsize = sk_and_U_size[clique_size-1];
    }


    void phased_local_search_adjlist::reinitialize()
    {
        //initialize U
        for(int i=0; i<n; i++)
        {
            sk_and_U_size[i]=0;
            U[i]=false;
        }

        // select a new vertex to add to clique
        int v_new = (int)( rand() * ((n-1) + 1.0) / (1.0 + RAND_MAX) );

        // increment v_new while v_new is already in clique
        while(num_of_adj_in_clique[v_new] > n)
        {
            v_new++;
            if(v_new>=n)
            {
                v_new=0;
            }
        }

        set_adjArray(v_new);

        //delete vertices nonadjacent to v_new from clique
        for(int i=0;i<clique_size;i++)
        {
            int v_in_clique=clique[i];
            if(adjArray[v_in_clique] != v_new) //if nonadjacent
            {
                drop_info info={i, -vertex_weight[v_in_clique] - ewsum[v_in_clique]};
                drop_vertex_from_clique(&info);
                i--;
            }
        }
        //add v_new to clique
        add_info info={v_new,vertex_weight[v_new]+ewsum[v_new]};
        add_vertex_to_clique(&info);

        c0size = sk_list_size[clique_size];
        c1size = sk_list_size[clique_size-1];
        c1andUsize = sk_and_U_size[clique_size-1];
    }

    inline void phased_local_search_adjlist::add_to_sk(int v, int k)
    {
        if(U[v])
        {
            //add v to last
            int u = sk_list_last[k];
            sk_list_next[u] = v;
            sk_list_prev[v] = u;
            sk_list_next[v] = n;
            sk_list_last[k] = v;
            sk_list_size[k]++;
            sk_and_U_size[k]++;
            if( u == n )
            {
                sk_list_top[k] = v;
            }
        }
        else
        {
            //add v to top
            int u = sk_list_top[k];
            sk_list_prev[u] = v;
            sk_list_prev[v] = n;
            sk_list_next[v] = u;
            sk_list_top[k] = v;
            sk_list_size[k]++;
            if( u == n )
            {
                sk_list_last[k] = v;
            }
        }
    }

    inline void phased_local_search_adjlist::remove_from_sk(int v, int k)
    {
        int prev = sk_list_prev[v];
        int next = sk_list_next[v];
        sk_list_prev[next] = prev;
        sk_list_next[prev] = next;
        if( prev == n )
        {
            sk_list_top[k] = next;
        }
        if( next == n )
        {
            sk_list_last[k] = prev;
        }
        sk_list_size[k]--;
        sk_list_prev[v] = n;
        sk_list_next[v] = n;
        if(U[v])
        {
            sk_and_U_size[k]--;
        }
    }

    void phased_local_search_adjlist::set_adjArray(int v)
    {
        int *adj=adjacency_list[v];
        int d=degree[v];
        for( int i=0; i<d; i++)
        {
            int u=adj[i];
            adjArray[u]=v;
        }
    }

    bool phased_local_search_adjlist::check_status()
    {
        if(!graph->is_clique(clique,clique_size,clique_weight)) return false;
        if(c0size != sk_list_size[clique_size]) return false;
        if(c1size != sk_list_size[clique_size-1]) return false;
        if(c1andUsize != sk_and_U_size[clique_size-1]) return false;

        {
            int *ewsum_test=new int[n];
            for(int i=0; i<n; i++)
            {
                ewsum_test[i]=0;
            }
            for(int i=0; i<clique_size; i++)
            {
                int v=clique[i];
                int *adj=adjacency_list[v];
                int *ew=edge_weight_list[v];
                int d=degree[v];
                for(int j=0; j<d; j++)
                {
                    int u=adj[j];
                    ewsum_test[u] += ew[j];
                }
            }
            for(int i=0; i<n; i++)
            {
                if(ewsum[i] != ewsum_test[i])
                {
                    delete[] ewsum_test;
                    return false;
                }
            }
            delete[] ewsum_test;
        }
        {
            bool *checked = new bool[n+1];
            for(int i=0; i<n; i++)
            {
                checked[i] = false;
            }

            for(int i=0; i<clique_size; i++)
            {
                int v=clique[i];
                if(sk_list_prev[v] != n)
                {
                    delete[] checked;
                    return false;
                }
                if(sk_list_next[v] != n)
                {
                    delete[] checked;
                    return false;
                }
                if(num_of_adj_in_clique[v] <= n)
                {
                    delete[] checked;
                    return false;
                }
                if( checked[v] == true )
                {
                    delete[] checked;
                    return false;
                }
                checked[v] = true;
            }
            for(int i=0; i<n; i++)
            {
                if(sk_list_prev[i] == i)
                {
                    delete[] checked;
                    return false;
                }
                if(sk_list_next[i] == i)
                {
                    delete[] checked;
                    return false;
                }

                int size=0;
                int v = sk_list_top[i];
                int prev = n; //dummy
                bool check_U=false;
                int Usize=0;
                while( v != n )
                {
                    if( checked[v] == true )
                    {
                        delete[] checked;
                        return false;
                    }
                    if(sk_list_prev[v] != prev)
                    {
                        delete[] checked;
                        return false;
                    }
                    if(num_of_adj_in_clique[v] != i)
                    {
                        delete[] checked;
                        return false;
                    }

                    checked[v] = true;
                    size++;
                    prev = v;

                    if( (check_U == true) && (U[v] == false) )
                    {
                        delete[] checked;
                        return false;
                    }
                    if(U[v])
                    {
                        check_U = true;
                        Usize++;
                    }

                    v = sk_list_next[v];
                }
                if( size != sk_list_size[i] )
                {
                    delete[] checked;
                    return false;
                }
                if( Usize != sk_and_U_size[i] )
                {
                    delete[] checked;
                    return false;
                }
            }
            for(int i=0; i<n; i++)
            {
                if(checked[i] == false)
                {
                    delete[] checked;
                    return false;
                }
            }
            delete[] checked;
        }
        return true;
    }

    phased_local_search_adjlist::~phased_local_search_adjlist()
    {
        delete[] num_of_adj_in_clique;
        delete[] ewsum;
        delete[] sk_list_size;
        delete[] sk_list_prev;
        delete[] sk_list_next;
        delete[] sk_list_top;
        delete[] sk_list_last;
        delete[] sk_and_U_size;
        delete[] adjArray;
    }
}
