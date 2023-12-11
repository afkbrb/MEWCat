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
#ifndef phased_local_search_adjlist_h
#define phased_local_search_adjlist_h

#include"phased_local_search.h"
#include"weighted_graph.h"
#include<functional>

namespace PHASED_LOCAL_SEARCH
{
    class phased_local_search_adjlist : public phased_local_search
    {
        public:
            phased_local_search_adjlist(weighted_graph *graph);
            ~phased_local_search_adjlist();

        private:
            int** adjacency_list;
            int** edge_weight_list;

            int* num_of_adj_in_clique;
            int* sk_list_size;
            int* sk_list_prev;
            int* sk_list_next;
            int* sk_list_top;
            int* sk_list_last;
            int* sk_and_U_size;

            int* ewsum;

            int *adjArray;

            void init_U();

            add_info  select_from_C0_random();
            swap_info select_from_C1_random();
            add_info  select_from_C0_degree();
            swap_info select_from_C1_degree();
            add_info  select_from_C0_penalty();
            swap_info select_from_C1_penalty();

            void reinitialize();
            void initialize();

            void add_vertex_to_clique(add_info *info);
            void drop_vertex_from_clique(drop_info *info);

            inline void add_to_sk(int v,int k);
            inline void remove_from_sk(int v,int k);

            void set_adjArray(int v);

            bool check_status();
    };
}
#endif
