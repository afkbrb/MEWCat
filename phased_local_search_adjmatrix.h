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
#ifndef phased_local_search_adjmatrix_h
#define phased_local_search_adjmatrix_h

#include"phased_local_search.h"
#include"weighted_graph.h"
#include<functional>

namespace PHASED_LOCAL_SEARCH
{
    class phased_local_search_adjmatrix : public phased_local_search
    {
        public:
            phased_local_search_adjmatrix(weighted_graph *graph);
            ~phased_local_search_adjmatrix();

        private:
            int** edge_weight;
            int** nonadjacency_list;

            int* num_of_nonadj_in_clique;
            int* c0list;
            int* c1list;
            int* c0index;
            int* c1index;

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

            inline void add_to_c0(int v);
            inline void remove_from_c0(int v);
            inline void add_to_c1(int v);
            inline void remove_from_c1(int v);

            int calc_ewsum(int v);

            bool check_status();
    };
}
#endif
