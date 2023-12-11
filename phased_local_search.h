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
#ifndef phased_local_search_h
#define phased_local_search_h

#include"weighted_graph.h"
#include<functional>
#include <ctime>

namespace PHASED_LOCAL_SEARCH
{
    typedef struct
    {
        int v_add;
        int weight_add;
    } add_info;

    typedef struct
    {
        int index_drop;
        int weight_drop;
    } drop_info;

    typedef struct
    {
        add_info ainfo;
        drop_info dinfo;
    } swap_info;

    class phased_local_search
    {
        public:
            phased_local_search(weighted_graph *graph);
            virtual ~phased_local_search();

            void search(int iterations);

            int* best_clique;
            int  best_clique_size;
            int  best_clique_weight;
            double elapsed_time_sec;

        protected:
            weighted_graph *graph;
            int n;
            int m;
            int* degree;
            int* vertex_weight;

            clock_t start_time;

            int* clique;
            int  clique_size;
            int  clique_weight;

            int c0size;
            int c1size;
            int c1andUsize;
            bool* U;

            int* penalty;
            int penalty_deley;
            int penalty_cycle_count;
            int num_of_penalised_vertices;

            void phase(int iterations,
                    std::function<add_info()> select_from_C0,
                    std::function<swap_info()> select_from_C1,
                    std::function<void()> perturb);
            void update_penalties();

            virtual void init_U();

            virtual add_info  select_from_C0_random() = 0;
            virtual swap_info select_from_C1_random() = 0;
            virtual add_info  select_from_C0_degree() = 0;
            virtual swap_info select_from_C1_degree() = 0;
            virtual add_info  select_from_C0_penalty() = 0;
            virtual swap_info select_from_C1_penalty() = 0;

            virtual void reinitialize() = 0;
            virtual void initialize() = 0;

            virtual void add_vertex_to_clique(add_info *info) = 0;
            virtual void drop_vertex_from_clique(drop_info *info) = 0;

            virtual bool check_status();
    };
}
#endif
