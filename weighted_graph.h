//
//  Copyright (C) 2018 Satoshi Shimizu
//
//
//  This file is part of MECQ shown in the following paper:
//  - Satoshi Shimizu, Kazuami Yamaguchi, Sumio Masuda,
//    ``A Maximum Edge-Weight Clique Extraction Algorithm Based on Branch-and-Bound,''
//    https://arxiv.org/abs/1810.10258, 2018.
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
#ifndef weighted_graph_h
#define weighted_graph_h

#include <limits>

#define NOT_ADJACENT (std::numeric_limits<int>::min()) /* vertices are not adjacent if edge_weight is this value */

class weighted_graph
{
    public:
        int n; /* number of vertices */
        int m; /* number of edges */
        int **edge_weight; /* This is also used as adjacency-matrix */
        int** nonadjacency_list;
        int** adjacency_list;
        int** edge_weight_list;
        int* degree;
        int* vertex_weight;

        weighted_graph(char *inFile);
        ~weighted_graph();
        bool is_clique(int* clique, int clique_size,int clique_weight);

        void create_nonadjlist();
        void delete_nonadjlist();
};

#endif
