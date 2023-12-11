#include <memory.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <algorithm>

#include "phased_local_search.h"
#include "phased_local_search_adjlist.h"
#include "phased_local_search_adjmatrix.h"
#include "weighted_graph.h"

using namespace std;

class Solver {
   public:
    int *maximum_weight_clique;
    int maximum_weight_clique_size;
    int maximum_weight_clique_weight;
    double elapsed_time_sec;

    unsigned long branch_count;
    unsigned long internal_count;

    int n;
    int m;
    int *degree;
    int *degree_original;
    int **adjlist;
    int **edge_weight;

    clock_t start_time;

    int *clique;
    int clique_size;
    int clique_weight;

    int *weight;
    int *weight_add;
    int *upper;
    int *upper_add;

    int average_degree = 0;

    vector<int> get_branches_dense(int *set, int set_size,
                                   int *pseudo_vertex_weight) {
        vector<int> B;

        // candidate set is the set of remaining vertices that can be add to
        // current independent set I
        int candidate_size = set_size;

        for (int i = 0; i < set_size; i++) {
            int vi = set[i];
            weight[vi] = pseudo_vertex_weight[i];
            weight_add[vi] = 0;
            upper[vi] = 0;
            upper_add[vi] = 0;
        }

        int v_index = -1;

        // find the vertex with max degree
        int max_degree = -1;
        for (int i = 0; i < set_size; i++) {
            int vi = set[i];
            if (degree[vi] > max_degree) {
                max_degree = degree[vi];
                v_index = i;
            }
        }

        // process vertices
        for (int i = 0; i < set_size; i++) {
            // [0, i) is processed
            // [i, i + candidate_size) are candidates for the current
            // independent set
            // [i + candidate_size, set_size) are not processed
            // and are not candidates.

            // swap i and v_index
            int v = set[v_index];
            int v_pseudo_weight = pseudo_vertex_weight[v_index];

            set[v_index] = set[i];
            pseudo_vertex_weight[v_index] = pseudo_vertex_weight[i];

            set[i] = v;
            pseudo_vertex_weight[i] = v_pseudo_weight;

            // before adding the pseudo weight
            // upper[v] = \sum_{i < \tau(v)} max{sigma(u) + w(u, v) | u \in I_i
            // \cap N(v)}
            upper[v] += pseudo_vertex_weight[i];

            // delete v from candidates
            candidate_size--;

            if (clique_weight + upper[v] > maximum_weight_clique_weight) {
                // then we do not add v to the independent set
                // add to B instead
                B.push_back(i);
            } else {
                // now we are adding v to the independent set
                int *ewv = edge_weight[v];

                // update weight_add and upper_add
                for (int j = i + 1; j < set_size; j++) {
                    int vj = set[j];
                    if (ewv[vj] != NOT_ADJACENT) {
                        weight_add[vj] = max(weight_add[vj], ewv[vj]);
                        upper_add[vj] = max(upper_add[vj], weight[v] + ewv[vj]);
                    }
                }

                // update candidate set
                int lo = i + 1;
                int hi = i + candidate_size;
                while (true) {
                    while (lo <= hi && ewv[set[lo]] == NOT_ADJACENT) lo++;
                    while (hi >= lo && ewv[set[hi]] != NOT_ADJACENT) hi--;

                    if (lo > hi) break;

                    // swap
                    int tmp_set = set[lo];
                    int tmp_pseudo_weight = pseudo_vertex_weight[lo];

                    set[lo] = set[hi];
                    pseudo_vertex_weight[lo] = pseudo_vertex_weight[hi];

                    set[hi] = tmp_set;
                    pseudo_vertex_weight[hi] = tmp_pseudo_weight;

                    lo++;
                    hi--;
                }
                candidate_size = hi - i;
            }

            // find from the candidate set the vertex with maximum degree
            int max_degree = -1;
            for (int j = i + 1; j < i + 1 + candidate_size; j++) {
                int vj = set[j];
                if (degree[vj] > max_degree) {
                    max_degree = degree[vj];
                    v_index = j;
                }
            }

            if (candidate_size == 0) {
                // now the independent set is maximal
                // open a new empty one

                candidate_size = set_size - i - 1;

                int max_degree = -1;
                for (int j = i + 1; j < set_size; j++) {
                    int vj = set[j];

                    // update weight and upper
                    weight[vj] += weight_add[vj];
                    upper[vj] += upper_add[vj];
                    weight_add[vj] = 0;
                    upper_add[vj] = 0;

                    // find the vertex with max degree
                    if (max_degree < degree[vj]) {
                        max_degree = degree[vj];
                        v_index = j;
                    }
                }
            }
        }

        return B;
    }

    vector<int> get_branches_sparse(int *set, int set_size,
                                    int *pseudo_vertex_weight) {
        vector<int> B;

        vector<int> order_in_set(n, set_size);
        // the vertices from set sorted by degree desc
        vector<int> degree_order(set_size);
        vector<int> order_in_degree(n, set_size);
        for (int i = 0; i < set_size; i++) {
            order_in_set[set[i]] = i;
            degree_order[i] = set[i];
        }
        // sort by degree desc
        sort(degree_order.begin(), degree_order.end(),
             [&](int u, int v) { return degree[u] > degree[v]; });
        for (int i = 0; i < set_size; i++) {
            order_in_degree[degree_order[i]] = i;
        }

        vector<vector<bool>> Is;
        // T[i] = true <=> degree_order[i] is a remaining vertex
        vector<bool> T(set_size + 1, true);
        T[set_size] = false;
        int T_index = 0;

        while (T_index < set_size) {
            int I_size = 0;
            // I[v] = true <=> v is in I
            vector<bool> I(n, false);
            vector<bool> X = T;
            int X_index = T_index;
            while (X_index < set_size) {
                // the first vertex in X
                // is the one with max degree
                int v = degree_order[X_index];
                int v_index = order_in_set[v];
                X[X_index] = false;
                T[X_index] = false;
                weight[v] = pseudo_vertex_weight[v_index];
                int upper = clique_weight + pseudo_vertex_weight[v_index];
                int degree_v = degree_original[v];
                int *adjlist_v = adjlist[v];
                int *ewv = edge_weight[v];
                for (auto &I : Is) {
                    int weight_add = 0;
                    int upper_add = 0;
                    for (int u_index = 0; u_index < degree_v; u_index++) {
                        int u = adjlist_v[u_index];
                        if (!I[u]) continue;
                        int w = ewv[u];
                        weight_add = max(weight_add, w);
                        upper_add = max(upper_add, weight[u] + w);
                    }
                    weight[v] += weight_add;
                    upper += upper_add;
                }

                if (upper <= maximum_weight_clique_weight) {
                    I[v] = true;
                    I_size++;
                    for (int u_index = 0; u_index < degree_v; u_index++) {
                        int u = adjlist_v[u_index];
                        X[order_in_degree[u]] = false;
                    }
                } else {
                    B.push_back(v_index);
                }

                while (X_index < set_size && !X[X_index]) X_index++;
            }

            while (T_index < set_size && !T[T_index]) T_index++;

            if (I_size == 0) {
                break;
            }

            Is.push_back(I);
        }

        return B;
    }

    vector<int> get_branches_degree(int *set, int set_size) {
        vector<int> B;

        // find the vertex with max degree
        int v_index = -1;
        int max_degree = -1;
        for (int i = 0; i < set_size; i++) {
            int vi = set[i];
            if (max_degree < degree[vi]) {
                max_degree = degree[vi];
                v_index = i;
            }
        }

        // include v and its non-neighbors in B
        int v = set[v_index];
        int *ewv = edge_weight[v];
        for (int i = 0; i < set_size; i++) {
            int vi = set[i];
            if (ewv[vi] == NOT_ADJACENT || vi == v) {
                B.push_back(i);
            }
        }

        return B;
    }

    void search(int *set, int set_size, int *pseudo_vertex_weight) {
        branch_count++;

        if (set_size == 0) {
            if (clique_weight > maximum_weight_clique_weight) {
                memcpy(maximum_weight_clique, clique,
                       sizeof(int) * clique_size);
                maximum_weight_clique_size = clique_size;
                maximum_weight_clique_weight = clique_weight;
                printf("node: %lu, size: %d, weight: %d\n", branch_count,
                       clique_size, clique_weight);
            }
            return;
        }

        // udpate degree
        if (branch_count == 1) {
            // no need to recompute at the root node
            for (int i = 0; i < set_size; i++) {
                int vi = set[i];
                degree[vi] = degree_original[vi];
            }
        } else {
            for (int i = 0; i < set_size; i++) {
                int vi = set[i];
                degree[vi] = 0;
            }
            for (int i = 0; i < set_size; i++) {
                int vi = set[i];
                for (int j = i + 1; j < set_size; j++) {
                    int vj = set[j];
                    if (edge_weight[vi][vj] != NOT_ADJACENT) {
                        degree[vi]++;
                        degree[vj]++;
                    }
                }
            }
        }

        // different implementations of the GetBranches function
        // described in the paper
        vector<int> B;
        if (set_size > 5 * average_degree) {
            // if the size of an independent set is likely to be larger
            // than the neighbor size of a vertex
            // we use this implementation

            // usually triggered at the first layers of very sparse graphs
            B = get_branches_sparse(set, set_size, pseudo_vertex_weight);
        } else {
            B = get_branches_dense(set, set_size, pseudo_vertex_weight);
        }

        vector<int> B1 = get_branches_degree(set, set_size);
        if (B.size() > B1.size()) {
            B = B1;
        }

        if (B.size()) {
            internal_count++;
        }

        // sort by degree asc
        sort(B.begin(), B.end(), [&](int i, int j) {
            int vi = set[i];
            int vj = set[j];
            return degree[vi] < degree[vj];
        });

        vector<bool> removed(set_size, false);

        int *set2 = new int[set_size];
        int set2_size = 0;
        int *pseudo_vertex_weight2 = new int[set_size];
        for (int i : B) {
            int vi = set[i];

            clique[clique_size++] = vi;
            clique_weight += pseudo_vertex_weight[i];

            set2_size = 0;

            int *ew_vi = edge_weight[vi];
            for (int j = 0; j < set_size; j++) {
                int vj = set[j];
                if (removed[j]) continue;
                if (ew_vi[vj] == NOT_ADJACENT) continue;
                set2[set2_size] = vj;
                pseudo_vertex_weight2[set2_size] =
                    pseudo_vertex_weight[j] + ew_vi[vj];
                set2_size++;
            }

            search(set2, set2_size, pseudo_vertex_weight2);

            removed[i] = true;

            clique_size--;
            clique_weight -= pseudo_vertex_weight[i];
        }

        delete[] set2;
        delete[] pseudo_vertex_weight2;
    }

    Solver(weighted_graph *graph, int *initial_clique = NULL,
           int initial_clique_size = 0, int initial_clique_weight = 0) {
        start_time = clock();

        n = graph->n;
        m = graph->m;
        degree_original = graph->degree;
        adjlist = graph->adjacency_list;

        edge_weight = graph->edge_weight;
        average_degree = (2 * m) / n;

        branch_count = 0;
        internal_count = 0;

        clique = new int[n];
        clique_size = 0;
        clique_weight = 0;

        maximum_weight_clique = new int[n];
        memcpy(maximum_weight_clique, initial_clique,
               sizeof(int) * initial_clique_size);
        maximum_weight_clique_size = initial_clique_size;
        maximum_weight_clique_weight = initial_clique_weight;

        clique_size = 0;
        clique_weight = 0;

        degree = new int[n];
        weight = new int[n];
        weight_add = new int[n];
        upper = new int[n];
        upper_add = new int[n];

        int *set = new int[n];
        for (int i = 0; i < n; i++) {
            set[i] = i;
        }
        int set_size = n;
        int *pseudo_vertex_weight = new int[n];
        for (int i = 0; i < n; i++) {
            pseudo_vertex_weight[i] = 0;
        }
        search(set, set_size, pseudo_vertex_weight);

        delete[] pseudo_vertex_weight;
        delete[] set;
        delete[] clique;
        delete[] degree;
        delete[] weight;
        delete[] weight_add;
        delete[] upper;
        delete[] upper_add;

        clock_t time2 = clock();
        elapsed_time_sec = ((double)(time2 - start_time) / CLOCKS_PER_SEC);
    }

    ~Solver() { delete[] maximum_weight_clique; }
};

int main(int argc, char **argv) {
    if (argc != 2) {
        cout << "usage: " << argv[0] << " <dataset>" << endl;
        return 0;
    }

    printf("dataset: %s\n", argv[1]);
    weighted_graph *graph = new weighted_graph((char *)argv[1]);

    int n = graph->n;
    int m = graph->m;
    double density = 2.0 * m / ((long)n * (n - 1));
    printf("n: %d, m: %d, density: %f\n", n, m, density);

    int *initial_clique = NULL;
    int initial_clique_size = 0;
    int initial_clique_weight = 0;
    double pls_time_sec = 0;

    int pls_iteration = 10;
    if (pls_iteration > 0) {
        using namespace PHASED_LOCAL_SEARCH;
        double edge_density = ((double)graph->m / (n * (n - 1) / 2));
        if (edge_density < 0.5) {
            phased_local_search *pls = new phased_local_search_adjlist(graph);
            pls->search(pls_iteration);

            initial_clique = new int[pls->best_clique_size];
            memcpy(initial_clique, pls->best_clique,
                   sizeof(int) * pls->best_clique_size);
            initial_clique_size = pls->best_clique_size;
            initial_clique_weight = pls->best_clique_weight;

            pls_time_sec = pls->elapsed_time_sec;

            delete pls;
        } else {
            graph->create_nonadjlist();
            phased_local_search *pls = new phased_local_search_adjmatrix(graph);
            pls->search(pls_iteration);

            initial_clique = new int[pls->best_clique_size];
            memcpy(initial_clique, pls->best_clique,
                   sizeof(int) * pls->best_clique_size);
            initial_clique_size = pls->best_clique_size;
            initial_clique_weight = pls->best_clique_weight;

            pls_time_sec = pls->elapsed_time_sec;

            delete pls;
            graph->delete_nonadjlist();
        }
    }

    printf("PLS time: %.2f\n", pls_time_sec);
    printf("lb: %d\n", initial_clique_weight);
    assert(graph->is_clique(initial_clique, initial_clique_size,
                            initial_clique_weight));

    fflush(stdout);
    Solver *clq = new Solver(graph, initial_clique, initial_clique_size,
                             initial_clique_weight);
    assert(graph->is_clique(clq->maximum_weight_clique,
                            clq->maximum_weight_clique_size,
                            clq->maximum_weight_clique_weight));

    printf("search time: %.2f\n", clq->elapsed_time_sec);
    printf("node explored: %ld\n", clq->branch_count);
    printf("internal node: %ld\n", clq->internal_count);
    printf("total time: %.2f\n", clq->elapsed_time_sec + pls_time_sec);
    printf("opt: %d\n", clq->maximum_weight_clique_weight);

    // printf("The maximum weight clique has %d vertices,\n [",
    //        clq->maximum_weight_clique_size);
    // for (int i = 0; i < clq->maximum_weight_clique_size; ++i) {
    //     printf(" %d", clq->maximum_weight_clique[i] + 1);
    // }
    // printf(" ]\n");

    delete[] initial_clique;
    delete clq;
    delete graph;
    return 0;
}
