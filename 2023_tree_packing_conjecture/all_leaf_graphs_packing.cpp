/*
 * File:   all_leaf_graphs_packing.cpp
 * Author: Nikolay Ulyanov
 *
 * Created on 26 September 2023
 */

#include <algorithm>
#include <random>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdio>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <map>
#include <string>

#include "main.h"


using namespace std;

const int maxn = 15; // numbers are slightly arbitrary

static int totalAmounts[maxn] = {0};
short graphsDatabaseEdges[maxn][30000][maxn][2];
short graphsDatabaseN[maxn][30000];

int n;

int label[maxn][maxn];
vector<vector<int>> modGracefulLabels;
vector<int> curModGracefulLabel;

// chosen graphs for packing
vector<int> parents[maxn][maxn];
int ns[maxn];
int degrees[maxn][maxn];
bool usedEdge[maxn][maxn];

bool usedEdgeInKn[maxn][maxn];
int unusedDegreeInKn[maxn];

int vertices[maxn][maxn];

bool res;

int cnt;

int depth;
const int maxDepth = 10000;

int profile[maxn];

clock_t beginClock, curClock, prevClock, endClock;


// graph6 specification:
// https://cs.anu.edu.au/~bdm/data/formats.txt
// Convert graph6 character sequence to 6-bit integers
vector<int> string_to_graph6(const string& s) {
    vector<int> data;
    for (const auto& c : s) {
        data.push_back(static_cast<int>(c) - 63);
    }
    int min_v = 64;
    int max_v = 0;
    for (const auto& v : data) {
        min_v = min(min_v, v);
        max_v = max(max_v, v);
    }
    if (data.size() > 0 && (min_v < 0 or max_v > 63)) {
        data.clear();
        return data;
    }
    return data;
}

// read initial one-, four- or eight-unit value from graph6 integer sequence
// return value and rest of sequence
int graph6_to_vertex_count(const vector<int>& data) {
    if (data[0] <= 62) {
        return data[0];
    } else {
        throw std::invalid_argument("expected small vertex count");
    }
    return 0;
}

vector<bool> graph6_to_edge_pairs(const vector<int>& data) {
    vector<bool> edge_pairs;
    for (int i = 1; i < data.size(); ++i) {
        for (int bit_index = 5; bit_index >= 0; --bit_index) {
            int bit = (data[i] >> bit_index);
            if ((bit & 1) == 1) {
                edge_pairs.push_back(true);
            } else {
                edge_pairs.push_back(false);
            }
        }
    }
    return edge_pairs;
}

bool decode_graph6(istream& input, vector<pair<int, int>>& graph_edges) {
    string s;
    getline(input, s);
    if (s.empty()) {
        return false;
    }
    vector<int> data = string_to_graph6(s);
    int vertex_count = graph6_to_vertex_count(data);
    int edge_pair_count = (vertex_count * (vertex_count - 1) / 2 + 5) / 6;
    if (data.size() - 1 != edge_pair_count) {
        throw std::invalid_argument("wrong edge pair count");
    }
    graph_edges.clear();
    vector<bool> edge_pairs = graph6_to_edge_pairs(data);
    int v1 = 0;
    int v2 = 1;
    for (const auto& ep : edge_pairs) {
        if (ep) {
            graph_edges.push_back(make_pair(v1, v2));
        }
        if (v2 - v1 > 1) {
            v1 += 1;
        } else {
            v1 = 0;
            v2 += 1;
        }
    }
    return true;
}

bool gen(int num, int v) {
    if (depth > maxDepth) {
        return false;
    }
    ++depth;
    if (v == ns[num]) {
        if (num == 2) {
            ++cnt;
            for (int k = 2; k <= n; ++k) {
                curModGracefulLabel.clear();
                for (int i = 0; i < ns[k]; ++i) {
                    curModGracefulLabel.push_back(label[k][i]);
                }
                modGracefulLabels.push_back(curModGracefulLabel);
            }
            return true;
        } else {
            return gen(num - 1, 0);
        }
    } else {
        int upper = n - v; // this is for speed purposes
        // optimization: no need to permutate the biggest graph
        if (num == n) {
            upper = 1;
        }
        // we hide used labels in the back of the vector
        // so, we can't use labels after the 'upper' index
        for (int ii = 0; ii < upper; ++ii) {
            int v1_label = vertices[num][ii];

            label[num][v] = v1_label;
            // optimization:
            // compare vertex degrees
            // one in graph, one in Kn; first one shouldn't be bigger
            if (degrees[num][v] > unusedDegreeInKn[v1_label]) {
                continue;
            }

            vector<int> v2_labels;
            bool fits = true;

            for (const auto& parent : parents[num][v]) {
                int v2_label = label[num][parent];
                if (usedEdgeInKn[v1_label][v2_label]) {
                    fits = false;
                    break;
                }
                v2_labels.push_back(v2_label);
            }
            if (!fits) {
                continue;
            }

            for (const auto& v2_label: v2_labels) {
                usedEdgeInKn[v1_label][v2_label] = true;
                usedEdgeInKn[v2_label][v1_label] = true;
                unusedDegreeInKn[v1_label]--;
                unusedDegreeInKn[v2_label]--;
            }
            swap(vertices[num][ii], vertices[num][n - v - 1]);

            if (gen(num, v + 1)) {
                return true;
            }

            swap(vertices[num][ii], vertices[num][n - v - 1]);
            for (const auto& v2_label: v2_labels) {
                usedEdgeInKn[v1_label][v2_label] = false;
                usedEdgeInKn[v2_label][v1_label] = false;
                unusedDegreeInKn[v1_label]++;
                unusedDegreeInKn[v2_label]++;
            }
        }
    }
    return false;
}


int get_vertex_count(const vector<pair<int, int>>& graph_edges) {
    int graph_n = 0;
    for (const auto& edge: graph_edges) {
        graph_n = max(max(graph_n, edge.first + 1), edge.second + 1);
    }
    return graph_n;
}


int get_min_degree(const vector<pair<int, int>>& graph_edges) {
    int graph_n = get_vertex_count(graph_edges);
    vector<int> graph_degrees;
    for (int v = 0; v < graph_n; ++v) {
        graph_degrees.push_back(0);
    }
    for (const auto& edge: graph_edges) {
        graph_degrees[edge.first]++;
        graph_degrees[edge.second]++;
    }
    sort(begin(graph_degrees), end(graph_degrees));
    return graph_degrees[0];
}

int get_second_min_degree(const vector<pair<int, int>>& graph_edges) {
    int graph_n = get_vertex_count(graph_edges);
    vector<int> graph_degrees;
    for (int v = 0; v < graph_n; ++v) {
        graph_degrees.push_back(0);
    }
    for (const auto& edge: graph_edges) {
        graph_degrees[edge.first]++;
        graph_degrees[edge.second]++;
    }
    sort(begin(graph_degrees), end(graph_degrees));
    return graph_degrees[1];
}

int get_max_degree(const vector<pair<int, int>>& graph_edges) {
    int graph_n = get_vertex_count(graph_edges);
    vector<int> graph_degrees;
    for (int v = 0; v < graph_n; ++v) {
        graph_degrees.push_back(0);
    }
    for (const auto& edge: graph_edges) {
        graph_degrees[edge.first]++;
        graph_degrees[edge.second]++;
    }
    sort(begin(graph_degrees), end(graph_degrees));
    return graph_degrees[graph_n - 1];
}


int main(int argc, char* argv[])
{
    beginClock = clock();
    prevClock = beginClock;
    srand(time(NULL));

    n = stoi(argv[1]);
    int start_num = 0;
    if (argc >= 3) {
        start_num = stoi(argv[2]);
    }
    int seed = 42;
    if (argc >= 4) {
        seed = stoi(argv[3]);
    }

    auto rng = default_random_engine{};
    rng.seed(seed);

    // so, we have Kn, we expect graphs with n-1 and less edges

    for (int edge_cnt = 1; edge_cnt <= n - 1; edge_cnt++) {
        ifstream ifs(string("./ge") + to_string(edge_cnt) + "c.g6");
        while (true) {
            vector<pair<int, int>> graph_edges;
            if (!decode_graph6(ifs, graph_edges)) {
                break;
            }
            int graph_n = get_vertex_count(graph_edges);

            // remove "smooth" graphs
            if (get_min_degree(graph_edges) > 1) {
                continue;
            }

            // some subsets of examples

            // if (edge_cnt >= 9) {
            //     if (graph_n > edge_cnt - 3) {
            //         continue;
            //     }
            //     // if (get_max_degree(graph_edges) > 4) {
            //     //     continue;
            //     // }
            // }
            // if (edge_cnt >= 11) {
            //     if (graph_n > 7) {
            //         continue;
            //     }
            // } else if (edge_cnt >= 8) {
            //     if (graph_n > 6) {
            //         continue;
            //     }
            // }

            // if ((edge_cnt >= 4) && get_second_min_degree(graph_edges) == 1) {
            //     continue;
            // }

            int graph_idx = totalAmounts[edge_cnt + 1];
            totalAmounts[edge_cnt + 1]++;
            graphsDatabaseN[edge_cnt + 1][graph_idx] = graph_n;
            // optimization: reversing vertex numbers increases bruteforce speed
            // TODO: possible optimization here
            // relabel the vertices even better
            // so that they have even better order
            // e. g., choose the order which prioritizes centroids
            // or heavy vertices
            for (int edge_idx = 0; edge_idx < graph_edges.size(); ++edge_idx) {
                graphsDatabaseEdges[edge_cnt + 1][graph_idx][edge_idx][0] =
                    graph_n - 1 - graph_edges[edge_idx].first;
                graphsDatabaseEdges[edge_cnt + 1][graph_idx][edge_idx][1] =
                    graph_n - 1 - graph_edges[edge_idx].second;
            }
        }
    }

    cout << "n=" << n << endl;
    for (int i = 1; i < n; i++) {
        cout << "# " << i << "-edged-graphs: " << totalAmounts[i + 1] << endl;
    }
    // return 0;

    for (int num = 0; num <= n; ++num) {
        profile[num] = 0;
    }

    int profile_cnt = -1;
    while (true) {
        profile_cnt += 1;
        if (profile_cnt >= start_num) {
        // if (profile[n] >= start_num) {
            for (int num = 2; num <= n; ++num) {
                // reset graphs
                int graph_idx = profile[num];
                ns[num] = graphsDatabaseN[num][graph_idx];
                for (int v = 0; v < ns[num]; ++v) {
                    parents[num][v].clear();
                    degrees[num][v] = 0;
                }
                for (int edge_idx = 0; edge_idx < num - 1; ++edge_idx) {
                    int v1 = graphsDatabaseEdges[num][graph_idx][edge_idx][0];
                    int v2 = graphsDatabaseEdges[num][graph_idx][edge_idx][1];
                    if (v1 < v2) {
                        parents[num][v2].push_back(v1);
                    } else {
                        parents[num][v1].push_back(v2);
                    }
                    degrees[num][v1]++;
                    degrees[num][v2]++;
                }
            }

            res = false;
            while (!res) {
                depth = 0;

                cnt = 0;
                modGracefulLabels.clear();

                for (int i = 0; i <= n; ++i) {
                    if (i < n) {
                        unusedDegreeInKn[i] = n - 1;
                    }
                    for (int k = 0; k < n; ++k) {
                        usedEdgeInKn[i][k] = false;
                    }
                    // optimization: shuffle vertices labels
                    vector<int> labels;
                    for (int k = 0; k < n; ++k) {
                        labels.push_back(k);
                    }
                    shuffle(begin(labels), end(labels), rng);
                    for (int k = 0; k < n; ++k) {
                        vertices[i][k] = labels[k];
                    }
                }

                for (int num = 1; num <= n; ++num) {
                    for (int i = 0; i <= n; ++i) {
                        usedEdge[num][i] = false;
                    }
                    usedEdge[num][0] = true;
                }

                res = gen(n, 0);
            }

            if (cnt == 0) {
                cout << "CNT==0" << endl;
                cerr << "CNT==0" << endl;
            }

            if (profile_cnt % 1000 == 0) {
                cerr << "profile: ";
                cerr << profile_cnt << "; ";
                for (int j = 2; j <= n; ++j)
                    cerr << j << ":" << profile[j] << " ";
                cerr << endl;
            }

            if (cnt == 0) {
            // if (true) {
                // print out all labelings
                cout << "profile: ";
                cout << profile_cnt << "; ";
                for (int j = 2; j <= n; ++j)
                    cout << j << ":" << profile[j] << " ";
                cout << endl;

                cout << "trees:" << endl;
                for (int num = n; num >= 2; --num) {
                    int graph_idx = profile[num];
                    for (int edge_idx = 0; edge_idx < num - 1; ++edge_idx) {
                        int v1 = graphsDatabaseEdges[num][graph_idx][edge_idx][0];
                        int v2 = graphsDatabaseEdges[num][graph_idx][edge_idx][1];
                        cout <<  v2 << "->" << v1 << "; ";
                    }
                    for (int j = 0; j < graphsDatabaseN[num][graph_idx]; ++j) {
                        cout << modGracefulLabels[num - 2][j] << " ";
                    }
                    cout << endl;
                }
                cout << "cnt=" << cnt << ";" << endl;

                cout << endl;
                cout << endl;
            }
        }
        int idx = 2;
        while ((idx <= n) && (profile[idx] == totalAmounts[idx] - 1)) {
            profile[idx] = 0;
            ++idx;
        }
        if (idx > n) {
            break;
        }
        ++profile[idx];

        if (profile[n] < start_num) {
            continue;
        }

        if (idx == n) {
            cerr << profile[n] << ",";
        }

    }

    endClock = clock();
    double elapsed_secs = double(endClock - beginClock) / CLOCKS_PER_SEC;

    cerr << "Time: " << elapsed_secs  << "s" << endl;
    cout << "Time: " << elapsed_secs  << "s" << endl;
    cerr << "the end" << endl;
    return 0;
}
