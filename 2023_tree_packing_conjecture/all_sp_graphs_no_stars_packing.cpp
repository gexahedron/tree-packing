/*
 * File:   all_sp_graphs_no_stars_packing.cpp
 * Author: Nikolay Ulyanov
 *
 * Created on 25 September 2023
 */

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

const int maxn = 10;  // numbers are slightly arbitrary

static int totalAmounts[maxn] = {0};
short graphsDatabaseEdges[maxn][10000][maxn][2];
short graphsDatabaseN[maxn][10000];

int n;

int label[maxn][maxn];
vector<vector<int>> modGracefulLabels;
vector<int> curModGracefulLabel;

// chosen graphs for packing
vector<int> parents[maxn][maxn];
int ns[maxn];
bool usedEdge[maxn][maxn];

bool usedEdgeInKn[maxn][maxn];

int vertices[maxn][maxn];

bool res;

int cnt;

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
    if (v == ns[num]) {
        if (num == 4) { // ignoring T2 and T1
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
            int i = vertices[num][ii];

            label[num][v] = i;
            // TODO: possible optimization here
            // compare vertex degrees
            // one in graph, one in Kn

            int v1_label = label[num][v];

            vector<pair<int, int>> knEdges;
            bool fits = true;

            for (const auto& parent : parents[num][v]) {
                int v2_label = label[num][parent];
                if (usedEdgeInKn[v1_label][v2_label]) {
                    fits = false;
                    break;
                }
                knEdges.push_back(make_pair(v1_label, v2_label));
            }
            if (!fits) {
                continue;
            }

            for (const auto& knEdge: knEdges) {
                usedEdgeInKn[knEdge.first][knEdge.second] = true;
                usedEdgeInKn[knEdge.second][knEdge.first] = true;
            }
            swap(vertices[num][ii], vertices[num][n - v - 1]);

            if (gen(num, v + 1)) {
                return true;
            }

            swap(vertices[num][ii], vertices[num][n - v - 1]);
            for (const auto& knEdge: knEdges) {
                usedEdgeInKn[knEdge.first][knEdge.second] = false;
                usedEdgeInKn[knEdge.second][knEdge.first] = false;
            }
        }
    }
    return false;
}

bool is_star(const vector<pair<int, int>>& graph_edges) {
    int center_vertex = -1;
    if (graph_edges[0].first == graph_edges[1].first) {
        center_vertex = graph_edges[0].first;
    }
    if (graph_edges[0].second == graph_edges[1].first) {
        center_vertex = graph_edges[0].second;
    }
    if (graph_edges[0].first == graph_edges[1].second) {
        center_vertex = graph_edges[0].first;
    }
    if (graph_edges[0].second == graph_edges[1].second) {
        center_vertex = graph_edges[0].second;
    }
    if (center_vertex == -1) {
        return false;
    }
    for (const auto& edge: graph_edges) {
        if ((edge.first != center_vertex) && (edge.second != center_vertex)) {
            return false;
        }
    }
    return true;
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

    // so, we have Kn, we expect graphs with n-1 and less edges

    for (int i = 2; i <= n + 1; i++) {
        ifstream ifs(string("./planar_conn.") + to_string(i) + ".g6");
        while (true) {
            vector<pair<int, int>> graph_edges;
            if (!decode_graph6(ifs, graph_edges)) {
                break;
            }
            if (graph_edges.size() >= n) {
                continue;
            }

            // count treewidth
            Light_Graph * gs = 0;
            Tree_Decomposition * TD = 0;
            gs = new Light_Graph;
            gs->Load_TW (i, graph_edges);
            TD = new Tie_Break_Min_Fill_Triangulation_Tree_Decomposition (gs, new Neighbor_Future_Degree_Tie_Break_Policy(gs->N()));
            if (TD->Get_W() > 2) {  // not a series-parallel graph
                continue;
            }
            delete TD;
		    delete gs;

            // check that it's not a star
            if ((graph_edges.size() > 2) && is_star(graph_edges)) {
                continue;
            }

            int edge_cnt = graph_edges.size();
            int graph_idx = totalAmounts[edge_cnt + 1];
            totalAmounts[edge_cnt + 1]++;
            graphsDatabaseN[edge_cnt + 1][graph_idx] = i;
            for (int edge_idx = 0; edge_idx < graph_edges.size(); ++edge_idx) {
                graphsDatabaseEdges[edge_cnt + 1][graph_idx][edge_idx][0] =
                    graph_edges[edge_idx].first;
                graphsDatabaseEdges[edge_cnt + 1][graph_idx][edge_idx][1] =
                    graph_edges[edge_idx].second;
            }
        }
    }

    cout << "n=" << n << endl;
    for (int i = 1; i < n; i++) {
        cout << "# " << i << "-edged-graphs: " << totalAmounts[i + 1] << endl;
    }

    for (int num = 0; num <= n; ++num) {
        profile[num] = 0;
    }

    while (true) {
        if (profile[n] >= start_num) {
            for (int num = 2; num <= n; ++num) {
                // reset graphs
                int graph_idx = profile[num];
                ns[num] = graphsDatabaseN[num][graph_idx];
                for (int v = 0; v < ns[num]; ++v) {
                    parents[num][v].clear();
                }
                for (int edge_idx = 0; edge_idx < num - 1; ++edge_idx) {
                    int v1 = graphsDatabaseEdges[num][graph_idx][edge_idx][0];
                    int v2 = graphsDatabaseEdges[num][graph_idx][edge_idx][1];
                    if (v1 < v2) {
                        parents[num][v2].push_back(v1);
                    } else {
                        parents[num][v1].push_back(v2);
                    }
                }
            }

            cnt = 0;
            modGracefulLabels.clear();

            for (int i = 0; i <= n; ++i) {
                for (int k = 0; k <= n; ++k) {
                    usedEdgeInKn[i][k] = false;
                    vertices[i][k] = k;
                }
            }

            for (int num = 1; num <= n; ++num) {
                for (int i = 0; i <= n; ++i) {
                    usedEdge[num][i] = false;
                }
                usedEdge[num][0] = true;
            }

            res = gen(n, 0);

            if (cnt == 0) {
                cout << "CNT==0" << endl;
                cerr << "CNT==0" << endl;
            }

            // print out all labelings
            cout << "profile: ";
            for (int j = 2; j <= n; ++j)
                cout << j << ":" << profile[j] << " ";
            cout << endl;

            cout << "trees:" << endl;
            for (int num = n; num > 3; --num) {
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
