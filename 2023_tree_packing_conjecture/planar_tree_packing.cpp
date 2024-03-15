/*
 * File:   planar_tree_packing.cpp
 * Author: Nikolay Ulyanov
 *
 * Created on 22 September 2023
 */

#include <iostream>
#include <cstdio>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <map>
#include <string>
#include "freetree.h"

using namespace std;

const int maxn = 30;
int n;
vector<int> gr[maxn][maxn];
int label[maxn][maxn];
vector<vector<int>> modGracefulLabels;
vector<int> curModGracefulLabel;
int parent[maxn][maxn];
int deg[maxn][maxn];
int link[maxn][maxn];
bool isCopy[maxn][maxn];
int revLinkCount[maxn][maxn];
bool usedEdge[maxn][maxn];

bool usedEdgeInKnn[maxn][maxn];

int diam[maxn];
int branch[maxn][maxn];
int distTo0[maxn][maxn];
int diams[maxn][maxn];

int vertices[maxn][maxn];

bool res;

int cnt;

clock_t beginClock, curClock, prevClock, endClock;

void countDiam(int num) {
    for (int i = 1; i < num; ++i) {
        if (parent[num][i] == 0) {
            branch[num][i] = i;
            diams[num][i] = 1;
        } else {
            branch[num][i] = branch[num][parent[num][i]];
        }
    }
    distTo0[num][0] = 0;
    for (int i = 1; i < num; ++i) {
        distTo0[num][i] = distTo0[num][parent[num][i]] + 1;
    }
    for (int i = 1; i < num; ++i) {
        diams[num][branch[num][i]] = max(diams[num][branch[num][i]], distTo0[num][i]);
    }
    diam[num] = diams[num][gr[num][0][0]] + diams[num][gr[num][0][1]];
}

bool checkTrees(int num, int v1, int v2) {
    if (deg[num][v1] != deg[num][v2]) {
        return false;
    }
    for (int i = 1; i < (int) gr[num][v1].size(); ++i) {
        if (!(checkTrees(num, gr[num][v1][i], gr[num][v2][i]))) {
            return false;
        }
    }
    return true;
}

void markAsCopy(int num, int v1) {
    isCopy[num][v1] = true;
    for (int i = 1; i < (int) gr[num][v1].size(); ++i) {
        markAsCopy(num, gr[num][v1][i]);
    }
}

void buildLinks(int num) {
    for (int i = 0; i < num; ++i) {
        link[num][i] = -1;
        isCopy[num][i] = false;
    }
    for (int center = 0; center < num; ++center) {
        for (int ii = 0; ii < (int) gr[num][center].size()-1; ++ii) {
            int neib1 = gr[num][center][ii];
            if (neib1 < center) {
                continue;
            }
            int neib2 = gr[num][center][ii + 1];
            if (checkTrees(num, neib1, neib2)) {
                link[num][neib2] = neib1;
                if (!isCopy[num][neib2]) {
                    markAsCopy(num, neib2);
                }
            }
        }
    }

    for (int i = 0; i < num; ++i) {
        revLinkCount[num][i] = 0;
    }
    for (int i = num - 1; i >= 0; --i) {
        if (link[num][i] != -1) {
            revLinkCount[num][link[num][i]] = revLinkCount[num][i] + 1;
        }
    }
}

bool gen(int num, int v) {
    if (v == num) {
        if (num == 2) {
            ++cnt;
            for (int k = n; k >= 2; --k) {
                curModGracefulLabel.clear();
                for (int i = 0; i < k; ++i) {
                    curModGracefulLabel.push_back(label[k][i]);
                }
                modGracefulLabels.push_back(curModGracefulLabel);
            }
            return true;
        } else {
            return gen(num - 1, 0);
        }
    } else {
        int upper = n - v;
        for (int ii = 0; ii < upper; ++ii) {
            int i = vertices[num][ii];

            label[num][v] = i;
            int minVal = -1;

            int diff;

            if (parent[num][v] != -1) {
                int v1 = label[num][v];
                int v2 = label[num][parent[num][v]];
                if (v2 < v1) {
                    swap(v1, v2);
                }
                bool isPlanar = true;
                for (int prev_v = 1; prev_v < v; ++prev_v) {
                    int v3 = label[num][prev_v];
                    int v4 = label[num][parent[num][prev_v]];
                    if (v4 < v3) {
                        swap(v3, v4);
                    }
                    if ((v1 == v3) || (v1 == v4)) {
                        continue;
                    }
                    if ((v2 == v3) || (v2 == v4)) {
                        continue;
                    }
                    if (((v1 < v3) && (v3 < v2) && (v2 < v4)) ||
                            ((v3 < v1) && (v1 < v4) && (v4 < v2))) {
                        isPlanar = false;
                        break;
                    }
                }
                if (!isPlanar) {
                    continue;
                }

                diff = abs(abs(label[num][v]) - abs(label[num][parent[num][v]]));
                minVal = min(abs(label[num][v]), abs(label[num][parent[num][v]]));
                if (usedEdgeInKnn[minVal][diff]) {
                    continue;
                }

                if (link[num][v] != -1) {
                    // code specific only for planar tree packing
                    label[num][v] = i;
                    if (label[num][v] > label[num][link[num][v]]) {
                        continue;
                    }
                }
                // if (edgeVal < revLinkCount[num][v]) {
                    // continue;
                // }

                usedEdgeInKnn[minVal][diff] = true;
            }
            swap(vertices[num][ii], vertices[num][n - v - 1]);

            if (gen(num, v + 1)) {
                return true;
            }

            swap(vertices[num][ii], vertices[num][n - v - 1]);
            if (parent[num][v] != -1) {
                usedEdgeInKnn[minVal][diff] = false;
            }
        }
    }
    return false;
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
    for (int i = 1; i <= n; ++i) {
        freeTree(i);
    }

    cout << "n=" << n << "; #trees: " << totalAmounts[n] << endl;
    for (int ttt = 0; ttt < totalAmounts[n]; ++ttt) {
        if (ttt < start_num) {
            continue;
        }

        if (ttt % 1 == 0) {
            cerr << ttt << ",";
        }

        // select one random biggest tree
        for (int v = 0; v < n; ++v) {
            parent[n][v] = treeDatabase[n][ttt][v + 1];
        }

        // other (smaller) trees are stars
        // for (int num = 1; num < n; ++num) {
        //     parent[num][0] = -1;
        //     for (int j = 1; j < num; ++j) {
        //         parent[num][j] = 0;
        //     }
        // }

        // other (smaller) trees are subtrees of the biggest tree
        for (int num = 1; num < n; ++num) {
            for (int v = 0; v < num; ++v) {
                parent[num][v] = treeDatabase[n][ttt][v + 1];
            }
        }

        for (int num = 1; num <= n; ++num) {
            for (int i = 0; i < n; ++i)
                gr[num][i].clear();

            for (int i = 1; i < num; ++i) {
                gr[num][parent[num][i]].push_back(i);
                gr[num][i].push_back(parent[num][i]);
            }

            for (int i = 0; i < num; ++i) {
                deg[num][i] = 0;
            }

            for (int i = 1; i < num; ++i) {
                ++deg[num][parent[num][i]];
                ++deg[num][i];
            }
            buildLinks(num);
        }


        for (int t = 0; t < 1 << (n - 1); ++t) {
            cnt = 0;
            modGracefulLabels.clear();

            for (int i = 0; i <= n; ++i) {
                for (int k = 0; k <= n; ++k) {
                    usedEdgeInKnn[i][k] = false;
                    bool isPlus = ((t >> (k - 1)) % 2 == 0);
                    vertices[i][k] = k;
                    if (!isPlus)
                        vertices[i][k] *= -1;
                }
            }

            for (int num = 1; num <= n; ++num) {
                for (int i = 0; i <= n; ++i) {
                    usedEdge[num][i] = false;
                }
                usedEdge[num][0] = true;
            }

            res = gen(n, 0);
            if (res) {
                break;
            }
        }

        if (cnt == 0) {
            cout << "CNT==0" << endl;
            cerr << "CNT==0" << endl;
        }
        // print out all labelings
        cout << "trees:" << endl;
        for (int num = n; num >= 2; --num) {
            for (int j = 1; j < num; ++j) {
                cout << parent[num][j] << "->" << j << "; ";
            }
            for (int j = 0; j < num; ++j) {
                cout << modGracefulLabels[n - num][j] << " ";
            }
            cout << endl;
        }
        // for (int i = 0; i < modGracefulLabels.size(); ++i) {
        //     cout << "modGraceful labels: ";
        //     for (int j = 0; j < n - i; ++j) {
        //         cout << modGracefulLabels[i][j] << " ";
        //     }
        //     cout << endl;

        //     cout << "|";
        //     cout << endl;
        // }

        // cout << "tree: ";
        // for (int j = 1; j < n; ++j) {
        //     cout << parent[n][j] << "->" << j << "; ";
        // }
        // cout << "trees:" << endl;
        // for (int num = n; num >= 2; --num) {
        //     for (int j = 1; j < num; ++j) {
        //         cout << parent[num][j] << "->" << j << "; ";
        //     }
        //     cout << endl;
        // }
        cout << endl;
        countDiam(n);
        int maxDeg = 0;
        for (int i = 0; i < n; ++i) {
            maxDeg = max(maxDeg, deg[n][i]);
        }
        cout << "diam=" << diam[n] << "; maxDeg=" << maxDeg << endl;
        cout << "cnt=" << cnt << ";" << endl;

        cout << endl;
        cout << endl;
    }

    endClock = clock();
    double elapsed_secs = double(endClock - beginClock) / CLOCKS_PER_SEC;

    cerr << "Time: " << elapsed_secs  << "s" << endl;
    cout << "Time: " << elapsed_secs  << "s" << endl;
    cerr << "the end" << endl;
    return 0;
}
