/*
 * File:   graceful_tree_packing.cpp
 * Author: Nikolay Ulyanov
 *
 * Created on 25 September 2016
 */

// labeling conventions: 
// V \in [0, 1, ..., n-1]
// E \in [1, ..., n-1]
// n < maxn
// Current set of restrictions:
// graceful: f(xy) = f(x)+f(y) (mod n)
// works for n <= ?
// doesn't work for P4

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
vector <int> gr[maxn][maxn];
int label[maxn][maxn];
vector <vector <int> > gracefulLabels;
vector <int> curGracefulLabel;
int parent[maxn][maxn];
int deg[maxn][maxn];
int link[maxn][maxn];
bool isCopy[maxn][maxn];
int revLinkCount[maxn][maxn];
bool usedVertex[maxn][maxn];
bool usedEdge[maxn][maxn];

bool usedEdgeInKnn[maxn][maxn];

int diam[maxn];
int branch[maxn][maxn];
int distTo0[maxn][maxn];
int diams[maxn][maxn];

bool res;

int kol;

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
            ++kol;
            for (int k = n; k >= 2; --k) {
                curGracefulLabel.clear();
                for (int i = 0; i < k; ++i) {
                    curGracefulLabel.push_back(label[k][i]);
                }
                gracefulLabels.push_back(curGracefulLabel);
            }
            return true;
        } else {
            for (int j = 0; j < n; ++j) {
                label[num - 1][0] = j;
                usedVertex[num - 1][j] = true;
                if (gen(num - 1, 1)) {
                    return true;
                }
                usedVertex[num - 1][j] = false;
            }
        }
	} else {
		int edgeVal;

		for (int i = 0; i < n; ++i) {
			if (usedVertex[num][i]) {
				continue;
            } else {
				label[num][v] = i;
                int minVal = -1;

				if (parent[num][v] != -1) {
					edgeVal = abs(label[num][v] - label[num][parent[num][v]]);
                    minVal = min(i, label[num][parent[num][v]]);
					if (edgeVal > num || usedEdge[num][edgeVal] || usedEdgeInKnn[minVal][edgeVal]) {
						continue;
                    }

					/*if (link[num][v] != -1) {
						if (edgeVal > abs(label[num][link[num][v]] - label[num][parent[num][link[num][v]]])) {
							continue;
                        }
                    }
					if (edgeVal < revLinkCount[num][v]) {
						continue;
                    }*/
					usedEdge[num][edgeVal] = true;
                    usedEdgeInKnn[minVal][edgeVal] = true;
				}
				usedVertex[num][label[num][v]] = true;

				if (gen(num, v + 1)) {
					return true;
                }

				usedVertex[num][label[num][v]] = false;
				if (parent[num][v] != -1) {
					usedEdge[num][edgeVal] = false;
                    usedEdgeInKnn[minVal][edgeVal] = false;
                }
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
    for (int i = 1; i <= n; ++i) {
    	freeTree(i);
    }

	cout << "n=" << n << "; #trees: " << totalAmounts[n] << endl;
	for (int ttt = 0; ttt < totalAmounts[n]; ++ttt) {
		int tt = ttt;
		
		if (ttt % 1 == 0) {
			cerr << tt << ",";
        }

		for (int v = 0; v < n; ++v) {
			parent[n][v] = treeDatabase[n][tt][v + 1];
        }

        for (int num = 1; num < n; ++num) {
            parent[num][0] = -1;
            for (int j = 1; j < num; ++j) {
                parent[num][j] = 0;
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

        kol = 0;
        gracefulLabels.clear();
       
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i <= n; ++i) {
                for (int k = 0; k <= n; ++k) {
                    usedEdgeInKnn[i][k] = false;
                }
            }

            for (int num = 1; num <= n; ++num) {
                for (int i = 0; i <= n; ++i) {
                    usedVertex[num][i] = false;
                    usedEdge[num][i] = false;
                }
                usedEdge[num][0] = true;
            }
            
            label[n][0] = j;
            usedVertex[n][j] = true;
            res = gen(n, 1);
            if (res) {
                break;
            }
            usedVertex[n][j] = false;
        }

        if (!kol) {
            cout << "OHSHI~" << endl;
            cerr << "OHSHI~" << endl;
        }
        // print out all labelings
        for (int i = 0; i < gracefulLabels.size(); ++i) {
            cout << "graceful labels: ";
            for (int j = 0; j < n - i; ++j) {
                cout << gracefulLabels[i][j] << " ";
            }
            cout << endl;
            
            cout << "|";
            cout << endl;
        }

        cout << "tree: ";
        for (int j = 1; j < n; ++j) {
            cout << parent[n][j] << "->" << j << "; ";
        }
        cout << endl;
        countDiam(n);
        int maxDeg = 0;
        for (int i = 0; i < n; ++i) {
            maxDeg = max(maxDeg, deg[n][i]);
        }
        cout << "diam=" << diam[n] << "; maxDeg=" << maxDeg << endl;
        cout << "kol=" << kol << ";" << endl;

        cout << endl;
        cout << endl;
	}

	endClock = clock();
	double elapsed_secs = double(endClock - beginClock) / CLOCKS_PER_SEC;

	cerr << "Time: " << elapsed_secs  << "s" << endl;
	cout << "Time: " << elapsed_secs  << "s" << endl;
	cerr << "the end" << endl;
	//system ("pause");
	return 0;
}

