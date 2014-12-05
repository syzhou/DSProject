#include "kmeansDNA.h"
#include <iostream>
#include <limits>

using namespace std;

static vector<int> kmeans(const vector<string> &dnas, int k);

int main(int argc, char **argv) {
	if (argc != 4) {
		cout<<"Usage: <input file> <output file> <k>"<<endl;
		exit(1);
	}

	vector<string> dnas;
	read_dnas(argv[1], dnas);
	int k = atoi(argv[3]);
	vector<int> belonging = kmeans(dnas, k);
	output_clusters(argv[2], dnas, belonging);

	return 0;
}

static vector<int> kmeans(const vector<string> &dnas, int k) {
	vector<string> centers = select_centers(dnas, k);
	vector<int> belonging(dnas.size(), -1);
	bool flag;

	while (true) {
		flag = false;
		for (int i = 0; i < dnas.size(); ++i) {
			int min_dis = numeric_limits<int>::max();
			int belong;
			for (int j = 0; j < k; ++j) {
				int dis = cal_dis(dnas[i], centers[j]);
				if (dis < min_dis) {
					min_dis = dis;
					belong = j;
				}
			}
			if (belonging[i] != belong) {
				flag = true;
				belonging[i] = belong;
			}
		}
		if (flag == false)
			break;
		cal_centers(dnas, belonging, centers);
	}

	return belonging;
}
