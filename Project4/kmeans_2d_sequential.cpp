#include "kmeans2D.h"
#include <numeric>
#include <iostream>
#include <limits>

using namespace std;


static vector<int> kmeans(const vector<Point> &points, int k);

int main(int argc, char **argv) {
	if (argc != 4) {
		cout<<"Usage: <input file> <output file> <k>"<<endl;
		exit(1);
	}

	vector<Point> points;
	read_points(argv[1], points);
	int k = atoi(argv[3]);
	vector<int> belonging = kmeans(points, k);
	output_clusters(argv[2], points, belonging);

	return 0;
}

static vector<int> kmeans(const vector<Point> &points, int k) {
	vector<Point> centers = select_centers(points, k);
	vector<int> belonging(points.size(), -1);
	bool flag;

	while (true) {
		flag = false;
		for (int i = 0; i < points.size(); ++i) {
			double min_dis = numeric_limits<double>::max();
			int belong;
			for (int j = 0; j < k; ++j) {
				double dis = cal_dis(points[i], centers[j]);
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
		cal_centers(points, belonging, centers);
	}

	return belonging;
}




