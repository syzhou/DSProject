/*
 * kmeans2D.cpp
 *
 *  Created on: Dec 3, 2014
 *      Author: Siyuan Zhou
 *      		Zichang Feng
 */

#include "kmeans2D.h"
#include <cmath>
#include <fstream>
#include <random>
#include <iostream>
#include <cstdlib>


using namespace std;

void read_points(char *input, vector<Point> &points) {
	ifstream fin(input);
	if (!fin) {
		cout<<"Can not open file "<<input<<endl;
		exit(2);
	}
	double x, y;
	char c;

	while (fin>>x>>c>>y){
		points.push_back(Point(x, y));
	}
	fin.close();
}

vector<Point> select_centers(const vector<Point> &points, int k) {
	vector<Point> centers(k);
	random_device rd;
	default_random_engine generator(rd());

	for (int i = 0; i < k; ++i)
		centers[i] = points[i];
	for (int i = k; i < points.size(); ++i) {
		uniform_int_distribution<int> distribution(0, i);
		int r = distribution(generator);
		if (r < k) {
			centers[r] = points[i];
		}
	}
	return centers;
}

void cal_centers(const vector<Point> &points, const vector<int> &belonging,
			vector<Point> &centers) {
	vector<int> cnt(centers.size(), 0);

	for (int i = 0; i < centers.size(); ++i)
		centers[i].x = centers[i].y = 0;
	for (int i = 0; i < belonging.size(); ++i) {
		centers[belonging[i]].x += points[i].x;
		centers[belonging[i]].y += points[i].y;
		cnt[belonging[i]]++;
	}

	for (int i = 0; i < centers.size(); ++i) {
		centers[i].x /= cnt[i];
		centers[i].y /= cnt[i];
	}
}

double cal_dis(const Point &point, const Point &center) {
	double dis = 0;

	dis += (point.x - center.x) * (point.x - center.x);
	dis += (point.y - center.y) * (point.y - center.y);
	return sqrt(dis);
}

void output_clusters(char *output, const vector<Point> &points,
				const vector<int> &belonging) {
	ofstream fout(output);

	for (int i = 0; i < points.size(); ++i) {
		fout<<points[i].x<<" "<<points[i].y<<" "<<belonging[i]<<endl;
	}
	fout.close();
}



