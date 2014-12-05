/*
 * kmeansDNA.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: siyuanzhou
 */

#include "kmeansDNA.h"
#include <iostream>
#include <fstream>
#include <random>
#include <numeric>
#include <cmath>
#include <limits>

using namespace std;


void read_dnas(char *input, vector<string> &dnas) {
	ifstream fin(input);
	if (!fin) {
		cout<<"Can not open file "<<input<<endl;
		exit(2);
	}
	string dna;
	while (fin>>dna){
		dnas.push_back(dna);
	}
	fin.close();
}

vector<string> select_centers(const vector<string> &dnas, int k) {
	vector<string> centers(k);
	random_device rd;
	default_random_engine generator(rd());

	for (int i = 0; i < k; ++i)
		centers[i] = dnas[i];
	for (int i = k; i < dnas.size(); ++i) {
		uniform_int_distribution<int> distribution(0, i);
		int r = distribution(generator);
		if (r < k) {
			centers[r] = dnas[i];
		}
	}
	return centers;
}

int get_idx(char base) {
	if (base == 'A')
		return 0;
	else if (base == 'C')
		return 1;
	else if (base == 'G')
		return 2;
	else
		return 3;
}

char get_base(int idx) {
	if (idx == 0)
		return 'A';
	else if (idx == 1)
		return 'C';
	else if (idx == 2)
		return 'G';
	else
		return 'T';
}


void cal_centers(const vector<string> &dnas, const vector<int> &belonging,
			vector<string> &centers) {
	vector<vector<Bases> > cnt(centers.size(),
					vector<Bases>(dnas[0].size()));
	for (int i = 0; i < dnas.size(); ++i) {
		for (int j = 0; j < dnas[i].size(); ++j) {
			int idx = get_idx(dnas[i][j]);
			cnt[belonging[i]][j].increase(idx);
		}
	}
	for (int i = 0; i < centers.size(); ++i) {
		string center;
		for (int j = 0; j < dnas[0].size(); ++j) {
			center += get_base(cnt[i][j].get_max());
		}
		centers[i] = center;
	}
}

int cal_dis(const string &dna, const string &center) {
	int diff = 0;

	for (int i = 0; i < dna.size(); ++i) {
		if (dna[i] != center[i])
			diff++;
	}
	return diff;
}

void output_clusters(char *output, const vector<string> &dnas,
				const vector<int> &belonging) {
	ofstream fout(output);

	for (int i = 0; i < dnas.size(); ++i) {
		fout<<dnas[i]<<" "<<belonging[i]<<endl;
	}
	fout.close();
}
