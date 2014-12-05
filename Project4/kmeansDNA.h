/*
 * kmeansDNA.h
 *
 *  Created on: Dec 4, 2014
 *      Author: siyuanzhou
 */

#ifndef KMEANSDNA_H_
#define KMEANSDNA_H_

#include <vector>
#include <string>

struct Bases {
private:
	int cnt[4];
public:
	Bases() {
		for (int i = 0; i < 4; ++i)
			cnt[i] = 0;
	}

	void increase(int idx) {
		cnt[idx]++;
	}

	int get_max() {
		int max = 0;
		int idx;

		for (int i = 0; i < 4; ++i) {
			if (cnt[i] > max) {
				max = cnt[i];
				idx = i;
			}
		}
		return idx;
	}

	Bases &operator += (const Bases &other) {
		for (int i = 0; i < 4; i++) {
			cnt[i] += other.cnt[i];
		}
		return *this;
	}
};

void read_dnas(char *input, std::vector<std::string> &dnas);

std::vector<std::string> select_centers(const std::vector<std::string> &dnas, int k);

int get_idx(char base);

char get_base(int idx);

void cal_centers(const std::vector<std::string> &dnas, const std::vector<int> &belonging,
			std::vector<std::string> &centers);

int cal_dis(const std::string &dna, const std::string &center);

void output_clusters(char *output, const std::vector<std::string> &dnas,
				const std::vector<int> &belonging);


#endif /* KMEANSDNA_H_ */
