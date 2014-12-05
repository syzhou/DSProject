/*
 * kmeansDNAParallel.cpp
 *
 *  Created on: Dec 4, 2014
 *      Author: siyuanzhou
 */

#include "kmeansDNA.h"
#include "mpi.h"
#include <limits>
#include <iostream>

using namespace std;

const int DEFAULT_TAG = 1;
const int MASTER_RANK = 0;


static vector<pair<int, int> > devideTasks(int sizeCluster, int numPoints);
static bool isMaster(int rank);
static vector<int> kmeans(const vector<string> &dnas, int k);


int main(int argc, char **argv) {
    int rank;
    //Initialize MPI and get task ID and size
    MPI::Init();
    rank = MPI::COMM_WORLD.Get_rank();
	if (argc != 4) {
		cout<<"Usage: <input file> <output file> <k>"<<endl;
		MPI::COMM_WORLD.Abort(1);
	}

	vector<string> dnas;
	read_dnas(argv[1], dnas);
	int k = atoi(argv[3]);
	vector<int> belonging = kmeans(dnas, k);


	if (isMaster(rank)) {
		output_clusters(argv[2], dnas, belonging);
	}
	MPI::Finalize();
	return 0;
}


static vector<pair<int, int> > devideTasks(int sizeCluster, int numPoints) {
	int numPointsPerTask = (int) (numPoints / sizeCluster);
	int currIndex = 0;
	vector<pair<int, int> > tasks(sizeCluster);
	for (int i = 0 ; i < sizeCluster; i++) {
		if (i == sizeCluster - 1) {
			tasks[i] = make_pair(currIndex, numPoints);
		} else {
			tasks[i] = make_pair(currIndex, currIndex + numPointsPerTask);
			currIndex += numPointsPerTask;
		}
	}
	return tasks;
}

static bool isMaster(int rank) {
	return rank == MASTER_RANK;
}

static vector<int> kmeans(const vector<string> &dnas, int k) {
    int rank = MPI::COMM_WORLD.Get_rank();
    int size = MPI::COMM_WORLD.Get_size();
    int dnaLength = dnas[0].length();

    //Define the string for dna as MPI type
    MPI::Datatype dnaType = MPI::CHAR.Create_contiguous(dnaLength + 1);
    dnaType.Commit();

    //Define the class Bases as MPI type
    MPI::Datatype basesType, oldTypes[1];
    MPI::Aint offsets[1];
    int blockCounts[1];
    oldTypes[0] = MPI::INT;
    offsets[0] = 0;
    blockCounts[0] = 4;

    basesType = MPI::Datatype::Create_struct(1, blockCounts, offsets, oldTypes);
    basesType.Commit();

	vector<pair<int, int> > tasks(devideTasks(size, dnas.size()));
	int startIndex = tasks[rank].first;
	int endIndex = tasks[rank].second;
	vector<string> prevCenters(k, string(dnaLength, 'A'));
	vector<string> centers(k);
    if (isMaster(rank)) {
    	centers = select_centers(dnas, k);
    }

	vector<int> belonging(dnas.size());
	char *buffer = new char[dnaLength + 1];
	while (true) {
		//Send the new centers to all tasks
		if (isMaster(rank)) {
			for (int i = 0; i < size; i++) {
				if (!isMaster(i)) {
					for (int j = 0; j < k; j++) {
						MPI::COMM_WORLD.Send(centers[j].c_str(), 1, dnaType, i, DEFAULT_TAG);
					}
				}
			}
		}
		//All tasks update their centers
		if (!isMaster(rank)) {
			for (int j = 0; j < k; j++) {
				MPI::COMM_WORLD.Recv(buffer, 1, dnaType, MASTER_RANK, DEFAULT_TAG);
				centers[j] = buffer;
			}
		}
		MPI::COMM_WORLD.Barrier();

		if (equal(centers.begin(), centers.end(), prevCenters.begin())) {
			break;
		}
		prevCenters = centers;
		vector<vector<Bases> > sumClusters(k, vector<Bases>(dnaLength));
		for (int i = startIndex; i < endIndex; i++) {
			int min_dis = numeric_limits<int>::max();
			int belong;
			for (int j = 0; j < k; ++j) {
				int dis = cal_dis(dnas[i], centers[j]);
				if (dis < min_dis) {
					min_dis = dis;
					belong = j;
				}
			}
			belonging[i] = belong;
			for (int j = 0; j < dnaLength; j++) {
				int idx = get_idx(dnas[i][j]);
				sumClusters[belong][j].increase(idx);
			}
		}
		if (!isMaster(rank)) {
			for (int i = 0; i < k; i++) {
				MPI::COMM_WORLD.Send(&(sumClusters[i][0]), dnaLength, basesType, MASTER_RANK, DEFAULT_TAG);
			}
		}
		if (isMaster(rank)) {
			vector<vector<Bases> > aggrSumClusters(k, vector<Bases>(dnaLength));
			for (int i = 0; i < k; i++) {
				for (int j = 0; j < dnaLength; j++) {
					aggrSumClusters[i][j] += sumClusters[i][j];
					//cout << "cluster " << i << " dna " << j << " count " << aggrSumClusters[i][j].get_max() << endl;
				}
			}
			for (int i = 0; i < size; i++) {
				if (isMaster(i)) {
					continue;
				}
				for (int j = 0; j < k; j++) {
					MPI::COMM_WORLD.Recv(&(sumClusters[j][0]), dnaLength, basesType, i, DEFAULT_TAG);
				}

				for (int j = 0; j < k; j++) {
					for (int n = 0; n < dnaLength; n++) {
						aggrSumClusters[j][n] += sumClusters[j][n];
					}
				}
			}
			for (int i = 0; i < centers.size(); ++i) {
				string center;
				for (int j = 0; j < dnas[0].size(); ++j) {
					center += get_base(aggrSumClusters[i][j].get_max());
				}
				centers[i] = center;
			}
		}
		MPI::COMM_WORLD.Barrier();
	}

	delete[] buffer;
	if (!isMaster(rank)) {
		MPI::COMM_WORLD.Send(&belonging[startIndex], endIndex - startIndex, MPI::INT, MASTER_RANK, DEFAULT_TAG);
	}
	if (isMaster(rank)) {
		for (int i = 0; i < size; i++) {
			if (isMaster(i)) {
				continue;
			}
			int taskStart = tasks[i].first;
			int taskEnd = tasks[i].second;
			MPI::COMM_WORLD.Recv(&belonging[taskStart], taskEnd - taskStart, MPI::INT, i, DEFAULT_TAG);
		}
	}

	return belonging;
}
