/*
 * kmeans2DParallel.cpp
 *
 *  Created on: Dec 3, 2014
 *      Author: siyuanzhou
 */

#include "kmeans2D.h"
#include "mpi.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <limits>


using namespace std;

const int DEFAULT_TAG = 1;
const int MASTER_RANK = 0;

/*
 * Divide points to sizeCluster parts.
 * Return the task division with each pair being starting(inclusive) and
 * end(exclusive) point index.
 */
static vector<pair<int, int> > devideTasks(int sizeCluster, int numPoints);
static bool isMaster(int rank);
static vector<int> kmeans(const vector<Point> &points, int k);

int main(int argc, char **argv) {
    int rank;
    //Initialize MPI and get task ID and size
    MPI::Init();
    rank = MPI::COMM_WORLD.Get_rank();
	if (argc != 4) {
		cout << "Usage: <input file> <output file> <k>" << endl;
		MPI::COMM_WORLD.Abort(1);
	}

	vector<Point> points;
	read_points(argv[1], points);
	int k = atoi(argv[3]);


	vector<int> belonging = kmeans(points, k);

	if (isMaster(rank)) {
		output_clusters(argv[2], points, belonging);
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

static vector<int> kmeans(const vector<Point> &points, int k) {
    int rank = MPI::COMM_WORLD.Get_rank();
    int size = MPI::COMM_WORLD.Get_size();
    MPI::Datatype pointType, oldTypes[1];
    MPI::Aint offsets[1];
    int blockCounts[1];
    oldTypes[0] = MPI::DOUBLE;
    offsets[0] = 0;
    blockCounts[0] = 2;

    pointType = MPI::Datatype::Create_struct(1, blockCounts, offsets, oldTypes);
    pointType.Commit();

	vector<pair<int, int> > tasks(devideTasks(size, points.size()));
	int startIndex = tasks[rank].first;
	int endIndex = tasks[rank].second;
	vector<Point> prevCenters(k, Point(numeric_limits<double>::max(), numeric_limits<double>::max()));
	vector<Point> centers(k);
	if (isMaster(rank)) {
		centers = select_centers(points, k);
	}
	vector<int> belonging(points.size());
	while (true) {
		/*
		//Tell the tasks whether to stop
		if (isMaster(rank)) {
			for (int i = 0; i < size; i++) {
				MPI::COMM_WORLD.Send(&stop, 1, MPI::BOOL, i, DEFAULT_TAG);
			}
		}
		MPI::COMM_WORLD.Recv(&stop, 1, MPI::BOOL, MASTER_RANK, DEFAULT_TAG);
		if (stop) {
			break;
		}
		*/
		//Send the new centers to all tasks
		if (isMaster(rank)) {
			for (int i = 0; i < size; i++) {
				if (!isMaster(i)) {
					MPI::COMM_WORLD.Send(&centers[0], k, pointType, i, DEFAULT_TAG);
				}
			}
		}
		//All tasks update their centers
		if (!isMaster(rank)) {
			MPI::COMM_WORLD.Recv(&centers[0], k, pointType, MASTER_RANK, DEFAULT_TAG);
		}

		MPI::COMM_WORLD.Barrier();

		if (equal(centers.begin(), centers.end(), prevCenters.begin())) {
			break;
		}
		prevCenters = centers;
		vector<Point> sumClusters(k, Point(0, 0));
		vector<int> pointsPerCluster(k, 0);
		for (int i = startIndex; i < endIndex; i++) {
			double min_dis = numeric_limits<double>::max();
			int belong;
			const Point &point = points[i];
			for (int j = 0; j < k; ++j) {
				double dis = cal_dis(point, centers[j]);
				if (dis < min_dis) {
					min_dis = dis;
					belong = j;
				}
			}
			belonging[i] = belong;
			sumClusters[belong] += point;
			pointsPerCluster[belong]++;
		}
		if (!isMaster(rank)) {
			MPI::COMM_WORLD.Send(&sumClusters[0], k, pointType, MASTER_RANK, DEFAULT_TAG);
			MPI::COMM_WORLD.Send(&pointsPerCluster[0], k, MPI::INT, MASTER_RANK, DEFAULT_TAG);
		}
		if (isMaster(rank)) {
			vector<Point> aggrSumClusters(k, Point(0, 0));
			vector<int> aggrPointsPerCluster(k, 0);
			for (int j = 0; j < k; j++) {
				aggrSumClusters[j] += sumClusters[j];
				aggrPointsPerCluster[j] += pointsPerCluster[j];
			}
			for (int i = 0; i < size; i++) {
				if (isMaster(i)) {
					continue;
				}
				MPI::COMM_WORLD.Recv(&sumClusters[0], k, pointType, i, DEFAULT_TAG);
				MPI::COMM_WORLD.Recv(&pointsPerCluster[0], k, MPI::INT, i, DEFAULT_TAG);

				for (int j = 0; j < k; j++) {
					aggrSumClusters[j] += sumClusters[j];
					aggrPointsPerCluster[j] += pointsPerCluster[j];
				}
			}
			for (int i = 0; i < k; i++) {
				if (aggrPointsPerCluster[i] != 0) {
					centers[i] = aggrSumClusters[i] / aggrPointsPerCluster[i];
				}
			}
		}
		MPI::COMM_WORLD.Barrier();
	}
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


