/*
 * kmeans2D.h
 *
 *  Created on: Dec 3, 2014
 *      Author: Siyuan Zhou
 *      		Zichang Feng
 */

#ifndef KMEANS2D_H_
#define KMEANS2D_H_

#include <vector>

struct Point {
	double x;
	double y;
	Point(double x = 0, double y = 0) {
		this->x = x;
		this->y = y;
	}

	Point &operator+=(const Point &other) {
		this->x += other.x;
		this->y += other.y;
		return *this;
	}

	Point operator / (const int & divisor) {
		return Point(this->x / divisor, this->y / divisor);
	}

	bool operator == (const Point &b) const {
		return (this->x == b.x && this->y == b.y);
	}

};

void read_points(char *input, std::vector<Point> &points);
std::vector<Point> select_centers(const std::vector<Point> &points, int k);
void cal_centers(const std::vector<Point> &points,
		const std::vector<int> &belonging, std::vector<Point> &centers);
double cal_dis(const Point &point, const Point &center);
void output_clusters(char *output, const std::vector<Point> &points,
				const std::vector<int> &belonging);


#endif /* KMEANS2D_H_ */
