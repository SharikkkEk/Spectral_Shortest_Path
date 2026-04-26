export module Algorithms;

import std;
import Matrix;
import Graph;
import Solver;
import BinaryHeap;
using namespace std;

export double Dijkstra(size_t beg, size_t dst, Graph graph, MinHeap& dist) {
	//MinHeap dist(graph.verticesCount());
	Vertex current{ beg, 0 };

    while (!dist.empty() && current.number != dst) {
		for (Vertex vertex : graph.adjacentVertices(current.number)) {
			if (dist.has(vertex.number) && dist.findVertex(vertex.number).price > vertex.price + current.price)
				dist.decreaseKey(vertex.number, vertex.price + current.price);
		}
		current = dist.pop();
    }

    return current.price;
}


/*
Approximates shortest path on not-oriented graph
 - Works, because component of vertex in eigenVector becomes bigger, if adjacent vertices have bigger price and/or bigger value in eigenVector

Invariant: eigenVector[dst] == 0, eigenVector[i] >= 0 for all i

Why removing row and column correspondending to destination vertex:
1. Removes trivial eigenVector of laplassian (1, 1, 1, 1...)
2. Vertices adjacent to dst loose one neighbor, so their value in eigenvector becomes smaller, thanks to property of laplassian eigenvector

Complexity: O(E*k), where E - number of edges, k - number of iterations in gradient descent
*/
export double spectralPath(size_t curr, size_t dst, Graph& graph) {
	CuttedLaplassian laplassian{ dst, graph };
	EigvectorSolver solver{};
	std::vector<double> eigen_vector = solver(laplassian);
	eigen_vector.insert(eigen_vector.begin() + dst, 0);
	std::unordered_map<size_t, bool> visited;
	for (size_t i = 0; i < graph.verticesCount(); ++i) {
		visited[i] = false;
		eigen_vector[i] = std::abs(eigen_vector[i]);
	}
	double price = 0;
	
	while (curr != dst) {
		visited[curr] = true;
		auto& adjacentVertices = graph.adjacentVertices(curr);
		Vertex next = { 0, std::numeric_limits<double>::infinity()};
		double minEigenValue = std::numeric_limits<double>::infinity();
		for (auto vertex : adjacentVertices) {
			if (!visited[vertex.number] && eigen_vector[vertex.number] <= minEigenValue) {
				minEigenValue = eigen_vector[vertex.number];
				next = vertex;
			}
		}
		price += next.price;
		curr = next.number;
	}
	
	return price;
}

export pathHistory dijkstraPathHistory(size_t beg, size_t dst, Graph& graph) {
	std::map<size_t, double> dist;
	std::map<size_t, size_t> prev;

	for (size_t i = 0; i < graph.verticesCount(); ++i) {
		dist[i] = std::numeric_limits<double>::infinity();
	}

	dist[beg] = 0;
	Vertex current{ beg, 0 };

	while (current.number != dst) {
		dist.erase(current.number);

		for (Vertex vertex : graph.adjacentVertices(current.number)) {
			if (dist.find(vertex.number) != dist.end()) {
				double new_dist = vertex.price + current.price;
				if (new_dist < dist[vertex.number]) {
					dist[vertex.number] = new_dist;
					prev[vertex.number] = current.number;  
				}
			}
		}

		current = { dist.begin()->first, dist.begin()->second };

	}

	pathHistory history{ graph };
	std::vector<size_t> path;
	std::vector<double> segment_times;

	size_t node = dst;
	while (node != beg) {
		path.push_back(node);
		node = prev[node];
	}
	path.push_back(beg);
	std::reverse(path.begin(), path.end());

	for (size_t i = 0; i < path.size() - 1; ++i) {
		size_t from = path[i];
		size_t to = path[i + 1];

		for (Vertex v : graph.adjacentVertices(from)) {
			if (v.number == to) {
				segment_times.push_back(v.price);
				break;
			}
		}
	}

	for (size_t i = 0; i < path.size(); ++i) {
		double time = (i == 0) ? 0 : segment_times[i - 1];
		history.add(path[i], time);
	}
	return history;
}

export pathHistory spectralPathHistory(size_t curr, size_t dst, Graph& graph) {
	CuttedLaplassian laplassian{ dst, graph };
	EigvectorSolver solver{};
	std::vector<double> eigen_vector = solver(laplassian);
	eigen_vector.insert(eigen_vector.begin() + dst, 0);
	std::unordered_map<size_t, bool> visited;
	for (size_t i = 0; i < graph.verticesCount(); ++i) {
		visited[i] = false;
		eigen_vector[i] = std::abs(eigen_vector[i]);
	}

	pathHistory history{ graph };
	double price = 0;
	auto start = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	double elapsed_seconds = std::chrono::duration<double, milli>(end - start).count();
	history.add(curr, elapsed_seconds);
	
	while (curr != dst) {
		start = std::chrono::high_resolution_clock::now();

		visited[curr] = true;
		auto& adjacentVertices = graph.adjacentVertices(curr);
		Vertex next = adjacentVertices[0];
		for (auto vertex : adjacentVertices) {
			if (eigen_vector[vertex.number] <= eigen_vector[next.number] && !visited[vertex.number])
				next = vertex;
		}
		price += next.price;
		curr = next.number;

		end = std::chrono::high_resolution_clock::now();
		double elapsed_seconds = std::chrono::duration<double, milli>(end - start).count();
		history.add(curr, elapsed_seconds);
	}
	
	return history;
}

