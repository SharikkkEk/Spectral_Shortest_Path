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
		for (int i = 0; i < graph.adjacentVerticesCount(current.number); ++i) {
			size_t vertexNumber = graph.adjacentVertexNumber(current.number, i);
			double vertexPrice = graph.adjacentVertexPrice(current.number, i);
			if (dist.has(vertexNumber) && dist.findVertex(vertexNumber).price > vertexPrice + current.price ) 
				dist.decreaseKey(vertexNumber, vertexPrice + current.price);
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
	std::vector<double> eigen_vector = EigvectorSolver{} (laplassian);
	eigen_vector.insert(eigen_vector.begin() + dst, 0);
	std::unordered_map<size_t, bool> visited;

	for (size_t i = 0; i < graph.verticesCount(); ++i) {
		visited[i] = false;
		eigen_vector[i] = std::abs(eigen_vector[i]);
	}
	double price = 0;
	
	while (curr != dst) {
		visited[curr] = true;
		Vertex next = { 0, std::numeric_limits<double>::infinity()};
		double minEigenValue = std::numeric_limits<double>::infinity();

		for (int i = 0; i < graph.adjacentVerticesCount(curr); ++i) {
			size_t vertexNumber = graph.adjacentVertexNumber(curr, i);
			double vertexPrice = graph.adjacentVertexPrice(curr, i);

			if (!visited[vertexNumber] && eigen_vector[vertexNumber] < minEigenValue) {
				minEigenValue = eigen_vector[vertexNumber];
				next.number = vertexNumber;
				next.price = vertexPrice;
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

		for (int i = 0; i < graph.adjacentVerticesCount(current.number); ++i) {
			size_t vertexNumber = graph.adjacentVertexNumber(current.number, i);
			double vertexPrice = graph.adjacentVertexPrice(current.number, i);

			if (dist.find(vertexNumber) != dist.end()) {
				double new_dist = vertexPrice + current.price;
				if (new_dist < dist[vertexNumber]) {
					dist[vertexNumber] = new_dist;
					prev[vertexNumber] = current.number;  
				}
			}
		}

		current = { dist.begin()->first, dist.begin()->second };
	}

	pathHistory history{ graph };
	std::vector<size_t> path;

	size_t node = dst;
	while (node != beg) {
		path.push_back(node);
		node = prev[node];
	}
	path.push_back(beg);
	std::reverse(path.begin(), path.end());

	for (size_t i = 0; i < path.size(); ++i) {
		history.add(path[i], 0);
	}
	return history;
}

export pathHistory spectralPathHistory(size_t beg, size_t dst, Graph& graph) {
	CuttedLaplassian laplassian{ dst, graph };
	std::vector<double> eigen_vector = EigvectorSolver{} (laplassian);
	std::unordered_map<size_t, bool> visited;
	pathHistory history{ graph };

	eigen_vector.insert(eigen_vector.begin() + dst, 0);
	for (size_t i = 0; i < graph.verticesCount(); ++i) {
		visited[i] = false;
		eigen_vector[i] = std::abs(eigen_vector[i]);
	}

	double price = 0;
	size_t curr = beg;
	history.add(curr, 0);

	while (curr != dst) {
		auto start = std::chrono::steady_clock::now();

		visited[curr] = true;

		Vertex next = { 0, std::numeric_limits<double>::infinity() };
		double minEigenValue = std::numeric_limits<double>::infinity();

		for (int i = 0; i < graph.adjacentVerticesCount(curr); ++i) {
			size_t vertexNumber = graph.adjacentVertexNumber(curr, i);
			double vertexPrice = graph.adjacentVertexPrice(curr, i);

			if (!visited[vertexNumber] && eigen_vector[vertexNumber] < minEigenValue) {
				minEigenValue = eigen_vector[vertexNumber];
				next.number = vertexNumber;
				next.price = vertexPrice;
			}
		}

		price += next.price;
		curr = next.number;

		auto end = std::chrono::steady_clock::now();
		double elapsed = std::chrono::duration<double, milli>(end - start).count();
		history.add(curr, elapsed);
	}

	return history;
}

