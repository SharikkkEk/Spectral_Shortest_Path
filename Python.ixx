module;

#ifdef PYTHON

#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 

#endif

export module PythonFunctions;

#ifdef PYTHON

import std;
import Graph;
import Matrix;
import Solver;
import Algorithms;


std::vector<std::vector<double>> matrix_random_sp_graph(int n) {
	return randomGraph(GraphType::Sparse, n).adjacencyMatrix();
}

std::vector<std::vector<double>> matrix_random_tree(int n) {
	return randomGraph(GraphType::Tree, n).adjacencyMatrix();
}

std::vector<std::vector<double>> matrix_random_d_graph(int n) {
	return randomGraph(GraphType::Dense, n).adjacencyMatrix();
}

double spectral(std::vector<std::vector<double>> matrix) {
	Graph graph{ matrix };
	return spectralPath(0, matrix.size() - 1, graph);
}

double dijkstra(std::vector<std::vector<double>> matrix) {
	Graph graph{ matrix };
	return Dijkstra(0, matrix.size() - 1, graph);
}

pyPathHistory spectral_history(std::vector<std::vector<double>> matrix) {
	Graph graph{ matrix };
	return spectralPathHistory(0, matrix.size() - 1, graph).to_list();
}

pyPathHistory dijkstra_history(std::vector<std::vector<double>> matrix) {
	Graph graph{ matrix };
	return dijkstraPathHistory(0, matrix.size() - 1, graph).to_list();
}

pyAdjacencyList get_py_adjacency_list(std::vector<std::vector<double>> matrix) {
	Graph graph{ matrix };
	return getAdjacencyList(graph);
}

std::vector<double> eigenvector_2D(std::vector<std::vector<double>> matrix_of_graph) {
	Graph graph{ matrix_of_graph };
	CuttedLaplassian laplassian{ matrix_of_graph.size() - 1, graph };
	std::vector<double> eigvector = EigvectorSolver{}(laplassian);
	return mapVectorTo2D(eigvector);
}

PyEigvectorHistory eigenvector_2D_history(std::vector<std::vector<double>> matrix_of_graph) {
	Graph graph{ matrix_of_graph };
	CuttedLaplassian laplassian { matrix_of_graph.size()-1, graph};
	EigvectorSolver solver{};
	return solver.solveWithHistory(laplassian);
}

PYBIND11_MODULE(SpectralPath, m) {
	m.def("rand_sparse_graph", matrix_random_sp_graph, "");
	m.def("rand_dense_graph", matrix_random_d_graph, "");
	m.def("rand_tree", matrix_random_tree, "");
	m.def("spectral_history", spectral_history, "");
	m.def("dijkstra_history", dijkstra_history, "");
	m.def("spectral", spectral, "");
	m.def("dijkstra", dijkstra, "");
	m.def("get_adjacency_list", get_py_adjacency_list, "");
	m.def("eigenvector_2D", eigenvector_2D, "");
	m.def("eigenvector_2D_history", eigenvector_2D_history, "");
}

#endif
