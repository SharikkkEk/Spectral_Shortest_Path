#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
import std;
import Graph;
import Matrix;
import algorithms;

std::vector<std::vector<double>> matrix_random_sp_graph(int n) {
	return random_graph(Graph_type::Sparse, n).adjacency_matrix();
}

std::vector<std::vector<double>> matrix_random_tree(int n) {
	return random_graph(Graph_type::Tree, n).adjacency_matrix();
}

std::vector<std::vector<double>> matrix_random_d_graph(int n) {
	return random_graph(Graph_type::Dense, n).adjacency_matrix();
}

py_path_history spectral_history(std::vector<std::vector<double>> matrix) {
	Graph graph{ matrix };
	return spectral_path_history(0, matrix.size()-1, graph).to_list();
}

py_path_history dijkstra_history(std::vector<std::vector<double>> matrix){
	Graph graph{ matrix };
	return dijkstra_path_history(0, matrix.size()-1, graph).to_list();
}

std::vector<double> eigenvector(std::vector<std::vector<double>> matrix_of_graph) {
	return Eigvector_solver{}(Standard_matrix{ matrix_of_graph });
	//Graph graph{ matrix_of_graph };
	//Cutted_laplassian laplassian { matrix_of_graph.size()-1, graph};
}

py_eigvector_history eigenvector_history(std::vector<std::vector<double>> matrix_of_graph) {
	Graph graph{ matrix_of_graph };
	Cutted_laplassian laplassian { matrix_of_graph.size()-1, graph};
	Eigvector_solver solver{};
	return solver.solve_with_history(laplassian);
}

std::vector<std::vector<double>> test(std::vector<std::vector<double>> v) {
	return v;
}

PYBIND11_MODULE(SpectralPath, m) { 

	m.def("rand_sparse_graph", matrix_random_sp_graph, "");
	m.def("rand_dense_graph", matrix_random_d_graph, "");
	m.def("rand_tree", matrix_random_tree, "");
	m.def("spectral_history", spectral_history, "");
	m.def("dijkstra_history", dijkstra_history, "");
	m.def("eigenvector", eigenvector, "");
	m.def("eigenvector_history", eigenvector_history, "");
	m.def("test", test, "");
}

//int main() {
//	std::vector<std::vector<double>> v = {
//		{0, 1, 2},
//		{1, 0, 3},
//		{2, 3, 0}
//	};
//	Standard_matrix m{ v };
//	Eigvector_solver solver{};
//	std::vector<double> res = solver(m);
//	for (double i : res)
//		std::cout << i << '\n';
//}
