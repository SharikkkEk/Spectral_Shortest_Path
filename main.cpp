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

py_adjacency_list get_py_adjacency_list(std::vector<std::vector<double>> matrix) {
	Graph graph{ matrix };
	return get_adjacency_list(graph);
}

std::vector<double> eigenvector_2D(std::vector<std::vector<double>> matrix_of_graph) {
	Graph graph{ matrix_of_graph };
	Cutted_laplassian laplassian { matrix_of_graph.size()-1, graph};
	std::vector<double> eigvector = Eigvector_solver{}(laplassian);
	return convert_normalized_v_to2D(eigvector);
}

py_eigvector_history eigenvector_2D_history(std::vector<std::vector<double>> matrix_of_graph) {
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
	m.def("get_adjacency_list", get_py_adjacency_list, "");
	m.def("eigenvector_2D", eigenvector_2D, "");
	m.def("eigenvector_2D_history", eigenvector_2D_history, "");
	m.def("test", test, "");
}

//int main() {
//	std::vector<std::vector<double>> v = {
//		{0, 1, 2},
//		{1, 0, 3},
//		{2, 3, 0}
//	};
//	Graph g{ v };
//	eigenvector_2D(v);
//	eigenvector_2D_history(v);
//}
