export module Graph;

import std;
using namespace std;

export enum class Graph_type { Tree, Sparse, Dense };

export struct Vertex {
	size_t number;
	double price;
};


export class Graph {
private:
	vector<vector<Vertex>> adjacency_list;
public:
	Graph() : adjacency_list() {}
	Graph(size_t Vertices_count) : adjacency_list(Vertices_count) {}
	Graph(vector<vector<double>> lst);

	vector<vector<double>> adjacency_matrix() const;
	size_t vertices_count() const { return adjacency_list.size(); }
	const vector<Vertex>& adjacent_vertices(size_t v) const { return adjacency_list[v]; }

	void add_vertex() { adjacency_list.push_back(vector<Vertex>()); }
	void add_edge(size_t beg, size_t dst, double price);
};

export using py_path_history = std::vector<std::pair<size_t, double>>;
export using py_adjacency_list = std::vector<std::vector<std::pair<size_t, double>>>;
export py_adjacency_list get_adjacency_list(const Graph& g);

export struct path_history {
	Graph& graph;
	py_path_history history;
	path_history(Graph& g) : graph{ g }, history() {}

	void add(size_t vertex, double time) {
		history.push_back({ vertex, time });
	}

	py_path_history to_list() const {
		return history;
	}
};


// ===== IMPL =====

void Graph::add_edge(size_t beg, size_t dst, double price) {
	if (beg == dst) return;
	auto dst_ptr = find_if(adjacency_list[beg].begin(), adjacency_list[beg].end(),
		[dst](Vertex v) {
			return v.number == dst;
		}
	);

	if (dst_ptr == adjacency_list[beg].end()) {
		adjacency_list[beg].push_back(Vertex{ dst, static_cast<double>(price) });
		adjacency_list[dst].push_back(Vertex{ beg, static_cast<double>(price) });
	}
}

double round_double(double value) {
	return round(value * 100.0) / 100.0;
}
Graph::Graph(vector<vector<double>> matrix) {
	adjacency_list.resize(matrix.size());
	for (size_t i = 0; i < matrix.size(); ++i) {
		for (size_t j = 0; j < matrix[i].size(); ++j) {
			if (matrix[i][j] != 0)
				add_edge(i, j, round_double(matrix[i][j]));
		}
	}
}

vector<vector<double>> Graph::adjacency_matrix() const {
	vector<vector<double>> res(
		vertices_count(), vector<double>(vertices_count(), 0)
	);

	for (int i = 0; i < vertices_count(); ++i) {
		for (Vertex vertex : adjacent_vertices(i))
			res[i][vertex.number] = vertex.price;
	}

	return res;
}

export std::ostream& operator<<(std::ostream& ost, Graph& g) {
	vector<vector<double>> matrix = g.adjacency_matrix();

	for (auto col : matrix) {
		for (double i : col)
			ost << i << '\t';
		ost << '\n';
	}
	
	return ost;
}

Graph create_tree(size_t n, int weight_min, int weight_max) {
	Graph result(n);

	uniform_int_distribution<int> rand_weight(weight_min, weight_max);
	default_random_engine eng;

	for (size_t i = 1; i < n; ++i) {
		result.add_edge(i, i - 1, rand_weight(eng));
	}

	return result;
}

void add_random_edges(Graph& g, double chance, int weight_min, int weight_max) {
	if (g.vertices_count() < 2) return;

	uniform_real_distribution<double> rand_weight(weight_min, weight_max);
	default_random_engine eng;
	uniform_real_distribution<double> rand_chance(0, 1);

	for (size_t i = 0; i < g.vertices_count() - 2; ++i) {
		if (rand_chance(eng) < chance) {
			uniform_int_distribution<int> rand_vertex(i+1, g.vertices_count()-1);
			g.add_edge(i, rand_vertex(eng), rand_weight(eng));
		}
	}
}

export Graph random_graph(Graph_type type, size_t n, int weight_min = 1, int weight_max = 10) {
	Graph result = create_tree(n, weight_min, weight_max);
	
	if (type == Graph_type::Sparse) {
		add_random_edges(result, 0.05, weight_min, weight_max);
	}
	else if (type == Graph_type::Dense) {
		add_random_edges(result, 0.1, weight_min, weight_max);
	}
	
	return result;
}

py_adjacency_list get_adjacency_list(const Graph& g) {
	py_adjacency_list result;
	for (size_t i = 0; i < g.vertices_count(); ++i) {
		std::vector<std::pair<size_t, double>> row;
		for (Vertex v : g.adjacent_vertices(i)) {
			row.push_back({v.number, v.price});
		}
		result.push_back(row);
	}
	return result;
}
