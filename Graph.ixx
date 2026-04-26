export module Graph;
import ConstantsAndUtils;
import std;
using namespace std;

// ====== Declarations =====

export enum class GraphType { Tree, Sparse, Dense };

export struct Vertex {
	size_t number;
	double price;
	bool operator<(const Vertex& another) const { return this->price < another.price; }
	bool operator>(const Vertex& another) const { return this->price > another.price; }
	bool operator==(const Vertex& another) const { return this->price == another.price && this->number == another.number; }
};

export class Graph {
private:
	vector<vector<Vertex>> adjacencyList;
public:
	Graph() : adjacencyList() {}
	Graph(size_t VerticesCount) : adjacencyList(VerticesCount) {}
	Graph(vector<vector<double>> lst);

	vector<vector<double>> adjacencyMatrix() const;
	size_t verticesCount() const { return adjacencyList.size(); }
	const vector<Vertex>& adjacentVertices(size_t v) const { return adjacencyList[v]; }

	void add_vertex() { adjacencyList.push_back(vector<Vertex>()); }
	void addEdge(size_t beg, size_t dst, double price);
};

export using pyPathHistory = std::vector<std::pair<size_t, double>>;
export using pyAdjacencyList = std::vector<std::vector<std::pair<size_t, double>>>;
export pyAdjacencyList getAdjacencyList(const Graph& g);
export Graph randomGraph(GraphType type, size_t n, int weight_min = 1, int weight_max = 10);
export std::ostream& operator<<(std::ostream& ost, Graph& g);

export struct pathHistory {
	Graph& graph;
	pyPathHistory history;
	pathHistory(Graph& g) : graph{ g }, history() {}

	void add(size_t vertex, double time) {
		history.push_back({ vertex, time });
	}

	pyPathHistory to_list() const {
		return history;
	}
};

Graph createTree(size_t n, int weight_min, int weight_max);
double normalizeChance(double chance, int verticesCount);
void addRandomEdges(Graph& g, double chance, int weight_min, int weight_max);

// ===== IMPL =====

void Graph::addEdge(size_t beg, size_t dst, double price) {
	if (beg == dst) return;
	auto dst_ptr = find_if(adjacencyList[beg].begin(), adjacencyList[beg].end(),
		[dst](Vertex v) {
			return v.number == dst;
		}
	);

	if (dst_ptr == adjacencyList[beg].end()) {
		adjacencyList[beg].push_back(Vertex{ dst, static_cast<double>(price) });
		adjacencyList[dst].push_back(Vertex{ beg, static_cast<double>(price) });
	}
}

Graph::Graph(vector<vector<double>> matrix) {
	adjacencyList.resize(matrix.size());
	for (size_t i = 0; i < matrix.size(); ++i) {
		for (size_t j = 0; j < matrix[i].size(); ++j) {
			if (matrix[i][j] != 0)
				addEdge(i, j, roundDouble(matrix[i][j]));
		}
	}
}

vector<vector<double>> Graph::adjacencyMatrix() const {
	vector<vector<double>> res(
		verticesCount(), vector<double>(verticesCount(), 0)
	);

	for (int i = 0; i < verticesCount(); ++i) {
		for (Vertex vertex : adjacentVertices(i))
			res[i][vertex.number] = vertex.price;
	}

	return res;
}

Graph randomGraph(GraphType type, size_t n, int weight_min, int weight_max) {
	Graph result = createTree(n, weight_min, weight_max);

	// Из-за алгоритма генерации рандомных графов, даже разреженные графы при n > 300 становятся слишком плотными, 
	// т.к. мы каждую вершину пытаемся проверить на добавление дополнительного ребра
	// Поэтому нам нужно нормализовывать шанс, cлегка уменьшая его с ростом количества вершин
	// Степенная функция a^n при a чуть меньшим 1 идеально подходит для этого

	if (type == GraphType::Sparse) {
		addRandomEdges(result, normalizeChance(sparseChance, n), weight_min, weight_max);
	}
	else if (type == GraphType::Dense) {
		addRandomEdges(result, normalizeChance(denseChance, n), weight_min, weight_max);
	}

	return result;
}

// ===== HELP =====

pyAdjacencyList getAdjacencyList(const Graph& g) {
	pyAdjacencyList result;
	for (size_t i = 0; i < g.verticesCount(); ++i) {
		std::vector<std::pair<size_t, double>> row;
		for (Vertex v : g.adjacentVertices(i)) {
			row.push_back({v.number, v.price});
		}
		result.push_back(row);
	}
	return result;
}

Graph createTree(size_t n, int weight_min, int weight_max) {
	Graph result(n);

	uniform_int_distribution<int> rand_weight(weight_min, weight_max);
	std::time_t seed = std::time(nullptr);
	default_random_engine eng(seed);

	for (size_t i = 1; i < n; ++i) {
		result.addEdge(i, i - 1, rand_weight(eng));
	}

	return result;
}

double normalizeChance(double chance, int verticesCount) {
	return std::pow(normalizationRatio, verticesCount) * chance;
}

void addRandomEdges(Graph& g, double chance, int weight_min, int weight_max) {
	if (g.verticesCount() < 2) return;

	uniform_real_distribution<double> rand_weight(weight_min, weight_max);
	std::time_t seed = std::time(nullptr);
	default_random_engine eng(seed);
	uniform_real_distribution<double> rand_chance(0, 1);

	for (size_t i = 0; i < g.verticesCount() - 2; ++i) {
		if (rand_chance(eng) < chance) {
			uniform_int_distribution<int> rand_vertex(i + 1, g.verticesCount() - 1);
			g.addEdge(i, rand_vertex(eng), rand_weight(eng));
		}
	}
}

std::ostream& operator<<(std::ostream& ost, Graph& g) {
	vector<vector<double>> matrix = g.adjacencyMatrix();

	for (auto col : matrix) {
		for (double i : col)
			ost << i << '\t';
		ost << '\n';
	}

	return ost;
}

