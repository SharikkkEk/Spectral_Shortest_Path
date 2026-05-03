export module Graph;
import Matrix;
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
	vector<int> adjacencyList;
	vector<double> priceList;
	vector<int> verticesPos; // Последний элемент = длина adjacencyList, для удобного вычисления количества смежных вершин
public:
	Graph() : adjacencyList{}, priceList{}, verticesPos{} {}
	Graph(size_t VerticesCount) : verticesPos(VerticesCount+1, 0) { }
	//Graph(vector<vector<double>> lst);

	vector<vector<double>> adjacencyMatrix() const;
	size_t verticesCount() const { return verticesPos.size() - 1; }
	
	size_t adjacentVerticesCount(size_t v) const { return verticesPos[v + 1] - verticesPos[v]; }
	size_t adjacentVertexNumber(size_t vertex, size_t i) const {
		return adjacencyList[verticesPos[vertex] + i]; };
	double adjacentVertexPrice(size_t vertex, size_t i) const { return priceList[verticesPos[vertex] + i]; }

	void add_vertex();
	void addEdge(size_t beg, size_t dst, double price);
};

export class CuttedLaplassian : public SquareMatrix {
private:
	size_t _cutted;
	Graph& _graph;
public:
	CuttedLaplassian(size_t dst, Graph& graph)
		: _cutted{ dst }, _graph{ graph } {
	}
	vector<double> operator*(const vector<double>&) const override;
	double calcRaleigh(const vector<double>&) const override;
	double productRow(const vector<double>& v, size_t row) const override;
	size_t size() const override { return _graph.verticesCount(); }
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

void Graph::add_vertex() {
	verticesPos.insert(verticesPos.begin() + adjacencyList.size(), adjacencyList.size());
}

void Graph::addEdge(size_t beg, size_t dst, double price) {
	adjacencyList.insert(adjacencyList.begin() + verticesPos[beg], dst);
	priceList.insert(priceList.begin() + verticesPos[beg], price);
	for (int shiftedPositions = beg + 1; shiftedPositions < verticesPos.size(); ++shiftedPositions)
		verticesPos[shiftedPositions]++;

	adjacencyList.insert(adjacencyList.begin() + verticesPos[dst], beg);
	priceList.insert(priceList.begin() + verticesPos[dst], price);
	for (int shiftedPositions = dst + 1; shiftedPositions < verticesPos.size(); ++shiftedPositions)
		verticesPos[shiftedPositions]++;

}

/*Graph::Graph(vector<vector<double>> matrix) {
	adjacencyList.resize(matrix.size());
	for (size_t i = 0; i < matrix.size(); ++i) {
		for (size_t j = 0; j < matrix[i].size(); ++j) {
			if (matrix[i][j] != 0)
				addEdge(i, j, roundDouble(matrix[i][j]));
		}
	}
}*/

vector<vector<double>> Graph::adjacencyMatrix() const {
	vector<vector<double>> res(
		verticesCount(), vector<double>(verticesCount(), 0)
	);

	for (int i = 0; i < verticesCount(); ++i) {
		for (int j = 0; j < adjacentVerticesCount(i); ++j)
			res[i][adjacentVertexNumber(i, j)] = adjacentVertexPrice(i, j);
	}

	return res;
}

// ===== Laplassian IMPL =====

double CuttedLaplassian::productRow(const vector<double>& v, size_t row) const {
	double sum = 0;
	for (int j = 0; j < _graph.adjacentVerticesCount(row); ++j) {
		int vertexNumber = _graph.adjacentVertexNumber(row, j);
		double vertexPrice = _graph.adjacentVertexPrice(row, j);
		if (vertexNumber != _cutted)
			sum += vertexPrice * (v[row] - v[vertexNumber]);
		else
			sum += vertexPrice * v[row];
	}
	return sum;
}

double CuttedLaplassian::calcRaleigh(const vector<double>& v) const {
	double sum = 0;
	for (size_t i = 0; i < _graph.verticesCount(); ++i) {
		for (size_t j = 0; j < _graph.adjacentVerticesCount(i); ++j) {
			size_t vertexNumber = _graph.adjacentVertexNumber(i, j);
			size_t vertexPrice = _graph.adjacentVertexPrice(i, j);
			if (vertexNumber >= i) { // Матрица симметрична, элементы ниже главной диагонали являются повторами соответствующших элементов выше главной диагонали
				if (vertexNumber != _cutted) {
					double diff = v[i] - v[vertexNumber];
					sum += vertexPrice * diff * diff;
				}
				else {
					sum += vertexPrice * v[i] * v[i];
				}
			}
		}
	}
	return sum;
}

vector<double> CuttedLaplassian::operator*(const vector<double>& mult_vector) const {
	vector<double> res(mult_vector.size(), 0);
	for (size_t vertex = 0; vertex != size(); ++vertex)
		res[vertex] = productRow(mult_vector, vertex);

	res[_cutted] = 0;

	return res;
}

// ===== Generation IMPL =====

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
		for (int j = 0; j < g.adjacentVerticesCount(i); ++j) {
			row.push_back({g.adjacentVertexNumber(i, j), g.adjacentVertexPrice(i, j)});
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

