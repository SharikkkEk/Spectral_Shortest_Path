export module Matrix;

import std;
import Graph;

using namespace std;

// ===== DECLARATIONS =====

export struct SquareMatrix {
	// В принципе удалить оператор произведения сейчас можно
	virtual vector<double> operator*(const vector<double>&) const = 0;
	virtual size_t size() const = 0;
	virtual double calcRaleigh(const vector<double>&) const = 0;
	virtual double productRow(const vector<double>& v, size_t row) const = 0;
};

export class CuttedLaplassian : public SquareMatrix {
private:
	size_t _cutted;
	Graph& _graph;
public:
	CuttedLaplassian(size_t dst, Graph& graph) 
		: _cutted{ dst }, _graph{ graph } {}
	vector<double> operator*(const vector<double>&) const override;
	double calcRaleigh(const vector<double>&) const override;
	double productRow(const vector<double>& v, size_t row) const override;
	size_t size() const override { return _graph.verticesCount(); }
};

export class StandardMatrix : public SquareMatrix {
private:
	vector<vector<double>> _data;
public:
	StandardMatrix(const vector<vector<double>>& Data) : _data{ Data } {}
	vector<double> operator*(const vector<double>& vector) const override;
	size_t size() const override { return _data.size(); }
};

export using PyEigvectorHistory = std::vector<std::pair<std::vector<double>, double>>;

export double eucledeNorm(const vector<double>& v);
export double norm(const vector<double>& v);
export void normalize(vector<double>& v);

// В данный момент это всё дерьмо лишнее
export double dot(const vector<double>& lvector, const vector<double>& rvector);
export vector<double> operator-(const vector<double>& v1, const vector<double>& v2);
export vector<double> operator*(double scale, const vector<double>& v);
export vector<double>& operator*=(vector<double>& v, double scale);

// ===== IMPL ======

double CuttedLaplassian::productRow(const vector<double>& v, size_t row) const {
	double sum = 0;
	const vector<Vertex>& adjacentVertices = _graph.adjacentVertices(row);
	for (Vertex vertex : adjacentVertices) {
		if (vertex.number != _cutted)
			sum += vertex.price * (v[row] - v[vertex.number]);
		else
			sum += vertex.price * v[row];
	}
	return sum;
}

double CuttedLaplassian::calcRaleigh(const vector<double>& v) const {
	double sum = 0;
	for (int i = 0; i < _graph.verticesCount(); ++i) {
		const vector<Vertex>& adjacentVertices = _graph.adjacentVertices(i);
		for (Vertex vertex : adjacentVertices) {
			if (vertex.number >= i) {
				if (vertex.number != _cutted) {
					double diff = v[i] - v[vertex.number];
					sum += vertex.price * diff * diff;
				}
				else {
					sum += vertex.price * v[i] * v[i];
				}
			}
		}
	}
	return sum;
}

vector<double> CuttedLaplassian::operator*(const vector<double>& mult_vector) const {
	vector<double> res(mult_vector.size(), 0);
	size_t vertex = 0;
	
	for (vertex = 0; vertex != size(); ++vertex)
		if (vertex != _cutted) {
			res[vertex] = 0;

			for (Vertex adjacentVertex : _graph.adjacentVertices(vertex))
				if (adjacentVertex.number != _cutted)
					res[vertex] += adjacentVertex.price * (mult_vector[vertex] - mult_vector[adjacentVertex.number]);
				else
					res[vertex] += adjacentVertex.price * mult_vector[vertex];
		}

	res[_cutted] = 0;

	return res;
}

vector<double> operator-(const vector<double>& v1, const vector<double>& v2){
	vector<double> res(v1.size());
	for (int i = 0; i < v1.size(); ++i)
		res[i] = v1[i] - v2[i];
	return res;
}

vector<double> operator*(double scale, const vector<double>& v) {
	vector<double> res(v.size());
	for (int i = 0; i < v.size(); ++i)
		res[i] = scale * v[i];
	return res;
}

double norm(const vector<double>& v) {
	double abs_maxi = abs(v[0]);
	for (double val : v)
		abs_maxi = abs(val) > abs_maxi ? abs(val) : abs_maxi;
	return abs_maxi;
}

/*double StandardMatrix::calcRaleigh(const vector<double>& v) const {
	double sum = 0;
	for (size_t i = 0; i < _data.size(); ++i) {
		for (size_t j = 0; j < _data.size(); ++j) {
		}
	}
}

double StandardMatrix::productRow(const vector<double>& v, size_t row) const {

}*/

vector<double> StandardMatrix::operator*(const vector<double>& v) const {
	vector<double> res(v.size(), 0);

	for (int row = 0; row < _data.size(); ++row) {
		for (int i = 0; i < _data[row].size(); ++i)
			res[row] += _data[row][i] * v[i];
	}

	return res;
}

double eucledeNorm(const vector<double>& v) {
	double sum = 0;
	for (double i : v) sum += i * i;
	return sqrt(sum);
}

void normalize(vector<double>& v) {
	double x_norm = eucledeNorm(v);
	for (double& val : v)
		val /= x_norm;
}

double dot(const vector<double>& lvector, const vector<double>& rvector) {
	double res = 0;
	for (int i = 0; i < lvector.size(); ++i)
		res += lvector[i] * rvector[i];
	return res;
}

vector<double>& operator*=(vector<double>& v, double scale) {
	for (double& x : v)
		x *= scale;
	return v;
}
