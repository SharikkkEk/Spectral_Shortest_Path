export module Matrix;

import std;
import Graph;

using namespace std;

double euclede_norm(const vector<double>& v) {
	double sum = 0;
	for (double i : v) sum += i * i;
	return sqrt(sum);
}

void normalize(vector<double>& v) {
	double x_norm = euclede_norm(v);
	for (double& val : v)
		val /= x_norm;
}

double dot(const vector<double>& lvector, const vector<double>& rvector) {
	double res = 0;
	for (int i = 0; i < lvector.size(); ++i)
		res += lvector[i] * rvector[i];
	return res;
}

export class Matrix {
public:
	virtual vector<double> operator*(const vector<double>&) const = 0;
	virtual size_t size() const = 0;
};

export class Cutted_laplassian : public Matrix {
private:
	size_t _cutted;
	Graph& _graph;

	double get_function_diffs(size_t vertex, const vector<double>& function) const;
public:
	Cutted_laplassian(size_t dst, Graph& graph) 
		: _cutted{ dst }, _graph{ graph } {}
	vector<double> operator*(const vector<double>&) const override;
	size_t size() const override { return _graph.vertices_count(); }
};

export using py_eigvector_history = std::vector<std::pair<std::vector<double>, double>>;

export class Eigvector_solver {
private:
	vector<double> _eigvector;
	vector<double> _deriv;
public:
	double error;
	std::vector<double> operator()(const Matrix&);
	Eigvector_solver(double Error = 0.1) : error{ Error } {}
	py_eigvector_history solve_with_history(const Matrix&);
	double armijo(const Matrix& A, const vector<double>& x);
};

export class Standard_matrix : public Matrix {
private:
	vector<vector<double>> _data;
public:
	Standard_matrix(const vector<vector<double>>& Data) : _data{ Data } {}
	vector<double> operator*(const vector<double>& vector) const override;
	size_t size() const override { return _data.size(); }
};


double operator*(const vector<double>& lvector, const vector<double>& rvector);
vector<double> operator-(const vector<double>& v1, const vector<double>& v2);
double norm(const vector<double>& v);
double euclede_norm(const vector<double>& v);
vector<double> operator*(double scale, const vector<double>& v);

// ===== IMPL ======




double Eigvector_solver::armijo(const Matrix& A, const vector<double>& x) {
    double alpha = 1.0;
    double c = 1e-4;

    vector<double> Ax = A * x;
    double lambda = dot(x, Ax);

	vector<double>& g = _deriv;
    vector<double> p = -1.0 * g;

    double f0 = dot(x, Ax);

    while (alpha > 1e-6) {

        vector<double> x_new = x;
        for (size_t i = 0; i < x_new.size(); ++i)
            x_new[i] += alpha * p[i];

        normalize(x_new);

        double f_new = dot(x_new, A * x_new);

        if (f_new <= f0 - c * alpha * dot(g, g))
            break;

        alpha *= 0.5;
    }

    return alpha;
}

vector<double> gradient(const Matrix& A, const vector<double>& x) {
	double nx = euclede_norm(x);
	vector<double> Ax = A * x;
	double Rx = dot(x, Ax) / nx;
	vector<double> grad = Ax - Rx * x;
	grad = 2 / (nx * nx) * grad;
	return grad;
}

vector<double> riemann_gradient(const Matrix& A, const vector<double>& x) {
	double nx = euclede_norm(x);
	vector<double> Ax = A * x;
	double Rx = dot(x, Ax) / nx;
	vector<double> grad = Ax - Rx * x;
	grad = 2 * grad;
	return grad;
}

vector<double> Eigvector_solver::operator()(const Matrix& matrix) {
	_eigvector = vector<double>(matrix.size(), 1);
	_deriv = vector<double>(matrix.size(), 1);
	double learning_rate = 1;
	
	while (norm(_deriv) > error) {
		_deriv = riemann_gradient(matrix, _eigvector);
		learning_rate = armijo(matrix, _eigvector);
		for (int i = 0; i < _eigvector.size(); ++i) {
			_eigvector[i] -= learning_rate * _deriv[i]; // ��� ������
		}
		normalize(_eigvector);
	}
#ifdef DEBUG
	auto check = matrix * _eigvector;
	for (size_t i = 0; i < _eigvector.size(); ++i)
		std::cout << _eigvector[i] << ' ' << check[i] << '\n';
#endif
	return _eigvector;
}

double get_cos_normalized_v(const std::vector<double>& v) {
	return std::accumulate(v.begin(), v.end(), 0.0) / sqrt(v.size());
}

export std::vector<double> convert_normalized_v_to2D(const std::vector<double>& v) {
	double cos_v = get_cos_normalized_v(v);
	double sin_v = sqrt(1 - cos_v * cos_v);
	std::vector<double> res(2);
	res[0] = cos_v;
	res[1] = sin_v;
	return res;
}

py_eigvector_history Eigvector_solver::solve_with_history(const Matrix& matrix) {
	py_eigvector_history history;
	_eigvector = vector<double>(matrix.size(), 1);
	_deriv = vector<double>(matrix.size(), 1);
	double learning_rate = 1;
	
	while (norm(_deriv) > error) {
		auto start = std::chrono::high_resolution_clock::now();

		_deriv = riemann_gradient(matrix, _eigvector);
		learning_rate = armijo(matrix, _eigvector);
		for (int i = 0; i < _eigvector.size(); ++i) {
			_eigvector[i] -= learning_rate * _deriv[i]; // ��� ������
		}
		normalize(_eigvector);

		auto end = std::chrono::high_resolution_clock::now();
		double elapsed_seconds = std::chrono::duration<double, milli>(end - start).count();
		history.push_back({ convert_normalized_v_to2D(_eigvector), elapsed_seconds });
	}
	
#ifdef DEBUG
	auto check = matrix * _eigvector;
	for (size_t i = 0; i < _eigvector.size(); ++i)
		std::cout << _eigvector[i] << ' ' << check[i] << '\n';
#endif

	return history;
}

vector<double> Cutted_laplassian::operator*(const vector<double>& mult_vector) const {
	vector<double> res(mult_vector.size(), 0);
	size_t vertex = 0;
	
	for (vertex = 0; vertex != size(); ++vertex)
		if (vertex != _cutted)
			res[vertex] = get_function_diffs(vertex, mult_vector);

	res[_cutted] = 0;

	return res;
}

double Cutted_laplassian::get_function_diffs(size_t vertex, const vector<double>& function) const {
	double res = 0;
	for (Vertex adjacent_vertex : _graph.adjacent_vertices(vertex))
		if (adjacent_vertex.number != _cutted)
			res += adjacent_vertex.price * (function[vertex] - function[adjacent_vertex.number]);
		else
			res += adjacent_vertex.price * function[vertex];
	
	return res;
}

// ===== HELP FUNCTIONS =====

double operator*(const vector<double>& lvector, const vector<double>& rvector) {
	double res = 0;
	for (int i = 0; i < lvector.size(); ++i)
		res += lvector[i] * rvector[i];
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

vector<double> Standard_matrix::operator*(const vector<double>& v) const {
	vector<double> res(v.size(), 0);

	for (int row = 0; row < _data.size(); ++row) {
		for (int i = 0; i < _data[row].size(); ++i)
			res[row] += _data[row][i] * v[i];
	}

	return res;
}
