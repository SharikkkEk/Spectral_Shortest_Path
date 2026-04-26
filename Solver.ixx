export module Solver;
import ConstantsAndUtils;
import Graph;
import Matrix;
import std;
using namespace std;

ostream& operator<<(ostream& ost, vector<double>& v) {
	for (double x : v)
		ost << x << '\n';
	return ost;
}

export class EigvectorSolver {
private:
	vector<double> _eigvector;
	vector<double> _gradient;
	vector<double> _auxiliaryVector;
	double _currentFunctionValue;
	double armijo(const SquareMatrix& A, const vector<double>& x);
	void riemannienGradient(const SquareMatrix& A, const vector<double>& x);
public:
	double error;
	std::vector<double> operator()(const SquareMatrix&);
	EigvectorSolver(double Error = constError) : error{ Error } {}
	PyEigvectorHistory solveWithHistory(const SquareMatrix&);
};

vector<double> EigvectorSolver::operator()(const SquareMatrix& matrix) {
	_eigvector = vector<double>(matrix.size(), 1);
	normalize(_eigvector);
	_gradient = vector<double>(matrix.size(), 1);
	_auxiliaryVector = vector<double>(matrix.size(), 1);
	double learning_rate = 1;
	_currentFunctionValue = matrix.calcRaleigh(_eigvector);
	double diff = 0;

	do {
		double oldFunctionValue = _currentFunctionValue;
		riemannienGradient(matrix, _eigvector);
		learning_rate = armijo(matrix, _eigvector);
		for (int i = 0; i < _eigvector.size(); ++i) {
			_eigvector[i] -= learning_rate * _gradient[i]; // ��� ������
		}
		normalize(_eigvector);
		_currentFunctionValue = matrix.calcRaleigh(_eigvector);
		diff = abs(_currentFunctionValue - oldFunctionValue);
	} while (diff > error);

	return _eigvector;
}

double EigvectorSolver::armijo(const SquareMatrix& A, const vector<double>& x) {
    double alpha = 1.0;
    double smallConstant = 1e-4;
	double p = -dot(_gradient, _gradient);
	
    while (alpha > armijoRatio) {
        for (size_t i = 0; i < _auxiliaryVector.size(); ++i)
            _auxiliaryVector[i] = x[i] - alpha * _gradient[i];

        normalize(_auxiliaryVector);

		double newFunctionValue = A.calcRaleigh(_auxiliaryVector);
        if (newFunctionValue <= _currentFunctionValue + smallConstant * alpha * p)
            break;

        alpha *= 0.5;
    }

    return alpha;
}

// Substract from gradent component that takes it off the sphere (dot product is projection of one vector on another)
void EigvectorSolver::riemannienGradient(const SquareMatrix& A, const vector<double>& x) {
	for (int i = 0; i < x.size(); ++i) {
		_gradient[i] = 2 * (A.productRow(x, i) - _currentFunctionValue * x[i]);
	}
}

/*
*	Косинус получаем по отношению к единичному вектору
*/
double getNormalizedVectorCos(const std::vector<double>& v) {
	double len = sqrt(v.size()); // У единичного вектора в n-ной размерности длина равна корню из количества его элементов
	return std::accumulate(v.begin(), v.end(), 0.0) / len;
}

export std::vector<double> mapVectorTo2D(const std::vector<double>& v) {
	double cos_v = getNormalizedVectorCos(v);
	double sin_v = sqrt(1 - cos_v * cos_v); // Основное тригонометрическое тождество
	std::vector<double> res(2);
	res[0] = cos_v;
	res[1] = sin_v;
	return res;
}

// Мне кажется это сейчас не работает
PyEigvectorHistory EigvectorSolver::solveWithHistory(const SquareMatrix& matrix) {
	PyEigvectorHistory history;
	_eigvector = vector<double>(matrix.size(), 1);
	_gradient = vector<double>(matrix.size(), 1);
	double learning_rate = 1;
	
	while (norm(_gradient) > error) {
		auto start = std::chrono::high_resolution_clock::now();

		riemannienGradient(matrix, _eigvector);
		learning_rate = armijo(matrix, _eigvector);
		for (int i = 0; i < _eigvector.size(); ++i) {
			_eigvector[i] -= learning_rate * _gradient[i]; // ��� ������
		}
		normalize(_eigvector);

		auto end = std::chrono::high_resolution_clock::now();
		double elapsed_seconds = std::chrono::duration<double, milli>(end - start).count();
		history.push_back({ mapVectorTo2D(_eigvector), elapsed_seconds });
	}
	
#ifdef DEBUG
	auto check = matrix * _eigvector;
	for (size_t i = 0; i < _eigvector.size(); ++i)
		std::cout << _eigvector[i] << ' ' << check[i] << '\n';
#endif

	return history;
}
