export module algorithms;
import std;
import Matrix;
import Graph;

using namespace std;

export double Dijkstra(size_t beg, size_t dst, Graph graph) {
	map<size_t, double> dist;
	for (size_t i = 0; i < graph.vertices_count(); ++i)
		dist[i] = numeric_limits<double>::infinity();
	Vertex current{ beg, 0 };

    while (current.number != dst) {
		dist.erase(current.number);
		
		for (Vertex vertex : graph.adjacent_vertices(current.number)) {
			if (dist.find(vertex.number) != dist.end())
				dist[vertex.number] = min(dist[vertex.number], vertex.price + current.price);
		}
		
		current = { dist.begin()->first, dist.begin()->second };

#ifdef DEBUG
		if (current.price == numeric_limits<double>::infinity())
			return numeric_limits<double>::infinity();
#endif
    }

    return current.price;
}

export double spectral_path(size_t curr, size_t dst, Graph& graph) {
	Cutted_laplassian laplassian{ dst, graph };
	Eigvector_solver solver{};
	std::vector<double> eigen_vector = solver(laplassian);
	eigen_vector.insert(eigen_vector.begin() + dst, 0);
	std::unordered_map<size_t, bool> visited;
	for (size_t i = 0; i < graph.vertices_count(); ++i) {
		visited[i] = false;
		eigen_vector[i] = std::abs(eigen_vector[i]);
	}
	double price = 0;
	
	while (curr != dst) {
		visited[curr] = true;
		auto& adjacent_vertices = graph.adjacent_vertices(curr);
		Vertex next = adjacent_vertices[0];
		for (auto vertex : adjacent_vertices) {
			if (eigen_vector[vertex.number] <= eigen_vector[next.number] && !visited[vertex.number])
				next = vertex;
		}
		price += next.price;
		curr = next.number;
	}
	
	return price;
}

export path_history dijkstra_path_history(size_t beg, size_t dst, Graph& graph) {
	path_history history{ graph };

	map<size_t, double> dist;
	for (size_t i = 0; i < graph.vertices_count(); ++i)
		dist[i] = numeric_limits<double>::infinity();
	Vertex current{ beg, 0 };

	auto end = std::chrono::high_resolution_clock::now();
	auto start = std::chrono::high_resolution_clock::now();

    while (current.number != dst) {
		double elapsed_seconds = std::chrono::duration<double, milli>(end - start).count();
		history.add(current.number, elapsed_seconds);
		start = std::chrono::high_resolution_clock::now();

		dist.erase(current.number);
		
		for (Vertex vertex : graph.adjacent_vertices(current.number)) {
			if (dist.find(vertex.number) != dist.end())
				dist[vertex.number] = min(dist[vertex.number], vertex.price + current.price);
		}
		
		current = { dist.begin()->first, dist.begin()->second };

		end = std::chrono::high_resolution_clock::now();
    }

	return history;
}

export path_history spectral_path_history(size_t curr, size_t dst, Graph& graph) {
	Cutted_laplassian laplassian{ dst, graph };
	Eigvector_solver solver{};
	std::vector<double> eigen_vector = solver(laplassian);
	eigen_vector.insert(eigen_vector.begin() + dst, 0);
	std::unordered_map<size_t, bool> visited;
	for (size_t i = 0; i < graph.vertices_count(); ++i) {
		visited[i] = false;
		eigen_vector[i] = std::abs(eigen_vector[i]);
	}

	path_history history{ graph };
	double price = 0;
	auto start = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	
	while (curr != dst) {
		double elapsed_seconds = std::chrono::duration<double, milli>(end - start).count();
		history.add(curr, elapsed_seconds);
		start = std::chrono::high_resolution_clock::now();

		visited[curr] = true;
		auto& adjacent_vertices = graph.adjacent_vertices(curr);
		Vertex next = adjacent_vertices[0];
		for (auto vertex : adjacent_vertices) {
			if (eigen_vector[vertex.number] <= eigen_vector[next.number] && !visited[vertex.number])
				next = vertex;
		}
		price += next.price;
		curr = next.number;

		end = std::chrono::high_resolution_clock::now();
	}
	
	return history;
}

