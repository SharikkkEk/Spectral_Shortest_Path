import std;
import Matrix;
import Graph;
import Algorithms;
import BinaryHeap;

using namespace std;

export module Benchmark;

unordered_map<GraphType, string> nameFileType{
 {GraphType::Dense, "DenseGraphTest.csv"},
{GraphType::Sparse, "SparseGraphTest.csv"},
{GraphType::Tree, "TreeGraphTest.csv"}
};

struct dataForSave {
	double spectralTime = 0;
	double spectralSum = 0;
	double dijkstraTime = 0;
	double dijkstraSum = 0;
};

export void testAlgorithms(GraphType type, int startCount = 100, int testsCount = 10, int step = 10) {
	ofstream ost{ nameFileType[type] };
	ost << "Dijkstra sum, dijkstra time, spectral sum, spectral time\n";
	vector<dataForSave> savedData;
	savedData.reserve(testsCount);

	for (int i = 0; i < testsCount; ++i) {
		Graph randomedGraph = randomGraph(type, startCount);
		dataForSave save;
		
		auto startTime = chrono::high_resolution_clock::now();
		save.spectralSum = spectralPath(0, randomedGraph.verticesCount() - 1, randomedGraph);
		auto endTime = chrono::high_resolution_clock::now();
		save.spectralTime = chrono::duration<double>(endTime - startTime).count();
		
		startTime = chrono::high_resolution_clock::now();
		MinHeap binheap(randomedGraph.verticesCount());
		save.dijkstraSum = Dijkstra(0, randomedGraph.verticesCount() - 1, randomedGraph, binheap);
		endTime = chrono::high_resolution_clock::now();
		save.dijkstraTime = chrono::duration<double>(endTime - startTime).count();

		savedData.push_back(save);
		startCount += step;
	}

	for (dataForSave data : savedData) {
		ost << data.dijkstraSum << ',' << data.dijkstraTime << ',' << data.spectralSum << ',' << data.spectralTime << '\n';
	}
}
