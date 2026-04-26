export module BinaryHeap;

import std;
import Graph;

constexpr size_t sizeMax = std::numeric_limits<size_t>::max();
const double infinity = std::numeric_limits<double>::infinity();

export struct MinHeap {
	bool empty() { return data.empty(); }
	void decreaseKey(size_t vertex, double price);
	Vertex pop();
	bool has(size_t vertex);
	size_t size() { return data.size(); }
	Vertex& findVertex(size_t vertex) { return data[pos[vertex]]; }
	MinHeap(size_t N);
private: 
	std::vector<Vertex> data;
	std::vector<size_t> pos;
	void heapify(size_t index);
	void swapVertices(size_t, size_t);
	Vertex& parent(size_t vertex) { return data[(pos[vertex] - 1) / 2]; }
};

void MinHeap::swapVertices(size_t v1, size_t v2) {
	std::swap(findVertex(v1), findVertex(v2));
	std::swap(pos[findVertex(v1).number], pos[findVertex(v2).number]);
}

MinHeap::MinHeap(size_t N)
	: data(N), pos(N)
{
	for (size_t i = 0; i < N; ++i) {
		data[i] = Vertex{ i, infinity };
		pos[i] = i;
	}
}

void MinHeap::decreaseKey(size_t vertex, double price) {
	findVertex(vertex).price = price;
	size_t index = pos[vertex];
	if (index == 0) return;
	size_t parentVertex = parent(vertex).number;
	while (findVertex(vertex).price < findVertex(parentVertex).price) {
		swapVertices(vertex, parentVertex);
		index = pos[vertex];
		if (index == 0) break;
		parentVertex = parent(vertex).number;
	}
}

bool MinHeap::has(size_t vertex) {
	if (pos[vertex] == sizeMax)
		return false;
	return true;
}

Vertex MinHeap::pop() {
	Vertex first = data[0];
	if (data.size() > 1) {
		data[0] = data.back();
		pos[data[0].number] = 0;
		data.pop_back();
		pos[first.number] = sizeMax;
		if (!data.empty())
			heapify(0);
	}
	else {
		pos[first.number] = sizeMax;
		data.pop_back();
	}
	return first;
}

void MinHeap::heapify(size_t index) {
	size_t smallestChild = 0;
	size_t lchild = 0;
	size_t rchild = 0;
	
	for (;;) {
		lchild = 2 * index + 1;
		rchild = 2 * index + 2;
		smallestChild = index;

		if (lchild < size() && data[lchild] < data[smallestChild])
			smallestChild = lchild;
		if (rchild < size() && data[rchild] < data[smallestChild])
			smallestChild = rchild;

		if (smallestChild == index)
			return;
		
		swapVertices(data[index].number, data[smallestChild].number);

		index = smallestChild;
	}
}
