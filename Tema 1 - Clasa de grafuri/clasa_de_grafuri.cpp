#include <fstream>
#include <algorithm>
#include <queue>
#include <iterator>
#include <vector>
#include <stack>

using namespace std;

ifstream inputFile("dfs.in");
ofstream outputFile("dfs.out");

constexpr int maxSize = 100001;

void countingSort(vector<int>&);
template <class T>
void printVector(vector<T>);

class Graph {

	vector<vector<int>> adjacencyLists;
	int nrNodes;
	bool isDirected;

public:
	Graph(int, bool);
	void initializeAdjacencyLists();

	void readEdges(int);
	void printDistancesToNode(int);
	void BFS(int, vector<int>&);
	int numberOfConnectedComponents();
	void DFS(int, vector<int>&);
	void DFS(int, vector<int>&, stack<int>&);
	bool havelHakimi(vector<int>);
	vector<int> topologicalSort();
};

Graph::Graph(int size = maxSize, bool isDirected = false) {
	nrNodes = size;
	this->isDirected = isDirected;
	initializeAdjacencyLists();
}


/// <summary>
/// initializes empty adjacency lists (vectors)
/// </summary>
void Graph::initializeAdjacencyLists() {
	for (int i = 0; i < nrNodes; i++) {
		vector<int> emptyAdjacencyVector;
		adjacencyLists.push_back(emptyAdjacencyVector);
	}
}

/// <summary>
/// reads a number of edges given as parameter
/// </summary>
void Graph::readEdges(int nrEdges) {
	int inNode, outNode;

	if (isDirected) {
		for (int i = 0; i < nrEdges; i++) {
			inputFile >> outNode >> inNode;
			outNode--;
			inNode--;
			adjacencyLists[outNode].push_back(inNode);
		}
	}
	else {
		for (int i = 0; i < nrEdges; i++) {
			inputFile >> outNode >> inNode;
			outNode--;
			inNode--;

			adjacencyLists[outNode].push_back(inNode);
			adjacencyLists[inNode].push_back(outNode);
		}
	}
}
/// <summary>
/// Prints the distance from a node given as parameter to all nodes in the graph, or -1 if they are inaccesible.
/// </summary>
/// <param name="startNode">the node for which distances are calculated</param>
void Graph::printDistancesToNode(int startNode) {

	vector<int> distances(nrNodes, -1);
	distances[startNode] = 0;          //distance from starting node to starting node is 0
	BFS(startNode, distances);

	printVector(distances);
}


/// <summary>
/// Does a breadth-first search and maps the distances to a vector of ints given as parameter.
/// </summary>
/// <param name="startNode"> starting node of search </param>
/// <param name="distances"> vector to map distances to </param>
void Graph::BFS(int startNode, vector<int>& distances) {
	int currentNode, currentDistance;
	queue<int> toVisit;
	toVisit.push(startNode);

	// while there still are accesible nodes that were not visited
	while (toVisit.empty() != true) {
		currentNode = toVisit.front();
		currentDistance = distances[currentNode];

		/// iterate through the current node's neighbors
		for (int neighboringNode : adjacencyLists[currentNode])
			if (distances[neighboringNode] == -1) {
				toVisit.push(neighboringNode);
				distances[neighboringNode] = currentDistance + 1;
			}
		toVisit.pop();
	}
}


/// <summary>
/// Computes number of connected components
/// </summary>
/// <returns></returns>
int Graph::numberOfConnectedComponents() {
	int nr = 0;
	vector<int> visited(nrNodes, 0); //all nodes ar unvisited at start

	//go through all nodes
	for (int i = 0; i < nrNodes; i++)
		//if we still have univisted nodes that means we need another DFS => another connected component
		if (visited[i] == 0) {
			nr++;
			DFS(i, visited);
		}
	return nr;
}

/// <summary>
/// Iterative depth-first traversal that maps visited nodes in a vector of ints given as parameter.
/// </summary>
void Graph::DFS(int startNode, vector<int>& visited){

	int currentNode;
	stack<int> toVisit;

	toVisit.push(startNode);

	// while there still are accesible nodes that were not visited
	while (toVisit.empty() != true) {
		currentNode = toVisit.top();

		toVisit.pop();

		// if current node was not visited before
		if (visited[currentNode] == 0) {
			// iterate through the current node's neighbors
			for (auto node : adjacencyLists[currentNode])
				if (visited[node] == 0)
					toVisit.push(node);

			visited[currentNode] = 1;
		}
	}
}

/// <summary>
/// Recursive DFS overload that pushes nodes to a stack when returning from them. Used in topological sort. 
/// </summary>
void Graph::DFS(int currentNode, vector<int>& visited, stack<int>& solution) {

	visited[currentNode] = 1;

	for (int node : adjacencyLists[currentNode]) 
		if (!visited[node]) 
			DFS(node, visited, solution);
		
	//add 1 that was subtracted on read
	solution.push(currentNode + 1);
}

/// <summary>
/// Checks if a collection of degrees can form a graph.
/// </summary>
bool Graph::havelHakimi(vector<int> degrees) {

	//check if sum of degrees is even or odd
	int sum = 0;
	for (int i = 0; i < nrNodes; i++)
		sum += degrees[i];

	//if odd, degrees can not form a graph
	if (sum % 2)
		return 0;

	//sort descending
	countingSort(degrees);

	//if biggest degree is larger than the number of nodes - 1 the degrees can't form a graph
	if (degrees[0] > nrNodes - 1) {
		return 0;
	}

	while (degrees[0] != 0) {
		outputFile << endl;
		//for the next degrees[0] nodes, subtract one because we connect this node to them
		for (int i = 1; i <= degrees[0]; i++)
			if (degrees[i] != 0)
				degrees[i]--;
			else // if degrees[i] = 0 there are no more nodes to connect to, degrees can't form a graph
				return 0;

		//the current node was connected the required number of times and is set to 0
		degrees[0] = 0;
		countingSort(degrees);
	}

	//at this point degrees vector is empty so degrees can form a graph
	return 1;
}

/// <summary>
/// Computes a topological sort of an oriented graph.
/// Note: assumes graph is acyclic.
/// </summary>
/// <returns>Vector of topological sort or empty vector if graph is not directed</returns>
vector<int> Graph::topologicalSort() {
	vector<int> sortedItems;
	stack<int> solution;

	//if graph is not directed returns empty vector
	if (!isDirected) {
		return sortedItems;
	}
	else {
		vector<int> visited(nrNodes, 0);
		for (int i = 0; i < nrNodes; i++) 
			if (!visited[i]) 
				DFS(i, visited, solution);
	}

	//move items from solution stack to sorted vector I.E. reverse order
	while (!solution.empty()) {
		sortedItems.push_back(solution.top());
		solution.pop();
	}
	return sortedItems;
}

int nrNoduri, nrMuchii;

int main()
{
	inputFile >> nrNoduri >> nrMuchii;

	Graph graph(nrNoduri, false);
	graph.readEdges(nrMuchii);
	outputFile << graph.numberOfConnectedComponents();

	/* MAIN HAVEL HAKIMI
	//files: havelhakimi.in, havelhakimi.out
	inputFile >> nrNoduri;

	Graph graph(nrNoduri, false);

	vector<int> degrees(nrNoduri);

	for (int i = 0; i < nrNoduri; i++) {
		inputFile >> degrees[i];
	}

	if (graph.havelHakimi(degrees)) {
		outputFile << "Degrees can form a graph.";
	}
	else {
		outputFile << "Degrees can't form a graph.";
	}
	*/
}


/// <summary>
/// Basic descending counting sort algorithm, only works for numbers up to 1000;
/// </summary>
void countingSort(vector<int>& toSort) {
	vector<int> count(1000, 0);
	int maxElem = -1;

	for (int i = 0; i < toSort.size(); i++) {
		count[toSort[i]]++;

		//check for maximum number to cut down execution time since we expect small numbers
		if (toSort[i] > maxElem) {
			maxElem = toSort[i];
		}
	}

	int currentIndex = 0;

	for (int i = maxElem; i >= 0; i--) {
		while (count[i] != 0) {
			toSort[currentIndex++] = i;
			count[i]--;
		}
	}
}

template <class T>
void printVector(vector<T> toPrint) {
	for (int i = 0; i < toPrint.size(); i++)
		outputFile << toPrint[i] << " ";
}