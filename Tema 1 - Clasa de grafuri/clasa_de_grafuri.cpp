#include <fstream>
#include <algorithm>
#include <queue>
#include <iterator>
#include <vector>
#include <stack>
#include <list>
#include <unordered_set>

using namespace std;

ifstream inputFile("ctc.in");
ofstream outputFile("ctc.out");

constexpr int maxSize = 100001;

//Auxiliary functions not related to graphs, defined after main.
/// <summary>
/// Basic descending counting sort algorithm.
/// </summary>
/// <param name="size">maximum size of sorted numbers; maximum value = 100.000.000</param>
void countingSort(vector<int>&, int = 1000);

template <class T>
void printVector(vector<T>);

//--------------------------------------------------------------


class Graph {

	vector<vector<int>> adjacencyLists;
	int nrNodes;
	bool isDirected;

public:
	Graph(int, bool);
	void readEdges(int);
	vector<int> getDistancesToNode(int);
	vector<int> getTopologicalSort();
	list<vector<int>> getStronglyConnectedComponents();
	int getNumberOfConnectedComponents();
	friend bool havelHakimi(vector<int>);

private:
	void initializeAdjacencyLists();
	void tarjanDFS(int, int&, vector<int>&, vector<int>&, list<vector<int>>&, stack<int>&, unordered_set<int>&);
	void BFS(int, vector<int>&);
	void DFS(int, vector<int>&);
	void DFS(int, vector<int>&, stack<int>&);
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
/// Returns a vector of distances from a node given as parameter to all nodes in the graph, or -1 if they are inaccesible from that node.
/// </summary>
/// <param name="startNode">the node for which distances are calculated</param>
vector<int> Graph::getDistancesToNode(int startNode) {

	vector<int> distances(nrNodes, -1);
	distances[startNode] = 0;          //distance from starting node to starting node is 0
	BFS(startNode, distances);

	return distances;
}


/// <summary>
/// Does an iterative breadth-first search and maps the distances from starting node to a vector of ints given as parameter.
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
/// Computes number of connected components.
/// </summary>
int Graph::getNumberOfConnectedComponents() {
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
void Graph::DFS(int startNode, vector<int>& visited) {

	int currentNode;
	stack<int> toVisit;

	toVisit.push(startNode);

	// while there still are accesible nodes that were not visited
	while (toVisit.empty() != true) {

		currentNode = toVisit.top();
		toVisit.pop();

		if (visited[currentNode] == 0) {
			visited[currentNode] = 1;

			// iterate through the current node's neighbors
			for (auto neighboringNode : adjacencyLists[currentNode])
				if (visited[neighboringNode] == 0)
					toVisit.push(neighboringNode);
		}
	}
}

/// <summary>
/// Recursive DFS that pushes nodes to a stack when returning from them. Used in topological sort. 
/// </summary>
void Graph::DFS(int currentNode, vector<int>& visited, stack<int>& solution) {

	visited[currentNode] = 1;

	for (int neighboringNode : adjacencyLists[currentNode])
		if (!visited[neighboringNode])
			DFS(neighboringNode, visited, solution);

	//add 1 that was subtracted on read
	solution.push(currentNode + 1);
}

/// <summary>
/// Checks if a collection of degrees can form a graph.
/// </summary>
bool havelHakimi(vector<int> degrees) {
	int nrNodes = degrees.size();

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
/// Computes a topological sort of a directed graph.
/// Note: assumes graph is acyclic.
/// </summary>
vector<int> Graph::getTopologicalSort() {
	vector<int> sortedItems;
	stack<int> solution;

	if (!isDirected) {
		outputFile << "Can't compute topological sort of undirected graph";
		return sortedItems;
	}
	else {
		vector<int> visited(nrNodes, 0);
		for (int node = 0; node < nrNodes; node++)
			if (!visited[node])
				DFS(node, visited, solution);
	}

	//move items from solution stack to sorted vector
	while (!solution.empty()) {
		sortedItems.push_back(solution.top());
		solution.pop();
	}
	return sortedItems;
}

/// <summary>
/// Tarjan's algorithm used in strongly connected components computing.
/// </summary>
/// <param name="counter">auxiliary for mapping traversal order to visitIndex</param>
/// <param name="visitIndex">order of traversal in DFS graph forest</param>
/// <param name="lowestAncestor">index of lowest ancestor reachable by back edges</param>
/// <param name="connectedComponents">list where connected components are added when found</param>
/// <param name="notCompleted">collection of nodes that have been visited but are not part of strongly complete connected components yet; used as a stack</param>
/// <param name="onStack">hash containing all nodes that are currently on notCompleted stack</param>
void Graph::tarjanDFS(int currentNode, int& counter, vector<int>& visitIndex, vector<int>& lowestAncestor, list<vector<int>>& connectedComponents, stack<int>& notCompleted, unordered_set<int>& onStack) {

	visitIndex[currentNode] = counter++;
	lowestAncestor[currentNode] = visitIndex[currentNode];

	notCompleted.push(currentNode);
	onStack.insert(currentNode);

	for (int neighboringNode : adjacencyLists[currentNode])
		if (visitIndex[neighboringNode] == -1) {
			tarjanDFS(neighboringNode, counter, visitIndex, lowestAncestor, connectedComponents, notCompleted, onStack);

			if (lowestAncestor[neighboringNode] < lowestAncestor[currentNode])
				lowestAncestor[currentNode] = lowestAncestor[neighboringNode];

		}
		else if (onStack.find(neighboringNode) != onStack.end() && lowestAncestor[neighboringNode] < lowestAncestor[currentNode])
			lowestAncestor[currentNode] = lowestAncestor[neighboringNode];

	if (lowestAncestor[currentNode] == visitIndex[currentNode]) {
		vector<int> connectedComponent;
		do {
			//add 1 subtracted on read
			connectedComponent.push_back(notCompleted.top() + 1);
			notCompleted.pop();
			onStack.erase(connectedComponent.back() - 1);
		} while (connectedComponent.back() - 1 != currentNode);

		connectedComponents.push_back(connectedComponent);
	}
}

/// <summary>
/// Prints number of strongly connected components and the members of each component on separate lines.
/// </summary>
list<vector<int>> Graph::getStronglyConnectedComponents() {

	list<vector<int>> connectedComponents;
	vector<int> visitIndex(nrNodes, -1);  //order of traversal in DFS graph forest
	vector<int> lowestAncestor(nrNodes, 0); //index of lowest ancestor reachable by back edges
	stack<int> notCompleted; //contains visited elements that are not part of completed strongly connected components
	unordered_set<int> onStack; //contains all elements that are currently on notCompleted stack
	int counter = 0;

	if (!isDirected) {
		outputFile << "The term \"strongly connected component\" applies only to directed graphs";
		return connectedComponents;
	}
	else for (int node = 0; node < nrNodes; node++)
		if (visitIndex[node] == -1) {
			tarjanDFS(node, counter, visitIndex, lowestAncestor, connectedComponents, notCompleted, onStack);
		}

	return connectedComponents;
}

int nodeNr, edgeNr;

int main()
{
	inputFile >> nodeNr >> edgeNr;

	Graph graph(nodeNr, true);
	graph.readEdges(edgeNr);

	list<vector<int>> stronglyConnectedComponents = graph.getStronglyConnectedComponents();
	outputFile << stronglyConnectedComponents.size() << "\n";
	for (auto elem : stronglyConnectedComponents)
		printVector(elem);
	/* MAIN HAVEL HAKIMI
	//files: havelhakimi.in, havelhakimi.out
	inputFile >> nrNoduri;

	Graph graph(nrNoduri, false);

	vector<int> degrees(nrNoduri);

	for (int i = 0; i < nrNoduri; i++) {
		inputFile >> degrees[i];
	}

	if (havelHakimi(degrees)) {
		outputFile << "Degrees can form a graph.";
	}
	else {
		outputFile << "Degrees can't form a graph.";
	}
	*/
}

void countingSort(vector<int>& toSort, int maxSize) {
	if (maxSize > 10000000)
		maxSize = 100000000;

	vector<int> count(maxSize, 0);
	int maxElem = -1;

	for (int i = 0; i < toSort.size(); i++) {
		count[toSort[i]]++;

		//check for maximum number to cut down execution time since mostly small numbers are expected
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
	outputFile << "\n";
}
