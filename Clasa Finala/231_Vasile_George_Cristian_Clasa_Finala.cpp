#include <fstream>
#include <algorithm>
#include <queue>
#include <iterator>
#include <vector>
#include <stack>
#include <list>
#include <unordered_set>

using namespace std;

ifstream inputFile;
ofstream outputFile;

constexpr int maxSize = 100001;

//Auxiliary functions not related to graphs, defined after main.

/// <summary>
/// Basic descending counting sort algorithm.
/// </summary>
/// <param name="size"> - maximum size of sorted numbers; maximum value = 100.000.000</param>
void countingSort(vector<int>& toSort, int maxSize = 1000);


/// <param name="startIndex"> - can be used to exclude first "startIndex" elements from print</param>
template <class T>
void printVector(ostream& outputFile, vector<T> toPrint, int startIndex = 0);


/// <summary>
/// Initializes .in and .out files with a given filename.
/// </summary>
void initializeFiles(const string filename);
//--------------------------------------------------------------

class Disjoint_sets {

	vector <int> parent;
	vector <int> height;

	void mergeSets(const int, const int);
public:
	//creates n disjoint sets
	Disjoint_sets(const int);
	//finds the root of a node
	int findRoot(int);
	//merges two sets containing given nodes
	void mergeContainingSets(const int, const int);
};

Disjoint_sets::Disjoint_sets(const int nrNodes) {
	parent = vector<int>(nrNodes, -1);
	height = vector<int>(nrNodes, 1);
}

/// <summary>
/// Finds root of given node.
/// </summary>
int Disjoint_sets::findRoot(int node) {

	int lastNode;
	int root = node;

	//go up the tree
	while (parent[root] != -1)
		root = parent[root];

	//path compression
	while (parent[node] != -1) {
		lastNode = node;
		node = parent[node];
		parent[lastNode] = root;
	}

	return root;
}

/// <summary>
/// Merges two disjoint sets.
/// </summary>
void Disjoint_sets::mergeSets(const int root1, const int root2) {

	if (height[root1] > height[root2])
		parent[root2] = root1;
	else {
		parent[root1] = root2;
		//if both sets had same height increment resulting set's height 
		if (height[root1] == height[root2])
			height[root2]++;
	}
}

/// <summary>
/// Searches for the roots of two given nodes and merges the sets the nodes are contained in
/// </summary>
void Disjoint_sets::mergeContainingSets(const int node1, const int node2) {
	mergeSets(findRoot(node1), findRoot(node2));
}

class Graph {

	int nrNodes, nrEdges;
	bool isDirected, isWeighted;
	vector<vector<int>> adjacencyLists;
	vector<vector<int>> adjacencyMatrix;
	struct WeightedEdge { int node, weight; WeightedEdge(int n, int w) : node(n), weight(w) {} };
	vector<vector<WeightedEdge>> weightedAdjacencyLists;
	friend Disjoint_sets;

	void initializeAdjacencyLists();
	void initializeAdjacencyMatrix();
	void addEdge(const int, const int);
	void addWeightedEdge(const int, const int, const int);
	void BFS(const int, vector<int>&);
	void DFS(const int, vector<bool>&);
	void DFS(const int, vector<bool>&, stack<int>&);
	void tarjanDFS_SCC(const int, int&, vector<int>&, vector<int>&, stack<int>&, unordered_set<int>&, list<vector<int>>&);
	void tarjanDFS_CC(const int, const int, int&, vector<int>&, vector<int>&, vector<vector<int>>&);
	void tarjanDFS_BC(const int, const int, int&, vector<int>&, vector<int>&, stack<int>&, list<vector<int>>&);

public:
	//creates a graph with nrNodes, isDirected and isWeighted
	Graph(const int = maxSize, const bool = false, const bool = false);
	//reads edges from given input stream
	void readEdges(istream&, const int, const bool);
	//reads edges from given input stream, used only for eulerian cycle computation
	void readAdjacencyMatrix(istream&);
	//returns a vector of distances from a node given as parameter to all nodes in the graph
	vector<int> getDistancesToNode(const int);
	//returns a topological sort if the graph is directed
	vector<int> getTopologicalSort();
	//returns all SCCs in the graph if it is directed
	list<vector<int>> getStronglyConnectedComponents();
	//returns all CCs (bridges) in the graph
	vector<vector<int>> getCriticalConnections();
	//returns all biconnected components in the graph if it is directed
	list<vector<int>> getBiconnectedComponents();
	//returns the number of connected components
	int getNumberOfConnectedComponents();
	//returns length of longest chain in the tree
	int getTreeDiameter();
	//checks if collection of degrees can form a graph
	friend bool havelHakimi(vector<int>);
	//returns a MST of the graph
	vector<vector<int>> getMinimumSpanningTree(int&);
	//returns minimum weight path from a node to all others if all edges have positive weights (Djikstra)
	vector<int> getLightestPathFromNode(const int);
	//checks if graph has negative weight cycles, and if not then returns minimum weight path (Bellman-Ford)
	vector<int> getLightestPathFromNode(const int, bool&);
	//returns a matrix of distances between all pairs of nodes
	vector<vector<int>> getDistancesBetweenNodes();
	//checks if graph is eulerian, and if it is it returns an eulerian cycle using the parameter
	bool isEulerian(vector<int>&);
};


/// <param name="nrNodes"> - number of vertexes</param>
/// <param name="isDirected"> - is the graph directed?</param>
/// <param name="isWeighted"> - is the graph weighted?</param>
Graph::Graph(const int nrNodes, const bool isDirected, const bool isWeighted) {
	this->nrNodes = nrNodes;
	this->isDirected = isDirected;
	this->isWeighted = isWeighted;
	this->nrEdges = 0;
}


/// <summary>
/// initializes empty adjacency lists (vectors)
/// </summary>
void Graph::initializeAdjacencyLists() {

	adjacencyLists = vector<vector<int>>(nrNodes, vector<int>());

	//!isDirected for eulerian cycle computation
	if (isWeighted || !isDirected)
		weightedAdjacencyLists = vector<vector<WeightedEdge>>(nrNodes, vector<WeightedEdge>());
}


/// <summary>
/// initializes an adjacency matrix filled with 0
/// </summary>
void Graph::initializeAdjacencyMatrix() {
	adjacencyMatrix = vector<vector<int>>(nrNodes, vector<int>(nrNodes, 0));
}

/// <summary>
/// Adds an edge to the adjacency list.
/// </summary>
void Graph::addEdge(const int outNode, const int inNode) {
	adjacencyLists[outNode].push_back(inNode);
}

/// <summary>
/// Adds a weighted edge to the weightedAdjacency list.
/// </summary>
void Graph::addWeightedEdge(const int outNode, const int inNode, const int weight) {
	weightedAdjacencyLists[outNode].push_back(WeightedEdge(inNode, weight));
}

/// <summary>
/// Reads a number of edges given as parameter.
/// </summary>
void Graph::readEdges(istream& inputFile, const int nrEdges, const bool forEulerian = false) {
	int inNode, outNode, weight;

	this->nrEdges = nrEdges;

	initializeAdjacencyLists();

	if (!forEulerian) {
		if (!isWeighted) {
			if (isDirected) {
				for (int i = 0; i < nrEdges; i++) {
					inputFile >> outNode >> inNode;
					outNode--;
					inNode--;
					addEdge(outNode, inNode);
				}
			}
			else {
				for (int i = 0; i < nrEdges; i++) {
					inputFile >> outNode >> inNode;
					outNode--;
					inNode--;

					addEdge(outNode, inNode);
					addEdge(inNode, outNode);
				}
			}
		}
		else {
			if (isDirected) {
				for (int i = 0; i < nrEdges; i++) {
					inputFile >> outNode >> inNode >> weight;
					outNode--;
					inNode--;

					addEdge(outNode, inNode);
					addWeightedEdge(outNode, inNode, weight);
				}
			}
			else {
				for (int i = 0; i < nrEdges; i++) {
					inputFile >> outNode >> inNode >> weight;
					outNode--;
					inNode--;

					addEdge(outNode, inNode);
					addWeightedEdge(outNode, inNode, weight);
					addEdge(inNode, outNode);
					addWeightedEdge(inNode, outNode, weight);
				}
			}
		}
	}
	else {
		for (int i = 0; i < nrEdges; i++) {
			inputFile >> outNode >> inNode;
			outNode--;
			inNode--;

			//if the graph is undirected and unweighted then we can use the unused "weight" property as an unique identifier
			//used for isEulerian function
			addWeightedEdge(inNode, outNode, i);
			addWeightedEdge(outNode, inNode, i);
		}
	}
}


/// <summary>
/// Reads the graph's adjacency matrix.
/// </summary>
void Graph::readAdjacencyMatrix(istream& inputFile) {

	initializeAdjacencyMatrix();
	for (int i = 0; i < nrNodes; i++)
		for (int j = 0; j < nrNodes; j++)
			inputFile >> adjacencyMatrix[i][j];

}

/// <summary>
/// Returns a vector of distances from a node given as parameter to all nodes in the graph, or -1 if they are inaccesible from that node.
/// </summary>
/// <param name="startNode">- the node for which distances are calculated</param>
vector<int> Graph::getDistancesToNode(const int startNode) {

	vector<int> distances(nrNodes, -1);
	distances[startNode] = 0;          //distance from starting node to starting node is 0
	BFS(startNode, distances);

	return distances;
}


/// <summary>
/// Does an iterative breadth-first search and maps the distances from starting node to a vector of ints given as parameter.
/// </summary>
/// <param name="startNode">- starting node of search </param>
/// <param name="distances">- vector to map distances to </param>
void Graph::BFS(const int startNode, vector<int>& distances) {
	int currentNode, currentDistance;
	queue<int> toVisit;
	toVisit.push(startNode);

	while (toVisit.empty() != true) {
		currentNode = toVisit.front();
		currentDistance = distances[currentNode];

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
	vector<bool> visited(nrNodes, 0); //all nodes ar unvisited at start

	//go through all nodes
	for (int i = 0; i < nrNodes; i++)
		//if there is an unvisited node do a DFS and increment counter
		if (visited[i] == 0) {
			nr++;
			DFS(i, visited);
		}
	return nr;
}

/// <summary>
/// Iterative depth-first traversal that maps visited nodes in a vector of ints given as parameter.
/// </summary>
void Graph::DFS(const int startNode, vector <bool>& visited) {

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
/// Recursive DFS that pushes nodes to a stack when returning from them. NOTE: it adds one to every node 
/// </summary>
void Graph::DFS(const int currentNode, vector<bool>& visited, stack<int>& solution) {

	visited[currentNode] = 1;

	for (int neighboringNode : adjacencyLists[currentNode])
		if (!visited[neighboringNode])
			DFS(neighboringNode, visited, solution);

	//add 1 that was subtracted on read
	solution.push(currentNode + 1);
}

/// <summary>
/// Checks if a sequence of degrees can form a graph.
/// </summary>
bool havelHakimi(vector<int> degrees) {
	int nrNodes = degrees.size();
	int currentIndex, nrConnectedNodes;

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

	currentIndex = 0;

	while (degrees[currentIndex] != 0) {
		nrConnectedNodes = 0;
		//for the next degrees[currentIndex] nodes, connect current node by subtracting one
		for (int i = currentIndex + 1; i < nrNodes && nrConnectedNodes < degrees[currentIndex]; i++)
			if (degrees[i] != 0) {
				degrees[i]--;
				nrConnectedNodes++;
			}

		if (degrees[currentIndex] != nrConnectedNodes) // if there were not enough nodes to connect to, degrees can't form a graph
			return 0;

		//the current node was connected the required number of times and is set to 0
		degrees[currentIndex++] = 0;

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
		//throw x;
		return sortedItems;
	}
	else {
		vector<bool> visited(nrNodes, 0);
		for (int node = 0; node < nrNodes; node++)
			if (!visited[node])
				DFS(node, visited, solution);
	}

	//reverse solution order by moving items from solution stack to vector 
	while (!solution.empty()) {
		sortedItems.push_back(solution.top());
		solution.pop();
	}
	return sortedItems;
}

/// <summary>
/// Tarjan's algorithm, used in finding strongly connected components.
/// </summary>
/// <param name="counter">- auxiliary used in mapping traversal order to visitIndex</param>
/// <param name="visitIndex">- order of traversal in DFS graph forest</param>
/// <param name="lowestAncestor">- index of lowest ancestor reachable by back edges</param>
/// <param name="stronglyConnectedComponents">- list where connected components are added when found</param>
/// <param name="notCompleted">- collection of nodes that have been visited but are not part of strongly complete connected components yet</param>
/// <param name="onStack">- hash containing all nodes that are currently on notCompleted stack</param>
void Graph::tarjanDFS_SCC(const int currentNode, int& counter, vector<int>& visitIndex, vector<int>& lowestAncestorReachable, stack<int>& notCompleted, unordered_set<int>& onStack, list<vector<int>>& stronglyConnectedComponents) {

	visitIndex[currentNode] = counter++;
	lowestAncestorReachable[currentNode] = visitIndex[currentNode];

	notCompleted.push(currentNode);
	onStack.insert(currentNode);

	for (int neighboringNode : adjacencyLists[currentNode])
		if (visitIndex[neighboringNode] == -1) {

			tarjanDFS_SCC(neighboringNode, counter, visitIndex, lowestAncestorReachable, notCompleted, onStack, stronglyConnectedComponents);

			if (lowestAncestorReachable[neighboringNode] < lowestAncestorReachable[currentNode])
				lowestAncestorReachable[currentNode] = lowestAncestorReachable[neighboringNode];
		}
	//must check if node is on notCompleted stack to not consider cross edges 
		else if (lowestAncestorReachable[neighboringNode] < lowestAncestorReachable[currentNode] && onStack.find(neighboringNode) != onStack.end())
			lowestAncestorReachable[currentNode] = lowestAncestorReachable[neighboringNode];

	//if current node is the root of a strongly connected component
	if (lowestAncestorReachable[currentNode] == visitIndex[currentNode]) {
		vector<int> newConnectedComponent;
		//remove the node and all of it's successors from stack and add them to a new connected component
		do {
			//add 1 subtracted on read
			newConnectedComponent.push_back(notCompleted.top() + 1);
			notCompleted.pop();
			onStack.erase(newConnectedComponent.back() - 1);
		} while (newConnectedComponent.back() - 1 != currentNode);

		stronglyConnectedComponents.push_back(newConnectedComponent);
	}
}

/// <summary>
/// Computes strongly connected components in the graph.
/// </summary>
list<vector<int>> Graph::getStronglyConnectedComponents() {

	list<vector<int>> stronglyConnectedComponents;
	vector<int> visitIndex(nrNodes, -1);  //order of traversal in DFS graph forest
	vector<int> lowestAncestorReachable(nrNodes, 0); //index of lowest ancestor reachable by back edges
	stack<int> notCompleted; //contains visited elements that are not part of completed strongly connected components
	unordered_set<int> onStack; //contains all elements that are currently on notCompleted stack
	int counter = 0; //auxiliary used in mapping traversal order to visitIndex

	if (!isDirected) {
		outputFile << "The term \"strongly connected component\" exists only in the context of directed graphs";
		//throw x;
		return stronglyConnectedComponents;
	}
	else for (int node = 0; node < nrNodes; node++)
		if (visitIndex[node] == -1) {
			tarjanDFS_SCC(node, counter, visitIndex, lowestAncestorReachable, notCompleted, onStack, stronglyConnectedComponents);
		}

	return stronglyConnectedComponents;
}

/// <summary>
/// Tarjan's DFS used in searching for critical connections.
/// </summary>
/// <param name="counter">- auxiliary used in mapping traversal order to visitIndex</param>
/// <param name="visitIndex">- order of traversal in DFS graph forest</param>
/// <param name="lowestAncestor">- index of lowest ancestor reachable by back edges</param>
/// <param name="parent">- parent node (in DFS tree) of current node</param>
/// <param name="criticalConnections">- critical connections container</param>
void Graph::tarjanDFS_CC(const int currentNode, const int parent, int& counter, vector<int>& visitIndex, vector<int>& lowestAncestorReachable, vector<vector<int>>& criticalConnections) {

	visitIndex[currentNode] = counter++;
	lowestAncestorReachable[currentNode] = visitIndex[currentNode];

	for (int neighboringNode : adjacencyLists[currentNode])
		if (visitIndex[neighboringNode] == -1) {

			tarjanDFS_CC(neighboringNode, currentNode, counter, visitIndex, lowestAncestorReachable, criticalConnections);

			if (lowestAncestorReachable[neighboringNode] < lowestAncestorReachable[currentNode])
				lowestAncestorReachable[currentNode] = lowestAncestorReachable[neighboringNode];

			//critical connection found
			if (visitIndex[currentNode] < lowestAncestorReachable[neighboringNode])
				//add one subtracted on read
				criticalConnections.push_back({ currentNode + 1, neighboringNode + 1 });
		}
		else if (parent != neighboringNode && visitIndex[neighboringNode] < lowestAncestorReachable[currentNode])
			lowestAncestorReachable[currentNode] = visitIndex[neighboringNode];
}

/// <summary>
/// Computes and returns all critical connections (bridges) in the graph
/// </summary>
vector<vector<int>> Graph::getCriticalConnections() {
	vector<vector<int>> criticalConnections;
	vector<int> visitIndex(nrNodes, -1);  //order of traversal in DFS graph forest
	vector<int> lowestAncestorReachable(nrNodes, 0); //index of lowest ancestor reachable by back edges
	int counter = 0; //auxiliary used in mapping traversal order to visitIndex

	for (int node = 0; node < nrNodes; node++)
		if (visitIndex[node] == -1)
			tarjanDFS_CC(node, -1, counter, visitIndex, lowestAncestorReachable, criticalConnections);
	//parent of first node in DFS tree is -1

	return criticalConnections;
}

/// <summary>
/// Tarjan's DFS used in finding biconnected components. NOTE: only works with directed graphs.
/// </summary>
/// <param name="counter">- auxiliary used in mapping traversal order to visitIndex</param>
/// <param name="parent">- parent node (in DFS tree) of current node</param>
/// <param name="visitIndex">- order of traversal in DFS graph forest</param>
/// <param name="lowestAncestor">- index of lowest ancestor reachable by back edges</param>
/// <param name="notCompleted">- all visited nodes that are not yet part of completed biconnected components</param>
/// <param name="biconnectedComponents">- found biconnected components are stored here</param>
void Graph::tarjanDFS_BC(int currentNode, const int parent, int& counter, vector<int>& visitIndex, vector<int>& lowestAncestorReachable, stack<int>& notCompleted, list<vector<int>>& biconnectedComponents) {

	visitIndex[currentNode] = counter++;
	lowestAncestorReachable[currentNode] = visitIndex[currentNode];

	notCompleted.push(currentNode);

	for (int neighboringNode : adjacencyLists[currentNode]) {
		bool firstChildFound = false;

		if (visitIndex[neighboringNode] == -1) {

			tarjanDFS_BC(neighboringNode, currentNode, counter, visitIndex, lowestAncestorReachable, notCompleted, biconnectedComponents);

			if (lowestAncestorReachable[neighboringNode] < lowestAncestorReachable[currentNode])
				lowestAncestorReachable[currentNode] = lowestAncestorReachable[neighboringNode];

			//articulation point cases:
			//1 - root node has multiple children
			//2 - non-root node has a child that has no back edge reaching to one of it's ancestors
			if ((parent == -1 && firstChildFound) || (parent != -1 && lowestAncestorReachable[neighboringNode] >= visitIndex[currentNode])) {
				vector<int> newBiconnectedComponent;

				//add articulation point at the end of stack
				notCompleted.push(currentNode);

				//remove articulation point's successors and add them to the stack
				do {
					//add 1 subtracted on read
					newBiconnectedComponent.push_back(notCompleted.top() + 1);
					notCompleted.pop();
				} while (newBiconnectedComponent.back() - 1 != neighboringNode); //until articulation point's child is reached

				biconnectedComponents.push_back(newBiconnectedComponent);
			}

			//check if root has multiple children
			if (parent == -1)
				firstChildFound = true;
		}
		else if (parent != neighboringNode && visitIndex[neighboringNode] < lowestAncestorReachable[currentNode]) //back edge found
			lowestAncestorReachable[currentNode] = visitIndex[neighboringNode];
	}

	//first biconnected component containing root node will be left on stack
	if (parent == -1 && notCompleted.size() != 0) {
		vector<int> newBiconnectedComponent;
		//remove all the node's successors from stack and add them to a new biconnected component
		do {
			//add 1 subtracted on read
			newBiconnectedComponent.push_back(notCompleted.top() + 1);
			notCompleted.pop();
		} while (notCompleted.size() != 0);

		biconnectedComponents.push_back(newBiconnectedComponent);
	}
}

/// <summary>
/// Computes and returns all biconnected components in the graph. NOTE: Only works with undirected graphs.
/// </summary>
list<vector<int>> Graph::getBiconnectedComponents() {
	list<vector<int>> biconnectedComponents;
	vector<int> visitIndex(nrNodes, -1);  //order of traversal in DFS graph forest
	vector<int> lowestAncestorReachable(nrNodes, 0); //index of lowest ancestor reachable by back edges
	stack<int> notCompleted; //contains all visited elements that are not yet part of a completed biconnected component
	int counter = 0; //auxiliary used in mapping traversal order to visitIndex

	if (isDirected) {
		outputFile << "This function only works with directed graphs.";
		//throw x
		return biconnectedComponents;
	}
	else {
		for (int i = 0; i < nrNodes; i++)
			if (visitIndex[i] == -1)
				tarjanDFS_BC(i, -1, counter, visitIndex, lowestAncestorReachable, notCompleted, biconnectedComponents); 		//parent of first node in DFS tree is -1

		return biconnectedComponents;
	}
}

/// <summary>
/// Prim's algorithm used in determining minimum spanning tree.
/// </summary>
/// <param name="totalWeight"> - total weight of MST edges returned in this parameter</param>
/// <returns>Vector of edges (stored in vectors) that make up the MST.</returns>
vector<vector<int>> Graph::getMinimumSpanningTree(int& totalWeight) {

	struct Edge {
		int outNode, inNode, weight;
		Edge(int o, int i, int w) : outNode(o), inNode(i), weight(w) {};
		bool operator<(const Edge& ob) const {
			return weight > ob.weight;
		}
	};

	vector<vector<int>> minimumSpanningTree;
	vector<bool> isAdded(nrNodes, false);
	priority_queue<Edge> lightestEdge;

	if (!isWeighted) {
		outputFile << "Graph must be weighted.";
		//throw x
		return minimumSpanningTree;
	}
	else {
		totalWeight = 0;
		int addedNodes = 1;
		int lastAddedNode = 0;
		isAdded[0] = true;

		while (addedNodes != nrNodes) {
			//add all edges of last added node that lead to undiscovered nodes to heap
			for (auto edge : weightedAdjacencyLists[lastAddedNode])
				if (!isAdded[edge.node])
					lightestEdge.push(Edge(lastAddedNode, edge.node, edge.weight));

			//while lightest edge leads to already added node, pop
			while (isAdded[lightestEdge.top().inNode])
				lightestEdge.pop();

			//add lightest edge to MST and advance to destination node
			vector<int> newEdge;
			newEdge.push_back(lightestEdge.top().outNode + 1);
			newEdge.push_back(lightestEdge.top().inNode + 1);
			minimumSpanningTree.push_back(newEdge);

			totalWeight += lightestEdge.top().weight;

			lastAddedNode = lightestEdge.top().inNode;
			isAdded[lastAddedNode] = true;
			addedNodes++;
		}
		return minimumSpanningTree;
	}
}

/// <summary>
/// Dijkstra's algorithm used in determining lightest path from starting node to all others. Note: only works if all edges have positive weights.
/// </summary>
vector<int> Graph::getLightestPathFromNode(const int startNode) {

	const int maxValue = 1000000000;

	struct Node {
		int node, weight;
		Node(int n, int w) : node(n), weight(w) {};
		bool operator<(const Node& ob) const {
			return weight > ob.weight;
		}
	}; //auxiliary node struct used in closestNode heap

	unordered_set <int> checkedNodes;
	vector<int> lightestPath(nrNodes, maxValue); //we assume all nodes are inaccessible on initialization
	priority_queue<Node> closestNode;
	int currentNode;

	currentNode = startNode;
	lightestPath[currentNode] = 0;
	closestNode.push(Node(currentNode, 0));

	while (closestNode.size() > 0)
		if (checkedNodes.find(closestNode.top().node) != checkedNodes.end()) //if closest node is already checked pop from stack
			closestNode.pop();
		else {
			currentNode = closestNode.top().node;
			for (auto neighboringNode : weightedAdjacencyLists[currentNode]) //check all neighbors and update distance
				if (lightestPath[neighboringNode.node] > lightestPath[currentNode] + neighboringNode.weight) {
					lightestPath[neighboringNode.node] = lightestPath[currentNode] + neighboringNode.weight;
					closestNode.push(Node(neighboringNode.node, lightestPath[neighboringNode.node]));
				}

			checkedNodes.insert(currentNode);
		}

	for (int i = 0; i < nrNodes; i++)
		if (lightestPath[i] == maxValue)
			lightestPath[i] = 0;

	return lightestPath;
}

/// <summary>
/// Bellman-Ford algorithm used in determining lightest path from starting node to all others.
/// </summary>
/// <param name="hasNegativeCycles"> - used to check if graph has negative weight cycles</param>
vector<int> Graph::getLightestPathFromNode(const int startNode, bool& hasNegativeCycles) {

	const int maxValue = 1000000000;

	struct Node {
		int node, weight;
		Node(int n, int w) : node(n), weight(w) {};
		bool operator<(const Node& ob) const {
			return weight > ob.weight;
		}
	}; //auxiliary node struct used in closestNode heap

	vector<int> lightestPath(nrNodes, maxValue); //we assume all nodes are inaccessible on initialization
	vector<int> checkedNodes(nrNodes, 0); //if a node is checked nrNodes times the graph has a negative-weight cycle
	queue<Node> nodeQueue;
	int currentNode, currentWeight;

	hasNegativeCycles = false;

	lightestPath[startNode] = 0;
	nodeQueue.push(Node(startNode, 0));

	while (nodeQueue.empty() == false) {

		currentNode = nodeQueue.front().node;
		currentWeight = nodeQueue.front().weight;
		nodeQueue.pop();

		if (lightestPath[currentNode] == currentWeight) {
			checkedNodes[currentNode] ++;

			if (checkedNodes[currentNode] == nrNodes) {
				hasNegativeCycles = true;
				break;
			}

			for (auto& edge : weightedAdjacencyLists[currentNode])
				if (lightestPath[edge.node] > currentWeight + edge.weight) {
					lightestPath[edge.node] = currentWeight + edge.weight;
					nodeQueue.push(Node(edge.node, lightestPath[edge.node]));
				}
		}

	}
	return lightestPath;
}

/// <summary>
/// Computes and returns the length of the longest chain in the tree.
/// </summary>
int Graph::getTreeDiameter() {

	vector<int> distances = getDistancesToNode(0);

	int max = -1, farthestNode = 0;
	for (int i = 0; i < nrNodes; i++)
		if (distances[i] > max) {
			max = distances[i];
			farthestNode = i;
		}

	distances = getDistancesToNode(farthestNode);

	max = -1;
	for (int i = 0; i < nrNodes; i++)
		if (distances[i] > max)
			max = distances[i];

	// if distance from A to B is 7, the chain length is 8
	return max + 1;
}

/// <summary>
/// Computes and returns a matrix containing distances between each pair of nodes in the graph. (Floyd-Warshall/Roy-Floyd)
/// </summary>
vector<vector<int>> Graph::getDistancesBetweenNodes() {
	vector<vector<int>> distanceMatrix;
	distanceMatrix = adjacencyMatrix;

	for (int k = 0; k < nrNodes; k++)
		for (int i = 0; i < nrNodes; i++)
			for (int j = 0; j < nrNodes; j++)
				if (distanceMatrix[i][k] != 0 && distanceMatrix[k][j] != 0
					&& (distanceMatrix[i][j] > distanceMatrix[i][k] + distanceMatrix[k][j]
						|| (i != j && distanceMatrix[i][j] == 0)))

					distanceMatrix[i][j] = distanceMatrix[i][k] + distanceMatrix[k][j];
	return distanceMatrix;
}

/// <summary>
/// Checks if a graph is eulerian and if it is, it returns one eulerian cycle using a parameter.
/// </summary>
/// <param name="eulerianCycle"> - if the graph is eulerian, one eulerian cycle will be returned using this parameter</param>
/// <returns> - bool representing if the graph is eulerian or not </returns>
bool Graph::isEulerian(vector<int>& eulerianCycle) {

	stack<int> stack;
	vector<bool> usedEdges(nrEdges, false);
	vector<int> startIndex(nrNodes, 0);
	int currentNode;

	for (int i = 0; i < nrNodes; i++)
		//if at least one vertex has an odd degree then the graph is not eulerian
		if (weightedAdjacencyLists[i].size() % 2 != 0)
			return false;

	//otherwise, the algorithm searches for an eulerian cycle
	stack.push(0);
	while (stack.size() != 0) {
		currentNode = stack.top();

		int i;
		for (i = startIndex[currentNode]; i < weightedAdjacencyLists[currentNode].size(); i++) {
			WeightedEdge edge = weightedAdjacencyLists[currentNode][i];
			//weight is not actually an edge's weight, but it's identifier
			if (!usedEdges[edge.weight]) {

				usedEdges[edge.weight] = true;
				stack.push(edge.node);
				//the first i edges have already been checked and will be ignored on next iteration
				startIndex[currentNode] = i + 1;

				break;
			}
		}

		//if there are no more edges to check, add current vertex to solution
		if (i == weightedAdjacencyLists[currentNode].size()) {
			stack.pop();
			eulerianCycle.push_back(currentNode + 1);
		}
	}

	//last vertex added will be the start node, which is not needed in our solution
	eulerianCycle.pop_back();

	return true;
}

int nodeNr, edgeNr;

int main()
{

	initializeFiles("dfs");

	inputFile >> nodeNr >> edgeNr;
	Graph graph(nodeNr, false, false);
	graph.readEdges(inputFile, edgeNr);
	outputFile << graph.getNumberOfConnectedComponents();


	/*
	//MAIN HAVEL HAKIMI
	initializeFiles("havelhakimi");
	inputFile >> nodeNr;

	Graph graph(nodeNr, false, false);

	vector<int> degrees(nodeNr);

	for (int i = 0; i < nodeNr; i++) {
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
	if (maxSize > 100000000)
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
void printVector(ostream& outputFile, vector<T> toPrint, int startIndex) {
	for (int i = startIndex; i < toPrint.size(); i++)
		outputFile << toPrint[i] << " ";
	outputFile << "\n";
}

void initializeFiles(const string filename) {
	inputFile = ifstream(filename + ".in");
	outputFile = ofstream(filename + ".out");
}