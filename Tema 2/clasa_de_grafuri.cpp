#include <fstream>
#include <algorithm>
#include <queue>
#include <iterator>
#include <vector>
#include <stack>
#include <list>
#include <unordered_set>

using namespace std;

ifstream inputFile("dijkstra.in");
ofstream outputFile("dijkstra.out");

constexpr int maxSize = 100001;

//Auxiliary functions not related to graphs, defined after main.
/// <summary>
/// Basic descending counting sort algorithm.
/// </summary>
/// <param name="size"> - maximum size of sorted numbers; maximum value = 100.000.000</param>
void countingSort(vector<int>& toSort, int maxSize = 1000);


/// <param name="startIndex"> - can be used to exclude first "startIndex" elements from print</param>
template <class T>
void printVector(vector<T> toPrint, int startIndex = 0);
//--------------------------------------------------------------

struct Edge {
	int outNode, inNode, weight;
	Edge(int o, int i, int w) : outNode(o), inNode(i), weight(w) {};
	bool operator<(const Edge& ob) const{
		return weight > ob.weight;
	}
};

class Graph {

	vector<vector<int>> adjacencyLists;
	struct WeightedEdge { int node, weight; WeightedEdge(int n, int w) : node(n), weight(w) {} };
	vector<vector<WeightedEdge>> weightedAdjacencyLists;
	int nrNodes;
	bool isDirected, isWeighted;

public:
	Graph(int, bool, bool);
	void readEdges(int);
	vector<int> getDistancesToNode(int);
	vector<int> getTopologicalSort();
	list<vector<int>> getStronglyConnectedComponents();
	vector<vector<int>> getCriticalConnections();
	list<vector<int>> getBiconnectedComponents();
	int getNumberOfConnectedComponents();
	friend bool havelHakimi(vector<int>);
	vector<vector<int>> getMinimumSpanningTree(int&);
	vector<int> getShortestPathFromNode(int);

private:
	void initializeAdjacencyLists();
	void BFS(int, vector<int>&);
	void DFS(int, vector<int>&);
	void DFS(int, vector<int>&, stack<int>&);
	void tarjanDFS_SCC(int, int&, vector<int>&, vector<int>&, list<vector<int>>&, stack<int>&, unordered_set<int>&);
	void tarjanDFS_CC(int, int&, vector<int>&, vector<int>&, int, vector<vector<int>>&);
	void tarjanDFS_BC(int, int, int&, vector<int>&, vector<int>&, stack<int>&, list<vector<int>>&);
};

Graph::Graph(int size = maxSize, bool isDirected = false, bool isWeighted = false) {
	nrNodes = size;
	this->isDirected = isDirected;
	this->isWeighted = isWeighted;
	initializeAdjacencyLists();
}


/// <summary>
/// initializes empty adjacency lists (vectors)
/// </summary>
void Graph::initializeAdjacencyLists() {
	if(!isWeighted)
		for (int i = 0; i < nrNodes; i++) {
			vector<int> emptyAdjacencyVector;
			adjacencyLists.push_back(emptyAdjacencyVector);
		}
	else
		for (int i = 0; i < nrNodes; i++) {
			vector<int> emptyAdjacencyVector;
			adjacencyLists.push_back(emptyAdjacencyVector);

			vector<WeightedEdge> emptyWeightVector;
			weightedAdjacencyLists.push_back(emptyWeightVector);
		}
}

/// <summary>
/// reads a number of edges given as parameter
/// </summary>
void Graph::readEdges(int nrEdges) {
	int inNode, outNode, weight;

	if (!isWeighted) {
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
	else {
		if (isDirected) {
			for (int i = 0; i < nrEdges; i++) {
				inputFile >> outNode >> inNode >> weight;
				outNode--;
				inNode--;
				adjacencyLists[outNode].push_back(inNode);
				weightedAdjacencyLists[outNode].push_back(WeightedEdge(inNode, weight));
			}
		}
		else {
			for (int i = 0; i < nrEdges; i++) {
				inputFile >> outNode >> inNode >> weight;
				outNode--;
				inNode--;

				adjacencyLists[outNode].push_back(inNode);
				weightedAdjacencyLists[outNode].push_back(WeightedEdge(inNode, weight));
				adjacencyLists[inNode].push_back(outNode);
				weightedAdjacencyLists[inNode].push_back(WeightedEdge(outNode, weight));
			}
		}
	}
}
/// <summary>
/// Returns a vector of distances from a node given as parameter to all nodes in the graph, or -1 if they are inaccesible from that node.
/// </summary>
/// <param name="startNode">- the node for which distances are calculated</param>
vector<int> Graph::getDistancesToNode(int startNode) {

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
void Graph::BFS(int startNode, vector<int>& distances) {
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
	vector<int> visited(nrNodes, 0); //all nodes ar unvisited at start

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
/// Recursive DFS that pushes nodes to a stack when returning from them.
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
/// Checks if a sequence of degrees can form a graph.
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
		//for the next degrees[0] nodes, connect current node by subtracting one
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
		//throw x;
		return sortedItems;
	}
	else {
		vector<int> visited(nrNodes, 0);
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
void Graph::tarjanDFS_SCC(int currentNode, int& counter, vector<int>& visitIndex, vector<int>& lowestAncestorReachable, list<vector<int>>& stronglyConnectedComponents, stack<int>& notCompleted, unordered_set<int>& onStack) {

	visitIndex[currentNode] = counter++;
	lowestAncestorReachable[currentNode] = visitIndex[currentNode];

	notCompleted.push(currentNode);
	onStack.insert(currentNode);

	for (int neighboringNode : adjacencyLists[currentNode])
		if (visitIndex[neighboringNode] == -1) {

			tarjanDFS_SCC(neighboringNode, counter, visitIndex, lowestAncestorReachable, stronglyConnectedComponents, notCompleted, onStack);

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
			tarjanDFS_SCC(node, counter, visitIndex, lowestAncestorReachable, stronglyConnectedComponents, notCompleted, onStack);
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
void Graph::tarjanDFS_CC(int currentNode, int& counter, vector<int>& visitIndex, vector<int>& lowestAncestorReachable, int parent, vector<vector<int>>& criticalConnections) {

	visitIndex[currentNode] = counter++;
	lowestAncestorReachable[currentNode] = visitIndex[currentNode];

	for (int neighboringNode : adjacencyLists[currentNode])
		if (visitIndex[neighboringNode] == -1) {

			tarjanDFS_CC(neighboringNode, counter, visitIndex, lowestAncestorReachable, currentNode, criticalConnections);

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
/// Computes and returns all critical connections (bridges) in graph
/// </summary>
vector<vector<int>> Graph::getCriticalConnections() {
	vector<vector<int>> criticalConnections;
	vector<int> visitIndex(nrNodes, -1);  //order of traversal in DFS graph forest
	vector<int> lowestAncestorReachable(nrNodes, 0); //index of lowest ancestor reachable by back edges
	int counter = 0; //auxiliary used in mapping traversal order to visitIndex

	for (int node = 0; node < nrNodes; node++)
		if (visitIndex[node] == -1)
			tarjanDFS_CC(node, counter, visitIndex, lowestAncestorReachable, -1, criticalConnections);
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
void Graph::tarjanDFS_BC(int currentNode, int parent, int& counter, vector<int>& visitIndex, vector<int>& lowestAncestorReachable, stack<int>& notCompleted, list<vector<int>>& biconnectedComponents) {

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

		while(addedNodes != nrNodes){
			//add all edges of last added node that lead to undiscovered nodes to heap
			for (auto edge : weightedAdjacencyLists[lastAddedNode])
				if(!isAdded[edge.node])
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
/// Dijkstra's algorithm used in determining shortest path to all nodes starting from a single point.
/// </summary>
vector<int> Graph::getShortestPathFromNode(int startNode) {

	const int maxValue = 1000000000;

	struct Node {
		int node, distance;
		Node(int n, int w) : node(n), distance(w) {};
		bool operator<(const Node& ob) const {
			return distance > ob.distance;
		}
	}; //auxiliary node struct used in closestNode heap

	unordered_set <int> checkedNodes;
	vector<int> minimumDistance(nrNodes, maxValue); //we assume all nodes are inaccessible on initialization
	priority_queue<Node> closestNode; 
	int currentNode;

	currentNode = startNode;
	minimumDistance[currentNode] = 0;
	closestNode.push(Node(currentNode, 0));

	while (closestNode.size() > 0)
		if (checkedNodes.find(closestNode.top().node) != checkedNodes.end()) //if closest node is already checked pop from stack
			closestNode.pop();
		else {
			currentNode = closestNode.top().node;
			for (auto neighboringNode : weightedAdjacencyLists[currentNode]) //check all neighbors and update distance
				if (minimumDistance[neighboringNode.node] > minimumDistance[currentNode] + neighboringNode.weight) {
					minimumDistance[neighboringNode.node] = minimumDistance[currentNode] + neighboringNode.weight;
					closestNode.push(Node(neighboringNode.node, minimumDistance[neighboringNode.node]));
				}
			
			checkedNodes.insert(currentNode);
		}

	for (int i = 0; i < nrNodes; i++)
		if (minimumDistance[i] == maxValue)
			minimumDistance[i] = 0;

	return minimumDistance;
}

int nodeNr, edgeNr;

int main()
{
	inputFile >> nodeNr >> edgeNr;

	Graph graph(nodeNr, true, true);
	graph.readEdges(edgeNr);

	printVector(graph.getShortestPathFromNode(0), 1);
	

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
void printVector(vector<T> toPrint, int startIndex) {
	for (int i = startIndex; i < toPrint.size(); i++)
		outputFile << toPrint[i] << " ";
	outputFile << "\n";
}
