
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

/// <summary>
/// Initializes .in and .out files with a given filename.
/// </summary>
void initializeFiles(const string filename);
//--------------------------------------------------------------

int nr, nodeNr, edgeNr, sursa, input;
bool check;
int predictedSolution[50000];
int actualSolution[50000];
struct Muchie {
	int inNode, weight;
	Muchie(int x, int y) : inNode(x), weight(y) {};
};
vector<vector<Muchie>> adjacencyLists;

void readEdges() {
	int outNode, inNode, weight;

	for (int i = 0; i < edgeNr; i++) {
		inputFile >> outNode >> inNode >> weight;
		outNode--;
		inNode--;

		adjacencyLists[outNode].push_back(Muchie(inNode, weight));
		adjacencyLists[inNode].push_back(Muchie(outNode, weight));
	}
}

void getLightestPathFromNode(const int startNode) {
	const int maxValue = 1000000000;

	vector<bool> checkedNodes(50000, false);

	struct Node {
		int node, weight;
		Node(int n, int w) : node(n), weight(w) {};
		bool operator<(const Node& ob) const {
			return weight > ob.weight;
		}
	}; //auxiliary node struct used in closestNode heap

	priority_queue<Node> closestNode;
	int currentNode;

	for (int i = 0; i < nodeNr; i++)
		actualSolution[i] = maxValue;

	currentNode = startNode;
	actualSolution[currentNode] = 0;
	closestNode.push(Node(currentNode, 0));

	while (closestNode.size() > 0)
		if (checkedNodes[closestNode.top().node]) //if closest node is already checked pop from stack
			closestNode.pop();
		else {
			currentNode = closestNode.top().node;
			for (auto neighboringNode : adjacencyLists[currentNode]) //check all neighbors and update distance
				if (actualSolution[neighboringNode.inNode] > actualSolution[currentNode] + neighboringNode.weight) {
					actualSolution[neighboringNode.inNode] = actualSolution[currentNode] + neighboringNode.weight;
					closestNode.push(Node(neighboringNode.inNode, actualSolution[neighboringNode.inNode]));
				}

			checkedNodes[currentNode] = true;
		}

	for (int i = 0; i < nodeNr; i++)
		if (actualSolution[i] == maxValue)
			actualSolution[i] = 0;

}

int main()
{

	initializeFiles("distante");

	adjacencyLists = vector<vector<Muchie>>(50000, vector<Muchie>());

	inputFile >> nr;
	for (int i = 0; i < nr; i++) {
		inputFile >> nodeNr >> edgeNr >> sursa;

		for (int j = 0; j < nodeNr; j++) {
			inputFile >> predictedSolution[j];
		}

		if (i > 0) {
			//clear adjacency lists
			for (int i = 0; i < nodeNr; i++)
				adjacencyLists[i].clear();
		}

		readEdges();

		getLightestPathFromNode(sursa - 1);

		check = true;
		for (int j = 0; j < nodeNr; j++)
			if (actualSolution[j] != predictedSolution[j]) {
				check = false;
				break;
			}

		if (check)
			outputFile << "DA\n";
		else
			outputFile << "NU\n";
	}

}

void initializeFiles(const string filename) {
	inputFile = ifstream(filename + ".in");
	outputFile = ofstream(filename + ".out");
}