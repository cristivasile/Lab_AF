#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

ifstream inputFile;
ofstream outputFile;

void initializeFiles(const string filename) {
    inputFile = ifstream(filename + ".in");
    outputFile = ofstream(filename + ".out");
}

struct Vertex {
    int index, price;
    Vertex(int x, int y) : index(x), price(y) {};
    bool operator< (Vertex const& obj) {
        return this->price < obj.price;
    }
};

int nrNodes, nrEdges, notVisited;
vector<Vertex> vertexes;
vector<vector<int>> adjacencyList;
vector<int> minimumPrice;

void initializeStorage() {
    minimumPrice = vector<int>(nrNodes, -1);
    adjacencyList = vector<vector<int>>(nrNodes, vector<int>());
    notVisited = nrNodes;
}

void read() {

    int aux, inNode, outNode;

    for (int i = 0; i < nrNodes; i++) {
        inputFile >> aux;
        vertexes.push_back(Vertex(i, aux));
    }

    for (int i = 0; i < nrEdges; i++) {
        inputFile >> inNode >> outNode;
        inNode--;
        outNode--;
        adjacencyList[outNode].push_back(inNode);
    }
}

template <class T>
void printVector(ostream& outputFile, vector<T> toPrint, const int startIndex = 0) {
    for (int i = startIndex; i < toPrint.size(); i++)
        outputFile << toPrint[i] << " ";
    outputFile << "\n";
}

void DFS(const int currentNode, const int price) {
    if (minimumPrice[currentNode] == -1) {

        notVisited--;
        minimumPrice[currentNode] = price;

        for (auto& neighboringNode : adjacencyList[currentNode])
            DFS(neighboringNode, price);
    }
}

int main()
{
    initializeFiles("srevni");

    inputFile >> nrNodes >> nrEdges;

    initializeStorage();
    read();
    sort(vertexes.begin(), vertexes.end());

    int index = 0;
    while (index < nrNodes && notVisited > 0) {

        DFS(vertexes[index].index, vertexes[index].price);
        index++;
    }

    printVector(outputFile, minimumPrice);
}