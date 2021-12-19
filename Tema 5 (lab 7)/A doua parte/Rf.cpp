#include <fstream>
#include <iostream>
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

int nodeNr;
int distanceMatrix[256][256], frequencyMatrix[256][256];

void readAdjacencyMatrix() {
	for (int i = 0; i < nodeNr; i++)
		for (int j = 0; j < nodeNr; j++) {
			inputFile >> distanceMatrix[i][j];
			if (distanceMatrix[i][j] != 0)
				frequencyMatrix[i][j] = 1;
		}
}

void printMatrix(int matrix[256][256]) {
	for (int i = 0; i < nodeNr; i++) {
		for (int j = 0; j < nodeNr; j++)
			outputFile << matrix[i][j] << " ";
		outputFile << "\n";
	}
}

void RoyFloyd() {
	for (int k = 0; k < nodeNr; k++)
		for (int i = 0; i < nodeNr; i++)
			for (int j = i + 1; j < nodeNr; j++)
				if (distanceMatrix[i][j] > distanceMatrix[i][k] + distanceMatrix[k][j]) {

					frequencyMatrix[i][j] = frequencyMatrix[i][k] + frequencyMatrix[k][j];
					distanceMatrix[i][j] = distanceMatrix[i][k] + distanceMatrix[k][j];

					frequencyMatrix[j][i] = frequencyMatrix[i][j];
					distanceMatrix[j][i] = distanceMatrix[i][j];
				}
				else if (frequencyMatrix[i][j] < frequencyMatrix[i][k] + frequencyMatrix[k][j]
					&& distanceMatrix[i][j] == distanceMatrix[i][k] + distanceMatrix[k][j]) {

					frequencyMatrix[i][j] = frequencyMatrix[i][k] + frequencyMatrix[k][j];
					frequencyMatrix[j][i] = frequencyMatrix[i][j];
				}

}

int main()
{

	initializeFiles("rf");
	inputFile >> nodeNr;

	readAdjacencyMatrix();

	RoyFloyd();

	printMatrix(distanceMatrix);
	printMatrix(frequencyMatrix);
}

void initializeFiles(const string filename) {
	inputFile = ifstream(filename + ".in");
	outputFile = ofstream(filename + ".out");
}