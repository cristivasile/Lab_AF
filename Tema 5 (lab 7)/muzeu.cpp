#include <fstream>
#include <vector>
#include <queue>

using namespace std;

ifstream inputFile;
ofstream outputFile;

void initializeFiles(const string filename) {
	inputFile = ifstream(filename + ".in");
	outputFile = ofstream(filename + ".out");
}

struct Position {
	int i, j;
	Position() { i = 0; j = 0; }
	Position(int x, int y) : i(x), j(y){};
};

bool isInMatrix(const int size, const int i, const int j) {
	if (i < 0 || j < 0 || i >= size || j >= size)
		return false;

	return true;
}

void Move(const int i, const int j, const int currentDistance, const int size, vector<vector<int>>& distanceMatrix, queue<Position>& toCheck) {
	if (isInMatrix(size, i, j) && distanceMatrix[i][j] == -1) {
		distanceMatrix[i][j] = currentDistance + 1;
		toCheck.push(Position(i, j));
	}
}

void Lee(const int size, vector<vector<int>> &distanceMatrix, queue<Position>& toCheck) {

	Position currentPosition;
	int i, j, currentDistance;

	while (toCheck.size() != 0) {
		currentPosition = toCheck.front();
		toCheck.pop();
		i = currentPosition.i;
		j = currentPosition.j;
		currentDistance = distanceMatrix[i][j];

		//go left
		Move(i, j - 1, currentDistance, size, distanceMatrix, toCheck);
		//go up
		Move(i - 1, j, currentDistance, size, distanceMatrix, toCheck);
		//go right
		Move(i, j + 1, currentDistance, size, distanceMatrix, toCheck);
		//go down
		Move(i + 1, j, currentDistance, size, distanceMatrix, toCheck);
	}
}

void Read(const int size, vector<vector<int>>& distanceMatrix, queue<Position>& toCheck) {

	char input;

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++) {
			inputFile >> input;
			if (input == '#')
				distanceMatrix[i][j] = -2;
			else if (input == 'P') {
				distanceMatrix[i][j] = 0;
				toCheck.push(Position(i, j));
			}
		}
}

void Write(const int size, vector<vector<int>>& distanceMatrix) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++)
			outputFile << distanceMatrix[i][j] << " ";
		outputFile << "\n";
	}
}

int main()
{
	int size;
	queue<Position> toCheck;

	initializeFiles("muzeu");
	inputFile >> size;
	vector<vector<int>> distanceMatrix(size, vector<int>(size, -1));

	Read(size, distanceMatrix, toCheck);
	Lee(size, distanceMatrix, toCheck);
	Write(size, distanceMatrix);

	return 0;
}

