//325763498
//michalshasha8@gmail.com

#include <iostream>
#include <vector>
#include <iostream>
#include <stdexcept>
#pragma once
using namespace std;
namespace ariel{

class Graph {

public:
    string printGraph(); 
    void loadGraph(vector<vector<int>>&);
    void setMatrix(vector<vector<int>>&);
    size_t getSize();
    const vector<vector<int>> getMatrix();
    bool getEdge();
    bool getDirected();
    bool isSymmetricMatrix(std::vector<std::vector<int>>&);
    bool isSquareMatrix(vector<vector<int>>& matrix);
    vector<int> getNeighbors(int u);
    int getWeight(int, int);

    friend ostream& operator<<(ostream& os, Graph& graph);
    Graph operator+(Graph& other);
    Graph& operator+=(int scalar);
    Graph operator+();
    Graph operator-(Graph& other);
    Graph& operator-=(int scalar);
    Graph operator-();

    Graph& operator*=(int scalar);
    Graph& operator/=(int scalar);
    Graph operator*(Graph& other);
    
    bool operator== (Graph& other); 
    bool operator!=(Graph& other);
    bool operator>(Graph& other); 
    bool operator<=(Graph& other);
    bool operator<(Graph& other);
    bool operator>=(Graph& other);

    Graph& operator++();
    Graph operator++(int);  
    Graph& operator--();
    Graph operator--(int);


    bool isGraphContained( Graph& smallerGraph, Graph& largerGraph);

    

    friend ostream& operator<<(ostream& os, Graph& graph);

private:
    vector<vector<int>> adjMatrix; 
    size_t numVertices; 
    int numEdges;
    bool edge;
    bool weight;
    bool directed;

};
}