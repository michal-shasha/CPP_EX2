//325763498
//michalshasha8@gmail.com

#include "Graph.hpp"
#include <sstream>
#include "numeric"
#include "algorithm"
using namespace std;

namespace ariel
{

    string Graph::printGraph()
    {
        stringstream ss;
        ss << *this;
        return ss.str();
    }

    bool Graph::isSquareMatrix(vector<vector<int>> &matrix)
    {
        size_t n = matrix.size();
        for (size_t i = 0; i < n; ++i)
        {
            if (matrix[i].size() != n)
            {
                return false;
            }
            if (matrix[i][i] != 0)
            {
                return false;
            }
        }
        return true;
    }

    bool Graph::isSymmetricMatrix(vector<vector<int>> &matrix)
    {

        for (size_t i = 0; i < matrix.size(); ++i)
        {
            for (size_t j = i + 1; j < matrix.size(); ++j)
            {
                if (matrix[i][j] != matrix[j][i])
                {
                    return false;
                }
            }
        }
        return true;
    }

    void Graph::loadGraph(vector<vector<int>> &g)
    {
        this->numEdges = 0;
        this->numVertices = 0;
        this->edge = false;
        this->weight = false;
        this->adjMatrix = g;

        // // Check if the matrix is square
        size_t n = g.size();
        if (n != g[0].size())
        {
            throw invalid_argument("Invalid graph: The graph is not a square matrix.");
        }

        // Check if the matrix contains edges or weights
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                if (this->adjMatrix[i][j] < 0)
                {
                    this->edge = true; // Negative value indicates an edge
                }
                if (this->adjMatrix[i][j] != 0 && this->adjMatrix[i][j] != 1)
                {
                    this->weight = true; // Value different from 0 and 1 indicates a weight
                }
                if (this->adjMatrix[i][j] != 0)
                {
                    this->numEdges++; // Counting the number of edges
                }
            }
        }

        this->numVertices = n;
        
        // Determine if the graph is directed or undirected
        this->directed = !isSymmetricMatrix(g);
        if(!this->directed){
            this->numEdges /= 2;
        }
    }

    void Graph::setMatrix(vector<vector<int>> &mat)
    {
        this->adjMatrix = mat;
    }

    size_t Graph::getSize()
    {
        return this->numVertices;
    }

    const vector<vector<int>> Graph::getMatrix()
    {
        return this->adjMatrix;
    }

    bool Graph::getEdge()
    {
        return this->edge;
    }

    bool Graph::getDirected()
    {
        return this->directed;
    }

    vector<int> Graph::getNeighbors(int u)
    {
        return this->adjMatrix[size_t(u)];
    }

    int Graph::getWeight(int u, int v)
    {
        return this->adjMatrix[static_cast<std::size_t>(u)][static_cast<std::size_t>(v)];
    }

    Graph Graph::operator+(Graph &other)
    {
       
        if (this->adjMatrix.size() != other.adjMatrix.size() || this->adjMatrix[0].size() != other.adjMatrix[0].size())
        {
        
            throw std::invalid_argument("Cannot add graphs of different sizes");
        }

        Graph result;

        vector<vector<int>> newMatrix;
        for (size_t i = 0; i < this->adjMatrix.size(); ++i)
        {
            vector<int> newRow;
            for (size_t j = 0; j < this->adjMatrix[0].size(); ++j)
            {
                newRow.push_back(this->adjMatrix[i][j] + other.adjMatrix[i][j]);
            }
            newMatrix.push_back(newRow);
        }
        // Set the main diagonal to zero
        for (size_t i = 0; i < newMatrix.size(); ++i)
        {
           newMatrix[i][i] = 0;
        }

        result.setMatrix(newMatrix);
        result.numVertices = this->numVertices;
        result.numEdges = this->numEdges + other.numEdges;

        return result;
    }

    Graph Graph::operator-(Graph &other)
    {
        // Check if the matrices are of the same size
        if (this->adjMatrix.size() != other.adjMatrix.size() || this->adjMatrix[0].size() != other.adjMatrix[0].size())
        {
            throw std::invalid_argument("Cannot add graphs of different sizes");
        }

        // Create a new graph object
        Graph result;

        // subtraction the matrices element-wise
        vector<vector<int>> newMatrix;
        for (size_t i = 0; i < this->adjMatrix.size(); ++i)
        {
            vector<int> newRow;
            for (size_t j = 0; j < this->adjMatrix[0].size(); ++j)
            {
                newRow.push_back(this->adjMatrix[i][j] - other.adjMatrix[i][j]);
            }
            newMatrix.push_back(newRow);
        }
         // Set the main diagonal to zero
        for (size_t i = 0; i < newMatrix.size(); ++i)
        {
           newMatrix[i][i] = 0;
        }

        // Set the new matrix for the result graph
        result.setMatrix(newMatrix);

        return result;
    }

    Graph Graph::operator+()
    {
        Graph result(*this);

        for (size_t i = 0; i < result.adjMatrix.size(); ++i)
        {
            for (size_t j = 0; j < result.adjMatrix[i].size(); ++j)
            {
                result.adjMatrix[i][j] *= 1;
            }
        }
        return result;
    }

    Graph Graph::operator-()
    {
        Graph result;

        vector<vector<int>> newMatrix = this->adjMatrix;
        for (size_t i = 0; i < newMatrix.size(); ++i)
        {
            for (size_t j = 0; j < newMatrix[i].size(); ++j)
            {
                newMatrix[i][j] = -newMatrix[i][j];
            }
        }

        result.setMatrix(newMatrix);
        result.numVertices = this->numVertices;
        result.numEdges = this->numEdges;
        result.edge = this->edge;
        result.weight = this->weight;
        result.directed = this->directed;

        return result;
    }


Graph& Graph::operator+=(int scalar)
{
    // Iterate over each row in the adjacency matrix
    for (auto& row : this->adjMatrix)
    {
        // Iterate over each element in the current row
        for (auto& elem : row)
        {
            // Add the scalar value to the current element
            elem += scalar;
        }
    }
    // Set the main diagonal to zero
    for (size_t i = 0; i < this->getSize(); ++i)
    {
        this->adjMatrix[i][i] = 0;
    }
    return *this;
}

Graph& Graph::operator-=(int scalar)
{
    // Iterate over each row in the adjacency matrix
    for (auto& row : this->adjMatrix)
    {
        // Iterate over each element in the current row
        for (auto& elem : row)
        {
            // Subtract the scalar value from the current element
            elem -= scalar;
        }
    }
    // Set the main diagonal to zero
    for (size_t i = 0; i < this->getSize(); ++i)
    {
        this->adjMatrix[i][i] = 0;
    }
    return *this;
}
    Graph &Graph::operator*=(int scalar)
    {
        for (size_t i = 0; i < adjMatrix.size(); ++i)
        {
            for (size_t j = 0; j < adjMatrix[i].size(); ++j)
            {
                adjMatrix[i][j] *= scalar;
            }
        }
        return *this;
    }

    Graph &Graph::operator/=(int scalar)
    {
        if (scalar == 0)
        {
            throw std::invalid_argument("Cannot divide by zero");
        }
        for (size_t i = 0; i < adjMatrix.size(); ++i)
        {
            for (size_t j = 0; j < adjMatrix[i].size(); ++j)
            {
                adjMatrix[i][j] /= scalar;
            }
        }

        return *this;
    }

    Graph Graph::operator*(Graph &other)
    {
        const vector<vector<int>> &matrix1 = this->getMatrix();
        const vector<vector<int>> &matrix2 = other.getMatrix();

        size_t size1 = matrix1.size();
        size_t size2 = matrix2.size();

        // Check if the number of columns in the first matrix is equal to the number of rows in the second matrix
        if (matrix1[0].size() != size2)
        {
            throw invalid_argument("The number of columns in the first matrix must be equal to the number of rows in the second matrix.");
        }

        // Create a new matrix to store the result of multiplication
        vector<vector<int>> result(size1, vector<int>(size2, 0));

        // Perform matrix multiplication
        for (size_t i = 0; i < size1; ++i)
        {
            for (size_t j = 0; j < size2; ++j)
            {
                for (size_t k = 0; k < matrix1[0].size(); ++k)
                {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
        // Set the main diagonal to zero
        for (size_t i = 0; i < size1; ++i)
        {
           result[i][i] = 0;
        }

        // Create a new graph and load the result matrix into it
        Graph multipliedGraph;
        multipliedGraph.loadGraph(result);

        return multipliedGraph;
    }



    bool Graph::operator==(Graph &other)
    {
      return this->adjMatrix == other.adjMatrix || (!(*this < other) && !(other < *this));
    }

    bool Graph::operator!=(Graph &other)
    {
        return !(*this == other);
    }

    bool Graph::isGraphContained( Graph &largerGraph)
    {
        
        if (this->getSize() > largerGraph.getSize())
        {
            return false;
        }

        // Create a vector for vertex mapping
        std::vector<size_t> mapping(largerGraph.getSize());
        std::iota(mapping.begin(), mapping.end(), 0);

        // Try all permutations of the vertex mapping
        do
        {
            bool is_subgraph = true;
            for (size_t i = 0; i < this->getSize() && is_subgraph; ++i)
            {
                for (size_t j = 0; j < this->getSize() && is_subgraph; ++j)
                {
                    if (this->adjMatrix[i][j] != 0 && largerGraph.adjMatrix[mapping[i]][mapping[j]] == 0)
                    {
                        is_subgraph = false;
                    }
                }
            }
            if (is_subgraph)
            {
                return true;
            }
        } while (std::next_permutation(mapping.begin(), mapping.end()));

        return false;
    }
    

    bool Graph::operator>(Graph &other)
    {
      return other<*this;
    }

    bool Graph::operator<=(Graph &other)
    {
        // Check if the graph is less than or equal to another graph based on containment, edge count, and order of the adjacency matrix
        return (*this < other) || (*this == other);
    }

    bool Graph::operator<(Graph &other)
    {
        if(this->adjMatrix==other.adjMatrix)
        {
            return false;
        }
         

        // Check if the graph is less than another graph based on containment, edge count, and order of the adjacency matrix
        if (isGraphContained(other)&& !(other.isGraphContained(*this)))
        {
            return true;
        }

        if (other.isGraphContained(*this)&&!(isGraphContained(other)))
        {
            return false; // If the other graph is contained in this, this graph is greater
        }

        // If neither graph is contained in the other, compare edge count
        if (this->numEdges != other.numEdges)
        {
            return this->numEdges < other.numEdges;
        }

        // If the number of edges is equal, compare the order of the adjacency matrix
        return this->adjMatrix.size() < other.adjMatrix.size();
    }

    bool Graph::operator>=(Graph &other)
    {
        // Check if the graph is greater than or equal to another graph based on containment, edge count, and order of the adjacency matrix
        return (*this == other)||(*this>other);
    }

    Graph &Graph::operator++()
    {
        for (size_t i = 0; i < adjMatrix.size(); ++i)
        {
            for (size_t j = 0; j < adjMatrix[i].size(); ++j)
            {
                adjMatrix[i][j]++;
            }
        }
        for (size_t i = 0; i < this->getSize(); ++i)
        {
        this->adjMatrix[i][i] = 0;
        }
        return *this;
    }

    Graph &Graph::operator--()
    {
        for (size_t i = 0; i < adjMatrix.size(); ++i)
        {
            for (size_t j = 0; j < adjMatrix[i].size(); ++j)
            {
                adjMatrix[i][j]--;
            }
        }
        for (size_t i = 0; i < this->getSize(); ++i)
        {
          this->adjMatrix[i][i] = 0;
        }

        return *this;
    }

    Graph Graph::operator++(int)
    {
        Graph temp(*this);
        for (size_t i = 0; i < adjMatrix.size(); ++i)
        {
            for (size_t j = 0; j < adjMatrix[i].size(); ++j)
            {
                adjMatrix[i][j]++;
            }
        }
        for (size_t i = 0; i < this->getSize(); ++i)
        {
          this->adjMatrix[i][i] = 0;
        }
        return temp;
    }

    Graph Graph::operator--(int)
    {
        Graph temp(*this);
        for (size_t i = 0; i < adjMatrix.size(); ++i)
        {
            for (size_t j = 0; j < adjMatrix[i].size(); ++j)
            {
                adjMatrix[i][j]--;
            }
        }
        for (size_t i = 0; i < this->getSize(); ++i)
        {
           this->adjMatrix[i][i] = 0;
        }
        return temp;
    }
    std::ostream &operator<<(std::ostream &out, Graph &graph)
    {
        vector<vector<int>> matrix = graph.getMatrix();
        for (size_t i = 0; i < matrix.size(); ++i)
        {
            out << "[";
            for (size_t j = 0; j < matrix[i].size(); ++j)
            {
                out << matrix[i][j];
                if (j != matrix[i].size() - 1)
                {
                    out << ", ";
                }
            }
            out << "]";
            if (i != matrix.size() - 1)
            {
                out << endl;
            }
        }
        return out;
    }
};
