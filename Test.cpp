//325763498
//michalshasha8@gmail.com

#include "doctest.h"
#include "Algorithms.hpp"
#include "Graph.hpp"

using namespace std;

TEST_CASE("Test graph multiplication operator (*)") {
    ariel::Graph g1, g2;
    vector<vector<int>> graph1 = {
        {0, 2},
        {3, 0}
    };
    vector<vector<int>> graph2 = {
        {0, 6},
        {7, 0}
    };
    g1.loadGraph(graph1);
    g2.loadGraph(graph2);

    // Perform multiplication operation
    ariel::Graph result = g1 * g2;

    // Check the result against the expected result
    CHECK(result.printGraph() == "[0, 0]\n[0, 0]");
    ariel::Graph result2 = (result*g1)+g2;
    CHECK(result2.printGraph() == "[0, 6]\n[7, 0]");

}


TEST_CASE("Test compound subtraction and division operators combined") {
    // Define specific graphs
    ariel::Graph g1, g2;
    vector<vector<int>> specificGraph1 = {
        {0, 20},
        {30, 0}
    };
    vector<vector<int>> specificGraph2 = {
        {0, 2},
        {2, 0}
    };
    g1.loadGraph(specificGraph1);
    g2.loadGraph(specificGraph2);

    // Perform compound subtraction and division operations combined
    g1 -= 5;
    g2 += 2;

    // Perform subtraction operation
    ariel::Graph result = g1 - g2;

    // Check the result against the expected result
    CHECK(result.printGraph() == "[0, 11]\n[21, 0]");
}

TEST_CASE("Test graph comparison operators") {
   ariel:: Graph g1, g2, g3;
    vector<vector<int>> Graph1 = {
        {0, 2},
        {3, 0}
    };
    vector<vector<int>> Graph2 = {
        {0, 2, 3},
        {4, 0, 6},
        {7, 8, 0}
    };
    vector<vector<int>> Graph3 = {
        {0, 2},
        {4, 0}
    };
    g1.loadGraph(Graph1);
    g2.loadGraph(Graph2);
    g3.loadGraph(Graph3);

    // Check equality
    CHECK(g1 == g1);
    CHECK_FALSE(g1 == g2);

    // Check inequality
    CHECK(g1 != g2);
    CHECK_FALSE(g1 != g1);

    // Check greater than
    CHECK_FALSE(g1 > g1);
    CHECK_FALSE(g1 > g2);

    // Check less than
    CHECK(g1 < g1);
    CHECK(g1 < g2);

    // Check greater than or equal to
    CHECK(g1 >= g1);
    CHECK_FALSE(g1 >= g2);

    // Check less than or equal to
    CHECK(g1 <= g1);
    CHECK(g1 <= g2);
}

TEST_CASE("Test compound addition and subtraction operators combined") {
    ariel::Graph g1, g2;
    vector<vector<int>> Graph1 = {
        {0, 20},
        {30, 0}
    };
    vector<vector<int>> Graph2 = {
        {0, 5},
        {5, 0}
    };
    g1.loadGraph(Graph1);
    g2.loadGraph(Graph2);

    // Perform compound addition and subtraction operations combined
    g1 += 7;
    g2 -= 3;

    // Perform addition and subtraction operations
    ariel::Graph result = g1 + g2;

    // Check the result against the expected result
    CHECK(result.printGraph() == "[0, 29]\n[39, 0]");
}

TEST_CASE("Test graph division operator with zero divisor") {
    // Define a specific graph
    ariel::Graph g;
    vector<vector<int>> Graph = {
        {0, 20},
        {30, 0}
    };
    g.loadGraph(Graph);

    // Attempt division by zero
    CHECK_THROWS_AS(g /= 0, std::invalid_argument);
}

TEST_CASE("Test graph addition")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, 1, 1},
        {1, 0, 2},
        {1, 2, 0}};
    g2.loadGraph(weightedGraph);
    ariel::Graph g3 = g1 + g2;
    vector<vector<int>> expectedGraph = {
        {0, 2, 1},
        {2, 0, 3},
        {1, 3, 0}};
    CHECK(g3.printGraph() == "[0, 2, 1]\n[2, 0, 3]\n[1, 3, 0]");
}

TEST_CASE("Test graph multiplication")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, 1, 1},
        {1, 0, 2},
        {1, 2, 0}};
    g2.loadGraph(weightedGraph);
    ariel::Graph g4 = g1 * g2;

    vector<vector<int>> expectedGraph = {
        {0, 0, 2},
        {1, 0, 1},
        {1, 0, 0}};
    CHECK(g4.printGraph() == "[0, 0, 2]\n[1, 0, 1]\n[1, 0, 0]");
}

TEST_CASE("Invalid operations")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, 1, 1, 1},
        {1, 0, 2, 1},
        {1, 2, 0, 1}};
    //g2.loadGraph(weightedGraph);
    ariel::Graph g5;
    vector<vector<int>> graph2 = {
        {0, 1, 0, 0, 1},
        {1, 0, 1, 0, 0},
        {0, 1, 0, 1, 0},
        {0, 0, 1, 0, 1},
        {1, 0, 0, 1, 0}};
    g5.loadGraph(graph2);
     CHECK_THROWS(g5 * g1);
    // CHECK_THROWS(g1 * g2);

    // Addition of two graphs with different dimensions
    ariel::Graph g6;
    vector<vector<int>> graph3 = {
        {0, 1, 0, 0, 1},
        {1, 0, 1, 0, 0},
        {0, 1, 0, 1, 0},
        {0, 0, 1, 0, 1},
        {1, 0, 0, 1, 0}};
    g6.loadGraph(graph3);
    CHECK_THROWS(g1 + g6);
}
TEST_CASE("Test unary minus operator")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);

    ariel:: Graph g2 = -g1;

    vector<vector<int>> expectedGraph = {
        {0, -1, 0},
        {-1, 0, -1},
        {0, -1, 0}};

    CHECK(g2.printGraph() == "[0, -1, 0]\n[-1, 0, -1]\n[0, -1, 0]");
}
