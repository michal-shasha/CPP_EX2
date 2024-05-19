//325763498
//michalshasha8@gmail.com

#include "Algorithms.hpp"
#include <limits>
#include <climits>
#include <sstream>
#include <algorithm>
#include <stack>
#include <set>

using namespace std;
namespace ariel{

Graph Algorithms::transposeGraph(Graph& graph) {
    const vector<vector<int>>& adjMat=graph.getMatrix();
    size_t size = graph.getSize();
    vector<vector<int>> transposed( size, vector<int>(size, 0));

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            transposed[j][i] = adjMat[i][j];
        }
    }
    graph.setMatrix(transposed);
    return graph;
}  

void Algorithms:: dfs( Graph& graph, int vertex, vector<bool>& visited) {
    // סימון הצומת הנוכחי כביקר
    visited[size_t(vertex)] = true;
    const vector<vector<int>>& adjMat=graph.getMatrix();

    // מעבר על כל הצמתים השכנים של הצומת הנוכחי
    for (size_t neighbor = 0; neighbor < adjMat.size(); ++neighbor) {
        // בדיקה האם יש קשת מהצומת הנוכחי לשכן זה
        if (adjMat[size_t(vertex)][size_t(neighbor)] && !visited[size_t(neighbor)]) {
            // נקרא לעצמנו רקורסיבית עבור השכן הנוכחי
            dfs(graph, neighbor, visited);
        }
    }
}

int Algorithms::isConnected(Graph &g) 
{
if(g.getSize()==0)   
    return 1;

vector<bool> visited(g.getSize(), false);
dfs(g, 0, visited);

// בדיקה האם כל הצמתים נגישים
bool allVisited = true;
for (bool isVisited : visited) {
    if (!isVisited) {
        allVisited = false;
        break;
    }
}

if (!allVisited) {
    return false; // לא כל הצמתים נגישים, הגרף אינו קשיר
}

// אם הגרף הוא קשיר, נבדוק את הגרף ההפוך
Graph gg= transposeGraph(g);
visited.assign(g.getSize(), false); // איפוס רשימת הביקור
dfs(gg, 0, visited);

// בדיקה האם כל הצמתים נגישים בגרף ההפוך
allVisited = true;
for (bool isVisited : visited) {
    if (!isVisited) {
        allVisited = false;
        break;
    }
}

if (!allVisited) {
    return false; // לא כל הצמתים נגישים בגרף ההפוך, הגרף אינו קשיר
}

// אם הגרף הוא קשיר וגם הגרף ההפוך הוא קשיר, הגרף כולו קשיר
return true;
}


   
string Algorithms::shortestPath(Graph& g, int start, int end) {
    int V = g.getSize();
    if(V==0)
       return "-1";
    vector<int> dist(size_t(V), numeric_limits<int>::max());
    vector<int> parent(size_t(V), -1);
    dist[size_t(start)] = 0;

    // Relax all edges V-1 times
    for (size_t  i = 0; i < V - 1; ++i) {
        for (size_t  u = 0; u < V; ++u) {
            for (size_t  v = 0; v < V; ++v) {
                int weight = g.getWeight(u, v);
                if (weight != 0 && dist[u] + weight < dist[v]) {
                    dist[v] = dist[size_t(u)] + weight;
                    parent[v] = u;
                    }
                }
            }
        }


    // Check for negative cycle
    for (size_t u = 0; u < V; ++u) {
        for (size_t v = 0; v < V; ++v) {
            int weight = g.getWeight(u, v);
            if (weight != 0 && dist[u] + weight < dist[v]) {
                // Negative cycle found
                return "-1";
            }
        }
    }


    // check for no path
    if( parent[(size_t)end] == -1){
        return "-1";
    }
    // Reconstruct path
    string path = to_string(end);
    int p = parent[(size_t)end];

    while(p!=-1){
        path = to_string(p) + "->" + path;
        p = parent[(size_t)p];
    }
    return path;
return "-1";
}


bool Algorithms::isContainsCycle(Graph& g) {
    if(g.getSize()==0)
      return false;
    vector<Color> colors(g.getSize(), WHITE);
    vector<int> parents(g.getSize(), -1);
    vector<int> path;

    for (size_t i = 0; i < g.getSize(); i++) {
        if (colors[i] == WHITE) {
            if (isContainsCycleUtil(g, i, &colors, &parents, &path)) {
                return true;
            }
        }
    }
    return false;
}

bool Algorithms::isContainsCycleUtil(Graph& g, size_t src, vector<Color>* colors, vector<int>* parents, vector<int>* path) 
{
    (*colors)[src] = GRAY;
    path->push_back(src);

    for (size_t v = 0; v < g.getSize(); v++) 
    {
        if (g.getMatrix()[src][v] != 0) 
        {
            if ((*colors)[v] == WHITE) 
            {
                (*parents)[v] = (int)src;
                if (isContainsCycleUtil(g, v, colors, parents, path)) 
                {
                    return true;
                }
            }
            else if ((*colors)[v] == GRAY) 
            {
                if (!g.getDirected() && (*parents)[src] == (int)v) 
                {
                    continue;
                }
                printCycle(*path, *parents, v);
                return true; // חזור עם ערך בוליאני true במקרה של התאמה לנקודת ההתחלה של המעגל
            }
        }
    }

    (*colors)[src] = BLACK;
    path->pop_back();
    return false;
}

void Algorithms::printCycle(vector<int>& path, vector<int>& parents, int start) {
    bool cycleFound = false;
    size_t idx = 0;
    for (idx = 0; idx < path.size(); ++idx) {
        if (path[idx] == start) {
            cycleFound = true;
            break;
        }
    }

    if (cycleFound) {
        cout << "The cycle is: ";
        for (size_t i = idx; i < path.size(); ++i) {
            cout << path[i];
            if (i != path.size() - 1) {
                cout << "->";
            }
        }
        cout << "->" << start << endl;
    }
}



string Algorithms::isBipartite( Graph &graph) 
{
    int n = graph.getSize();
    if (n == 0)  // Check if the graph has no vertices
      return "The graph is bipartite: A={}, B={}";
    vector<int> colors(size_t(n), -1); // Initialize all colors to -1 (unvisited)

    for (int i = 0; i < n; ++i) {
        if (colors[size_t(i)] == -1) {
            if ((!d(i, graph, colors, 0))) {
                return "0"; // Graph is not bipartite
            }
        }
    }

    // Construct sets for both partitions
    std::set<int> A, B;
    for (int i = 0; i < n; ++i) {
        if (colors[size_t(i)] == 0) {
            A.insert(i);
        } else {
            B.insert(i);
        }
    }

    // Construct the output string
    string result = "The graph is bipartite: A={";
    for (int vertex : A) {
        result += to_string(vertex) + ", ";
    }
    result.pop_back(); // Remove trailing comma
    result.pop_back(); // Remove space
    result += "}, B={";
    for (int vertex : B) {
        result += to_string(vertex) + ", ";
    }
    result.pop_back(); // Remove trailing comma
    result.pop_back(); // Remove space
    result += "}";
    return result;
}


bool Algorithms::d(int v, Graph &graph, vector<int> &colors, int color) {
    colors[size_t(v)] = color; // Assign color to the current vertex

    for (int u = 0; u < graph.getSize(); ++u) {
        if (graph.getWeight(v, u) != 0) { // If there's an edge between v and u
            if (colors[size_t(u)] == -1) { // If u is not colored yet
                if (!d(u, graph, colors, 1 - color)) { // Recursively call d with the opposite color
                    return false;
                }
            } else if (colors[size_t(u)] == color) { // If u is already colored and has the same color as v
                return false; // Graph is not bipartite
            }
        }
    }
    return true; // Graph is bipartite
}




string Algorithms::negativeCycle(Graph& g) 
{
    size_t V = g.getSize();
    vector<int> dist(V, 0);
    vector<int> parent(V, -1);

    // Step 1: Relax all edges V-1 times
    for (size_t i = 0; i < V - 1; ++i) 
    {
        for (size_t u = 0; u < V; ++u) 
        {
            for (size_t v = 0; v < V; ++v) 
            {
                int weight = g.getWeight(u, v);
                if (weight != 0 && dist[u] + weight < dist[v]) 
                {
                    dist[v] = dist[u] + weight;
                    parent[v] = u;
                }
            }
        }
    }

    // Step 2: Check for negative cycle
    vector<int> cycle;
    for (size_t u = 0; u < V; ++u) 
    {
        for (size_t v = 0; v < V; ++v) 
        {
            int weight = g.getWeight(u, v);
            if (weight != 0 && dist[u] + weight < dist[v]) 
            {
                // Negative cycle found, trace back to find the cycle
                int current = v;
                do {
                    cycle.push_back(current);
                    current = parent[size_t(current)];
                } while (current != v);

                // Add the starting vertex of the cycle
                cycle.push_back(v);

                // Reverse the cycle to get the correct order
                reverse(cycle.begin(), cycle.end());

                // Construct the cycle string
                string cycleStr = "Negative Cycle: ";
                for (size_t i = 0; i < cycle.size(); ++i) 
                {
                    cycleStr += to_string(cycle[i]);
                    if (i < cycle.size() - 1) cycleStr += " -> ";
                }

                // Return the cycle string
                return cycleStr;
            }
        }
    }

    // No negative cycle found
    return "No negative cycle found";
}
}