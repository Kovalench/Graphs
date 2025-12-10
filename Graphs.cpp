#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <string>

using namespace std;

const int INF = numeric_limits<int>::max();

const double PI = 3.141592653589793238;

class Graph_without_weight {
private:
    int numVertices; 
    vector<vector<int>> adjacencyList; 
    
public:
    Graph_without_weight(int vertices) : numVertices(vertices) {
        adjacencyList.resize(vertices, vector<int>(vertices, 0));
    }

    void addEdge(int u, int v) {
        adjacencyList[u][v] = 1;
        adjacencyList[v][u] = 1;
    }

    void dfs(int v, vector<bool>& visited);
    void dfsTraversal(int start);
    void bfs(int start);
};

void Graph_without_weight::dfs(int v, vector<bool>& visited) {
    cout << v << " ";
    visited[v] = true;

    for (int i = 0; i < numVertices; i++) {
        if (adjacencyList[v][i] == 1 && !visited[i]) {
            dfs(i, visited);
        }
    }
}

void Graph_without_weight::dfsTraversal(int start) {
    vector<bool> visited(numVertices, false);
    cout << "Обхiд в глибину, починаючи з вершини " << start << ": ";
    dfs(start, visited);
    cout << endl;
}

void Graph_without_weight::bfs(int start) {
    vector<bool> visited(numVertices, false);
    queue<int> q;

    q.push(start);
    visited[start] = true;

    cout << "Обхiд в ширину, починаючи з вершини " << start << ": ";

    while (!q.empty()) {
        int v = q.front();
        q.pop();
        cout << v << " ";

        for (int i = 0; i < numVertices; i++) {
            if (adjacencyList[v][i] == 1 && !visited[i]) {
                q.push(i);
                visited[i] = true;
            }
        }
    }
    cout << endl;
}

class Graph_with_weight {
private:
    int numVertices;
    vector<vector<pair<int, int>>> adjacencyList;

public:
    Graph_with_weight(int vertices) : numVertices(vertices), adjacencyList(vertices) {}

    void addEdge(int u, int v, int weight) {
        adjacencyList[u].emplace_back(v, weight);
        adjacencyList[v].emplace_back(u, weight);
    }

    vector<vector<int>> Kruskal();
    vector<vector<int>> Prim();
    vector<vector<int>> Dijkstra(int start);
    void generateDOTFile(const string& filename, const vector<vector<int>>& mst);
};

struct Edge {
    int weight, u, v;
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

struct DSU {
    vector<int> parent, rank;

    DSU(int n) {
        parent.resize(n);
        rank.resize(n, 1);
        for (int i = 0; i < n; i++) parent[i] = i;
    }

    int find(int x) {
        return (parent[x] == x) ? x : (parent[x] = find(parent[x]));
    }

    void unite(int x, int y) {
        int rx = find(x), ry = find(y);
        if (rx != ry) {
            if (rank[rx] > rank[ry]) swap(rx, ry);
            parent[rx] = ry;
            if (rank[rx] == rank[ry]) rank[ry]++;
        }
    }
};

vector<vector<int>> Graph_with_weight::Kruskal() {
    vector<Edge> edges;
    for (int i = 0; i < numVertices; i++) {
        for (const auto& neighbor : adjacencyList[i]) {
            int j = neighbor.first;
            int weight = neighbor.second;
            if (i < j) {
                edges.push_back({ weight, i, j });
            }
        }
    }

    sort(edges.begin(), edges.end());
    DSU dsu(numVertices);
    vector<vector<int>> tree(numVertices, vector<int>(numVertices, INF));
    int addedEdges = 0;

    for (const auto& e : edges) {
        if (dsu.find(e.u) != dsu.find(e.v)) {
            dsu.unite(e.u, e.v);
            tree[e.u][e.v] = e.weight;
            tree[e.v][e.u] = e.weight;
            if (++addedEdges == numVertices - 1) break;
        }
    }

    return tree;
}

vector<vector<int>> Graph_with_weight::Prim() {
    vector<vector<int>> tree(numVertices, vector<int>(numVertices, INF));
    vector<bool> inMST(numVertices, false);
    priority_queue<pair<int, pair<int, int>>, vector<pair<int, pair<int, int>>>, greater<>> pq;

    inMST[0] = true;
    for (const auto& neighbor : adjacencyList[0]) {
        pq.push({ neighbor.second, {0, neighbor.first} });
    }

    int edgesAdded = 0;
    while (!pq.empty() && edgesAdded < numVertices - 1) {
        auto current = pq.top();
        pq.pop();
        int weight = current.first;
        int u = current.second.first;
        int v = current.second.second;

        if (inMST[v]) continue;

        tree[u][v] = weight;
        tree[v][u] = weight;
        inMST[v] = true;
        edgesAdded++;

        for (const auto& neighbor : adjacencyList[v]) {
            if (!inMST[neighbor.first]) {
                pq.push({ neighbor.second, {v, neighbor.first} });
            }
        }
    }

    return tree;
}

vector<vector<int>> Graph_with_weight::Dijkstra(int start) {
    vector<int> distance(numVertices, INF);
    vector<int> predecessor(numVertices, -1);
    distance[start] = 0;

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
    pq.emplace(0, start);

    while (!pq.empty()) {
        int dist_u = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        if (dist_u > distance[u]) continue;

        for (const auto& [v, weight] : adjacencyList[u]) {
            if (distance[u] + weight < distance[v]) {
                distance[v] = distance[u] + weight;
                predecessor[v] = u;
                pq.emplace(distance[v], v);
            }
        }
    }

    cout << "Найкоротшi вiдстанi вiд вершини " << start << ":\n";
    for (int i = 0; i < numVertices; ++i) {
        cout << "До вершини " << i << ": " << (distance[i] == INF ? "INF" : to_string(distance[i])) << endl;
    }

    vector<vector<int>> result(numVertices, vector<int>(numVertices, INF));

    for (int v = 0; v < numVertices; ++v) {
        if (predecessor[v] != -1) {
            int u = predecessor[v];
            result[u][v] = distance[v] - distance[u];
            result[v][u] = distance[v] - distance[u];
        }
    }

    return result;
}

class diGraph_with_weight {
private:
    int numVertices;
    vector<vector<pair<int, int>>> adjacencyList;

public:
    diGraph_with_weight(int vertices) : numVertices(vertices), adjacencyList(vertices) {}

    void addEdge(int u, int v, int weight) {
        adjacencyList[u].emplace_back(v, weight);
    }

    vector<vector<int>> BellmanFord(int start); 
    vector<vector<int>> FloydWarshall();
    vector<vector<int>> Johnson();
    int FordFulkerson(int source, int sink);
    int EdmondsKarp(int source, int sink);
    void generateDOTFileDierected(const string& filename, const vector<vector<int>>& mst, int painted);
};

int dfs(vector<vector<int>>& residual, vector<bool>& visited, int u, int sink, int flow, vector<int>& path) {
    if (u == sink) {
        path.push_back(u);
        return flow;
    }
    visited[u] = true;

    for (int v = 0; v < residual.size(); ++v) {
        if (!visited[v] && residual[u][v] > 0) {
            int newFlow = dfs(residual, visited, v, sink, min(flow, residual[u][v]), path);
            if (newFlow > 0) {
                residual[u][v] -= newFlow;
                residual[v][u] += newFlow;
                path.push_back(u);
                return newFlow;
            }
        }
    }
    return 0;
}

bool bfs(const vector<vector<int>>& residual, int source, int sink, vector<int>& parent) {
    int n = residual.size();
    vector<bool> visited(n, false);
    queue<int> q;
    q.push(source);
    visited[source] = true;
    parent[source] = -1;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int v = 0; v < n; v++) {
            if (!visited[v] && residual[u][v] > 0) {
                parent[v] = u;
                if (v == sink) return true;
                visited[v] = true;
                q.push(v);
            }
        }
    }
    return false;
}

vector<vector<int>> diGraph_with_weight::BellmanFord(int start) {
    vector<int> distance(numVertices, INF);
    vector<int> predecessor(numVertices, -1);
    distance[start] = 0;

    for (int i = 0; i < numVertices - 1; ++i) {
        for (int u = 0; u < numVertices; ++u) {
            for (const auto& [v, weight] : adjacencyList[u]) {
                if (distance[u] != INF && distance[u] + weight < distance[v]) {
                    distance[v] = distance[u] + weight;
                    predecessor[v] = u;
                }
            }
        }
    }

    for (int u = 0; u < numVertices; ++u) {
        for (const auto& [v, weight] : adjacencyList[u]) {
            if (distance[u] != INF && distance[u] + weight < distance[v]) {
                cerr << "Граф містить від'ємний цикл!" << endl;
                return {};
            }
        }
    }

    vector<vector<int>> result(numVertices, vector<int>(numVertices, INF));
    for (int v = 0; v < numVertices; ++v) {
        if (predecessor[v] != -1) {
            result[predecessor[v]][v] = distance[v] - distance[predecessor[v]];
        }
    }

    cout << "Найкоротшi вiдстанi вiд вершини " << start << ":\n";
    for (int i = 0; i < numVertices; ++i) {
        cout << "До вершини " << i << ": " << (distance[i] == INF ? "INF" : to_string(distance[i])) << endl;
    }

    return result;
}

vector<vector<int>> diGraph_with_weight::FloydWarshall() {
    vector<vector<int>> result(numVertices, vector<int>(numVertices, INF));

    for (int i = 0; i < numVertices; ++i) {
        result[i][i] = 0;
    }

    for (int u = 0; u < numVertices; ++u) {
        for (const auto& [v, weight] : adjacencyList[u]) {
            result[u][v] = weight;
        }
    }

    for (int k = 0; k < numVertices; ++k) {
        for (int i = 0; i < numVertices; ++i) {
            for (int j = 0; j < numVertices; ++j) {
                if (result[i][k] != INF && result[k][j] != INF) {
                    result[i][j] = min(result[i][j], result[i][k] + result[k][j]);
                }
            }
        }
    }

    for (int i = 0; i < numVertices; ++i) {
        if (result[i][i] < 0) {
            cout << "Граф мiстить вiд'ємний цикл!" << endl;
            return {};
        }
    }


    cout << "Найкоротшi вiдстанi мiж усiма вершинами:\n";
    for (int i = 0; i < numVertices; ++i) {
        for (int j = 0; j < numVertices; ++j) {
            cout << "Вiд " << i << " до " << j << ": "
                << (result[i][j] == INF ? "INF" : to_string(result[i][j])) << endl;
        }
        cout << endl;
    }

    return result;
}

vector<vector<int>> diGraph_with_weight::Johnson() {
    int newVertex = numVertices;
    adjacencyList.push_back({});

    for (int i = 0; i < numVertices; ++i) {
        adjacencyList[newVertex].emplace_back(i, 0);
    }

    vector<int> h(numVertices + 1, INF);
    h[newVertex] = 0;

    for (int i = 0; i < numVertices; ++i) {
        for (int u = 0; u <= numVertices; ++u) {
            for (const auto& [v, weight] : adjacencyList[u]) {
                if (h[u] != INF && h[u] + weight < h[v]) {
                    h[v] = h[u] + weight;
                }
            }
        }
    }

    for (int u = 0; u <= numVertices; ++u) {
        for (const auto& [v, weight] : adjacencyList[u]) {
            if (h[u] != INF && h[u] + weight < h[v]) {
                cerr << "Граф мiстить вiд’ємний цикл!" << endl;
                return {};
            }
        }
    }

    adjacencyList.pop_back();

    vector<vector<int>> newWeights(numVertices, vector<int>(numVertices, INF));
    for (int u = 0; u < numVertices; ++u) {
        for (const auto& [v, weight] : adjacencyList[u]) {
            newWeights[u][v] = weight + h[u] - h[v];
        }
    }

    vector<vector<int>> result(numVertices, vector<int>(numVertices, INF));

    for (int src = 0; src < numVertices; ++src) {
        vector<int> dist(numVertices, INF);
        dist[src] = 0;
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
        pq.emplace(0, src);

        while (!pq.empty()) {
            int d = pq.top().first;
            int u = pq.top().second;
            pq.pop();

            if (d > dist[u]) continue;

            for (const auto& [v, weight] : adjacencyList[u]) {
                int newDist = dist[u] + newWeights[u][v];
                if (newDist < dist[v]) {
                    dist[v] = newDist;
                    pq.emplace(newDist, v);
                }
            }
        }

        for (int v = 0; v < numVertices; ++v) {
            if (dist[v] != INF) {
                result[src][v] = dist[v] + h[v] - h[src];
            }
        }
    }

    cout << "Найкоротшi вiдстанi мiж усiма вершинами:\n";
    for (int i = 0; i < numVertices; ++i) {
        for (int j = 0; j < numVertices; ++j) {
            cout << "Вiд " << i << " до " << j << ": "
                << (result[i][j] == INF ? "INF" : to_string(result[i][j])) << endl;
        }
        cout << endl;
    }

    return result;
}


int diGraph_with_weight::FordFulkerson(int source, int sink) {
    vector<vector<int>> capacity(numVertices, vector<int>(numVertices, 0));

    for (int u = 0; u < numVertices; ++u) {
        for (const auto& [v, w] : adjacencyList[u]) {
            capacity[u][v] = w;
        }
    }

    vector<vector<int>> residual = capacity;
    vector<bool> visited(numVertices);
    vector<pair<vector<int>, int>> paths;
    int maxFlow = 0, flow;

    do {
        fill(visited.begin(), visited.end(), false);
        vector<int> path;
        flow = dfs(residual, visited, source, sink, INT_MAX, path);

        if (flow > 0) {
            maxFlow += flow;
            reverse(path.begin(), path.end());
            paths.emplace_back(path, flow);
        }
    } while (flow > 0);

    cout << "Потоки по кожному шляху:\n";
    for (const auto& [path, flow] : paths) {
        cout << "Шлях: ";
        for (int v : path) {
            cout << v << " ";
        }
        cout << "| Потiк: " << flow << "\n";
    }

    return maxFlow;
}

int diGraph_with_weight::EdmondsKarp(int source, int sink) {
    vector<vector<int>> capacity(numVertices, vector<int>(numVertices, 0));

    for (int u = 0; u < numVertices; u++) {
        for (auto [v, w] : adjacencyList[u]) {
            capacity[u][v] = w;
        }
    }

    vector<vector<int>> residual = capacity;
    vector<int> parent(numVertices);
    int maxFlow = 0;

    cout << "Потоки по кожному шляху:\n";

    while (bfs(residual, source, sink, parent)) {
        int pathFlow = INF;

        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            if (residual[u][v] > 0) {
                pathFlow = min(pathFlow, residual[u][v]);
            }
        }

        vector<int> path;
        for (int v = sink; v != source; v = parent[v]) {
            path.push_back(v);
        }
        path.push_back(source);
        reverse(path.begin(), path.end());

        cout << "Шлях: ";
        for (int v : path) {
            cout << v << " ";
        }
        cout << "| Потiк: " << pathFlow << "\n";

        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            residual[u][v] -= pathFlow;
            residual[v][u] += pathFlow;
        }

        maxFlow += pathFlow;
    }

    return maxFlow;
}


void Graph_with_weight::generateDOTFile(const string& filename, const vector<vector<int>>& mst) {
    ofstream dotFile(filename);

    dotFile << "graph G {" << endl;
    dotFile << "  layout=neato;" << endl;
    dotFile << "  node [shape=circle];" << endl;

    double angle_step = 2 * PI / numVertices;
    for (int i = 0; i < numVertices; ++i) {
        double x = 10 * cos(i * angle_step);
        double y = 10 * sin(i * angle_step);
        dotFile << "  " << i << " [pos=\"" << x << "," << y << "!\"];" << endl;
    }

    for (int i = 0; i < numVertices; ++i) {
        for (const auto& neighbor : adjacencyList[i]) {
            int j = neighbor.first;
            int weight = neighbor.second;
            if (i < j) {
                string color = (mst[i][j] != INF && mst[i][j] == weight) ? "green" : "red";
                dotFile << "  " << i << " -- " << j << " [label=\"" << weight << "\", color=" << color << "]" << endl;
            }
        }
    }
    dotFile << "}" << endl;
}

void diGraph_with_weight::generateDOTFileDierected(const string& filename, const vector<vector<int>>& mst, int painted) {
    ofstream dotFile(filename);

    dotFile << "digraph G {" << endl;
    dotFile << "  layout=neato;" << endl;
    dotFile << "  node [shape=circle];" << endl;

    double angle_step = 2 * PI / numVertices;
    for (int i = 0; i < numVertices; ++i) {
        double x = 10 * cos(i * angle_step);
        double y = 10 * sin(i * angle_step);
        dotFile << "  " << i << " [pos=\"" << x << "," << y << "!\"];" << endl;
    }

    for (int i = 0; i < numVertices; ++i) {
        for (const auto& neighbor : adjacencyList[i]) {
            int j = neighbor.first;
            int weight = neighbor.second;

            string color;
            if (painted == false) {
                color = "black";
            }
            else {
                color = (mst[i][j] != INF && mst[i][j] == weight) ? "green" : "red";
            }

            dotFile << "  " << i << " -> " << j << " [label=\"" << weight << "\", color=" << color << "]" << endl;
        }
    }
    dotFile << "}" << endl;
}

int main() {
    setlocale(LC_ALL, "UKR");
    bool painted;
    int Working = 1;
    int Choice, Algorithm;
    cout << "------------------------------------------------------------------------------------------------------------------------" << std::endl;
    while (Working != 0) {
        cout << "Введiть який алгоритм роботи з графами вас цiкавить: " << endl;
        cout << "1 - Елементарнi алгоритми" << endl;
        cout << "2 - Мiнiмальнi остовнi дерева" << endl;
        cout << "3 - Найкоротшi шляхи" << endl;
        cout << "4 - Максимальний потiк" << endl;
        cout << "0 - Вихiд з програми" << endl;
        cin >> Choice;
        switch (Choice) {
        case 1:
            cout << "Оберiть алгоритм: " << endl;
            cout << "1 - Пошук в ширину " << endl;
            cout << "2 - Пошук в глибину " << endl;
            cin >> Algorithm;
            switch (Algorithm) {
            case 1:
            {
                Graph_without_weight g(7);
                g.addEdge(0, 1);
                g.addEdge(0, 3);
                g.addEdge(1, 2);
                g.addEdge(2, 6);
                g.addEdge(3, 4);
                g.addEdge(4, 5);
                g.bfs(0);
            }
            break;
            case 2:
            {
                Graph_without_weight g(7);
                g.addEdge(0, 1);
                g.addEdge(0, 3);
                g.addEdge(1, 2);
                g.addEdge(2, 6);
                g.addEdge(3, 4);
                g.addEdge(4, 5);
                g.dfsTraversal(0);
            }
            break;
            default:
                cout << "Введiть число, якi написанi в списку" << endl;
            }
        break;
        case 2:
            cout << "Оберiть алгоритм: " << endl;
            cout << "1 - Алгоритм Крускала " << endl;
            cout << "2 - Алгоритм Прiма " << endl;
            cin >> Algorithm;
            switch (Algorithm) {
            case 1:
            {
                Graph_with_weight g(7);
                g.addEdge(0, 1, 5);
                g.addEdge(0, 3, 4);
                g.addEdge(0, 4, 6);
                g.addEdge(1, 2, 2);
                g.addEdge(1, 6, 3);
                g.addEdge(2, 6, 5);
                g.addEdge(3, 4, 3);
                g.addEdge(4, 5, 1);

                vector<vector<int>> mst = g.Kruskal();
                g.generateDOTFile("Kruskal.dot", mst);

                system("circo -Tpng -Goverlap=prism Kruskal.dot -o Kruskal.png");
            }
            break;
            case 2:
            {
                Graph_with_weight g(7);
                g.addEdge(0, 1, 5);
                g.addEdge(0, 3, 4);
                g.addEdge(0, 4, 6);
                g.addEdge(1, 2, 2);
                g.addEdge(1, 6, 3);
                g.addEdge(2, 6, 5);
                g.addEdge(3, 4, 3);
                g.addEdge(4, 5, 1);

                vector<vector<int>> mst = g.Prim();
                g.generateDOTFile("Prim.dot", mst);

                system("circo -Tpng -Goverlap=prism Prim.dot -o Prim.png");
            }
            break;
            default:
                cout << "Введiть число, якi написанi в списку" << endl;
            }
        cout << "Алгоритм усiпшно виконав роботу i зображення було збережено" << endl;
        break;
        case 3:
            cout << "Оберiть алгоритм: " << endl;
            cout << "1 - Алгоритм Беллмана-Форда " << endl;
            cout << "2 - Алгоритм Дейкстри " << endl;
            cout << "3 - Алгоритм Флойда-Варшалла " << endl;
            cout << "4 - Алгоритм Джонсона " << endl;
            cin >> Algorithm;
            painted = true;
            switch (Algorithm) {
            case 1:
            {
                diGraph_with_weight g(5);
                g.addEdge(0, 1, -1);
                g.addEdge(0, 2, 4);
                g.addEdge(1, 2, 3);
                g.addEdge(1, 3, 2);
                g.addEdge(1, 4, 2);
                g.addEdge(3, 2, 5);
                g.addEdge(3, 1, 1);
                g.addEdge(4, 3, -3);

                vector<vector<int>> mst = g.BellmanFord(0);
                g.generateDOTFileDierected("Bellman-Ford.dot", mst, painted);

                system("circo -Tpng -Goverlap=prism Bellman-Ford.dot -o Bellman-Ford.png");
            }
            break;
            case 2:
            {
                Graph_with_weight g(6);
                g.addEdge(0, 1, 7);
                g.addEdge(0, 2, 9);
                g.addEdge(0, 5, 14);
                g.addEdge(1, 2, 10);
                g.addEdge(1, 3, 15);
                g.addEdge(2, 5, 2);
                g.addEdge(2, 3, 11);
                g.addEdge(3, 4, 6);
                g.addEdge(4, 5, 9);

                vector<vector<int>> mst = g.Dijkstra(0);
                g.generateDOTFile("Dijkstra.dot", mst);

                system("circo -Tpng -Goverlap=prism Dijkstra.dot -o Dijkstra.png");
            }
            break;
            case 3:
            {
                diGraph_with_weight g(4);
                g.addEdge(0, 1, 1);
                g.addEdge(0, 2, 6);
                g.addEdge(1, 2, 4);
                g.addEdge(1, 3, 1);
                g.addEdge(3, 2, 1);

                vector<vector<int>> mst = g.FloydWarshall();
                g.generateDOTFileDierected("FloydWarshall.dot", mst, painted);

                system("circo -Tpng -Goverlap=prism FloydWarshall.dot -o FloydWarshall.png");
            }
            break;
            case 4:
            {
                diGraph_with_weight g(5);
                g.addEdge(0, 1, -1);
                g.addEdge(0, 2, 4);
                g.addEdge(1, 2, 3);
                g.addEdge(1, 3, 2);
                g.addEdge(1, 4, 2);
                g.addEdge(3, 2, 5);
                g.addEdge(3, 1, 1);
                g.addEdge(4, 3, -3);

                vector<vector<int>> mst = g.Johnson();
                g.generateDOTFileDierected("Johnson.dot", mst, painted);

                system("circo -Tpng -Goverlap=prism Johnson.dot -o Johnson.png");
            }
            break;
            default:
                cout << "Введiть число, якi написанi в списку" << endl;
            }
        cout << "Алгоритм усiпшно виконав роботу i зображення було збережено" << endl;
        break;
        case 4:
            cout << "Оберiть алгоритм: " << endl;
            cout << "1 - Алгоритм Форда-Фалкерсона" << endl;
            cout << "2 - Алгоритм Едмондса-Карпа" << endl;
            cin >> Algorithm;
            painted = false;
            switch (Algorithm) {
            case 1:
            {
                int vertices = 6;
                diGraph_with_weight graph(vertices);

                graph.addEdge(0, 1, 7);
                graph.addEdge(0, 2, 4);
                graph.addEdge(1, 2, 4);
                graph.addEdge(1, 4, 2);
                graph.addEdge(2, 4, 8);
                graph.addEdge(2, 3, 4);
                graph.addEdge(3, 5, 12);
                graph.addEdge(4, 3, 4);
                graph.addEdge(4, 5, 5);

                int maxFlow = graph.FordFulkerson(0, 5);
                cout << "Максимальний потiк: " << maxFlow << endl;

                vector<vector<int>> dummyMST(vertices, vector<int>(vertices, INF));
                graph.generateDOTFileDierected("FordFulkerson.dot", dummyMST, 0);

                system("circo -Tpng -Goverlap=prism FordFulkerson.dot -o FordFulkerson.png");
            }
            break;
            case 2:
            {
                int vertices = 6;
                diGraph_with_weight graph(vertices);

                graph.addEdge(0, 1, 7);
                graph.addEdge(0, 2, 4);
                graph.addEdge(1, 2, 4);
                graph.addEdge(1, 4, 2);
                graph.addEdge(2, 4, 8);
                graph.addEdge(2, 3, 4);
                graph.addEdge(3, 5, 12);
                graph.addEdge(4, 3, 4);
                graph.addEdge(4, 5, 5);

                int maxFlow = graph.EdmondsKarp(0, 5);
                cout << "Максимальний потiк: " << maxFlow << endl;

                vector<vector<int>> dummyMST(vertices, vector<int>(vertices, INF));
                graph.generateDOTFileDierected("EdmondsKarp.dot", dummyMST, 0);

                system("circo -Tpng -Goverlap=prism EdmondsKarp.dot -o EdmondsKarp.png");
            }
            break;
            default:
                cout << "Введiть число, якi написанi в списку" << endl;
            }
        cout << "Алгоритм усiпшно виконав роботу i зображення було збережено" << endl;
        break;
        case 0:
        {
            Working = 0;
        }
        break;
        default:
            cout << "Введiть число, якi написанi в списку" << endl;
        }
        cout << "------------------------------------------------------------------------------------------------------------------------" << endl;
    }

    cout << "\n" << endl;
    system("pause");
    return 0;
}
