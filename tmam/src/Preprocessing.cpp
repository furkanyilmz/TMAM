#define _SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING

#include "Preprocessing.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <limits>

bool Preprocessing::Cut(const std::string& inputFilePath, const std::string& outputFilePath) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if (inputFilePath.size() > 4 && inputFilePath.compare(inputFilePath.size() - 4, 4, ".obj") == 0) {
        igl::readOBJ(inputFilePath, V, F);
    }
    else {
        igl::readOFF(inputFilePath, V, F);
    }
    if (isWatertight(V, F)) {
        auto [topIdx, bottomIdx] = findTopBottomVertices(V);

        // Compute the shortest path using Dijkstra's algorithm
        Eigen::VectorXi path = computeShortestPath(V, F, topIdx, bottomIdx);

        // Convert the path to cuts format required by igl::cut_mesh
        Eigen::MatrixXi cuts = Eigen::MatrixXi::Zero(F.rows(), 3);
        for (int i = 0; i < path.size() - 1; ++i) {
            int v0 = path(i);
            int v1 = path(i + 1);

            // Iterate over each face
            for (int j = 0; j < F.rows(); ++j) {
                // Check each edge of the face
                for (int k = 0; k < F.cols(); ++k) {
                    int fk0 = F(j, k);
                    int fk1 = F(j, (k + 1) % F.cols());
                    if ((fk0 == v0 && fk1 == v1) || (fk0 == v1 && fk1 == v0)) {
                        cuts(j, k) = 1; // Mark this edge to be cut
                    }
                }
            }
        }
        Eigen::MatrixXd Vcut(V.rows(), 3);
        Eigen::MatrixXi Fcut(F.rows(), 3);
        igl::cut_mesh(V, F, cuts, Vcut, Fcut);
        igl::writeOBJ(outputFilePath, Vcut, Fcut);
    }
    else {
        igl::writeOBJ(outputFilePath, V, F);
    }
    return true;
}

std::pair<int, int> Preprocessing::findTopBottomVertices(const Eigen::MatrixXd& V) {
    int topIdx = 0, bottomIdx = 0;
    double maxY = V(0, 1), minY = V(0, 1);

    for (int i = 1; i < V.rows(); ++i) {
        if (V(i, 1) > maxY) {
            maxY = V(i, 1);
            topIdx = i;
        }
        if (V(i, 1) < minY) {
            minY = V(i, 1);
            bottomIdx = i;
        }
    }

    return { topIdx, bottomIdx };
}



Eigen::VectorXi Preprocessing::computeShortestPath(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int startIdx, int endIdx) {
    int numVertices = V.rows();
    std::vector<std::vector<std::pair<int, double>>> adjList(numVertices);

    // Create an adjacency list for the mesh
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < F.cols(); ++j) {
            int v1 = F(i, j);
            int v2 = F(i, (j + 1) % F.cols());
            double weight = euclideanDistance(V.row(v1), V.row(v2));
            adjList[v1].emplace_back(v2, weight);
            adjList[v2].emplace_back(v1, weight); // For undirected graph
        }
    }

    // Dijkstra's algorithm
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> minHeap;
    std::vector<double> distances(numVertices, std::numeric_limits<float>::infinity());
    std::vector<int> previous(numVertices, -1);

    minHeap.push({ 0.0, startIdx });
    distances[startIdx] = 0.0;

    while (!minHeap.empty()) {
        int u = minHeap.top().second;
        minHeap.pop();

        if (u == endIdx) {
            break; // Shortest path to endIdx found
        }

        for (const auto& [v, weight] : adjList[u]) {
            double distThroughU = distances[u] + weight;
            if (distThroughU < distances[v]) {
                distances[v] = distThroughU;
                previous[v] = u;
                minHeap.push({ distances[v], v });
            }
        }
    }

    // Reconstruct the shortest path
    std::vector<int> path;
    for (int at = endIdx; at != -1; at = previous[at]) {
        path.push_back(at);
    }
    std::reverse(path.begin(), path.end());

    // Convert path to Eigen::VectorXi
    Eigen::VectorXi pathVec(path.size());
    for (size_t i = 0; i < path.size(); ++i) {
        pathVec(i) = path[i];
    }

    return pathVec;
}

bool Preprocessing::isWatertight(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) {
    Eigen::MatrixXi B;
        igl::boundary_facets(F, B);
        // Check for boundary edges and whether the mesh is both edge and vertex manifold
        if (B.rows() != 0 || !igl::is_edge_manifold(F) || !igl::is_vertex_manifold(F)) {
            return false;
        }
        std::cout << "This object is watertight" << std::endl;
        return true;
}


double Preprocessing::euclideanDistance(const Eigen::RowVector3d& a, const Eigen::RowVector3d& b) {
    return (a - b).norm();
}


bool Preprocessing::CheckifOffFileisWatertightandConvert(const std::string& inputFilePath, const std::string& outputFilePath) {
    if (!ReadOffFile(inputFilePath)) {
        std::cerr << "Failed to read OFF file.\n";
        return false;
    }
    if (isWatertight(vertices, triangles)) {
        std::cout << "This object is watertight.\n";
        Preprocessing::splitVerticesAlongSeam(findSeamPath());

    }
    if (!WriteObjFile(outputFilePath)) {
        std::cerr << "Failed to write OBJ file.\n";
        return false;
    }
    return true;
}

bool Preprocessing::CheckifObjFileisWatertight(const std::string& inputFilePath, const std::string& outputFilePath){
    if (!ReadObjFile(inputFilePath)) {
        std::cerr << "Failed to read OFF file.\n";
        return false;
    }
    if (isWatertight(vertices, triangles)) {
        std::cout << "This object is watertight.\n";
        Preprocessing::splitVerticesAlongSeam(findSeamPath());
    }
    if (!WriteObjFile(outputFilePath)) {
        std::cerr << "Failed to write OBJ file.\n";
        return false;
    }
    return true;
}


bool Preprocessing::ReadOffFile(const std::string& filePath) {
    std::ifstream offFile(filePath);
    if (!offFile.is_open()) {
        std::cerr << "Failed to open .off file\n";
        return false;
    }

    std::string line;
    getline(offFile, line); // Read the header
    if (line != "OFF") {
        std::cerr << "Not a valid OFF file\n";
        return false;
    }

    int numVertices, numFaces, numEdges;
    offFile >> numVertices >> numFaces >> numEdges;

    vertices.resize(numVertices);
    for (int i = 0; i < numVertices; ++i) {
        offFile >> vertices[i].x >> vertices[i].y >> vertices[i].z;
    }

    triangles.resize(numFaces);
    for (int i = 0; i < numFaces; ++i) {
        int vertsInFace;
        offFile >> vertsInFace;
        if (vertsInFace != 3) {
            std::cerr << "Non-triangle face encountered.\n";
            return false;
        }
        offFile >> triangles[i].v1 >> triangles[i].v2 >> triangles[i].v3;
    }

    return true;
}
bool Preprocessing::ReadObjFile(const std::string& filePath) {
    std::ifstream objFile(filePath);
    if (!objFile.is_open()) {
        std::cerr << "Failed to open .obj file\n";
        return false;
    }
    std::string line;
    while (getline(objFile, line)) {
        std::stringstream ss(line);
        std::string lineHeader;
        ss >> lineHeader;
        if (lineHeader == "vt") {
            // Vertex definition
            TextureCoord coord;
            ss >> coord.u >> coord.v;
            textureCoords.push_back(coord);
        }
        else if (lineHeader == "vn") {
            // Vertex definition
            Normal normal;
            ss >> normal.x >> normal.y >> normal.z;
            normals.push_back(normal);
        }
        else if (lineHeader == "v") {
            // Vertex definition
            Vertex vertex;
            ss >> vertex.x >> vertex.y >> vertex.z;
            vertices.push_back(vertex);
        }
        else if (lineHeader == "f") {
            Triangle triangle;
            std::string v1, v2, v3;
            ss >> v1 >> v2 >> v3;
            std::istringstream v1ss(v1), v2ss(v2), v3ss(v3);
            v1ss >> triangle.v1; v1ss.ignore(); v1ss >> triangle.vt1; v1ss.ignore(); v1ss >> triangle.vn1;
            v2ss >> triangle.v2; v2ss.ignore(); v2ss >> triangle.vt2; v2ss.ignore(); v2ss >> triangle.vn2;
            v3ss >> triangle.v3; v3ss.ignore(); v3ss >> triangle.vt3; v3ss.ignore(); v3ss >> triangle.vn3;

            // Adjust indices from 1-based to 0-based
            triangle.v1--; triangle.v2--; triangle.v3--;
            triangle.vt1--; triangle.vt2--; triangle.vt3--;
            triangle.vn1--; triangle.vn2--; triangle.vn3--;

            triangles.push_back(triangle);
        }
        // Only vertices and triangular faces are processed
    }

    return true;
}

void Preprocessing::addEdge(std::map<std::pair<int, int>, int>& edgeCount, int v1, int v2) {
    auto edge = v1 < v2 ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
    edgeCount[edge]++;
}

bool Preprocessing::isWatertight(const std::vector<Vertex>& vertices, const std::vector<Triangle>& triangles) {
    std::map<std::pair<int, int>, int> edgeCount;

    for (const auto& triangle : triangles) {
        addEdge(edgeCount, triangle.v1, triangle.v2);
        addEdge(edgeCount, triangle.v2, triangle.v3);
        addEdge(edgeCount, triangle.v3, triangle.v1);
    }

    for (const auto& edge : edgeCount) {
        if (edge.second != 2) {  // If any edge is not shared exactly twice
            return false;
        }
    }
    return true;  // All edges are shared exactly twice
}


bool Preprocessing::WriteObjFile(const std::string& filePath) {
    std::ofstream objFile(filePath);
    if (!objFile.is_open()) {
        std::cerr << "Failed to open .obj file\n";
        return false;
    }

    for (const auto& vertex : vertices) {
        objFile << "v " << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
    }

    for (const auto& textureCoord : textureCoords) {
        objFile << "vt " << textureCoord.u << " " << textureCoord.v << "\n";
    }

    for (const auto& normal : normals) {
        objFile << "vn " << normal.x << " " << normal.y << " " << normal.z << "\n";
    }

    for (const auto& triangle : triangles) {
        objFile << "f ";
        if (!vertices.empty()) {
            objFile << triangle.v1+1;
            if (!textureCoords.empty() || !normals.empty()) objFile << "/";
        }
        if (!textureCoords.empty()) {
            objFile << triangle.vt1+1;
            if (!normals.empty()) objFile << "/";
        }
        if (!normals.empty()) {
            objFile << triangle.vn1+1;
        }

        objFile << " ";

        if (!vertices.empty()) {
            objFile << triangle.v2+1;
            if (!textureCoords.empty() || !normals.empty()) objFile << "/";
        }
        if (!textureCoords.empty()) {
            objFile << triangle.vt2+1;
            if (!normals.empty()) objFile << "/";
        }
        if (!normals.empty()) {
            objFile << triangle.vn2+1;
        }

        objFile << " ";

        if (!vertices.empty()) {
            objFile << triangle.v3+1;
            if (!textureCoords.empty() || !normals.empty()) objFile << "/";
        }
        if (!textureCoords.empty()) {
            objFile << triangle.vt3+1;
            if (!normals.empty()) objFile << "/";
        }
        if (!normals.empty()) {
            objFile << triangle.vn3+1;
        }

        objFile << "\n";
    }

    return true;
}

int Preprocessing::findHighestVertexIndex() {
    int highestIndex = 0;
    float highestY = this->vertices[0].y;
    for (size_t i = 1; i < vertices.size(); ++i) {
        if (vertices[i].y > highestY) {
            highestY = vertices[i].y;
            highestIndex = i;
        }
    }
    return highestIndex;
}

int Preprocessing::findLowestVertexIndex() {
    int lowestIndex = 0;
    float lowestY = vertices[0].y;
    for (size_t i = 1; i < vertices.size(); ++i) {
        if (vertices[i].y < lowestY) {
            lowestY = vertices[i].y;
            lowestIndex = i;
        }
    }
    return lowestIndex;
}


std::vector<int> Preprocessing::dijkstra(const Graph& graph, int sourceIndex, int destinationIndex) {
    std::map<int, float> distances;
    std::map<int, int> predecessors;
    std::set<std::pair<float, int>> vertexQueue;
    auto start = std::chrono::high_resolution_clock::now();
    // Initialize distances and queue
    for (int i = 0; i < vertices.size(); i++) {
        distances[i] = std::numeric_limits<float>::infinity();
        predecessors[i] = -1;
        vertexQueue.insert({ distances[i], i });
    }
    vertexQueue.erase({ distances[sourceIndex], sourceIndex }); // Remove the initial entry
    distances[sourceIndex] = 0.0f;
    vertexQueue.insert({ distances[sourceIndex], sourceIndex }); // Insert with updated distance

    while (!vertexQueue.empty()) {
        int currentVertex = vertexQueue.begin()->second;
        vertexQueue.erase(vertexQueue.begin());
        //std::cout << "Current vertex: " << currentVertex << ", Distance: " << distances[currentVertex] << std::endl;
        if (currentVertex == destinationIndex) {
            break; // Destination reached
        }

        // Explore neighbors
        for (auto& edge : graph.adjacencyList.at(currentVertex)) {
            int neighbor = edge.vertexIndex2;
            //std::cout << "Current neighbor: " << neighbor << std::endl;

            float altDistance = distances[currentVertex] + edge.weight;

            if (altDistance < distances[neighbor]) {
                vertexQueue.erase({ distances[neighbor], neighbor });
                distances[neighbor] = altDistance;
                predecessors[neighbor] = currentVertex;
                vertexQueue.insert({ distances[neighbor], neighbor });
            }
        }
    }

    // Reconstruct the shortest path
    std::vector<int> path;
    if (predecessors[destinationIndex] != -1 || destinationIndex == sourceIndex) {
        for (int at = destinationIndex; at != -1; at = predecessors[at]) {
            path.push_back(at);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9; // convert to seconds
    std::cout << "Time taken by program is : " << std::fixed << time_taken << std::setprecision(9);
    std::cout << " sec" << std::endl;
    std::reverse(path.begin(), path.end());
    return path;
}

void Preprocessing::splitVerticesAlongSeam(const std::vector<int>& seamPath) {
    std::map<int, int> originalToDuplicateMap;  // Map from original vertex index to new duplicated vertex index

    // Step 1: Duplicate vertices along the seam
    for (int i = 0; i < seamPath.size(); i++) {
        Vertex originalVertex = vertices[seamPath[i]];
        vertices.push_back(originalVertex);  // Duplicate vertex

        int duplicatedIndex = vertices.size() - 1;
        originalToDuplicateMap[seamPath[i]] = duplicatedIndex;
    }

    for (const auto& pair : originalToDuplicateMap) {
        std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
    }


    // Step 2: Update triangle indices
    for (Triangle& triangle : triangles) {
        if (originalToDuplicateMap.count(triangle.v1)) {
            if (determineSide(triangle, seamPath, vertices)) {
                triangle.v1 = originalToDuplicateMap[triangle.v1];
            }
        }
        if (originalToDuplicateMap.count(triangle.v2)) {
            if (determineSide(triangle, seamPath, vertices)) {
                triangle.v2 = originalToDuplicateMap[triangle.v2];
            }
        }
        if (originalToDuplicateMap.count(triangle.v3)) {
            if (determineSide(triangle, seamPath, vertices)) {
                triangle.v3 = originalToDuplicateMap[triangle.v3];
            }
        }
    }
}

std::vector<int> Preprocessing::findSeamPath() {
    int highestVertexIndex = Preprocessing::findHighestVertexIndex();
    int lowestVertexIndex = Preprocessing::findLowestVertexIndex();
    Preprocessing::Graph graph = Preprocessing::createGraph();
    return Preprocessing::dijkstra(graph, highestVertexIndex, lowestVertexIndex);
}

bool Preprocessing::determineSide(const Triangle& triangle, const std::vector<int>& seamPath, const std::vector<Vertex>& vertices) {
    Vertex seamStart(0.0f, 0.0f, 0.0f);
    Vertex seamEnd(0.0f, 0.0f, 0.0f);
    int vertexnotonpathindex=-1;
    // Calculate the seam vector (from the first to second vertex of the seam)
    for (int vertice : {triangle.v1, triangle.v2, triangle.v3}){
        int found = 0;
        for (int i = 0;i < seamPath.size();i++) {
            if (vertice == seamPath[i]) {
                if (i == seamPath.size() - 1) {
                    seamStart = vertices[seamPath[i-1]];
                    seamEnd = vertices[seamPath[i]];
                    //std::cout << seamPath[i - 1] << " " << seamPath[i] << std::endl;
                }
                else {
                    seamStart = vertices[seamPath[i]];
                    seamEnd = vertices[seamPath[i+1]];
                    //std::cout << seamPath[i] << " " << seamPath[i+1] << std::endl;
                }
            found = 1;
            }
        }
        if (!found) {
            vertexnotonpathindex = vertice;
            if (std::find(seamPath.begin(), seamPath.end(), vertexnotonpathindex) != seamPath.end()) {
                std::cout << "error lololololol" << std::endl;
            }
        }
    }
    Vertex seamVector = { seamEnd.x - seamStart.x, seamEnd.y - seamStart.y, seamEnd.z - seamStart.z };

    // Calculate a normal to the seam (assuming Y-up coordinate system)
    Vertex upVector = { 0, 1, 0 };
    Vertex normalToSeam = crossProduct(seamVector, upVector);
    Vertex vertexnotonpath = vertices[vertexnotonpathindex];

    // Create a vector from seam start to this triangle vertex
    Vertex vectorToTriangle = { vertexnotonpath.x - seamStart.x, vertexnotonpath.y - seamStart.y, vertexnotonpath.z - seamStart.z };

    // Cross product to determine the side
    Vertex crossProd = crossProduct(seamVector, vectorToTriangle);
    return dotProduct(crossProd, normalToSeam) > 0; // Positive if on one side, negative if on the other
}

// Utility function to calculate cross product
Preprocessing::Vertex Preprocessing::crossProduct(const Vertex& a, const Vertex& b) {
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}

// Utility function to calculate dot product
float Preprocessing::dotProduct(const Vertex& a, const Vertex& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}