#pragma once

#ifndef OFFTOOBJCONVERTER_H
#define OFFTOOBJCONVERTER_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <queue>
#include <set>
#include <limits>
#include <algorithm>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/cut_mesh.h>
#include <igl/writeOBJ.h>
#include <igl/is_vertex_manifold.h>
#include <igl/is_edge_manifold.h>
#include <igl/boundary_facets.h>



class Preprocessing {
public:
    bool CheckifOffFileisWatertightandConvert(const std::string& inputFilePath, const std::string& outputFilePath);
    bool CheckifObjFileisWatertight(const std::string& inputFilePath, const std::string& outputFilePath);
    bool Cut(const std::string& inputFilePath, const std::string& outputFilePath);

private:
    struct Vertex {
        float x, y, z;
        Vertex() {}
        Vertex(float x, float y, float z) : x(x), y(y), z(z) {}
    };

    struct TextureCoord {
        float u, v;
    };

    struct Normal {
        float x, y, z;
    };

    struct Triangle {
        int v1, v2, v3;
        int vt1, vt2, vt3;
        int vn1, vn2, vn3;
    };

    struct Edge {
        int vertexIndex1;
        int vertexIndex2;
        float weight;
    };

    struct Graph {
        std::map<int, std::vector<Edge>> adjacencyList;

        void addEdge(int v1, int v2, float weight) {
            Edge edge = { v1, v2, weight };
            adjacencyList[v1].push_back(edge);

            Edge reverseEdge = { v2, v1, weight };
            adjacencyList[v2].push_back(reverseEdge);
        }

        std::vector<Edge> getEdges(int vertexIndex) {
            return adjacencyList[vertexIndex];
        }
    };

    Graph createGraph() {
        Graph graph;
        for (const auto& triangle : triangles) {
            float weight1 = calculateDistance(vertices[triangle.v1], vertices[triangle.v2]);
            float weight2 = calculateDistance(vertices[triangle.v2], vertices[triangle.v3]);
            float weight3 = calculateDistance(vertices[triangle.v3], vertices[triangle.v1]);

            graph.addEdge(triangle.v1, triangle.v2, weight1);
            graph.addEdge(triangle.v2, triangle.v3, weight2);
            graph.addEdge(triangle.v3, triangle.v1, weight3);
        }
        std::cout << "no errors I think?\n";
        return graph;
    }

    float calculateDistance(const Vertex& v1, const Vertex& v2) {
        return sqrt(pow(v1.x - v2.x, 2) + pow(v1.y - v2.y, 2) + pow(v1.z - v2.z, 2));
    }

    std::vector<Vertex> vertices;
    std::vector<Normal> normals;
    std::vector<TextureCoord> textureCoords;
    std::vector<Triangle> triangles;
    
    std::pair<int, int> findTopBottomVertices(const Eigen::MatrixXd& V);
    double euclideanDistance(const Eigen::RowVector3d& a, const Eigen::RowVector3d& b);
    Eigen::VectorXi computeShortestPath(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int startIdx, int endIdx);
    bool isWatertight(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

    bool ReadOffFile(const std::string& filePath);
    bool ReadObjFile(const std::string& filePath);
    bool WriteObjFile(const std::string& filePath);
    bool isWatertight(const std::vector<Vertex>& vertices, const std::vector<Triangle>& triangles);
    void addEdge(std::map<std::pair<int, int>, int>& edgeCount, int v1, int v2);
    int findLowestVertexIndex();
    int findHighestVertexIndex();
    std::vector<int> dijkstra(const Graph& graph,int sourceIndex, int destinationIndex);
    std::vector<int> findSeamPath();
    void splitVerticesAlongSeam(const std::vector<int>& seamPath);
    bool determineSide(const Triangle& triangle, const std::vector<int>& seamPath, const std::vector<Vertex>& vertices);
    Vertex crossProduct(const Vertex& a, const Vertex& b);
    float dotProduct(const Vertex& a, const Vertex& b);

};

#endif // OFFTOOBJCONVERTER_H
