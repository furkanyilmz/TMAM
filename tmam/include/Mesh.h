#pragma once
#pragma once

#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <sstream>
#include <unordered_set>
#include <fstream>
#include <Eigen/Sparse>

using namespace std;

struct UVCoordinates
{
	float u;
	float v;
};

struct Vertex
{
	float* coords;
	glm::vec3 normal;
	int idx; //who am i; verts[idx]
	bool isBoundary;
	bool visited; //to be used in finding the longest boundary chain
	UVCoordinates uvCoordinates;
	UVCoordinates bffuvCoordinates;
	vector< int > vertList; //adj vvertices;
	vector< int > triList;
	vector< int > edgeList;

	Vertex(int i, float* c) : idx(i), coords(c), normal(0.0f, 0.0f, 0.0f), isBoundary(0), uvCoordinates({1.1f, 1.1f}), visited(0) {};
};

struct Edge
{
	int idx; //edges[idx]
	int v1i, v2i; //endpnts
	float length;
	bool isBoundary;
	Edge(int id, int v1, int v2) : idx(id), v1i(v1), v2i(v2), length(-1), isBoundary(1) {};
};


struct Triangle
{
	int idx; //tris[idx]
	int v1i, v2i, v3i;
	glm::vec3 normal;
	Triangle(int id, int v1, int v2, int v3, glm::vec3 normal) : idx(id), v1i(v1), v2i(v2), v3i(v3), normal(normal) {};
};

struct DistortionResult {
	float averageDistortionMethod1;
	float averageDistortionMethod2;
};

struct AngleDistortionResult {
	float averageDistortionUV;
	float averageDistortionBFFUV;
};


class Mesh
{
private:
	void addTriangle(int v1, int v2, int v3);
	void addEdge(int v1, int v2);
	void addVertex(float x, float y, float z);
	bool makeVertsNeighbor(int v1i, int v2i);
	float computeVertexDistance(int v1i, int v2i);
	glm::vec3 computeTriangleNormal(int v1i, int v2i, int v3i);
	void computeVertexNormals();
	void addNormal(float x, float y, float z);
	void makeEdgeNonBoundary(int v1i, int v2i);
	void findBoundaryVertices();
	void findLongestBoundaryChain();
	std::vector<int> findBoundaryChain(Vertex* v);

public:
	vector< Vertex* > verts;
	vector< Triangle* > tris;
	vector< Edge* > edges;
	vector< glm::vec3 > normals;
	vector<int> boundaryChain;
	Eigen::SparseMatrix<double>* laplacianMatrix;
	Eigen::VectorXd vectorX;
	Eigen::VectorXd vectorY;
	int boundaryChainStartVertex;
	int boundaryChainLength;
	int linearUnwrapType; //0 for uniform, 1 for mean value, 2 for harmonic  
	Mesh() {};
	void createCube(float side);
	void loadOff(const char* name);
	void loadObj(const char* name);
	std::pair<int, int> findMiddleVertex(int v1i, int v2i);
	float angleBetweenVectors(const glm::vec3& v1, const glm::vec3& v2);
	DistortionResult calculateTriangleAreaDistortion();
	AngleDistortionResult calculateAngleBasedDistortion();
	void exportMesh(std::string& directory, const bool isTextureLoaded, bool useBffForTexture);
	float calculateTriangleArea3D(int triIndex);
	float calculateTriangleAreaUVMethod1(int triIndex);
	float calculateTriangleAreaUVMethod2(int triIndex);
};
