#pragma once

#include <Eigen/Sparse>
#include "mesh.h"
#include <cmath>
#include "Utils.h"
#include <random>
#include <chrono>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void computeLaplacianMatrix(Mesh* mesh, LinearUnwrappingType linearUnwrappingType);

void computeUVCoordinates(Mesh* mesh, std::vector<glm::vec2>& arbitraryPolygonVertices, LinearUnwrappingType linearUnwrappingType, MappingType mappingType);

void computeBoundaryUVCoordinates(std::vector< int >& boundaryChain, std::vector< Vertex* >& vertices, Eigen::VectorXd& vectorX, Eigen::VectorXd& vectorY, 
	std::vector<glm::vec2>& arbitraryPolygonVertices, MappingType mappingType);

void mapToDisc(std::vector< int >& boundaryChain, std::vector< Vertex* >& vertices, Eigen::VectorXd& vectorX, Eigen::VectorXd& vectorY);

void mapToSquare(std::vector< int >& boundaryChain, std::vector< Vertex* >& vertices, Eigen::VectorXd& vectorX, Eigen::VectorXd& vectorY);

void mapToConvexHull(std::vector< int >& boundaryChain, std::vector< Vertex* >& vertices, Eigen::VectorXd& vectorX, Eigen::VectorXd& vectorY);

void mapToArbitraryPolygon(std::vector< int >& boundaryChain, std::vector< Vertex* >& vertices, Eigen::VectorXd& vectorX, Eigen::VectorXd& vectorY,
	std::vector<glm::vec2>& arbitraryPolygonVertices);

double computeMeanValueWeight(std::vector<Vertex*>& vertices, int v1i, int v2i, std::pair<int, int> neighbor_indices);

double computeMeanValueWeight(std::vector<Vertex*>& vertices, int v1i, int v2i, int neighbor_index);

double computeHarmonicValueWeight(std::vector<Vertex*>& vertices, int v1i, int v2i, std::pair<int, int> neighbor_indices);

double computeHarmonicValueWeight(std::vector<Vertex*>& vertices, int v1i, int v2i, int neighbor_indices);

std::vector<int> findConvexHull(std::vector<int> boundaryChain, std::vector<Vertex*>& vertices);

double computePolarAngle(std::vector<Vertex*>& vertices, int v1i, int v2i);

bool comparePolarAngle(std::vector<Vertex*>& vertices, int refVertex, int v1i, int v2i);

bool isCounterClockwise(std::vector<Vertex*>& vertices, int v1i, int v2i, int v3i);

void normalizeUVCoordinates(std::vector< int >& boundaryChain, std::vector< Vertex* >& vertices, Eigen::VectorXd& vectorX, Eigen::VectorXd& vectorY);

void improveUnwrapping(Mesh* mesh);

void solveLaplacianMatrix(Mesh* mesh);
