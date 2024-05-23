#include "Unwrap.h"
#include <iostream>

void computeLaplacianMatrix(Mesh* mesh, LinearUnwrappingType linearUnwrappingType)
{
	std::random_device rd;
	std::mt19937 gen(rd());

	// Define a distribution for floating-point values between 0 and 10
	std::uniform_real_distribution<double> dis(0.0, 10.0);

	std::vector<Eigen::Triplet<double>> tripletList;

	for (int v_index = 0; v_index < mesh->verts.size(); v_index++)
	{
		Vertex* v1 = mesh->verts[v_index];

		if (v1->uvCoordinates.u != 1.1f && v1->uvCoordinates.v != 1.1f)
		{
			tripletList.push_back(Eigen::Triplet<double>(v_index, v_index, 1));
		}
		else
		{
			double weight_sum = 0;
			for (int neighbor_v_index = 0; neighbor_v_index < v1->vertList.size(); neighbor_v_index++)
			{
				double weight;
				std::pair<int, int> middles;

				switch (linearUnwrappingType)
				{
				case LinearUnwrappingType::UNIFORM_WEIGHTS:
					weight = 1;
					break;

				case LinearUnwrappingType::MEAN_VALUE_WEIGHTS:
					middles = mesh->findMiddleVertex(v1->idx, v1->vertList[neighbor_v_index]);
					weight = computeMeanValueWeight(mesh->verts, v1->idx, v1->vertList[neighbor_v_index], middles);
					break;

				case LinearUnwrappingType::HARMONIC_WEIGHTS:
					middles = mesh->findMiddleVertex(v1->idx, v1->vertList[neighbor_v_index]);
					weight = computeHarmonicValueWeight(mesh->verts, v1->idx, v1->vertList[neighbor_v_index], middles);

					break;

				case LinearUnwrappingType::RANDOM_WEIGHTS:
					weight = dis(gen);
					break;
				default:
					break;
				}
				tripletList.push_back(Eigen::Triplet<double>(v_index, v1->vertList[neighbor_v_index], -weight));

				weight_sum += weight;
			}
			tripletList.push_back(Eigen::Triplet<double>(v_index, v_index, weight_sum));
		}
	}

	mesh->laplacianMatrix->setFromTriplets(tripletList.begin(), tripletList.end());
}

void computeUVCoordinates(Mesh* mesh, std::vector<glm::vec2>& arbitraryPolygonVertices, LinearUnwrappingType linearUnwrappingType, MappingType mappingType)
{
	for (int i = 0; i < mesh->verts.size(); i++)
	{
		Vertex* v = mesh->verts[i];
		v->uvCoordinates = { 1.1f, 1.1f };
	}

	computeBoundaryUVCoordinates(mesh->boundaryChain, mesh->verts, mesh->vectorX, mesh->vectorY, arbitraryPolygonVertices, mappingType);

	computeLaplacianMatrix(mesh, linearUnwrappingType);

	solveLaplacianMatrix(mesh);

}

void solveLaplacianMatrix(Mesh* mesh)
{
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.analyzePattern(*(mesh->laplacianMatrix));
	solver.factorize(*(mesh->laplacianMatrix));

	if (solver.info() != Eigen::Success)
	{
		std::cout << "Error in uv unwrapping" << endl;
	}

	Eigen::VectorXd uVector = solver.solve(mesh->vectorX);

	Eigen::VectorXd vVector = solver.solve(mesh->vectorY);

	for (int i = 0; i < mesh->verts.size(); i++)
	{ 
		uVector[i] = uVector[i];
		vVector[i] = vVector[i];

		Vertex* v = mesh->verts[i];
		v->uvCoordinates.u = uVector[i];
		v->uvCoordinates.v = vVector[i];
	}
}

void computeBoundaryUVCoordinates(std::vector<int>& boundaryChain, std::vector<Vertex*>& vertices, Eigen::VectorXd& vectorX, Eigen::VectorXd& vectorY, 
	std::vector<glm::vec2>& arbitraryPolygonVertices, MappingType mappingType)
{
	switch (mappingType)
	{
	case MappingType::DISC:
		mapToDisc(boundaryChain, vertices, vectorX, vectorY);
		break;

	case MappingType::SQUARE:
		mapToSquare(boundaryChain, vertices, vectorX, vectorY);
		break;

	case MappingType::CONVEX_HULL:
		mapToConvexHull(boundaryChain, vertices, vectorX, vectorY);
		break;

	case MappingType::ARBITRARY:
		mapToArbitraryPolygon(boundaryChain, vertices, vectorX, vectorY, arbitraryPolygonVertices);
		break;

	default:
		break;
	}

}

void mapToDisc(std::vector< int >& boundaryChain, std::vector< Vertex* >& vertices, Eigen::VectorXd& vectorX, Eigen::VectorXd& vectorY)
{
	double angle = 2.0 * M_PI / boundaryChain.size();
	double angle_delta = angle;

	for (int i = 0; i < boundaryChain.size(); i++)
	{
		float x = 0.5 * (1.0 + cos(angle));
		float y = 0.5 * (1.0 + sin(angle));

		int v_i = boundaryChain[i];

		vertices[v_i]->uvCoordinates = { x, y };

		vectorX[v_i] = x;
		vectorY[v_i] = y;

		angle += angle_delta;
	}
}

void mapToSquare(std::vector< int >& boundaryChain, std::vector< Vertex* >& vertices, Eigen::VectorXd& vectorX, Eigen::VectorXd& vectorY)
{
	int top_right_corner = boundaryChain.size() / 4;
	int bottom_right_corner = boundaryChain.size() / 2;
	int bottom_left_corner = boundaryChain.size() * 3 / 4;

	float x = 0, y = 0;

	double delta_movement = 1 / (float)top_right_corner;

	for (int i = 0; i < boundaryChain.size(); i++)
	{
		int v_i = boundaryChain[i];

		vertices[v_i]->uvCoordinates = { x, y };

		vectorX[v_i] = x;
		vectorY[v_i] = y;

		if (i < top_right_corner)
		{
			x += delta_movement;
		}
		else if (i < bottom_right_corner)
		{
			y += delta_movement;
		}
		else if (i < bottom_left_corner)
		{
			x -= delta_movement;
		}
		else
		{
			y -= delta_movement;
		}
	}
}

void mapToConvexHull(std::vector< int >& boundaryChain, std::vector< Vertex* >& vertices, Eigen::VectorXd& vectorX, Eigen::VectorXd& vectorY)
{
	std::vector<int> convexHull = findConvexHull(boundaryChain, vertices);

	auto elem0 = std::find(boundaryChain.begin(), boundaryChain.end(), convexHull[0]);
	auto elem1 = std::find(boundaryChain.begin(), boundaryChain.end(), convexHull[1]);
	int direction = elem1 > elem0 ? 1 : -1;
	if (direction == -1)
	{
		std::reverse(convexHull.begin(), convexHull.end());
		direction = 1;
	}
	int hull_iterator = 0;
	int chain_start_index = elem0 - boundaryChain.begin();
	int chain_iterator = chain_start_index;

	while (hull_iterator <= convexHull.size())
	{
		if (convexHull[hull_iterator % convexHull.size()] != boundaryChain[chain_iterator])
		{
			int vertices_in_between = 1;
			while (boundaryChain[(chain_iterator + (vertices_in_between)*direction) % boundaryChain.size()] != convexHull[hull_iterator % convexHull.size()])
			{
				vertices_in_between++;
			}

			Vertex* start_v = vertices[convexHull[(hull_iterator - 1) % convexHull.size()]];
			Vertex* end_v = vertices[convexHull[hull_iterator % convexHull.size()]];

			float x = start_v->coords[0];
			float y = start_v->coords[1];
			float step_size_x = (end_v->coords[0] - start_v->coords[0]) / ((vertices_in_between + 1) * direction);
			float step_size_y = (end_v->coords[1] - start_v->coords[1]) / ((vertices_in_between + 1) * direction);

			for (int i = 0; i < vertices_in_between + 1; i++)
			{
				int v_i = boundaryChain[(chain_iterator + i * direction) % boundaryChain.size()];
				Vertex* middle_v = vertices[v_i];

				x += step_size_x;
				y += step_size_y;

				middle_v->uvCoordinates = { x, y };

			}
			chain_iterator = (chain_iterator + vertices_in_between * direction) % boundaryChain.size();
		}
		else
		{
			Vertex* v = vertices[boundaryChain[chain_iterator]];
			v->uvCoordinates = { v->coords[0], v->coords[1] };
		}

		hull_iterator++;
		chain_iterator = chain_iterator >= boundaryChain.size() - 1 ? 0 : chain_iterator + direction;

	}

	normalizeUVCoordinates(boundaryChain, vertices, vectorX, vectorY);
}

void mapToArbitraryPolygon(std::vector< int >& boundaryChain, std::vector< Vertex* >& vertices, Eigen::VectorXd& vectorX, Eigen::VectorXd& vectorY,
	std::vector<glm::vec2>& arbitraryPolygonVertices)
{	
	// calculate the edge length for the given arbitrary polygon
	std::vector<double> edge_lengths;

	for (int i = 0; i < arbitraryPolygonVertices.size(); i++)
	{
		glm::vec2 first_p = arbitraryPolygonVertices[i], second_p = arbitraryPolygonVertices[(i + 1) % arbitraryPolygonVertices.size()];

		double edge_len = sqrt((first_p.x - second_p.x) * (first_p.x - second_p.x) +
								(first_p.y - second_p.y) * (first_p.y - second_p.y));

		edge_lengths.push_back(edge_len);
	}

	// find the relative lengths (for finding how many points we will map to each edge)
	double sum = 0;
	for (int i = 0; i < edge_lengths.size(); i++) sum += edge_lengths[i];

	float sum_test = 0;
	for (int i = 0; i < edge_lengths.size(); i++) {
		edge_lengths[i] = edge_lengths[i] / sum * boundaryChain.size();
		sum_test += edge_lengths[i];
	}

	double accumulator = 0;
	for (int i = 0; i < edge_lengths.size(); i++)
	{
		glm::vec2 start_point = arbitraryPolygonVertices[i];
		glm::vec2 end_point = arbitraryPolygonVertices[(i + 1) % edge_lengths.size()];

		int vertices_to_map = ((int)(accumulator + edge_lengths[i])) - accumulator;

		glm::vec2 step = (end_point - start_point) / glm::vec2{vertices_to_map-1, vertices_to_map-1};

		for (int j = accumulator; j < accumulator + vertices_to_map; j++)
		{
			int v_i = boundaryChain[j];

			vertices[v_i]->uvCoordinates.u = start_point.x;
			vertices[v_i]->uvCoordinates.v = start_point.y;

			vectorX[v_i] = start_point.x;
			vectorY[v_i] = start_point.y;

			start_point += step;
		}
		accumulator += edge_lengths[i];
	}

}

void improveUnwrapping(Mesh* mesh)
{
	std::vector<float> tri_areas;
	std::vector<float> tri_sigmas;
	for (int i = 0; i < mesh->tris.size(); i++)
	{
		Triangle* triangle = mesh->tris[i];

		Vertex* v1 = mesh->verts[triangle->v1i];
		Vertex* v2 = mesh->verts[triangle->v2i];
		Vertex* v3 = mesh->verts[triangle->v3i];

		glm::vec3 q1 = { v1->coords[0], v1->coords[1], v1->coords[2] };
		glm::vec3 q2 = { v2->coords[0], v2->coords[1], v2->coords[2] };
		glm::vec3 q3 = { v3->coords[0], v3->coords[1], v3->coords[2] };

		UVCoordinates uv1 = v1->uvCoordinates;
		UVCoordinates uv2 = v2->uvCoordinates;
		UVCoordinates uv3 = v3->uvCoordinates;
		float area = mesh->calculateTriangleAreaUVMethod1(i);

		glm::vec3 s_s = (q1 * (uv2.v - uv3.v) + q2 * (uv3.v - uv1.v) + q3 * (uv1.v - uv2.v)) / (2 * area);
		glm::vec3 s_t = (q1 * (uv3.u - uv2.u) + q2 * (uv1.u - uv3.u) + q3 * (uv2.u - uv1.u)) / (2 * area);

		float a = glm::dot(s_s, s_s);
		float b = glm::dot(s_s, s_t);
		float c = glm::dot(s_t, s_t);

		float max_val = sqrt(0.5 * (a + c + sqrt(std::pow(a - c, 2) + 4 * b * b)));
		float min_val = sqrt(0.5 * (a + c - sqrt(std::pow(a - c, 2) + 4 * b * b)));

		float sigma = sqrt((a + c) / 2);

		tri_areas.push_back(area);
		tri_sigmas.push_back(sigma);
	}

	std::vector<float> stretches;
	for (int i = 0; i < mesh->verts.size(); i++)
	{
		Vertex* v = mesh->verts[i];

		std::vector<int>& tris = v->triList;

		float numerator = 0;
		float denominator = 0;
		for (int j = 0; j < tris.size(); j++)
		{
			numerator += tri_areas[tris[j]] * tri_sigmas[tris[j]] * tri_sigmas[tris[j]];
			denominator += tri_areas[tris[j]];
		}

		float stretch = sqrt(numerator / denominator);
		stretches.push_back(stretch);
	}

	for (int i = 0; i < mesh->verts.size(); i++)
	{
		Vertex* v = mesh->verts[i];
		if (v->isBoundary) continue;
		std::vector<int>& neighbors = v->vertList;

		float sum = 0;
		for (int j = 0; j < neighbors.size(); j++)
		{
			if(! std::isnan(mesh->laplacianMatrix->coeff(i, neighbors[j]) / stretches[neighbors[j]]))
				mesh->laplacianMatrix->coeffRef(i, neighbors[j]) /= stretches[neighbors[j]];
			sum += mesh->laplacianMatrix->coeff(i, neighbors[j]);
		}
		mesh->laplacianMatrix->coeffRef(i, i) = -sum;
	}

	solveLaplacianMatrix(mesh);
}

void normalizeUVCoordinates(std::vector< int >& boundaryChain, std::vector< Vertex* >& vertices, Eigen::VectorXd& vectorX, Eigen::VectorXd& vectorY)
{
	Vertex* first_v = vertices[boundaryChain[0]];

	float min_y = first_v->coords[1], max_y = first_v->coords[1];
	float min_x = first_v->coords[0], max_x = first_v->coords[0];

	for (int i = 0; i < boundaryChain.size(); i++)
	{
		UVCoordinates coords = vertices[boundaryChain[i]]->uvCoordinates;

		if (coords.u < min_x) min_x = coords.u;
		else if (coords.u > max_x) max_x = coords.u;

		if (coords.v < min_y) min_y = coords.v;
		else if (coords.v > max_y) max_y = coords.v;
	}

	float normalizer_x = 1 / (max_x - min_x);
	float normalizer_y = 1 / (max_y - min_y);

	for (int i = 0; i < boundaryChain.size(); i++)
	{
		UVCoordinates* coords = &(vertices[boundaryChain[i]]->uvCoordinates);

		coords->u = (coords->u - min_x) * normalizer_x;
		coords->v = (coords->v - min_y) * normalizer_y;

		vectorX[boundaryChain[i]] = coords->u;
		vectorY[boundaryChain[i]] = coords->v;

	}
}

double computeMeanValueWeight(std::vector<Vertex*>& vertices, int v1i, int v2i, std::pair<int, int> neighbor_indices)
{
	if (neighbor_indices.first == -1)
	{
		return computeMeanValueWeight(vertices, v1i, v2i, neighbor_indices.second);
	}

	else if (neighbor_indices.second == -1)
	{
		return computeMeanValueWeight(vertices, v1i, v2i, neighbor_indices.first);
	}

	Vertex* v1 = vertices[v1i];
	Vertex* v2 = vertices[v2i];
	std::pair<Vertex*, Vertex*> neighbor_vertices(vertices[neighbor_indices.first], vertices[neighbor_indices.second]);


	glm::vec3 vector1(v2->coords[0] - v1->coords[0], v2->coords[1] - v1->coords[1], v2->coords[2] - v1->coords[2]);
	glm::vec3 vector2(neighbor_vertices.first->coords[0] - v1->coords[0],
		neighbor_vertices.first->coords[1] - v1->coords[1],
		neighbor_vertices.first->coords[2] - v1->coords[2]);

	glm::vec3 vector3(neighbor_vertices.second->coords[0] - v1->coords[0],
		neighbor_vertices.second->coords[1] - v1->coords[1],
		neighbor_vertices.second->coords[2] - v1->coords[2]);

	double first_angle = std::acos(glm::dot(vector1, vector2) / (glm::length(vector1) * glm::length(vector2)));
	double first_angle_tan = std::atan(first_angle / 2);

	double second_angle = std::acos(glm::dot(vector1, vector3) / (glm::length(vector1) * glm::length(vector3)));
	double second_angle_tan = std::atan(second_angle / 2);

	return (first_angle_tan + second_angle_tan) / (2 * glm::length(vector1));
}

double computeMeanValueWeight(std::vector<Vertex*>& vertices, int v1i, int v2i, int neighbor_index)
{
	Vertex* v1 = vertices[v1i];
	Vertex* v2 = vertices[v2i];
	Vertex* neighbor_vertex = vertices[neighbor_index];

	glm::vec3 vector1(v2->coords[0] - v1->coords[0], v2->coords[1] - v1->coords[1], v2->coords[2] - v1->coords[2]);

	glm::vec3 vector2(neighbor_vertex->coords[0] - v1->coords[0],
		neighbor_vertex->coords[1] - v1->coords[1],
		neighbor_vertex->coords[2] - v1->coords[2]);

	double angle = std::acos(glm::dot(vector1, vector2) / (glm::length(vector1) * glm::length(vector2)));
	double angle_tan = std::atan(angle / 2);

	return angle_tan / (2 * glm::length(vector1));
}

double computeHarmonicValueWeight(std::vector<Vertex*>& vertices, int v1i, int v2i, std::pair<int, int> neighbor_indices)
{
	if (neighbor_indices.first == -1)
	{
		return computeHarmonicValueWeight(vertices, v1i, v2i, neighbor_indices.second);
	}

	else if (neighbor_indices.second == -1)
	{
		return computeHarmonicValueWeight(vertices, v1i, v2i, neighbor_indices.first);
	}

	Vertex* v1 = vertices[v1i];
	Vertex* v2 = vertices[v2i];
	std::pair<Vertex*, Vertex*> neighbor_vertices(vertices[neighbor_indices.first], vertices[neighbor_indices.second]);

	glm::vec3 vector_31(v1->coords[0] - neighbor_vertices.first->coords[0],
		v1->coords[1] - neighbor_vertices.first->coords[1],
		v1->coords[2] - neighbor_vertices.first->coords[2]);

	glm::vec3 vector_32(v2->coords[0] - neighbor_vertices.first->coords[0],
		v2->coords[1] - neighbor_vertices.first->coords[1],
		v2->coords[2] - neighbor_vertices.first->coords[2]);

	glm::vec3 vector_41(v1->coords[0] - neighbor_vertices.second->coords[0],
		v1->coords[1] - neighbor_vertices.second->coords[1],
		v1->coords[2] - neighbor_vertices.second->coords[2]);

	glm::vec3 vector_42(v2->coords[0] - neighbor_vertices.second->coords[0],
		v2->coords[1] - neighbor_vertices.second->coords[1],
		v2->coords[2] - neighbor_vertices.second->coords[2]);

	double first_angle = std::acos(glm::dot(vector_31, vector_32) / (glm::length(vector_31) * glm::length(vector_32)));
	double first_angle_cot = 1 / std::atan(first_angle);

	double second_angle = std::acos(glm::dot(vector_41, vector_42) / (glm::length(vector_41) * glm::length(vector_42)));
	double second_angle_cot = 1 / std::atan(second_angle);

	return (first_angle_cot + second_angle_cot) / 2;
}

double computeHarmonicValueWeight(std::vector<Vertex*>& vertices, int v1i, int v2i, int neighbor_index)
{
	Vertex* v1 = vertices[v1i];
	Vertex* v2 = vertices[v2i];
	Vertex* neighbor_vertex = vertices[neighbor_index];

	glm::vec3 vector1(v1->coords[0] - neighbor_vertex->coords[0],
		v1->coords[1] - neighbor_vertex->coords[1],
		v1->coords[2] - neighbor_vertex->coords[2]);

	glm::vec3 vector2(v2->coords[0] - neighbor_vertex->coords[0],
		v2->coords[1] - neighbor_vertex->coords[1],
		v2->coords[2] - neighbor_vertex->coords[2]);

	double angle = std::acos(glm::dot(vector1, vector2) / (glm::length(vector1) * glm::length(vector2)));
	double angle_cot = 1 / std::atan(angle);

	return angle_cot / 2;
}

std::vector<int> findConvexHull(std::vector<int> boundaryChain, std::vector<Vertex*>& vertices)
{

	// find the bottom-most points
	int lowest_point_i = boundaryChain[0];
	float min_y = vertices[boundaryChain[0]]->coords[1];
	for (int i = 0; i < boundaryChain.size(); i++)
	{

		int cur_i = boundaryChain[i];
		float cur_y = vertices[cur_i]->coords[1];

		if (cur_y < min_y)
		{
			lowest_point_i = cur_i;
			min_y = cur_y;
		}
	}

	// sort the chain based on the polar angle between the point and the bottom-most point
	std::sort(boundaryChain.begin(), boundaryChain.end(), [&lowest_point_i, &vertices](int v1i, int v2i) {
		return comparePolarAngle(vertices, lowest_point_i, v1i, v2i);
		});

	std::vector<int> convexHull = { boundaryChain[0], boundaryChain[1], boundaryChain[2] };

	int hull_size = convexHull.size();
	for (int i = 3; i < boundaryChain.size(); i++)
	{
		int current_point = boundaryChain[i];

		while (!(isCounterClockwise(vertices, convexHull[hull_size - 2], convexHull[hull_size - 1], current_point)) && hull_size >= 2)
		{
			convexHull.pop_back();
			hull_size--;
		}

		convexHull.push_back(current_point);
		hull_size++;
	}

	return convexHull;
}

double computePolarAngle(std::vector<Vertex*>& vertices, int v1i, int v2i)
{
	float x = vertices[v2i]->coords[0] - vertices[v1i]->coords[0];
	float y = vertices[v2i]->coords[1] - vertices[v1i]->coords[1];

	return std::atan2(y, x);
}

bool comparePolarAngle(std::vector<Vertex*>& vertices, int refVertex, int v1i, int v2i)
{
	return computePolarAngle(vertices, refVertex, v1i) < computePolarAngle(vertices, refVertex, v2i);
}

bool isCounterClockwise(std::vector<Vertex*>& vertices, int v1i, int v2i, int v3i)
{
	float* v1_coords = vertices[v1i]->coords;
	float* v2_coords = vertices[v2i]->coords;
	float* v3_coords = vertices[v3i]->coords;

	float first = (v3_coords[1] - v2_coords[1]) * (v2_coords[0] - v1_coords[0]);
	float second = (v2_coords[1] - v1_coords[1]) * (v3_coords[0] - v2_coords[0]);

	return first > second;
}