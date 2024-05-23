#include "Mesh.h"
#include "Unwrap.h"


#pragma warning(disable : 4996) // Disable warning for deprecated functions
void Mesh::loadOff(const char* name)
{
    FILE* fPtr = fopen(name, "r");
    char str[334];

    fscanf(fPtr, "%s", str);

    int nVerts, nFaces, nEdges, i = 0;
    float x, y, z;

	fscanf(fPtr, "%d %d %d\n", &nVerts, &nFaces, &nEdges);
	while (i++ < nVerts)
	{
		fscanf(fPtr, "%f %f %f", &x, &y, &z);
		addVertex(x/12, y/12, z/12);
	}

    for (i = 0; i < nFaces; i++)
    {
        int nVertsInFace;
        fscanf(fPtr, "%d", &nVertsInFace);
        std::vector<int> faceVertices(nVertsInFace);
        for (int j = 0; j < nVertsInFace; j++)
        {
            fscanf(fPtr, "%d", &faceVertices[j]);
        }

        for (int j = 0; j < nVertsInFace - 2; j++)
        {
            addTriangle(faceVertices[0], faceVertices[j + 1], faceVertices[j + 2]);
        }
    }
	laplacianMatrix = new Eigen::SparseMatrix<double>(this->verts.size(), this->verts.size());
	vectorX = Eigen::VectorXd::Zero(this->verts.size());
	vectorY = Eigen::VectorXd::Zero(this->verts.size());
	findBoundaryVertices();
	computeVertexNormals();
	findLongestBoundaryChain();
}

void Mesh::loadObj(const char* name) {
	glm::vec3 minExtent = glm::vec3(0.0, 0.0, 0.0), maxExtent = glm::vec3(0.0, 0.0, 0.0);
	glm::vec3 tempVertex;
	vector<UVCoordinates> temporaryuv;
	FILE* file = fopen(name, "r");
	if (file == NULL) {
		fprintf(stderr, "Unable to open OBJ file\n");
		return;
	}
	char lineHeader[128];
	while (fscanf(file, "%s", lineHeader) != EOF) {
		if (strcmp(lineHeader, "vn") == 0) {
			float nx, ny, nz;
			fscanf(file, "%f %f %f\n", &nx, &ny, &nz);
			addNormal(nx, ny, nz);
		}
		else if (strcmp(lineHeader, "vt") == 0) {
			float u, v;
			fscanf(file, "%f %f\n", &u, &v);
			UVCoordinates coords = { u, v };
			temporaryuv.push_back(coords);
		}
		else if (strcmp(lineHeader, "v") == 0) {
			float x, y, z;
			fscanf(file, "%f %f %f\n", &x, &y, &z);

			tempVertex.x = x; tempVertex.y = y; tempVertex.z = z;
			minExtent = glm::min(minExtent, tempVertex);
			maxExtent = glm::max(maxExtent, tempVertex);
			addVertex(x, y, z);
		}
		else if (strcmp(lineHeader, "f") == 0) {
			std::vector<int> vertices;
			std::vector<int> bffuvvertices;
			char line[256];
			fgets(line, sizeof(line), file); // read the entire line

			std::istringstream iss(line);
			std::string vertexStr;
			while (iss >> vertexStr) {
				int vertexIndex, textureIndex, normalIndex, bffuvvertexIndex;
				if (vertexStr.find('/') != std::string::npos) {
					// format with slashes
					std::istringstream iss(vertexStr);
					std::string token;
					std::vector<std::string> tokens;
					while (std::getline(iss, token, '/')) {
						tokens.push_back(token);
					}
					vertexIndex = std::stoi(tokens[0]);
					if (tokens.size() > 1 && !tokens[1].empty()) {
						bffuvvertexIndex = std::stoi(tokens[1]);
						bffuvvertices.push_back(bffuvvertexIndex - 1);
					}
				} else {
					// format without slashes
					sscanf(vertexStr.c_str(), "%d", &vertexIndex);
				}
				vertices.push_back(vertexIndex - 1);
				// subtract 1 to convert to 0-based indexing
				this->verts[vertexIndex - 1]->bffuvCoordinates.u = temporaryuv[bffuvvertexIndex - 1].u;
				this->verts[vertexIndex - 1]->bffuvCoordinates.v = temporaryuv[bffuvvertexIndex - 1].v;
			}
			for (int i = 0; i < vertices.size() - 2; ++i) {
				addTriangle(vertices[0], vertices[i + 1], vertices[i + 2]);
			}
		}
	}

	float sizeX = maxExtent.x - minExtent.x;
	float sizeY = maxExtent.y - minExtent.y;
	float sizeZ = maxExtent.z - minExtent.z;

	// Calculate a scaling factor based on the largest dimension
	float maxSize = std::max(sizeX, std::max(sizeY, sizeZ));
	if (maxSize == 0.0) maxSize = 1.0;
	for (int vertexIndex = 0; vertexIndex < verts.size(); ++vertexIndex) {
		verts[vertexIndex]->coords[0] /= maxSize;
		verts[vertexIndex]->coords[1] /= maxSize;
		verts[vertexIndex]->coords[2] /= maxSize;
	}

	fclose(file);
	
	laplacianMatrix = new Eigen::SparseMatrix<double>(this->verts.size(), this->verts.size());
	vectorX = Eigen::VectorXd::Zero(this->verts.size());
	vectorY = Eigen::VectorXd::Zero(this->verts.size());
	findBoundaryVertices();
	computeVertexNormals();
	findLongestBoundaryChain();
	
}

void Mesh::addNormal(float x, float y, float z) {
	glm::vec3 normal(x, y, z);
	normals.push_back(normal);
}


void Mesh::createCube(float sideLen)
{
	//coordinates
	float flbc[3] = { 0, 0, 0 }, deltaX = 0, deltaY = 0, deltaZ = 0;
	for (int v = 0; v < 8; v++)
	{
		switch (v)
		{
		case 1:
			deltaX = sideLen;
			break;
		case 2:
			deltaZ = -sideLen;
			break;
		case 3:
			deltaX = 0;
			break;
		case 4:
			deltaZ = 0;
			deltaY = sideLen;
			break;
		case 5:
			deltaX = sideLen;
			break;
		case 6:
			deltaZ = -sideLen;
			break;
		default:
			deltaX = 0;;
			break;
		}
		addVertex(flbc[0] + deltaX, flbc[1] + deltaY, flbc[2] + deltaZ);
	}

	addTriangle(0, 2, 1);
	addTriangle(0, 3, 2);

	addTriangle(1, 2, 5);
	addTriangle(2, 6, 5);

	addTriangle(2, 3, 6);
	addTriangle(3, 7, 6);

	addTriangle(3, 4, 7);
	addTriangle(3, 0, 4);

	addTriangle(4, 5, 6);
	addTriangle(4, 6, 7);

	addTriangle(0, 1, 5);
	addTriangle(0, 5, 4);
}

void Mesh::addTriangle(int v1, int v2, int v3)
{
	int idx = tris.size();
	
	glm::vec3 normal = computeTriangleNormal(v1, v2, v3);

	tris.push_back(new Triangle(idx, v1, v2, v3, normal));

	//set up structure

	verts[v1]->triList.push_back(idx);
	verts[v2]->triList.push_back(idx);
	verts[v3]->triList.push_back(idx);

	if (!makeVertsNeighbor(v1, v2))
		addEdge(v1, v2);
	else {
		makeEdgeNonBoundary(v1, v2);
	}

	if (!makeVertsNeighbor(v1, v3))
		addEdge(v1, v3);
		
	else {
		makeEdgeNonBoundary(v1, v3);
	}

	if (!makeVertsNeighbor(v2, v3))
		addEdge(v2, v3);
	else {
		makeEdgeNonBoundary(v2, v3);
	}

}

bool Mesh::makeVertsNeighbor(int v1i, int v2i)
{
	//returns true if v1i already neighbor w/ v2i; false o/w

	for (int i = 0; i < verts[v1i]->vertList.size(); i++)
		if (verts[v1i]->vertList[i] == v2i)
			return true;


	verts[v1i]->vertList.push_back(v2i);
	verts[v2i]->vertList.push_back(v1i);
	return false;
}

void Mesh::addVertex(float x, float y, float z)
{
	int idx = verts.size();
	float* c = new float[3];
	c[0] = x;
	c[1] = y;
	c[2] = z;

	verts.push_back(new Vertex(idx, c));
}

void Mesh::addEdge(int v1, int v2)
{
	int idx = edges.size();

	Edge* edge = new Edge(idx, v1, v2);
	edge->length = computeVertexDistance(v1, v2);

	edges.push_back(edge);
	
	verts[v1]->edgeList.push_back(idx);
	verts[v2]->edgeList.push_back(idx);
}

float Mesh::computeVertexDistance(int v1i, int v2i)
{
	float* v1_coords = verts[v1i]->coords;
	float* v2_coords = verts[v2i]->coords;

	float d_x = v1_coords[0] - v2_coords[0];
	float d_y = v1_coords[1] - v2_coords[1];
	float d_z = v1_coords[2] - v2_coords[2];

	return std::sqrt(d_x * d_x + d_y * d_y + d_z * d_z);
}

void Mesh::makeEdgeNonBoundary(int v1i, int v2i)
{
	vector< int > edge_indices = verts[v1i]->edgeList;

	for (int i = 0; i < edge_indices.size(); i++)
	{
		Edge* edge = edges[edge_indices[i]]; 

		if (edge->v1i == v2i || edge->v2i == v2i)
		{
			edge->isBoundary = false;
			verts[v1i]->isBoundary = false;
			verts[v2i]->isBoundary = false;
			break;
		}
	}
}

void Mesh::findBoundaryVertices()
{
	for (int i = 0; i < this->edges.size(); i++)
	{
		Edge* edge = this->edges[i];

		if (edge->isBoundary)
		{
			this->verts[edge->v1i]->isBoundary = 1;
			this->verts[edge->v2i]->isBoundary = 1;
		}
	}
}

glm::vec3 Mesh::computeTriangleNormal(int v1i, int v2i, int v3i)
{
	Vertex* v1 = verts[v1i];
	Vertex* v2 = verts[v2i];
	Vertex* v3 = verts[v3i];

	glm::vec3 first(v1->coords[0] - v2->coords[0], v1->coords[1] - v2->coords[1], v1->coords[2] - v2->coords[2]);
	glm::vec3 second(v1->coords[0] - v3->coords[0], v1->coords[1] - v3->coords[1], v1->coords[2] - v3->coords[2]);

	glm::vec3 normal = glm::normalize(glm::cross(first, second));

	return normal;
}

void Mesh::computeVertexNormals()
{
	for (int i = 0; i < verts.size(); i++)
	{
		Vertex* vertex = verts[i];

      	// If the vertex's normal is the zero vector, skip the computation for this vertex
        if (vertex->normal != glm::vec3(0.0f, 0.0f, 0.0f))
            continue;
		
		glm::vec3 avgNormal(0.0f);

		for (int j = 0; j < vertex->triList.size(); j++)
		{
			avgNormal += tris[vertex->triList[j]]->normal;
		}

		if (vertex->triList.size())
		{
			avgNormal /= static_cast<float>(vertex->triList.size());
			avgNormal = glm::normalize(avgNormal);
		}

		vertex->normal = avgNormal;
	}
}

void Mesh::findLongestBoundaryChain()
{
	int longest_chain = this->boundaryChain.size();


	for (int i = 0; i < verts.size(); i++)
	{
		Vertex* vertex_begin = verts[i];

		if (vertex_begin->isBoundary && !vertex_begin->visited)
		{
			std::vector<int> boundaryChain = findBoundaryChain(vertex_begin);

			if (boundaryChain.size() > this->boundaryChain.size()) {
				
				this->boundaryChain = boundaryChain;
			}
		}
	}
}

std::vector<int> Mesh::findBoundaryChain(Vertex* v)
{
	if (v->isBoundary)
	{
		v->visited = 1;

		std::vector<int> boundaryChain;
		boundaryChain.push_back(v->idx);

		Vertex* v_next = v;
		Vertex* v_prev = v;

		int count = 0;
		do
		{
			for (int i = 0; i < v_next->edgeList.size(); i++)
			{
				Edge* current_edge = edges[v_next->edgeList[i]];

				if (current_edge->isBoundary && ((current_edge->v1i != v_prev->idx && current_edge->v2i != v_prev ->idx) || !count))
				{
					Vertex* current_neighbor = current_edge->v1i == v_next->idx ? verts[current_edge->v2i] : verts[current_edge->v1i];

					if (!current_neighbor->visited || current_neighbor == v)
					{
						boundaryChain.push_back(current_neighbor->idx);

						v_prev = v_next;
						v_next = current_neighbor;
						v_next->visited = true;
						
						
						count++;
						break;
					}
					current_neighbor->visited = true;
				}
				
			}
		} while (v_next != v);

		boundaryChain.pop_back();
		return boundaryChain;
	}
}

std::pair<int, int> Mesh::findMiddleVertex(int v1i, int v2i)
{
	Vertex* v1 = verts[v1i];
	Vertex* v2 = verts[v2i];

	std::pair<int, int> middleVertexIndex = { -1, -1 };
	//TODO: use pair for multiple values (there can be 0, 1, 2)

	int count = 0;
	for (int i = 0; i < v1->triList.size(); i++)
	{
		Triangle* triangle = tris[v1->triList[i]];

		if (triangle->v1i == v2i)
		{
			if(count == 0)
				middleVertexIndex.first = triangle->v2i == v1i ? triangle->v3i : triangle->v2i;
			else 
				middleVertexIndex.second = triangle->v2i == v1i ? triangle->v3i : triangle->v2i;

			count++;
		}
		else if (triangle->v2i == v2i)
		{
			if(count == 0)
				middleVertexIndex.first = triangle->v1i == v1i ? triangle->v3i : triangle->v1i;
			else
				middleVertexIndex.second = triangle->v1i == v1i ? triangle->v3i : triangle->v1i;

			count++;
		}
		else if (triangle->v3i == v2i)
		{
			if(count == 0)
				middleVertexIndex.first = triangle->v1i == v1i ? triangle->v2i : triangle->v1i;
			else
				middleVertexIndex.second = triangle->v1i == v1i ? triangle->v2i : triangle->v1i;

			count++;
		}
	}

	return middleVertexIndex;
}

float Mesh::calculateTriangleAreaUVMethod1(int triIndex) {
	const Triangle* tri = tris[triIndex];
	const UVCoordinates& uv1 = verts[tri->v1i]->uvCoordinates;
	const UVCoordinates& uv2 = verts[tri->v2i]->uvCoordinates;
	const UVCoordinates& uv3 = verts[tri->v3i]->uvCoordinates;

	float area = 0.5f * abs((uv1.u * (uv2.v - uv3.v) + uv2.u * (uv3.v - uv1.v) + uv3.u * (uv1.v - uv2.v)));
	return area;
}

float Mesh::calculateTriangleAreaUVMethod2(int triIndex) {
	const Triangle* tri = tris[triIndex];
	const UVCoordinates& uv1 = verts[tri->v1i]->bffuvCoordinates;
	const UVCoordinates& uv2 = verts[tri->v2i]->bffuvCoordinates;
	const UVCoordinates& uv3 = verts[tri->v3i]->bffuvCoordinates;

	float area = 0.5f * abs((uv1.u * (uv2.v - uv3.v) + uv2.u * (uv3.v - uv1.v) + uv3.u * (uv1.v - uv2.v)));
	return area;
}

float Mesh::calculateTriangleArea3D(int triIndex) {
	const Triangle* tri = tris[triIndex];
	const float* v1 = verts[tri->v1i]->coords;
	const float* v2 = verts[tri->v2i]->coords;
	const float* v3 = verts[tri->v3i]->coords;

	glm::vec3 edge1 = glm::vec3(v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]);
	glm::vec3 edge2 = glm::vec3(v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]);
	glm::vec3 crossProduct = glm::cross(edge1, edge2);

	float area = 0.5f * glm::length(crossProduct);
	return area;
}

DistortionResult Mesh::calculateTriangleAreaDistortion() {
	float totalDistortionMethod1 = 0.0f;
	float totalDistortionMethod2 = 0.0f;

	for (int i = 0; i < tris.size(); ++i) {
		float area3D = calculateTriangleArea3D(i);
		float areaUV1 = calculateTriangleAreaUVMethod1(i);
		float areaUV2 = calculateTriangleAreaUVMethod2(i);

		float distortionMethod1 = abs(area3D - areaUV1) / area3D;
		float distortionMethod2 = abs(area3D - areaUV2) / area3D;

		totalDistortionMethod1 += distortionMethod1;
		totalDistortionMethod2 += distortionMethod2;
	}

	DistortionResult result;
	result.averageDistortionMethod1 = totalDistortionMethod1 / tris.size();
	result.averageDistortionMethod2 = totalDistortionMethod2 / tris.size();

	return result;
}

float Mesh::angleBetweenVectors(const glm::vec3& v1, const glm::vec3& v2) {
	float dot = glm::dot(v1, v2);
	float lenSq1 = glm::dot(v1, v1);
	float lenSq2 = glm::dot(v2, v2);
	float x = acos(dot / sqrt(lenSq1 * lenSq2));
	if (isnan(x)) {
		return 0;
	}
	else return x;
}

AngleDistortionResult Mesh::calculateAngleBasedDistortion() {
	float totalDistortionUV = 0.0f;
	float totalDistortionBFFUV = 0.0f;

	for (auto& tri : tris) {
		// Get 3D and UV vertices
		glm::vec3 v1_3D(verts[tri->v1i]->coords[0], verts[tri->v1i]->coords[1], verts[tri->v1i]->coords[2]);
		glm::vec3 v2_3D(verts[tri->v2i]->coords[0], verts[tri->v2i]->coords[1], verts[tri->v2i]->coords[2]);
		glm::vec3 v3_3D(verts[tri->v3i]->coords[0], verts[tri->v3i]->coords[1], verts[tri->v3i]->coords[2]);

		glm::vec2 v1_UV(verts[tri->v1i]->uvCoordinates.u, verts[tri->v1i]->uvCoordinates.v);
		glm::vec2 v2_UV(verts[tri->v2i]->uvCoordinates.u, verts[tri->v2i]->uvCoordinates.v);
		glm::vec2 v3_UV(verts[tri->v3i]->uvCoordinates.u, verts[tri->v3i]->uvCoordinates.v);

		// Calculate angles in 3D space
		float angle1_3D = angleBetweenVectors(v2_3D - v1_3D, v3_3D - v1_3D);
		float angle2_3D = angleBetweenVectors(v1_3D - v2_3D, v3_3D - v2_3D);
		float angle3_3D = angleBetweenVectors(v1_3D - v3_3D, v2_3D - v3_3D);

		// Calculate angles in UV space
		float angle1_UV = angleBetweenVectors(glm::vec3(v2_UV - v1_UV, 0.0f), glm::vec3(v3_UV - v1_UV, 0.0f));
		float angle2_UV = angleBetweenVectors(glm::vec3(v1_UV - v2_UV, 0.0f), glm::vec3(v3_UV - v2_UV, 0.0f));
		float angle3_UV = angleBetweenVectors(glm::vec3(v1_UV - v3_UV, 0.0f), glm::vec3(v2_UV - v3_UV, 0.0f));

		// Calculate distortion for each angle and sum them up
		float distortion = abs(angle1_3D - angle1_UV) + abs(angle2_3D - angle2_UV) + abs(angle3_3D - angle3_UV);
		totalDistortionUV += distortion;


		glm::vec2 v1_BFFUV(verts[tri->v1i]->bffuvCoordinates.u, verts[tri->v1i]->bffuvCoordinates.v);
		glm::vec2 v2_BFFUV(verts[tri->v2i]->bffuvCoordinates.u, verts[tri->v2i]->bffuvCoordinates.v);
		glm::vec2 v3_BFFUV(verts[tri->v3i]->bffuvCoordinates.u, verts[tri->v3i]->bffuvCoordinates.v);

		// Calculate angles in BFFUV space
		float angle1_BFFUV = angleBetweenVectors(glm::vec3(v2_BFFUV - v1_BFFUV, 0.0f), glm::vec3(v3_BFFUV - v1_BFFUV, 0.0f));
		float angle2_BFFUV = angleBetweenVectors(glm::vec3(v1_BFFUV - v2_BFFUV, 0.0f), glm::vec3(v3_BFFUV - v2_BFFUV, 0.0f));
		float angle3_BFFUV = angleBetweenVectors(glm::vec3(v1_BFFUV - v3_BFFUV, 0.0f), glm::vec3(v2_BFFUV - v3_BFFUV, 0.0f));

		// Calculate distortion for each angle and sum them up for BFFUV
		float distortionBFFUV = abs(angle1_3D - angle1_BFFUV) + abs(angle2_3D - angle2_BFFUV) + abs(angle3_3D - angle3_BFFUV);
		totalDistortionBFFUV += distortionBFFUV;

	}

	AngleDistortionResult result;
	result.averageDistortionUV = totalDistortionUV / (3.0f * tris.size()); // Average per angle for UV
	result.averageDistortionBFFUV = totalDistortionBFFUV / (3.0f * tris.size()); // Average per angle for BFFUV

	return result;

}


void Mesh::exportMesh(std::string& directory, const bool isTextureLoaded, bool useBffForTexture) {
	// exmaple: directory = "Y:\tmam_new\tmam\tmam\textures\test"
	// example directory2 = "Y:\tmam_new\tmam\tmam\textures\test.obj"

	// Find last backslash and get the filename
	size_t lastBackslash = directory.find_last_of("\\");
	std::string fileName = directory.substr(lastBackslash + 1);

	// If filename dont have .obj extension, add it to directory
	if (fileName.find(".obj") == std::string::npos) {
		directory = directory + ".obj";
		std::cerr << "Warning: Filename does not have .obj extension. Adding it to the save directory." << std::endl;
		std::cerr << "New directory: " << directory << std::endl;
	}

	// If filename has obj extension, remove it
	if (fileName.find(".obj") != std::string::npos) {
		fileName = fileName.substr(0, fileName.size() - 4);
		std::cerr << "Warning: Filename has .obj extension. Removing it from the filename." << std::endl;
		std::cerr << "New filename: " << fileName << std::endl;
	}

	std::ofstream objFile(directory);

	if (!objFile.is_open()) {
		std::cerr << "Error opening file: " << fileName << std::endl;
		return;
	}

    objFile << "mtllib ./" << fileName << ".obj.mtl" << std::endl;

	for (const auto& vertex : verts) {
		objFile << "v " << vertex->coords[0] << " " << vertex->coords[1] << " " << vertex->coords[2] << std::endl;
	}

	for (const auto& normal : normals) {
		objFile << "vn " << normal.x << " " << normal.y << " " << normal.z << std::endl;
	}
	if (useBffForTexture) {
	for (const auto& vertex : verts) {
        objFile << "vt " << vertex->bffuvCoordinates.u << " " << vertex->bffuvCoordinates.v << std::endl;
	}
	}
	else {
		for (const auto& vertex : verts) {
			objFile << "vt " << vertex->uvCoordinates.u << " " << vertex->uvCoordinates.v << std::endl;
		}
	}

    objFile << "usemtl material_0" << std::endl;

	for (const auto& triangle : tris) {
        objFile << "f "
                << triangle->v1i + 1 << "/" << triangle->v1i + 1 << " "
                << triangle->v2i + 1 << "/" << triangle->v2i + 1 << " "
                << triangle->v3i + 1 << "/" << triangle->v3i + 1 << std::endl;
	}

	objFile.close();

    if (isTextureLoaded) {
		std::ofstream mtlFile(directory + ".mtl");

		if (!mtlFile.is_open()) {
			std::cerr << "Error opening file: " << fileName << std::endl;
			return;
		}

		mtlFile << "newmtl material_0" << std::endl;
		mtlFile << "Ka 0.200000 0.200000 0.200000\n"
                << "Kd 0.752941 0.752941 0.752941\n"
                << "Ks 1.000000 1.000000 1.000000\n"
                << "Tr 0.000000\n"
                << "illum 2\n"
                << "Ns 0.000000" << std::endl;
        mtlFile << "map_Kd " << fileName << "_texture.jpg" << std::endl;
		mtlFile.close();
	}
}
#pragma warning(default : 4996) // Restore the warning level

