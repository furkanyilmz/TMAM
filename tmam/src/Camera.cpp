#include "Camera.h"

Camera::Camera(float camRadius, bool check3D, float camTheta, float camPhi,
	glm::vec3 camPos, glm::vec3 camTarget, glm::vec3 camWorldUp) :
	radius(camRadius),
	theta(camTheta),
	phi(camPhi),
	pos(camPos),
	target(camTarget),
	worldUp(camWorldUp),
	viewNeedsUpdate(true),
	initialRadius(camRadius),
	is3D(check3D),
	view(glm::mat4(1.0f)),
	rot2D(glm::mat4(1.0f))
{

	if (!is3D)
		target = glm::vec3(0.5, 0.5, 0.0);
}

void Camera::reset()
{
	this->radius = this->initialRadius;
	theta = 0.0;
	phi = M_PI_2;
	pos = glm::vec3(0.0, 0.0, 0.0);

	if(is3D)
		target = glm::vec3(0.0, 0.0, 0.0);
	else
		target = glm::vec3(0.5, 0.5, 0.0);

	worldUp = glm::vec3(0.0, 1.0, 0.0);
	viewNeedsUpdate = true;
}

void Camera::rotate(float dTheta, float dPhi) 
{
	if (is3D)
	{
		if (worldUp.y > 0.0)
			theta += dTheta;
		else theta -= dTheta;

		phi += dPhi;
		if (phi > 2 * M_PI)
			phi -= 2 * M_PI;
		else if (phi < -2 * M_PI)
			phi += 2 * M_PI;

		if (dPhi != 0) {
			if ((phi > 0 && phi < M_PI) || (phi < -M_PI && phi > -2 * M_PI))
				worldUp.y = 1.0;
			else worldUp.y = -1.0;
		}
	}
		
	
		else
		{

			glm::vec2 uvCenter(0.5, 0.5); 

			rot2D = glm::translate(rot2D, glm::vec3(uvCenter.x,uvCenter.y, 0.0f));

			float radian = 1.0f;

			// Perform the rotation
			rot2D = glm::rotate(rot2D, glm::radians(-radian), glm::vec3(0.0f, 0.0f, 1.0f));

			rot2D = glm::translate(rot2D, glm::vec3(-uvCenter.x, -uvCenter.y, 0.0f));

		}
	

	viewNeedsUpdate = true;
}

void Camera::zoom(float dist) {
	const float minRadius = 1.0;

	if (radius - dist > minRadius) {
		radius -= dist;

	}
	else {
		radius = minRadius;
	}
	viewNeedsUpdate = true;
}

void Camera::pan(float dx, float dy, float minX, float maxX, float minY, float maxY) {
	const float padding = 0.5;
	glm::vec3 look = glm::normalize(toCartesian());

	glm::vec3 right = glm::cross(worldUp, look);
	glm::vec3 up = glm::cross(right, look);

	glm::vec3 newTarget = target + (right * dx) + (up * dy);

	if (newTarget.x > minX+padding && newTarget.x < maxX-padding && newTarget.y > minY+0.8 && newTarget.y < maxY-padding) {
		target = newTarget;
		viewNeedsUpdate = true;
	}

}

glm::vec3 Camera::toCartesian() const {
	float x = radius * sin(phi) * sin(theta);
	float y = radius * cos(phi);
	float z = radius * sin(phi) * cos(theta);
	return glm::vec3(x, y, z);
}

 glm::vec3 Camera::position()
{
	if (viewNeedsUpdate) {
		pos = target + toCartesian();
	}

	return pos;
}


glm::mat4 Camera::orthoMatrix(float left, float right, float bottom, float top,
	 float near, float far) const
 {
	 return glm::ortho(left, right, bottom, top, near, far);
 }

 glm::mat4 Camera::projectionMatrix(float width, float height,
	 float near, float far) const
 {
	 return glm::perspective((float)M_PI_4, width / height, near, far);
 }

glm::mat4 Camera::viewMatrix()
 {
	 if (viewNeedsUpdate) {
		 view = glm::lookAt(position(), target, worldUp);
		 if (!is3D)
			 view *= rot2D;

		 viewNeedsUpdate = false;
		 
	 }
	 return view;
 }