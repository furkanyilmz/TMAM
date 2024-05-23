#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in float IsBoundary;
layout (location = 3) in vec2 aTexCoords;

out float isBoundary;
out vec3 FragPos;
out vec3 Normal;
out vec2 TexCoords1;
out vec2 TexCoords2;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 viewMesh;
uniform mat4 model;
// Transformations for the first texture
uniform vec2 textureOffset;
uniform float textureRotationAngle;
uniform float textureScale;
uniform vec2 textureStretch;
// Transformations for the second texture
uniform vec2 textureOffset2;
uniform float textureRotationAngle2;
uniform float textureScale2;
uniform vec2 textureStretch2;

vec2 transformTexCoords(vec2 aTexCoords, vec2 textureOffset, float textureRotationAngle, float textureScale, vec2 textureStretch) {
    vec2 textureCenter = vec2(0.5, 0.5);
    vec2 scaledTexCoordsAtOrigin = (aTexCoords - textureCenter) * textureScale * textureStretch;
    float cosA = cos(textureRotationAngle);
    float sinA = sin(textureRotationAngle);
    vec2 rotatedTexCoords = vec2(
        scaledTexCoordsAtOrigin.x * cosA - scaledTexCoordsAtOrigin.y * sinA,
        scaledTexCoordsAtOrigin.x * sinA + scaledTexCoordsAtOrigin.y * cosA
    );
    return rotatedTexCoords + textureCenter + textureOffset;
}

void main() {
    FragPos = vec3(model * vec4(aPos, 1.0f));
    Normal = normalize(vec3(model * viewMesh * vec4(aNormal, 0.0)));
    gl_Position = projection * view * vec4(FragPos, 1.0f);
    isBoundary = IsBoundary;

    TexCoords1 = transformTexCoords(aTexCoords, textureOffset, textureRotationAngle, textureScale, textureStretch);
    TexCoords2 = transformTexCoords(aTexCoords, textureOffset2, textureRotationAngle2, textureScale2, textureStretch2);
}
