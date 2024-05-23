#version 330 core
in vec3 FragPos;
in vec3 Normal;
in float isBoundary;
in vec2 TexCoords1;
in vec2 TexCoords2;
out vec4 FragColor;

uniform vec3 lightPos;
uniform vec3 lightColor;
uniform vec3 objectColor;
uniform vec3 boundaryColor;
uniform sampler2D textureSampler; // First texture sampler
uniform sampler2D textureSampler2; // Second texture sampler
uniform float blendFactor; // Blend factor between the two textures
uniform bool useTexture;
uniform bool isBoundaryColor;

void main() {
    vec3 lightDir = normalize(lightPos - FragPos);
    float ambientStrength = 0.15;
    vec3 ambient = ambientStrength * lightColor;
    float diff = max(dot(Normal, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    vec4 textureColor1 = texture(textureSampler, TexCoords1);
    vec4 textureColor2 = texture(textureSampler2, TexCoords2);

    vec4 blendedTextureColor = mix(textureColor1, textureColor2, blendFactor);

    vec3 finalColor = (diffuse + ambient) * objectColor;
    if (useTexture) finalColor *= blendedTextureColor.rgb;
    vec3 boundaryFinalColor = (diffuse + ambient) * boundaryColor;
    if(isBoundaryColor){
        if (isBoundary > 0.5) FragColor = vec4(boundaryColor, 1.0f);
        else FragColor = vec4(finalColor, 1.0f);
    }
    else FragColor = vec4(finalColor, 1.0f);
    
}
