#if defined(NANOGUI_GLAD)
#if defined(NANOGUI_SHARED) && !defined(GLAD_GLAPI_EXPORT)
#define GLAD_GLAPI_EXPORT
#endif

#include <glad/glad.h>
#else
#if defined(__APPLE__)
#define GLFW_INCLUDE_GLCOREARB
#else
#define GL_GLEXT_PROTOTYPES
#endif
#endif

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#include "stb/stb_image.h"

#include <GLFW/glfw3.h>

#include <nanogui/nanogui.h>
#include <iostream>
#include "Mesh.h"
#include "Camera.h"
#include "Unwrap.h"
#include "ViewArea.h"
#include "Utils.h"
#include "Preprocessing.h"

#include <fstream>
#include <iostream>
#include <filesystem>
#include <iomanip>

#include <opencv2/opencv.hpp>
#include <cmath>

#define eigen_assert(X) do { if(!(X)) throw std::runtime_error(#X); } while(false);
namespace fs = std::filesystem;

unsigned int shaderProgram;
Mesh* mesh = new Mesh();

//names for import and export
std::string loadedObjName;
std::string textureFileName;

std::string mainpath = fs::current_path().string();

glm::vec3 objectColor(0.7f, 0.7f, 0.7f);
glm::vec3 boundaryColor(1.0f, 0.0f, 0.0f);
bool isBoundaryColor = false;

nanogui::Screen* screen = nullptr;
nanogui::ref<nanogui::Window> linearDistWindow = nullptr;
nanogui::ref<nanogui::Window> nonLinearDistWindow = nullptr;
nanogui::Label* linearAreaDist = nullptr;
nanogui::Label* linearAngleDist = nullptr;
nanogui::Button* toggleUVsButton = nullptr;
nanogui::TextBox* rotateBox = nullptr;
nanogui::TextBox* scaleBox = nullptr;
nanogui::TextBox* offset_xBox = nullptr;
nanogui::TextBox* offset_yBox = nullptr;
nanogui::TextBox* stretch_xBox = nullptr;
nanogui::TextBox* stretch_yBox = nullptr;
nanogui::TextBox* rotateBox2 = nullptr;
nanogui::TextBox* scaleBox2 = nullptr;
nanogui::TextBox* offset_xBox2 = nullptr;
nanogui::TextBox* offset_yBox2 = nullptr;
nanogui::TextBox* stretch_xBox2 = nullptr;
nanogui::TextBox* stretch_yBox2 = nullptr;
ViewArea emptyArea; ViewArea meshArea;
ViewArea uvArea; ViewArea bffuvArea;
ViewArea fullArea;

glm::vec3 mouse = glm::vec3(0.0, 0.0, 0.0);

std::vector<Camera> cameras;

//boundaries for camera panning
glm::vec3 minExtent, maxExtent;
glm::vec3 minExtentUV, maxExtentUV;

//check if the mouse has been pressed
bool mouseDown = false;

//check if alt has been pressed
bool altDown = false;

//check if there is a model currently loaded
bool meshLoaded = false;

// Transformations for the first texture
glm::vec2 texOffset = glm::vec2(0.0, 0.0);
float textureRotationAngle = 0.0;
float textureScale = 1.0;
float textureMoveSpeed = 0.01;
glm::vec2 textureStretch = glm::vec2(1.0, 1.0);

// Transformations for the second texture
glm::vec2 texOffset2 = glm::vec2(0.0, 0.0);
float textureRotationAngle2 = 0.0;
float textureScale2 = 1.0;
float textureMoveSpeed2 = 0.01;
glm::vec2 textureStretch2 = glm::vec2(1.0, 1.0);

std::vector<float> vertices;
std::vector<unsigned int> indices;

unsigned int VBO, VAO, EBO;

// UV
std::vector<float> verticesUV;
std::vector<unsigned int> indicesUV;
std::vector<float> verticesbffUV;
std::vector<unsigned int> indicesbffUV;

unsigned int VBO_UV, VAO_UV, EBO_UV;
unsigned int bVBO_UV, bVAO_UV, bEBO_UV;

LinearUnwrappingType linearUnwrappingType = UNIFORM_WEIGHTS;
MappingType mappingType = DISC;

std::vector<glm::vec2> userPoints; 
int isMousePressed = 0;
int isWarningMessage = 0;
int isPointAdded = 0;
int isDrawn = 0;
nanogui::Screen* popupGui = nullptr;

//globals regarding texture
bool isTexture1 = true;
unsigned int textureID;
unsigned int textureID2;
float blendFactor = 0.0f;
bool useTexture = false;
bool useBffForTexture = false;
std::string textureName1 = "";
std::string textureName2 = "";

bool enabled = true;
int captureFlag = 0;


unsigned int initializeShaders(const char* vertexShaderPath, const char* fragmentShaderPath);
void initCameras(std::vector<Camera>& cameras);
void initGui(GLFWwindow* window);
nanogui::Widget* createLabelWidget(nanogui::ref<nanogui::Window> nanoguiWindow, string message, int size);
nanogui::Label* createTooltipLabel(nanogui::Widget* parent, string message, int labelWidth);
void setupAllMeshes(Mesh* mesh);
void setupMeshesExceptBff(Mesh* mesh);
void setupMeshesExceptBffWithoutComputingUVs(Mesh* mesh);
void setupMesh(Mesh* mesh, std::vector<float>& vertices, std::vector<unsigned int>& indices, unsigned int& VAO, unsigned int& VBO, unsigned int& EBO);
void setupUVMesh(Mesh* mesh, std::vector<float>& vertices, std::vector<unsigned int>& indices, unsigned int& VAO, unsigned int& VBO, unsigned int& EBO);
void setupbffUVMesh(Mesh* mesh, std::vector<float>& vertices, std::vector<unsigned int>& indices, unsigned int& VAO, unsigned int& VBO, unsigned int& EBO);

void loadTexture(std::string filename);
void loadTexture2(std::string filename);

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void exportImage(const ViewArea& viewArea);
void cursorPosCallback(double x, double y, glm::vec3 minExtent, glm::vec3 maxExtent, glm::vec3 minExtentUV, glm::vec3 maxExtentUV);
void mouseButtonCallback(GLFWwindow* window, int button, int action, int modifiers);
void keyButtonCallback(int key, int scancode, int action, int mods);
void scrollCallback(double x, double y);
void distGuiCallback();

void updateViewArea(int width, int height);
void clearViewArea(ViewArea& area);

void render(GLFWwindow* window, unsigned int& shaderProgram);
void renderMesh(ViewArea& area, Camera& camera, unsigned int& shaderProgram, unsigned int& VAO);
std::string readShaderFile(const std::string& shaderPath);

void createPolygonPopUp();
void drawPolygon(GLFWwindow* window, nanogui::Screen* popupScreen, int& isMousePressed, int& isWarningMessage);
bool isConvex(float x, float y);
void polygonMouseButtonCallback(GLFWwindow* popupWindow, int button, int action, int modifiers);
void initPopupGui(nanogui::Screen* popupGui);

cv::Mat transformTexture(const cv::Mat& image, const cv::Vec2f& offset, float rotationAngle, float scale);
cv::Mat loadAndProcessTexture(const std::string& filename, const cv::Vec2f& textureOffset, float textureRotationAngle, float textureScale);
void blendAndSaveTextures(std::string& directory);

const char* vertexShaderPath = "src/shaders/vertexShader.glsl";
const char* fragmentShaderPath = "src/shaders/fragmentShader.glsl";
const char* polygonVertexShaderPath = "src/shaders/polygonVertexShader.glsl";
const char* polygonFragmentShaderPath = "src/shaders/polygonFragmentShader.glsl";

int main(int /* argc */, char** /* argv */) {

    int width = 1600, height = 800;

    glfwInit();

    glfwSetTime(0);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    glfwWindowHint(GLFW_SAMPLES, 0);
    glfwWindowHint(GLFW_RED_BITS, 8);
    glfwWindowHint(GLFW_GREEN_BITS, 8);
    glfwWindowHint(GLFW_BLUE_BITS, 8);
    glfwWindowHint(GLFW_ALPHA_BITS, 8);
    glfwWindowHint(GLFW_STENCIL_BITS, 8);
    glfwWindowHint(GLFW_DEPTH_BITS, 24);
    glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

    GLFWwindow* window = glfwCreateWindow(width, height, "TMAM", nullptr, nullptr);
    if (window == nullptr) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

#if defined(NANOGUI_GLAD)
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
        throw std::runtime_error("Could not initialize GLAD!");
    glGetError();
#endif

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    updateViewArea(width, height);
    clearViewArea(fullArea); clearViewArea(meshArea); clearViewArea(uvArea); clearViewArea(bffuvArea);

    std::string name = "inputs/facem.off";
    loadedObjName = name.substr(0, name.length()-4);

    if (name.size() > 4 && name.compare(name.size() - 4, 4, ".off") == 0) {
        Preprocessing converter;
        if (converter.CheckifOffFileisWatertightandConvert(name, "output.obj")) {
            std::cout << "Conversion completed successfully.\n";
            name = "output.obj";
        }
        else {
            std::cerr << "Conversion failed.\n";
        }
    }

    if (name.size() > 4 && name.compare(name.size() - 4, 4, ".obj") == 0) {
        std::string command = mainpath + "\\bff\\bff.exe " + name + " outputBff.obj --normalizeUVs";

        system(command.c_str());
        mesh->loadObj("outputBff.obj");
        meshLoaded = true;
        remove((mainpath + "\\output.obj").c_str());
        remove((mainpath + "\\outputBff.obj").c_str());
    }

    setupAllMeshes(mesh);

    initCameras(cameras);

    screen = new nanogui::Screen();
    screen->initialize(window, true);

    glfwSetCursorPosCallback(window, [](GLFWwindow* window, double x, double y) {
        cursorPosCallback(x, y, minExtent, maxExtent, minExtentUV, maxExtentUV);
        });

    glfwSetMouseButtonCallback(window, [](GLFWwindow* window, int button, int action, int modifiers) {
        mouseButtonCallback(window, button, action, modifiers);
        });

    glfwSetKeyCallback(window, [](GLFWwindow* window, int key, int scancode, int action, int mods) {
        keyButtonCallback(key, scancode, action, mods);
        });

    glfwSetCharCallback(window, [](GLFWwindow* window, unsigned int codepoint) {
        screen->charCallbackEvent(codepoint);
        });

    glfwSetDropCallback(window, [](GLFWwindow* window, int count, const char** filenames) {
        screen->dropCallbackEvent(count, filenames);
        });

    glfwSetScrollCallback(window, [](GLFWwindow* window, double x, double y) {
        screen->scrollCallbackEvent(x, y);
        scrollCallback(x, y);
        });

    glfwSetFramebufferSizeCallback(window, [](GLFWwindow* window, int width, int height) {
        screen->resizeCallbackEvent(width, height);
        linearDistWindow->setPosition(nanogui::Vector2i(width * 0.5, height * 0.85));
        nonLinearDistWindow->setPosition(nanogui::Vector2i(width * 0.8, height * 0.85));
        glfwSetWindowSize(window, width, height);
        updateViewArea(width, height);
        });


    glfwGetFramebufferSize(window, &width, &height);
    glViewport(0, 0, width, height);
    glfwSwapInterval(0);
    glfwSwapBuffers(window);

    initGui(window);

    shaderProgram = initializeShaders(vertexShaderPath, fragmentShaderPath);

    render(window, shaderProgram);

    glfwTerminate();

    return 0;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, true);
    }
}

unsigned int initializeShaders(const char* vertexPath, const char* fragmentPath)
{
    std::string vertexCode = readShaderFile(vertexPath);
    std::string fragmentCode = readShaderFile(fragmentPath);
    const char* vertexShaderSource = vertexCode.c_str();
    const char* fragmentShaderSource = fragmentCode.c_str();

    unsigned int vertexShader;
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    int  success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);

    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    }

    unsigned int fragmentShader;
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    unsigned int shaderProgram;
    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return shaderProgram;
}


void initCameras(std::vector<Camera>& cameras){
    cameras.push_back(Camera(3.5, true));
    cameras.push_back(Camera(3.5, false));
    cameras.push_back(Camera(3.5, false));
}

void initGui(GLFWwindow* window) {
    
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    nanogui::FormHelper* gui = new nanogui::FormHelper(screen);
    nanogui::ref<nanogui::Window> nanoguiWindow = gui->addWindow(Eigen::Vector2i(10, 10), "");
    nanoguiWindow->setLayout(new nanogui::GroupLayout());

    gui->addButton("Load Mesh", [&]() {
        std::string source = nanogui::file_dialog({ {"off", "OFF Mesh File"}, {"obj", "OBJ Mesh File"}, {"ply2", "PLY2 Mesh File"} }, false);
        if (!source.empty()) {

            try
            {
                std::string filename = source.substr(source.find_last_of("/\\") + 1);
                loadedObjName = filename.substr(0, filename.length() - 4);

                if (filename.substr(filename.length() - 4) == "ply2")
                {
                    std::ifstream inp(source);
                    std::ofstream outp(mainpath + "\\" + "test.off");

                    outp << "OFF" << std::endl;
                    std::string comp_line;
                    std::string line;

                    for (int i = 0; i < 2; i++)
                    {
                        std::getline(inp, line);
                        comp_line.append(line + " ");
                    }
                    comp_line.append("10");
                    outp << comp_line << endl;

                    while (std::getline(inp, line))
                    {
                        outp << line << endl;
                    }
                    source = mainpath + "\\ply2tooff.off";
                }

                try {
                    fs::copy(source, mainpath + "\\" + filename, fs::copy_options::overwrite_existing);
                    std::cout << "File copied successfully to " << mainpath + "\\" + filename << std::endl;
                }
                catch (const fs::filesystem_error& e) {
                    std::cerr << e.what() << std::endl;
                }

                Preprocessing converter;
                if (!converter.Cut(mainpath + "\\" + filename, mainpath + "\\output.obj")) {
                    cout << "Error when cutting\n";
                }
                std::string command = mainpath + "\\bff\\bff.exe " + mainpath + "\\output.obj " + mainpath + "\\" + "outputBff.obj" + " --flattenToDisk" + " --normalizeUVs";
                Mesh* tempMesh = new Mesh();

                system(command.c_str());
                tempMesh->loadObj((mainpath + "\\outputBff.obj").c_str());
                remove((mainpath + "\\output.obj").c_str());
                remove((mainpath + "\\outputBff.obj").c_str());

                delete(mesh);
                mesh = tempMesh;
                vertices.clear();
                indices.clear();
                verticesUV.clear();
                indicesUV.clear();
                verticesbffUV.clear();
                indicesbffUV.clear();
                setupAllMeshes(mesh);
            }
            catch (const std::runtime_error& e)
            {
                 nanogui::MessageDialog* message = new nanogui::MessageDialog(screen, nanogui::MessageDialog::Type::Warning, "Importing Error",
                            "An error occurred while importing the selected file. Please ensure that the file is not corrupted.");
            }
        }
        });

    nanogui::Button* exportMeshButton = gui->addButton("Export 3D mesh", [&, width, height]() {
        captureFlag = 1;
        });

    exportMeshButton->setCallback([]() { captureFlag = 1; });

    nanogui::Widget* colorWheelWidget = new nanogui::Widget(nanoguiWindow);
    colorWheelWidget->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical,
        nanogui::Alignment::Middle, 0, 0));

    colorWheelWidget->add<nanogui::Label>("Change Mesh Color", "sans-bold", 18);
    nanogui::ColorWheel* colorWheel = colorWheelWidget->add<nanogui::ColorWheel>();

    colorWheel->setCallback([&](const nanogui::Color& color) {
        objectColor = glm::vec3(color[0], color[1], color[2]);
        });

    nanogui::CheckBox* isBoundaryColorCheckbox = colorWheelWidget->add<nanogui::CheckBox>();
    isBoundaryColorCheckbox->setCaption("Color Boundary");
    isBoundaryColorCheckbox->setCallback([](bool state) {
        isBoundaryColor = state;
        });


    nanogui::Widget* unwrappingLabels = createLabelWidget(nanoguiWindow, "Unwrapping Weights", 20);

    nanogui::ComboBox* unwrapComboBox = new nanogui::ComboBox(nanoguiWindow, { "Uniform" , "Harmonic", "Mean Value", "Random Value" });
    unwrapComboBox->setSelectedIndex(0);
    unwrapComboBox->setFixedHeight(25);

    unwrapComboBox->setCallback([&, width, height](int box) {
        switch (box) {
        case 0:
            linearUnwrappingType = UNIFORM_WEIGHTS;
            setupMeshesExceptBff(mesh);
            break;
        case 1:
            linearUnwrappingType = HARMONIC_WEIGHTS;
            setupMeshesExceptBff(mesh);
            break;
        case 2:
            linearUnwrappingType = MEAN_VALUE_WEIGHTS;
            setupMeshesExceptBff(mesh);
            break;
        case 3:
            linearUnwrappingType = RANDOM_WEIGHTS;
            setupMeshesExceptBff(mesh);
            break;
        }
        distGuiCallback(); 
        });

    nanogui::Button* improveButton = gui->addButton("Improve Weights", [&, width, height]() {
        improveUnwrapping(mesh);
        setupMeshesExceptBffWithoutComputingUVs(mesh);
        distGuiCallback();
        });

    nanogui::Widget* mappingLabels = createLabelWidget(nanoguiWindow, "Mapping Shape", 20);

    nanogui::ComboBox* mappingComboBox = new nanogui::ComboBox(nanoguiWindow, { "Disc" , "Square", "Convex Hull", "Arbitrary Polygon" });
    mappingComboBox->setSelectedIndex(0);
    mappingComboBox->setFixedHeight(25);

    mappingComboBox->setCallback([&](int box) {
        switch (box) {
        case 0:
            mappingType = DISC;
            userPoints.clear();
            setupMeshesExceptBff(mesh);
            break;
        case 1:
            mappingType = SQUARE;
            userPoints.clear();
            setupMeshesExceptBff(mesh);
            break;
        case 2:
            mappingType = CONVEX_HULL;
            userPoints.clear();
            setupMeshesExceptBff(mesh);
            break;
        case 3:
            MappingType temp = mappingType;
            mappingType = ARBITRARY;
            createPolygonPopUp();
            glfwMakeContextCurrent(window);
            if (userPoints.size() >= 3) {
                setupMeshesExceptBff(mesh);
            }
            else {
                mappingType = temp;
                nanogui::MessageDialog* message = new nanogui::MessageDialog(screen, nanogui::MessageDialog::Type::Warning, "Insufficent Points",
                    "You have to draw at least 3 points for arbitrary mapping.");
            }
            break;
        }
        distGuiCallback(); 
        });

    nanogui::Widget* loadTextureLabels = createLabelWidget(nanoguiWindow, "Load Textures", 20);

    nanogui::Widget* textureWidget = new nanogui::Widget(nanoguiWindow);

    textureWidget->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
        nanogui::Alignment::Middle, 2, 5));

    nanogui::Button* textureButton1 = textureWidget->add<nanogui::Button>("Texture1");
    textureButton1->setFixedHeight(25);
    textureButton1->setCallback([&]() {
        std::string filename = nanogui::file_dialog({ {"jpg", "JPEG Image"}, {"jpeg", "JPEG Image"}, {"png", "PNG Image"} }, false);
        if (!filename.empty()) {

            loadTexture(filename);
        }
        });

   nanogui::Button* textureButton2 = textureWidget->add<nanogui::Button>("Texture2");
   textureButton2->setFixedHeight(25);
   textureButton2->setCallback([&]() {
       std::string filename = nanogui::file_dialog({ {"jpg", "JPEG Image"}, {"jpeg", "JPEG Image"}, {"png", "PNG Image"} }, false);
       if (!filename.empty()) {

           loadTexture2(filename);
       }
       });
   nanogui::Widget* textureEditWidget = new nanogui::Widget(nanoguiWindow);
   textureEditWidget->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
       nanogui::Alignment::Middle, 2, 5));
   nanogui::PopupButton* textureEditButton = new  nanogui::PopupButton(nanoguiWindow, "Texture1 Detailed");
   nanogui::PopupButton* textureEditButton2 = new  nanogui::PopupButton(nanoguiWindow, "Texture2 Detailed");

   nanogui::Popup* textureEditPopup = textureEditButton->popup();
   textureEditPopup->setLayout(new nanogui::GroupLayout());
   textureEditButton->setFixedHeight(25);

   nanogui::Widget* scaleTextureWidget = new nanogui::Widget(textureEditPopup);
   scaleTextureWidget->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
       nanogui::Alignment::Middle, 0, 5));

   nanogui::Widget* rotateTextureWidget = new nanogui::Widget(textureEditPopup);
   rotateTextureWidget->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
       nanogui::Alignment::Middle, 0, 5));

   textureEditPopup->add<nanogui::Label>("Stretch", "sans-bold", 20);
   nanogui::Widget* stretchTextureWidget = new nanogui::Widget(textureEditPopup);
   stretchTextureWidget->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
       nanogui::Alignment::Middle, 0, 5));

   textureEditPopup->add<nanogui::Label>("Offset", "sans-bold", 20);
   nanogui::Widget* offsetTextureWidget = new nanogui::Widget(textureEditPopup);
   offsetTextureWidget->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
       nanogui::Alignment::Middle, 0, 5));

   scaleTextureWidget->add<nanogui::Label>("Scale", "sans-bold", 20);
   scaleBox = scaleTextureWidget->add<nanogui::TextBox>("");
   scaleBox->setEditable(true);
   scaleBox->setFixedSize(nanogui::Vector2i(100, 20));
   scaleBox->setValue("50");
   scaleBox->setDefaultValue("0.0");
   scaleBox->setFontSize(16);
   scaleBox->setFormat("[-]?[0-9]*\\.?[0-9]+");

   scaleBox->setCallback([](const std::string& value) {
       textureScale = std::stof(value);
       return true;
       });

   rotateTextureWidget->add<nanogui::Label>("Rotate", "sans-bold", 20);
   rotateBox = rotateTextureWidget->add<nanogui::TextBox>("");
   rotateBox->setEditable(true);
   rotateBox->setFixedSize(nanogui::Vector2i(100, 20));
   rotateBox->setValue("50");
   rotateBox->setDefaultValue("0.0");
   rotateBox->setFontSize(16);
   rotateBox->setFormat("[-]?[0-9]*\\.?[0-9]+");

   rotateBox->setCallback([](const std::string& value) {
       textureRotationAngle = std::stof(value);
       return true;
       });

   offsetTextureWidget->add<nanogui::Label>("x:", "sans-bold", 20);
   offset_xBox = offsetTextureWidget->add<nanogui::TextBox>("");
   offset_xBox->setEditable(true);
   offset_xBox->setFixedSize(nanogui::Vector2i(100, 20));
   offset_xBox->setValue("50");
   offset_xBox->setDefaultValue("0.0");
   offset_xBox->setFontSize(16);
   offset_xBox->setFormat("[-]?[0-9]*\\.?[0-9]+");

   offset_xBox->setCallback([](const std::string& value) {
       texOffset.x = std::stof(value);
       return true;
       });

   offsetTextureWidget->add<nanogui::Label>("y:", "sans-bold", 20);
   offset_yBox = offsetTextureWidget->add<nanogui::TextBox>("");
   offset_yBox->setEditable(true);
   offset_yBox->setFixedSize(nanogui::Vector2i(100, 20));
   offset_yBox->setValue("50");
   offset_yBox->setDefaultValue("0.0");
   offset_yBox->setFontSize(16);
   offset_yBox->setFormat("[-]?[0-9]*\\.?[0-9]+");

   offset_yBox->setCallback([](const std::string& value) {
       texOffset.y = std::stof(value);
       return true;
       });

   stretchTextureWidget->add<nanogui::Label>("x:", "sans-bold", 20);
   stretch_xBox = stretchTextureWidget->add<nanogui::TextBox>("");
   stretch_xBox->setEditable(true);
   stretch_xBox->setFixedSize(nanogui::Vector2i(100, 20));
   stretch_xBox->setValue("50");
   stretch_xBox->setDefaultValue("0.0");
   stretch_xBox->setFontSize(16);
   stretch_xBox->setFormat("[-]?[0-9]*\\.?[0-9]+");

   stretch_xBox->setCallback([](const std::string& value) {
       textureStretch[0] = std::stof(value);
       return true;
       });

   stretchTextureWidget->add<nanogui::Label>("y:", "sans-bold", 20);
   stretch_yBox = stretchTextureWidget->add<nanogui::TextBox>("");
   stretch_yBox->setEditable(true);
   stretch_yBox->setFixedSize(nanogui::Vector2i(100, 20));
   stretch_yBox->setValue("50");
   stretch_yBox->setDefaultValue("0.0");
   stretch_yBox->setFontSize(16);
   stretch_yBox->setFormat("[-]?[0-9]*\\.?[0-9]+");

   stretch_yBox->setCallback([](const std::string& value) {
       textureStretch[1] = std::stof(value);
       return true;
       });

   offset_xBox->setValue(std::to_string(texOffset.x));
   offset_yBox->setValue(std::to_string(texOffset.y));
   rotateBox->setValue(std::to_string(textureRotationAngle));
   scaleBox->setValue(std::to_string(textureScale));
   stretch_xBox->setValue(std::to_string(textureStretch[0]));
   stretch_yBox->setValue(std::to_string(textureStretch[1]));

   nanogui::Popup* textureEditPopup2 = textureEditButton2->popup();
   textureEditPopup2->setLayout(new nanogui::GroupLayout());
   textureEditButton2->setFixedHeight(25);

   nanogui::Widget* scaleTextureWidget2 = new nanogui::Widget(textureEditPopup2);
   scaleTextureWidget2->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
       nanogui::Alignment::Middle, 0, 5));

   nanogui::Widget* rotateTextureWidget2 = new nanogui::Widget(textureEditPopup2);
   rotateTextureWidget2->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
       nanogui::Alignment::Middle, 0, 5));

   textureEditPopup2->add<nanogui::Label>("Stretch", "sans-bold", 20);
   nanogui::Widget* stretchTextureWidget2 = new nanogui::Widget(textureEditPopup2);
   stretchTextureWidget2->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
       nanogui::Alignment::Middle, 0, 5));

   textureEditPopup2->add<nanogui::Label>("Offset", "sans-bold", 20);
   nanogui::Widget* offsetTextureWidget2 = new nanogui::Widget(textureEditPopup2);
   offsetTextureWidget2->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
       nanogui::Alignment::Middle, 0, 5));

   scaleTextureWidget2->add<nanogui::Label>("Scale", "sans-bold", 20);
   scaleBox2 = scaleTextureWidget2->add<nanogui::TextBox>("");
   scaleBox2->setEditable(true);
   scaleBox2->setFixedSize(nanogui::Vector2i(100, 20));
   scaleBox2->setValue("50");
   scaleBox2->setDefaultValue("0.0");
   scaleBox2->setFontSize(16);
   scaleBox2->setFormat("[-]?[0-9]*\\.?[0-9]+");

   scaleBox2->setCallback([](const std::string& value) {
       textureScale2 = std::stof(value);
       return true;
       });

   rotateTextureWidget2->add<nanogui::Label>("Rotate", "sans-bold", 20);
   rotateBox2 = rotateTextureWidget2->add<nanogui::TextBox>("");
   rotateBox2->setEditable(true);
   rotateBox2->setFixedSize(nanogui::Vector2i(100, 20));
   rotateBox2->setValue("50");
   rotateBox2->setDefaultValue("0.0");
   rotateBox2->setFontSize(16);
   rotateBox2->setFormat("[-]?[0-9]*\\.?[0-9]+");

   rotateBox2->setCallback([](const std::string& value) {
       textureRotationAngle2 = std::stof(value);
       return true;
       });

   offsetTextureWidget2->add<nanogui::Label>("x:", "sans-bold", 20);
   offset_xBox2 = offsetTextureWidget2->add<nanogui::TextBox>("");
   offset_xBox2->setEditable(true);
   offset_xBox2->setFixedSize(nanogui::Vector2i(100, 20));
   offset_xBox2->setValue("50");
   offset_xBox2->setDefaultValue("0.0");
   offset_xBox2->setFontSize(16);
   offset_xBox2->setFormat("[-]?[0-9]*\\.?[0-9]+");

   offset_xBox2->setCallback([](const std::string& value) {
       texOffset2.x = std::stof(value);
       return true;
       });

   offsetTextureWidget2->add<nanogui::Label>("y:", "sans-bold", 20);
   offset_yBox2 = offsetTextureWidget2->add<nanogui::TextBox>("");
   offset_yBox2->setEditable(true);
   offset_yBox2->setFixedSize(nanogui::Vector2i(100, 20));
   offset_yBox2->setValue("50");
   offset_yBox2->setDefaultValue("0.0");
   offset_yBox2->setFontSize(16);
   offset_yBox2->setFormat("[-]?[0-9]*\\.?[0-9]+");

   offset_yBox2->setCallback([](const std::string& value) {
       texOffset2.y = std::stof(value);
       return true;
       });

   stretchTextureWidget2->add<nanogui::Label>("x:", "sans-bold", 20);
   stretch_xBox2 = stretchTextureWidget2->add<nanogui::TextBox>("");
   stretch_xBox2->setEditable(true);
   stretch_xBox2->setFixedSize(nanogui::Vector2i(100, 20));
   stretch_xBox2->setValue("50");
   stretch_xBox2->setDefaultValue("0.0");
   stretch_xBox2->setFontSize(16);
   stretch_xBox2->setFormat("[-]?[0-9]*\\.?[0-9]+");

   stretch_xBox2->setCallback([](const std::string& value) {
       textureStretch2[0] = std::stof(value);
       return true;
       });

   stretchTextureWidget2->add<nanogui::Label>("y:", "sans-bold", 20);
   stretch_yBox2 = stretchTextureWidget2->add<nanogui::TextBox>("");
   stretch_yBox2->setEditable(true);
   stretch_yBox2->setFixedSize(nanogui::Vector2i(100, 20));
   stretch_yBox2->setValue("50");
   stretch_yBox2->setDefaultValue("0.0");
   stretch_yBox2->setFontSize(16);
   stretch_yBox2->setFormat("[-]?[0-9]*\\.?[0-9]+");

   stretch_yBox2->setCallback([](const std::string& value) {
       textureStretch2[1] = std::stof(value);
       return true;
       });

   offset_xBox->setValue(std::to_string(texOffset.x));
   offset_yBox->setValue(std::to_string(texOffset.y));
   rotateBox->setValue(std::to_string(textureRotationAngle));
   scaleBox->setValue(std::to_string(textureScale));
   stretch_xBox->setValue(std::to_string(textureStretch[0]));
   stretch_yBox->setValue(std::to_string(textureStretch[1]));

   offset_xBox2->setValue(std::to_string(texOffset2.x));
   offset_yBox2->setValue(std::to_string(texOffset2.y));
   rotateBox2->setValue(std::to_string(textureRotationAngle2));
   scaleBox2->setValue(std::to_string(textureScale2));
   stretch_xBox2->setValue(std::to_string(textureStretch2[0]));
   stretch_yBox2->setValue(std::to_string(textureStretch2[1]));

    toggleUVsButton = gui->addButton("Use BFF Texture", [&]() {
        useBffForTexture = useBffForTexture ? 0 : 1;
        toggleUVsButton->setCaption(useBffForTexture ? "Use Linear Texture" : "Use BFF Texture");
        setupMesh(mesh, vertices, indices, VAO, VBO, EBO);
        });

    nanogui::Widget* renderModeLabels = createLabelWidget(nanoguiWindow, "Rendering Modes", 20);

    nanogui::Widget* renderTools = new nanogui::Widget(nanoguiWindow);
    renderTools->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Horizontal,
        nanogui::Alignment::Middle, 2, 4));
    
    nanogui::CheckBox* useTextureCheckbox = renderTools->add<nanogui::CheckBox>();
    useTextureCheckbox->setCaption("Textured");
    useTextureCheckbox->setCallback([](bool state) {
        useTexture = state;
        });


    nanogui::CheckBox* isWireframeModeCheckbox = renderTools->add<nanogui::CheckBox>();
    isWireframeModeCheckbox->setCaption("Wireframe");
    isWireframeModeCheckbox->setCallback([](bool state) {
        isWireframeMode = state;
        });

    nanogui::Widget* textureOpLabels = createLabelWidget(nanoguiWindow, "Texture Operations", 20);
    
    nanogui::Widget* blendWidget = createLabelWidget(nanoguiWindow, "Blend Textures", 18);
    nanogui::Slider* blendSlider = new nanogui::Slider(blendWidget);
    blendSlider->setValue(0.0f);
    blendSlider->setFixedWidth(120);

    blendSlider->setCallback([](float value) {
        blendFactor = value;
        });

    nanogui::Label* moveTextureText = createTooltipLabel(nanoguiWindow, "Arrow keys to move the textures.", 180);
    nanogui::Label* rotateTextureText = createTooltipLabel(nanoguiWindow, "'Z' and 'X' to scale the textures.", 180);
    nanogui::Label* scaleTextureText = createTooltipLabel(nanoguiWindow, "'A' and 'D' to rotate the textures.", 180);
    nanogui::Label* switchTextureText = createTooltipLabel(nanoguiWindow, "'1' or '2' to switch between textures.", 180);

    DistortionResult areaDist = mesh->calculateTriangleAreaDistortion();
    AngleDistortionResult angleDist = mesh->calculateAngleBasedDistortion();

    string areaText = "Area Distortion: ";
    string angleText = "Angle Distortion: ";

    linearDistWindow = gui->addWindow(Eigen::Vector2i(10, 10), "Linear Method Distortion");
    linearDistWindow->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical,
        nanogui::Alignment::Middle, 2, 6));
    linearDistWindow->setPosition(Eigen::Vector2i(width * 0.5, height * 0.85));

    linearAreaDist = linearDistWindow->add<nanogui::Label>(areaText + to_string(areaDist.averageDistortionMethod1), "sans-bold", 15);
    linearAngleDist = linearDistWindow->add<nanogui::Label>(angleText + to_string(angleDist.averageDistortionUV), "sans-bold", 15);


    nonLinearDistWindow = gui->addWindow(Eigen::Vector2i(10, 10), "Non-Linear Method Distortion");
    nonLinearDistWindow->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical,
        nanogui::Alignment::Middle, 2, 6));
    nonLinearDistWindow->setPosition(nanogui::Vector2i(width * 0.8, height * 0.85));

    nanogui::Label* nonLinearAreaDist = nonLinearDistWindow->add<nanogui::Label>(areaText + to_string(areaDist.averageDistortionMethod2), "sans-bold", 15);
    nanogui::Label* nonLinearAngleDist = nonLinearDistWindow->add<nanogui::Label>(angleText + to_string(angleDist.averageDistortionBFFUV), "sans-bold", 15);


    screen->setVisible(true);
    screen->performLayout();
    
}

nanogui::Widget* createLabelWidget(nanogui::ref<nanogui::Window> nanoguiWindow, string message, int size) {

    nanogui::Widget* labels = new nanogui::Widget(nanoguiWindow);
    labels->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical,
        nanogui::Alignment::Middle, 2, 0));

    labels->add<nanogui::Label>(message, "sans-bold", size);
    return labels;
}

nanogui::Label* createTooltipLabel(nanogui::Widget* parent, string message, int labelWidth) {
    nanogui::Label* label = new nanogui::Label(parent, message, "sans-bold", 14);
    label->setColor(nanogui::Color(Eigen::Vector3f(0.6f, 0.812f, 1.0f)));
    label->setFixedWidth(labelWidth);
    return label;
}

void setupAllMeshes(Mesh* mesh)
{
    computeUVCoordinates(mesh, userPoints, linearUnwrappingType, mappingType);

    setupMesh(mesh, vertices, indices, VAO, VBO, EBO);
    setupUVMesh(mesh, verticesUV, indicesUV, VAO_UV, VBO_UV, EBO_UV);
    setupbffUVMesh(mesh, verticesbffUV, indicesbffUV, bVAO_UV, bVBO_UV, bEBO_UV);

    mesh->calculateTriangleAreaDistortion();
    mesh->calculateAngleBasedDistortion();
}

void setupMeshesExceptBff(Mesh* mesh)
{
    computeUVCoordinates(mesh, userPoints, linearUnwrappingType, mappingType);

    setupMesh(mesh, vertices, indices, VAO, VBO, EBO);
    setupUVMesh(mesh, verticesUV, indicesUV, VAO_UV, VBO_UV, EBO_UV);

    mesh->calculateTriangleAreaDistortion();
    mesh->calculateAngleBasedDistortion();
}

void setupMeshesExceptBffWithoutComputingUVs(Mesh* mesh)
{
    setupMesh(mesh, vertices, indices, VAO, VBO, EBO);
    setupUVMesh(mesh, verticesUV, indicesUV, VAO_UV, VBO_UV, EBO_UV);

    mesh->calculateTriangleAreaDistortion();
    mesh->calculateAngleBasedDistortion();
}

void setupMesh(Mesh* mesh, std::vector<float>& vertices, std::vector<unsigned int>& indices, unsigned int& VAO, unsigned int& VBO, unsigned int& EBO)
{
    vertices.clear();
    indices.clear();

    for (int i = 0; i < mesh->verts.size(); i++)
    {
        Vertex* vertex = mesh->verts[i];
           
        vertices.push_back(vertex->coords[0]);
        vertices.push_back(vertex->coords[1]);
        vertices.push_back(vertex->coords[2]);

        vertices.push_back(vertex->normal.x);
        vertices.push_back(vertex->normal.y);
        vertices.push_back(vertex->normal.z);

        vertices.push_back(mesh->verts[i]->isBoundary ? 1.0f : 0.0f);

        if (useBffForTexture)
        {
            vertices.push_back(vertex->bffuvCoordinates.u);
            vertices.push_back(vertex->bffuvCoordinates.v);
        }

        else
        {
            vertices.push_back(vertex->uvCoordinates.u);
            vertices.push_back(vertex->uvCoordinates.v);
        }
        
    }

    for (const auto& vertex : vertices) {
        minExtent = glm::min(minExtent, vertex);
        maxExtent = glm::max(maxExtent, vertex);
    }

    const void* vertices_ptr = vertices.data();

    for (int i = 0; i < mesh->tris.size(); i++)
    {
        indices.push_back(mesh->tris[i]->v1i);
        indices.push_back(mesh->tris[i]->v2i);
        indices.push_back(mesh->tris[i]->v3i);
    }
    const void* indices_ptr = indices.data();

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertices.size(), vertices_ptr, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices.size(), indices_ptr, GL_STATIC_DRAW);

    //Position
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    //Normal
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void*)(7 * sizeof(float)));
    glEnableVertexAttribArray(3);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindVertexArray(0);
}

void setupUVMesh(Mesh* mesh, std::vector<float>& vertices, std::vector<unsigned int>& indices, unsigned int& VAO, unsigned int& VBO, unsigned int& EBO)
{
    vertices.clear();
    indices.clear();
    for (int i = 0; i < mesh->verts.size(); i++)
    {
        Vertex* vertex = mesh->verts[i];

        vertices.push_back(vertex->uvCoordinates.u);
        vertices.push_back(vertex->uvCoordinates.v);

        vertices.push_back(vertex->normal.x);
        vertices.push_back(vertex->normal.y);
        vertices.push_back(vertex->normal.z);

        vertices.push_back(mesh->verts[i]->isBoundary ? 1.0f : 0.0f);

        vertices.push_back(vertex->uvCoordinates.u);
        vertices.push_back(vertex->uvCoordinates.v);

    }
        

    for (const auto& vertex : vertices) {
        minExtentUV = glm::min(minExtent, vertex);
        maxExtentUV = glm::max(maxExtent, vertex);
     }

    const void* vertices_ptr = vertices.data();

    for (int i = 0; i < mesh->tris.size(); i++)
    {
        indices.push_back(mesh->tris[i]->v1i);
        indices.push_back(mesh->tris[i]->v2i);
        indices.push_back(mesh->tris[i]->v3i);
    }
    const void* indices_ptr = indices.data();

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertices.size(), vertices_ptr, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices.size(), indices_ptr, GL_STATIC_DRAW);

    //Position
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    //Normal
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(5 * sizeof(float)));
    glEnableVertexAttribArray(2);

    glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(3);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindVertexArray(0);
 
}

void setupbffUVMesh(Mesh* mesh, std::vector<float>& vertices, std::vector<unsigned int>& indices, unsigned int& VAO, unsigned int& VBO, unsigned int& EBO)
{
    vertices.clear();
    indices.clear();
    for (int i = 0; i < mesh->verts.size(); i++)
    {
        Vertex* vertex = mesh->verts[i];

        vertices.push_back(vertex->bffuvCoordinates.u);
        vertices.push_back(vertex->bffuvCoordinates.v);

        vertices.push_back(vertex->normal.x);
        vertices.push_back(vertex->normal.y);
        vertices.push_back(vertex->normal.z);

        vertices.push_back(mesh->verts[i]->isBoundary ? 1.0f : 0.0f);

        vertices.push_back(vertex->uvCoordinates.u);
        vertices.push_back(vertex->uvCoordinates.v);
        }

    for (const auto& vertex : vertices) {
        minExtentUV = glm::min(minExtent, vertex);
        maxExtentUV = glm::max(maxExtent, vertex);
    }

    const void* vertices_ptr = vertices.data();

    for (int i = 0; i < mesh->tris.size(); i++)
    {
        indices.push_back(mesh->tris[i]->v1i);
        indices.push_back(mesh->tris[i]->v2i);
        indices.push_back(mesh->tris[i]->v3i);
    }
    const void* indices_ptr = indices.data();

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertices.size(), vertices_ptr, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices.size(), indices_ptr, GL_STATIC_DRAW);

    //Position
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    //Normal
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(5 * sizeof(float)));
    glEnableVertexAttribArray(2);

    glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(3);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindVertexArray(0);

}

void exportImage(const ViewArea& viewArea) {

    std::vector<std::pair<std::string, std::string>> fileTypes = {
        {"png", "Portable Network Graphics"},
        {"jpg", "JPEG Image"}
    };

    std::string savePath = nanogui::file_dialog(fileTypes, true);

    if (!savePath.empty()) {

        int screenWidth = static_cast<int>(viewArea.width);
        int screenHeight = static_cast<int>(viewArea.height);
        int screenX = static_cast<int>(viewArea.left);
        int screenY = static_cast<int>(viewArea.bottom);


        GLubyte* pixels = new GLubyte[3 * 1600 * 800];
        glReadPixels(0, 0, 1600, 800, GL_RGB, GL_UNSIGNED_BYTE, pixels);

        GLubyte* flippedPixels = new GLubyte[3 * screenWidth * screenHeight];
        for (int y = 0; y < screenHeight; ++y) {
            for (int x = 0; x < screenWidth; ++x) {
                // Calculate the index for the original and flipped pixel.
                // Add screenX and screenY to account for the offset in the larger image.
                int originalIndex = ((799 - y + screenY) * 1600 + (x + screenX)) * 3;
                int flippedIndex = (y * screenWidth + x) * 3;

                // Copy the pixel (RGB).
                flippedPixels[flippedIndex] = pixels[originalIndex];       // Red component
                flippedPixels[flippedIndex + 1] = pixels[originalIndex + 1]; // Green component
                flippedPixels[flippedIndex + 2] = pixels[originalIndex + 2]; // Blue component
            }
        }

        if (savePath.substr(savePath.find_last_of(".") + 1) == "png") {
            stbi_write_png(savePath.c_str(), screenWidth, screenHeight, 3, flippedPixels, screenWidth * 3);
        }
        else if (savePath.substr(savePath.find_last_of(".") + 1) == "jpg") {
            stbi_write_jpg(savePath.c_str(), screenWidth, screenHeight, 3, flippedPixels, 100); // 100 for quality
        }

        delete[] pixels;
        delete[] flippedPixels;
    }
}


void cursorPosCallback(double x, double y, glm::vec3 minExtent, glm::vec3 maxExtent, glm::vec3 minExtentUV, glm::vec3 maxExtentUV)
{
    
    if (!screen->cursorPosCallbackEvent(x, y) && meshLoaded) {
        
        if (mouseDown)
        {
            
            double dx = (mouse[0] - x);
            double dy =  (mouse[1] - y);

            if (altDown)
            {
                double cameraSpeed = 0.0015;
                dx *= cameraSpeed;
                dy *= cameraSpeed;

                if (mouse[0] < uvArea.left) {
                    cameras[0].pan(dx, dy, minExtent.x, maxExtent.y, minExtent.y, maxExtent.y);
                }
                else if (mouse[0] < bffuvArea.left) {
                    cameras[1].pan(dx, dy, minExtentUV.x, maxExtentUV.x, minExtentUV.y, maxExtentUV.y);
                }
                else {
                    cameras[2].pan(dx, dy, minExtentUV.x, maxExtentUV.x, minExtentUV.y, maxExtentUV.y);
                }
             
            }
            else
            {
                double dTheta = dx / 360;
                double dPhi = dy / 360;

                if (mouse[0] < uvArea.left) {
                    cameras[0].rotate(dTheta, dPhi);
                }
                else if (mouse[0] < bffuvArea.left) {
                    cameras[1].rotate(dTheta, dPhi);
                }
                else {
                    cameras[2].rotate(dTheta, dPhi);
                }
            }
        }

        mouse[0] = x; mouse[1] = y;
    }
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int modifiers)
{
    if (!screen->mouseButtonCallbackEvent(button, action, modifiers) && meshLoaded) {
        double x, y;
        glfwGetCursorPos(window, &x, &y);
        mouse[0] = x; mouse[1] = y;

        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
            mouseDown = true;

        }
        else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
            mouseDown = false;
        }
    }
}

void keyButtonCallback(int key, int scancode, int action, int mods)
{
    if (!screen->keyCallbackEvent(key, scancode, action, mods) && meshLoaded) {
        switch(key)
        {
            case GLFW_KEY_1:
                isTexture1 = true;
                break;
            case GLFW_KEY_2:
                isTexture1 = false;
                break;
            case GLFW_KEY_LEFT_ALT:
			    if (action == GLFW_PRESS) {
                    altDown = true;
                }
			    else if (action == GLFW_RELEASE) {
                    altDown = false;
                }
			    break;
        // move the texture
			case GLFW_KEY_RIGHT:
               (isTexture1 ? texOffset.x -= textureMoveSpeed : texOffset2.x -= textureMoveSpeed);
               offset_xBox->setValue(std::to_string(texOffset.x));
               offset_xBox2->setValue(std::to_string(texOffset2.x));
               break;
			case GLFW_KEY_LEFT:
                (isTexture1 ? texOffset.x += textureMoveSpeed : texOffset2.x += textureMoveSpeed);
                offset_xBox->setValue(std::to_string(texOffset.x));
                offset_xBox2->setValue(std::to_string(texOffset2.x));

				break;
			case GLFW_KEY_UP:
                (isTexture1 ? texOffset.y -= textureMoveSpeed : texOffset2.y -= textureMoveSpeed);
                offset_yBox->setValue(std::to_string(texOffset.y));
                offset_yBox2->setValue(std::to_string(texOffset2.y));
				break;
			case GLFW_KEY_DOWN:
                (isTexture1 ? texOffset.y += textureMoveSpeed : texOffset2.y += textureMoveSpeed);
                offset_yBox->setValue(std::to_string(texOffset.y));
                offset_yBox2->setValue(std::to_string(texOffset2.y));

				break;
        // rotate the texture
			case GLFW_KEY_D:
                (isTexture1 ? textureRotationAngle += 0.025 : textureRotationAngle2 += 0.025);
                rotateBox->setValue(std::to_string(textureRotationAngle));
                rotateBox2->setValue(std::to_string(textureRotationAngle2));

				break;
			case GLFW_KEY_A:
                (isTexture1 ? textureRotationAngle -= 0.025 : textureRotationAngle2 -= 0.025);
                rotateBox->setValue(std::to_string(textureRotationAngle));
                rotateBox2->setValue(std::to_string(textureRotationAngle2));

				break;
            // reset the texture
			case GLFW_KEY_R:
                if (isTexture1) {
                    texOffset = glm::vec2(0.0f, 0.0f);
                    textureRotationAngle = 0.0f;
                    textureScale = 1.0f;
                    textureMoveSpeed = 0.01f;
                    textureStretch = glm::vec2(0.0f, 0.0f);
                }
                else {
                    texOffset2 = glm::vec2(0.0f, 0.0f);
                    textureRotationAngle2 = 0.0f;
                    textureScale2 = 1.0f;
                    textureMoveSpeed2 = 0.01f;
                    textureStretch2 = glm::vec2(0.0f, 0.0f);

                }
                offset_xBox->setValue(std::to_string(texOffset.x));
                offset_yBox->setValue(std::to_string(texOffset.y));
                rotateBox->setValue(std::to_string(textureRotationAngle));
                scaleBox->setValue(std::to_string(textureScale));
                offset_xBox2->setValue(std::to_string(texOffset.x));
                offset_yBox2->setValue(std::to_string(texOffset.y));
                rotateBox2->setValue(std::to_string(textureRotationAngle));
                scaleBox2->setValue(std::to_string(textureScale));

				break;
        // scale the texture
			case GLFW_KEY_Z:
                if (isTexture1){
                    textureScale *= 1.1;
                    textureMoveSpeed *= 1.1;
                }
                else {
                    textureScale2 *= 1.1;
                    textureMoveSpeed2 *= 1.1;
                }
                scaleBox->setValue(std::to_string(textureScale));
                scaleBox2->setValue(std::to_string(textureScale2));

				break;
			case GLFW_KEY_X:
                if (isTexture1) {
                    textureScale *= 0.9;
                    textureMoveSpeed *= 0.9;
                }
                else {
                    textureScale2 *= 0.9;
                    textureMoveSpeed2 *= 0.9;
                }
                scaleBox->setValue(std::to_string(textureScale));
                scaleBox2->setValue(std::to_string(textureScale2));

				break;
        //stretch
            case GLFW_KEY_O:
                if (isTexture1)
                {
                    textureStretch[0] *= 1.1;
                }
                else
                {
                    textureStretch2[0] *= 1.1;
                }
                break;
            case GLFW_KEY_P:
                if (isTexture1)
                {
                    textureStretch[0] *= 0.9;
                }
                else
                {
                    textureStretch2[0] *= 0.9;
                }
                break;
            default:
				break;
        }
    }
}

void scrollCallback(double x, double y)
{
    if (mouse[0] < uvArea.left) {
        cameras[0].zoom(0.1 * y);
    }
    else if (mouse[0] < bffuvArea.left) {
        cameras[1].zoom(0.1 * y);
    }
    else{ 
        cameras[2].zoom(0.1 * y); 
    }
}

void distGuiCallback() {
    string areaText = "Area Distortion: ";
    string angleText = "Angle Distortion: ";

    DistortionResult areaDist = mesh->calculateTriangleAreaDistortion();
    AngleDistortionResult angleDist = mesh->calculateAngleBasedDistortion();

    linearAreaDist->setCaption(areaText + to_string(areaDist.averageDistortionMethod1));
    linearAngleDist->setCaption(angleText + to_string(angleDist.averageDistortionUV));

}

void updateViewArea(int width, int height)
{
    emptyArea = ViewArea(0, 0, 200, height);
    meshArea = ViewArea(200, 0, (width - 200) / 3, height);
    uvArea = ViewArea((width - 200) / 3 + 200, 0, (width - 200) / 3 + 1, height);
    bffuvArea = ViewArea((width - 200) * 2 / 3 + 200, 0, (width - 200) / 3 + 1, height);
    fullArea = ViewArea(0, 0, width, height);

}

void clearViewArea(ViewArea &area)
{
    glViewport(area.left, area.bottom, area.width, area.height);
    glScissor(area.left, area.bottom, area.width, area.height);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_SCISSOR_TEST);


    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void render(GLFWwindow* window, unsigned int& shaderProgram)
{
    while (!glfwWindowShouldClose(window))
    {
        glClearColor(255.0f, 255.0f, 255.0f, 1.0f); 
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

        clearViewArea(meshArea);

        // Render mesh
        renderMesh(meshArea, cameras[0], shaderProgram, VAO);

        clearViewArea(uvArea);
        // Render UV mesh
        renderMesh(uvArea, cameras[1], shaderProgram, VAO_UV);

        clearViewArea(bffuvArea);
        // Render UV mesh
        renderMesh(bffuvArea, cameras[2], shaderProgram, bVAO_UV);

        if (captureFlag == 1) {
            // Opem a nanogui file dialog to save obj file
            std::vector<std::pair<std::string, std::string>> fileTypes = {
                {"obj", "Wavefront .obj file"}
            };

            std::string savePath = nanogui::file_dialog(fileTypes, true);

            if (!savePath.empty()) {
                if (textureName1 == "" && textureName2 == "") {
                    mesh->exportMesh(savePath, false, useBffForTexture);
                }
                else{
                    mesh->exportMesh(savePath, true, useBffForTexture);
                    blendAndSaveTextures(savePath);
                }
            }
            captureFlag = 0;
        }

        screen->drawContents();
        screen->drawWidgets();

        glBindBuffer(GL_UNIFORM_BUFFER, 0);

        glfwSwapBuffers(window);
        if (isPointAdded) break;
    }
}

void renderMesh(ViewArea &area, Camera& camera, unsigned int &shaderProgram, unsigned int& VAO)
{
        glfwPollEvents();

        glm::mat4 projection = camera.projectionMatrix(area.width, area.height);
        glm::mat4 view = camera.viewMatrix();
        glm::mat4 viewMesh = cameras[0].viewMatrix();
        glm::mat4 model = glm::mat4(1.0f);

        glUseProgram(shaderProgram);

        GLint projLocation = glGetUniformLocation(shaderProgram, "projection");
        GLint viewLocation = glGetUniformLocation(shaderProgram, "view");
        GLint viewMeshLocation = glGetUniformLocation(shaderProgram, "viewMesh");
        GLint modelLocation = glGetUniformLocation(shaderProgram, "model");

        glUniformMatrix4fv(viewLocation, 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(viewMeshLocation, 1, GL_FALSE, glm::value_ptr(viewMesh));
        glUniformMatrix4fv(projLocation, 1, GL_FALSE, glm::value_ptr(projection));
        glUniformMatrix4fv(modelLocation, 1, GL_FALSE, glm::value_ptr(model));

        GLint objectColorLocation = glGetUniformLocation(shaderProgram, "objectColor");
        GLint lightPosLocation = glGetUniformLocation(shaderProgram, "lightPos");
        GLint lightColorLocation = glGetUniformLocation(shaderProgram, "lightColor");
        GLint boundaryColorLocation = glGetUniformLocation(shaderProgram, "boundaryColor");
        GLint useTextureLocation = glGetUniformLocation(shaderProgram, "useTexture");
        GLint isBoundaryColorLocation = glGetUniformLocation(shaderProgram, "isBoundaryColor");

        GLuint textureOffsetLocation = glGetUniformLocation(shaderProgram, "textureOffset");
        GLuint textureRotationAngleLocation = glGetUniformLocation(shaderProgram, "textureRotationAngle");
        GLint textureScaleLocation = glGetUniformLocation(shaderProgram, "textureScale");
        GLuint textureStretchLocation = glGetUniformLocation(shaderProgram, "textureStretch");

        GLuint textureOffset2Location = glGetUniformLocation(shaderProgram, "textureOffset2");
        GLuint textureRotationAngle2Location = glGetUniformLocation(shaderProgram, "textureRotationAngle2");
        GLint textureScale2Location = glGetUniformLocation(shaderProgram, "textureScale2");
        GLuint textureStretch2Location = glGetUniformLocation(shaderProgram, "textureStretch2");

        glUniform3fv(objectColorLocation, 1, glm::value_ptr(objectColor));
        glUniform3fv(lightPosLocation, 1, glm::value_ptr(glm::vec3(0.0, 2.0, 10.0)));
        glUniform3fv(lightColorLocation, 1, glm::value_ptr(glm::vec3(1.0, 1.0, 1.0)));
        glUniform3fv(boundaryColorLocation, 1, glm::value_ptr(boundaryColor));

        glUniform1i(useTextureLocation, useTexture);
        glUniform1i(isBoundaryColorLocation, isBoundaryColor);
        
        glUniform2fv(textureOffsetLocation, 1, glm::value_ptr(texOffset));
        glUniform1f(textureRotationAngleLocation, textureRotationAngle);
        glUniform1f(textureScaleLocation, textureScale);
        glUniform2fv(textureStretchLocation, 1, glm::value_ptr(textureStretch));

        glUniform2fv(textureOffset2Location, 1, glm::value_ptr(texOffset2));
        glUniform1f(textureRotationAngle2Location, textureRotationAngle2);
        glUniform1f(textureScale2Location, textureScale2);
        glUniform2fv(textureStretch2Location, 1, glm::value_ptr(textureStretch2));

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, textureID);
        glUniform1i(glGetUniformLocation(shaderProgram, "textureSampler"), 0);

        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, textureID2);
        glUniform1i(glGetUniformLocation(shaderProgram, "textureSampler2"), 1);

        glUniform1f(glGetUniformLocation(shaderProgram, "blendFactor"), blendFactor);

        glBindVertexArray(VAO);

        if (isWireframeMode)
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        }
        
        glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        GLenum error = glGetError();
        if (error != GL_NO_ERROR) {
            std::cout << "OpenGL error:" << error << std::endl;
        }

        glBindVertexArray(0);
  
}

void loadTexture(std::string filename)
{
    int width, height, nrChannels;

    if (textureID != 0) {
        glDeleteTextures(1, &textureID);
    }

    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);

    const char* newfilename = filename.c_str();

    unsigned char* data = stbi_load(newfilename, &width, &height, &nrChannels, 0);

    if (data)
    {

        if (nrChannels == 1)
        {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, width, height, 0, GL_RED, GL_UNSIGNED_BYTE, data);
        }
        else if (nrChannels == 3)
        {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
        }
        else if(nrChannels == 4)
        {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
        }
        glGenerateMipmap(GL_TEXTURE_2D);

        GLenum error = glGetError();
        if (error != GL_NO_ERROR) {
            std::cout << "OpenGL error: " << error << std::endl;
        }

    }
    else
    {
        std::cout << "Failed to load texture" << std::endl;
    }
    textureName1 = filename;

    stbi_image_free(data);
}

void loadTexture2(std::string filename)
{
    int width, height, nrChannels;

    if (textureID2 != 0) {
        glDeleteTextures(1, &textureID2);
    }

    glGenTextures(1, &textureID2);
    glBindTexture(GL_TEXTURE_2D, textureID2);

    const char* newfilename = filename.c_str();

    unsigned char* data = stbi_load(newfilename, &width, &height, &nrChannels, 0);

    if (data)
    {

        if (nrChannels == 1)
        {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, width, height, 0, GL_RED, GL_UNSIGNED_BYTE, data);
        }
        else if (nrChannels == 3)
        {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
        }
        else if(nrChannels == 4)
        {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
        }
        glGenerateMipmap(GL_TEXTURE_2D);

        GLenum error = glGetError();
        if (error != GL_NO_ERROR) {
            std::cout << "OpenGL error: " << error << std::endl;
        }

    }
    else
    {
        std::cout << "Failed to load texture" << std::endl;
    }
    textureName2 = filename;

    stbi_image_free(data);
}

std::string readShaderFile(const std::string& shaderPath) {
    std::string shaderCode;
    std::ifstream shaderFile;
    shaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    try {
        shaderFile.open(shaderPath);
        std::stringstream shaderStream;
        shaderStream << shaderFile.rdbuf();
        shaderFile.close();
        shaderCode = shaderStream.str();
    }
    catch (std::ifstream::failure& e) {
        std::cerr << "ERROR::SHADER::FILE_NOT_SUCCESSFULLY_READ: " << shaderPath << std::endl;
    }
    return shaderCode;
}

void initPopupGui(nanogui::Screen* popupGui) {
    
    nanogui::FormHelper* gui = new nanogui::FormHelper(popupGui);
    nanogui::ref<nanogui::Window> nanoguiPopupWindow = gui->addWindow(Eigen::Vector2i(8, 8), "Polygon Drawing Panel");
    nanoguiPopupWindow->setLayout(new nanogui::GroupLayout());

    
    nanogui::Widget* tools = new nanogui::Widget(nanoguiPopupWindow);
    tools->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical,
        nanogui::Alignment::Middle, 0, 6));

    nanogui::Label* drawInst = tools->add<nanogui::Label>("Draw a convex polygon by clicking on the screen.", "sans-bold", 15);
    nanogui::Label* clearInst = tools->add<nanogui::Label>("Press 'C' to clear the drawing.", "sans-bold", 15);
    nanogui::Label* escInst = tools->add<nanogui::Label>("Press 'ESC' to exit the canvas.", "sans-bold", 15);

    drawInst->setFixedSize(Eigen::Vector2i(150, 40));
    clearInst->setFixedSize(Eigen::Vector2i(150, 30));
    escInst->setFixedSize(Eigen::Vector2i(150, 30));

    popupGui->setVisible(true);
    popupGui->performLayout();
    nanoguiPopupWindow->setFixedHeight(25);
}


bool isConvex(float x, float y) {
    if (userPoints.size() < 3) return true;
    
    glm::vec2 endPoint = userPoints[userPoints.size() - 1];
    glm::vec2 beforeEndPoint = userPoints[userPoints.size() - 2];

    glm::vec2 afterStartPoint = userPoints[1];
    glm::vec2 startPoint = userPoints[0];

    float endConvexCheck = (beforeEndPoint.y - endPoint.y) * (endPoint.x - x) - (endPoint.y - y) * (beforeEndPoint.x - endPoint.x);
    float startConvexCheck = (startPoint.y - afterStartPoint.y) * (afterStartPoint.x - x) - (afterStartPoint.y - y) * (startPoint.x - afterStartPoint.x);
    float middleConvexCheck = (endPoint.y - y) * (x - startPoint.x) - (y - startPoint.y) * (endPoint.x - x);

    if (endConvexCheck < 0 || startConvexCheck < 0 || middleConvexCheck < 0) return false;

    return true;
}

void polygonMouseButtonCallback(GLFWwindow* popupWindow, int button, int action, int modifiers)
{
    if (!popupGui->mouseButtonCallbackEvent(button, action, modifiers)) {
        if (glfwGetMouseButton(popupWindow, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS && !isMousePressed && !isWarningMessage) {
            isMousePressed = 1;
        }

        if (glfwGetMouseButton(popupWindow, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE && isMousePressed && !isWarningMessage) {
            double mouseX, mouseY;
            glfwGetCursorPos(popupWindow, &mouseX, &mouseY);

            int width, height;
            glfwGetWindowSize(popupWindow, &width, &height);

            float normalizedX = static_cast<float>(mouseX) / width * 2 - 1;
            float normalizedY = -(static_cast<float>(mouseY) / height * 2 - 1);

            bool isConvexAfterAddition = true;

            if (isConvex(normalizedX, normalizedY)) {
                userPoints.emplace_back(normalizedX, normalizedY);
                if (userPoints.size() >= 3) isPointAdded = 1;
            }
            else {
                isWarningMessage = 1;
            }
            isMousePressed = 0;
        }
    }
}

void createPolygonPopUp()
{
    int popupWidth = 800, popupHeight = 600;

    GLFWwindow* popupWindow = glfwCreateWindow(popupWidth, popupHeight, "Polygon Canvas", NULL, screen->glfwWindow());
    if (!popupWindow) {
        std::cerr << "Failed to create GLFW popup window" << std::endl;
        return;
    }

    glfwMakeContextCurrent(popupWindow);
    glfwFocusWindow(popupWindow);
    isMousePressed = 0;
    isWarningMessage = 0;
    isDrawn = 0;
   
    popupGui = new nanogui::Screen();
    popupGui->initialize(popupWindow, false);


    glfwSetMouseButtonCallback(popupWindow, [](GLFWwindow* popupWindow, int button, int action, int modifiers) {
        polygonMouseButtonCallback(popupWindow, button, action, modifiers);
        });


    initPopupGui(popupGui);
    
    nanogui::MessageDialog* message = NULL;
    unsigned int popupShaderProgram = initializeShaders(polygonVertexShaderPath, polygonFragmentShaderPath);
    
    
    while (!glfwWindowShouldClose(popupWindow)) {

        glClearColor(255.0f, 255.0f, 255.0f, 1.0f); 
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
        glUseProgram(popupShaderProgram);
        
        drawPolygon(popupWindow, popupGui, isMousePressed, isWarningMessage);
        popupGui->drawContents();
        popupGui->drawWidgets();

        if (isWarningMessage) {
            if (isWarningMessage == 1) {
                message = new nanogui::MessageDialog(popupGui, nanogui::MessageDialog::Type::Warning, "Convexity Violation",
                    "Adding this point would break convexity!", "");
            } 
            isWarningMessage++;

            if (isWarningMessage > 180) {
                message->dispose();
                isWarningMessage = 0;
            }
        }

        glfwPollEvents();
        glfwSwapBuffers(popupWindow);

        if (isPointAdded && isDrawn) {
            glfwMakeContextCurrent(screen->glfwWindow());
            setupMeshesExceptBff(mesh);
            distGuiCallback();
            render(screen->glfwWindow(), shaderProgram);
            glfwMakeContextCurrent(popupWindow);
            isPointAdded = 0; isDrawn = 0;
        }
    }
    
    glfwDestroyWindow(popupWindow);
}

void drawPolygon(GLFWwindow* popupWindow, nanogui::Screen* popupScreen,  int& isMousePressed, int& isWarningMessage) {

    if (userPoints.size() > 1) {
        int width, height;
        glfwGetWindowSize(popupWindow, &width, &height);
        glViewport(0, 0, width, height);
        
        GLuint VAO, VBO;
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);

        glBindVertexArray(VAO);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, userPoints.size() * sizeof(glm::vec2), &userPoints[0], GL_STATIC_DRAW);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glDrawArrays(GL_LINE_LOOP, 0, userPoints.size());

        glBindVertexArray(0);
        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &VBO);

        if (isPointAdded) isDrawn = 1;
    }

    if (glfwGetKey(popupWindow, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(popupWindow, GLFW_TRUE);
    }

    if (glfwGetKey(popupWindow, GLFW_KEY_C) == GLFW_PRESS) {
        userPoints.clear();
    }
}

cv::Mat transformTexture(const cv::Mat& image, const cv::Vec2f& offset, float rotationAngle, float scale) {
    cv::Point2f center(image.cols / 2.0f, image.rows / 2.0f);

    float adjustedScale = 1.0f / scale;

    cv::Mat transformMat = cv::getRotationMatrix2D(center, (rotationAngle + CV_PI) * 180.0 / CV_PI, adjustedScale);

    transformMat.at<double>(0, 2) += offset[0] * image.cols;
    transformMat.at<double>(1, 2) += offset[1] * image.rows;

    cv::Mat transformedImage;
    cv::warpAffine(image, transformedImage, transformMat, image.size(), cv::INTER_LINEAR, cv::BORDER_WRAP);

    cv::Mat flippedImage;
    cv::flip(transformedImage, flippedImage, 1);

    return flippedImage;
}

cv::Mat loadAndProcessTexture(const std::string& filename, const cv::Vec2f& textureOffset, float textureRotationAngle, float textureScale) {
    int width, height, nrChannels;

    unsigned char* data = stbi_load(filename.c_str(), &width, &height, &nrChannels, 0);
    if (!data) {
        std::cerr << "Failed to load texture: " << filename << std::endl;
        return cv::Mat();
    }

    int cvType;
    if (nrChannels == 1) {
        cvType = CV_8UC1;
    } else if (nrChannels == 3) {
        cvType = CV_8UC3;
    } else if (nrChannels == 4) {
        cvType = CV_8UC4;
    } else {
        std::cerr << "Unsupported number of channels: " << nrChannels << std::endl;
        stbi_image_free(data);
        return cv::Mat();
    }

    cv::Mat image(height, width, cvType, data);

    if (nrChannels == 3 || nrChannels == 4) {
        cv::cvtColor(image, image, nrChannels == 3 ? cv::COLOR_RGB2BGR : cv::COLOR_RGBA2BGRA);
    }

    cv::Mat transformedImage = transformTexture(image, textureOffset, textureRotationAngle, textureScale);

    stbi_image_free(data);

    return transformedImage;
}

void blendAndSaveTextures(std::string& directory) {
	if (directory.find(".obj") != std::string::npos) {
        std::cerr << "Directory: " << directory << std::endl;
		directory = directory.substr(0, directory.size() - 4);
		directory += "_texture.jpg";
	}

    if (textureName1.empty() && textureName2.empty()) {
        std::cerr << "There are no textures to blend." << std::endl;
        return;
    }

    cv::Mat texture1, texture2;

    if (!textureName1.empty()) {
        texture1 = loadAndProcessTexture(textureName1, cv::Vec2f(texOffset.x, texOffset.y), textureRotationAngle, textureScale);
        if (texture1.empty()) {
            std::cerr << "Failed to load or process the first texture." << std::endl;
            return;
        }
    }

    if (!textureName2.empty()) {
        texture2 = loadAndProcessTexture(textureName2, cv::Vec2f(texOffset2.x, texOffset2.y), textureRotationAngle2, textureScale2);
        if (texture2.empty()) {
            std::cerr << "Failed to load or process the second texture." << std::endl;
            return;
        }
    }

    if (!textureName1.empty() && textureName2.empty()) {
        if (!cv::imwrite(directory, texture1)) {
            std::cerr << "Error: Could not save the processed image." << std::endl;
        } else {
            std::cout << "Processed image saved successfully." << std::endl;
        }
    } else if (textureName1.empty() && !textureName2.empty()) {
        if (!cv::imwrite(directory, texture2)) {
            std::cerr << "Error: Could not save the processed image." << std::endl;
        } else {
            std::cout << "Processed image saved successfully." << std::endl;
        }
    } else {
        bool resized = false;
        
        if (texture1.channels() != texture2.channels()) {
            if (texture1.channels() < texture2.channels()) {
                cv::cvtColor(texture1, texture1, (texture2.channels() == 4) ? cv::COLOR_BGR2BGRA : cv::COLOR_BGR2RGB);
            } else {
                cv::cvtColor(texture2, texture2, (texture1.channels() == 4) ? cv::COLOR_BGR2BGRA : cv::COLOR_BGR2RGB);
            }
            std::cerr << "Converted textures to have the same number of channels." << std::endl;
        }
        
        if (texture1.size() != texture2.size()) {
            if (texture1.total() > texture2.total()) {
                cv::resize(texture2, texture2, texture1.size());
            } else {
                cv::resize(texture1, texture1, texture2.size());
            }
        }

        if (texture1.type() != texture2.type()) {
            texture1.convertTo(texture1, CV_32F);
            texture2.convertTo(texture2, CV_32F);

            resized = true;
        }
 
        if (texture1.size() != texture2.size() || texture1.type() != texture2.type()) {
            std::cerr << "Textures must have the same size and type for blending." << std::endl;
            return;
        }

        cv::Mat blendedTexture;
        cv::addWeighted(texture1, 1.0 - blendFactor, texture2, blendFactor, 0.0, blendedTexture);

        if (resized) {
            int originalType = (texture1.depth() > texture2.depth()) ? texture1.type() : texture2.type();

            blendedTexture.convertTo(blendedTexture, originalType);
        }

        if (!cv::imwrite(directory, blendedTexture)) {
            std::cerr << "Error: Could not save the blended image." << std::endl;
        } else {
            std::cout << "Blended image saved successfully." << std::endl;
        }
    }

    texture1.release();
    texture2.release();
}