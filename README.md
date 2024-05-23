# TMAM

For running OpenGL we need to download and include some headers etc. This guide is for Windows only.

Go to https://cmake.org/download/ and download CMake
Go to https://www.glfw.org/download.html and download the source package of GLFW, and then extract
Open CMake, select the location of the extracted source package of GLFW as source code (eg. C:/Users/user/Downloads/glfw-3.3.8)
Create a folder named "build" under glfw-3.3.8, and then select it as the build location (eg. C:/Users/user/Downloads/glfw-3.3.8/build)
Click Configure, and from there select your Visual Studio version, leave other fields as default and click Finish
Click Configure, and then Generate, after the generation completes, you can close CMake.
Go to glfw-3.3.8\build and open GLFW.sln with Visual Studio and build
Go to https://glad.dav1d.de/, set language to C/C++, specification to OpenGL, profile to Core, select API gl as "Version 3.3" and click generate and download glad.zip
We now need to create 2 folders under C:, named Libs and Includes (the path is important, since the project is set up to look at those specific paths, namely C:\Libs, C:\Includes)
Go to your glfw-3.3.8\build\src\Debug and copy glfw3.lib to C:\Libs folder
Go to glfw-3.3.8\include and copy the folder named GLFW to C:\Includes
Unzip glad.zip, go to include and copy both folders to C:\Includes

After following these steps, you should be able to run the project successfully.
For NanoGUI Setup:

First, actually we will need to reverse some of the initialization we did due to conflicts regarding glad. (NanoGUI already includes glad). Delete the glad.h header from your Includes directory.
Clone the NanoGUI from https://github.com/wjakob/nanogui/tree/master to C:\ (Make sure that the folder name is nanogui and there are no sub-folders)
Create a folder named build in nanogui.
In build directory first run git submodule update --init --recursive and then run cmake -G "Visual Studio 17 2022" -A x64 ..
Now open the NanoGUI.sln in the build directory using Visual Studio.
Batch build from build tab, select all. (It may be unnecessary to build all but better safe than sorry)
Copy the nanogui.lib file inside the nanogui\build\Debug directory to C:\Libs
Copy the nanogu.dll file inside the nanogui\build\Debug directory to your working directory.

You should be able to run the program and see a NanoGUI screen embedding a glfw window now.
For writing images:

Go to https://github.com/nothings/stb/blob/master/stb_image_write.h and download the file
Put the header file into C:\Includes\stb

For including GLM (opengl mathematics):

Go to https://github.com/g-truc/glm/releases/tag/0.9.9.8 and download the zip.
Extract the files, go to the root file of the headers (glm-0.9.9.8\glm\glm) and copy it to your includes folder.

For libigl

Go to https://github.com/libigl/libigl and download the zip.
Extract the file, copy the igl folder inside the include folder to C:\Includes

To fix a bug:

Go to Includes/igl/adjacency_list.cpp and replace line 144 with:
maxCoeff = (std::max)(coeff, maxCoeff);
Go to Includes/igl/principal_curvature.cpp and replace line 228 with:
max_sp = - (std::numeric_limits::max)();
This ensures compiler properly recognizes the functions. (No idea why it would not work without this)

## Acknowledgements

- For the non-linear unwrapping, we directly used https://github.com/GeometryCollective/boundary-first-flattening
- We implemented "S. Yoshizawa, A. Belyaev and H. P. Seidel, "A fast and simple stretch-minimizing mesh parameterization," Proceedings Shape Modeling Applications, 2004., Genova, Italy, 2004, pp. 200-208, doi: 10.1109/SMI.2004.1314507"

