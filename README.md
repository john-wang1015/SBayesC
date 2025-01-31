# SBayesC

This repo contains the implementation of SBayesC model with different sampling algorithm. It depends on Eigen library and boost library. 

***
## How to run the code

The code has been tested in both Windows and MacOS enviroment. 

For windows user, add the path of Eigen and boost to the project path, then run Visual Studio. VS should compile and run it directly without any modifications.

For MacOS user, download the Eigen and boost library and then run the following code in terminal:
```
[Your Compiler] -std=c++[C++ Version] -I [Path to Eigen] [Path to boost] main.cpp data.cpp gibbs.cpp model.cpp stat.cpp inference.cpp -o main
```
make sure your MacOS has complier installer.

The make file is coming soon.
