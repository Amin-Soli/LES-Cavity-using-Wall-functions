
#include <iostream>
#include <math.h>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <string>

using namespace std;

#ifdef _WIN32
#include <windows.h>
#include <sys/stat.h>
#include <dirent.h>

#define DIV 1048576 
#define WIDTH 7
#endif

#ifdef linux
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
#endif

#include "../libs/flowModelType/flowModelType/flowModelType.C"
#include "../libs/simpleSolver/simpleSolver.C"
#include "../libs/wallFunctions/wallFunction/wallFunction.C"



