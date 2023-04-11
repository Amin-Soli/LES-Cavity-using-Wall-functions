
const char *path="../inputData/data.txt";

int numberOfNodes_yDirection, numberOfNodes_xDirection;
    
double Length, start_time;

ifstream file(path);

string searchForLength = "Length";
string searchForStartTime = "Startime";
string searchForNumberOfNodesInXDirection = "NumberOfNodesInXDirection";
string searchForNumberOfNodesInYDirection = "NumberOfNodesInYDirection";


string lineOfText;

for(;;)
{
	getline(file, lineOfText);

	if (file.eof()) break;

	if (lineOfText.find(searchForLength, 0) != string::npos)
	  file >> Length;
	  
	if (lineOfText.find(searchForStartTime, 0) != string::npos)
		  file >> start_time;
		  
	if (lineOfText.find(searchForNumberOfNodesInXDirection, 0) != string::npos)
	  file >> numberOfNodes_xDirection;

	if (lineOfText.find(searchForNumberOfNodesInYDirection, 0) != string::npos)
	  file >> numberOfNodes_yDirection;
}

file.close();

double dx, dy;
    
dx = Length/numberOfNodes_xDirection;
dy = Length/numberOfNodes_yDirection;

double **u, **v, **p, **mu_turbulence;

mu_turbulence = new double *[numberOfNodes_yDirection];    
u = new double *[numberOfNodes_yDirection];
v = new double *[numberOfNodes_yDirection+1];
p = new double *[numberOfNodes_yDirection];

int i, j;

for(j = 0; j <numberOfNodes_yDirection+1; j++)
	v[j] = new double[numberOfNodes_xDirection];
	
for(j = 0; j <numberOfNodes_yDirection; j++)
{
	u[j] = new double[numberOfNodes_xDirection+1];
	p[j] = new double[numberOfNodes_xDirection];
	mu_turbulence[j] = new double[numberOfNodes_xDirection];
}

if (start_time == 0)       // initialize field variables from t=0
{
	// initialize scalar field variables:
	for(j=0;j<numberOfNodes_yDirection;j++)
		for(i=0;i<numberOfNodes_xDirection;i++)
			p[j][i] = 0;

	// initialize u:
	for(j=0;j<numberOfNodes_yDirection;j++)
		for(i=0;i<numberOfNodes_xDirection+1;i++)
			u[j][i] = 0;

	// initialize v:
	for(j=0;j<numberOfNodes_yDirection+1;j++)
		for(i=0;i<numberOfNodes_xDirection;i++)
			v[j][i] = 0;
}

else           // initialize field variables from last results file
{ 
	#include "initializeFromResultsFile.H"
}
		
simpleSolver SIMPLE_Object(path, u, v, p, mu_turbulence);

flowModelType *FlowModel_Object;

FlowModel_Object = FlowModel_Object -> flowModelType::New(path, u, v, mu_turbulence);

wallFunction *wallFunction_Object;

wallFunction_Object = wallFunction_Object -> wallFunction::New(path, u, v);
