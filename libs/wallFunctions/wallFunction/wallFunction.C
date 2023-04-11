#include "wallFunction.H"
#include "../standardWallFunction/standardWallFunction.C"
#include "../enhancedWallFunction/enhancedWallFunction.C"
#include "../notWallFunction/notWallFunction.C"

//////////////// constructor ///////////////////////

wallFunction::wallFunction (const char *path, double **U, double **V)
{
	
	readInputData(path);
	
	u = U;
	v = V;
	
	time = start_time;
	
	count_for_writingSolutionField = 0;
	
	dx = Length/numberOfNodes_xDirection;
	dy = Length/numberOfNodes_yDirection;
	
	rho = new double *[numberOfNodes_yDirection];
    mu = new double *[numberOfNodes_yDirection];
    nu = new double *[numberOfNodes_yDirection];
    
    for(int j = 0; j <numberOfNodes_yDirection; j++)
    {
		rho[j] = new double[numberOfNodes_xDirection];
		mu[j] = new double[numberOfNodes_xDirection];
		nu[j] = new double[numberOfNodes_xDirection];
	}
	
	for(int j=0;j<numberOfNodes_yDirection;j++)
        for(int i=0;i<numberOfNodes_xDirection;i++)
        {
            rho[j][i] = rho_;
            mu[j][i] = mu_;
            nu[j][i] = mu[j][i]/rho[j][i];
        }
	
	Yplus_bottomWall = new double [numberOfNodes_xDirection];
	Yplus_leftWall = new double [numberOfNodes_yDirection];
	Yplus_rightWall = new double [numberOfNodes_yDirection];
	
	tau_bottomWall = new double [numberOfNodes_xDirection];
    tau_leftWall = new double [numberOfNodes_yDirection];
    tau_rightWall = new double [numberOfNodes_yDirection];
	
	nu_bottomWall = new double [numberOfNodes_xDirection];
    nu_leftWall = new double [numberOfNodes_yDirection];
    nu_rightWall = new double [numberOfNodes_yDirection];
	
}

//////////////// destructor  ///////////////////////

wallFunction::~wallFunction ()
{
	
	for(int j = 0; j < numberOfNodes_yDirection; ++j)
	{
		 delete[] mu[j];
		 delete[] nu[j];
         delete[] rho[j];
	}
	
	delete[] rho;
    delete[] mu;
    delete[] nu;
    
    delete[] Yplus_bottomWall;
    delete[] Yplus_leftWall;
    delete[] Yplus_rightWall;
    
    delete[] tau_bottomWall;   
    delete[] tau_leftWall;
    delete[] tau_rightWall;
    
    delete[] nu_bottomWall;   
    delete[] nu_leftWall;
    delete[] nu_rightWall;
	
}
	
//////////// public functions /////////////////////////////////////////////////////////////////////////

void wallFunction::update_time()
{
	check_writing_solutionField();
	time = time + dt;
	count_for_writingSolutionField++;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void wallFunction::check_writing_solutionField()
{
	if(count_for_writingSolutionField % numberOfDt_for_writingSolutionField == 0)
		write_solutionField();	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void wallFunction::write_solutionField()
{
	int i,j;
	
	stringstream ss;
	ss << time;
	string str = ss.str();
	string myName = "results/t=" + str ;	

    ofstream file1;
    string tau_bottomWallPath = myName + "/ShearStress_bottomWall.txt";
    const char* tau_bottomWallPathName = tau_bottomWallPath.c_str();	
    file1.open(tau_bottomWallPathName);

    file1 << fixed << setprecision(6);
    
    for (i=0;i<numberOfNodes_xDirection;i++)
		file1 << i*dx + dx/2.0 << ", " << tau_bottomWall[i] << endl;

    file1.close();
    
    ofstream file2;
    string tau_leftWallPath = myName + "/ShearStress_leftWall.txt";
    const char* tau_leftWallPathName = tau_leftWallPath.c_str();	
    file2.open(tau_leftWallPathName);

    file2 << fixed << setprecision(6);
    
    for (j=0;j<numberOfNodes_yDirection;j++)
		file2 << j*dy + dy/2.0 << ", " << tau_leftWall[j] << endl;

    file2.close();
    
    ofstream file3;
    string tau_rightWallPath = myName + "/ShearStress_rightWall.txt";
    const char* tau_rightWallPathName = tau_rightWallPath.c_str();	
    file3.open(tau_rightWallPathName);

    file3 << fixed << setprecision(6);
    
    for (j=0;j<numberOfNodes_yDirection;j++)
		file3 << j*dy + dy/2.0 << ", " << tau_rightWall[j] << endl;

    file3.close();
    
    ofstream file4;
    string skinFrictionPath = myName + "/WallsShearStressAverage.txt";
    const char* skinFrictionPathName = skinFrictionPath.c_str();	
    file4.open(skinFrictionPathName);

    file4 << fixed << setprecision(6);
    
	file4 << "Shear Stress Average in X direction for bottom wall is = " << average_tau_bottomWall << endl << endl;
	      
	file4 << "Shear Stress Average in Y direction for left wall is = " << average_tau_leftWall << endl << endl;
	      
	file4 << "Shear Stress Average in Y direction for right wall is = " << average_tau_rightWall << endl << endl;

    file4.close();

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void wallFunction::readInputData(const char *path)
{
	ifstream file(path);
	
	string searchForFlowModelType_Name = "FlowModelType";
	string searchForUsingWallFunctionFlag = "UsingWallFunction";
	string searchForLength = "Length";
	string searchForStartTime = "Startime";
	string searchForDelta_t = "Delta_t";
	string searchForNumberOfDt_for_writingSolutionFields = "NumberOfDtForWritingSolutionFields";
	string searchForMovingWallVelocity = "MovingWallVelocity";
	string searchForDensity = "Density";
	string searchForDynamicViscosity = "DynamicViscosity";
	string searchForNumberOfNodesInXDirection = "NumberOfNodesInXDirection";
	string searchForNumberOfNodesInYDirection = "NumberOfNodesInYDirection";
	string searchForKappa = "Kappa";
	string searchForE = "E_";

	string lineOfText;

	for(;;)
	{
		getline(file, lineOfText);

		if (file.eof()) break;

		if (lineOfText.find(searchForFlowModelType_Name, 0) != string::npos)
		  file >> flowModel;
		  
		if (lineOfText.find(searchForUsingWallFunctionFlag, 0) != string::npos)
		  file >> usingWallFunction_flag;
		
		if (lineOfText.find(searchForLength, 0) != string::npos)
		  file >> Length;
		  
		if (lineOfText.find(searchForStartTime, 0) != string::npos)
		  file >> start_time;
		  
		if (lineOfText.find(searchForMovingWallVelocity, 0) != string::npos)
		  file >> free_velocity;
		  
		if (lineOfText.find(searchForNumberOfDt_for_writingSolutionFields, 0) != string::npos)
		  file >> numberOfDt_for_writingSolutionField;
		  
		if (lineOfText.find(searchForDelta_t, 0) != string::npos)
		  file >> dt;
		
		if (lineOfText.find(searchForDensity, 0) != string::npos)
		  file >> rho_;

		if (lineOfText.find(searchForDynamicViscosity, 0) != string::npos)
		  file >> mu_;
		
		if (lineOfText.find(searchForNumberOfNodesInXDirection, 0) != string::npos)
		  file >> numberOfNodes_xDirection;

		if (lineOfText.find(searchForNumberOfNodesInYDirection, 0) != string::npos)
		  file >> numberOfNodes_yDirection;
		  
		if (lineOfText.find(searchForKappa, 0) != string::npos)
		  file >> kappa;
		  
		if (lineOfText.find(searchForE, 0) != string::npos)
		  file >> E;

    }

    file.close();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double wallFunction::calculate_Yplus_forEachFace_usingNewtonRaphson(double U, double d, double rho, double mu)
{
	
	double nu = mu/rho;
	double kappaRe = kappa*U*d/nu;
	
	double yp = 11.25;
	double ryPlusLam = 1.0/yp;
    
	int iter = 0;
	double yPlusLast = 0.0;
	do
	{
		yPlusLast = yp;
		yp = (kappaRe + yp)/(1.0 + log(E*yp));
	} while (fabs(ryPlusLam*(yp - yPlusLast)) > 0.01 && ++iter < 10 );
	
	return max(0.0, yp);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void wallFunction::calculate_Yplus_for_allBoundaryWalls()
{
	
	double velocity_x, velocity_y, velocity;
	
	for(int i=0;i<numberOfNodes_xDirection;i++)
	{
		velocity_x = (u[0][i] + u[0][i+1])/2.0;
		velocity_y = (v[0][i] + v[1][i])/2.0;
		velocity = sqrt(velocity_x*velocity_x + velocity_y*velocity_y);
		Yplus_bottomWall[i] = calculate_Yplus_forEachFace_usingNewtonRaphson(velocity, dy/2, rho[0][i], mu[0][i]);
	}
			
	for(int j=0;j<numberOfNodes_yDirection;j++)
	{
		velocity_x = (u[j][0] + u[j][1])/2.0;
		velocity_y = (v[j][0] + v[j+1][0])/2.0;
		velocity = sqrt(velocity_x*velocity_x + velocity_y*velocity_y);
		Yplus_leftWall[j] = calculate_Yplus_forEachFace_usingNewtonRaphson(velocity, dx/2, rho[j][0], mu[j][0]);	
	}
	
	for(int j=0;j<numberOfNodes_yDirection;j++)
	{
		velocity_x = (u[j][numberOfNodes_xDirection-1] + u[j][numberOfNodes_xDirection])/2.0;
		velocity_y = (v[j][numberOfNodes_xDirection-1] + v[j+1][numberOfNodes_xDirection-1])/2.0;
		velocity = sqrt(velocity_x*velocity_x + velocity_y*velocity_y);
		Yplus_rightWall[j] = calculate_Yplus_forEachFace_usingNewtonRaphson(velocity, dx/2, rho[j][numberOfNodes_xDirection-1], mu[j][numberOfNodes_xDirection-1]);
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void wallFunction::calculate_shearStress_for_allBoundaryWalls()
{
	int i,j; 
	
	calculate_wallsKinematicViscosity();

	for (i=0;i<numberOfNodes_xDirection;i++)
		tau_bottomWall[i] = rho[0][i]*nu_bottomWall[i]*(u[0][i] + u[0][i+1])/dy; 
		
	for (j=0;j<numberOfNodes_yDirection;j++)
	{
		tau_leftWall[j] = rho[j][0]*nu_leftWall[j]*(v[j][0] + v[j+1][0])/dx; 
		tau_rightWall[j] = rho[j][numberOfNodes_xDirection-1]*nu_rightWall[j]*(v[j][numberOfNodes_xDirection-1] + v[j+1][numberOfNodes_xDirection-1])/dx; 
	}
	
	average_tau_bottomWall = 0;
	
	for (i=0;i<numberOfNodes_xDirection;i++)
		average_tau_bottomWall = average_tau_bottomWall + tau_bottomWall[i]*dx;
		
	average_tau_bottomWall = average_tau_bottomWall/Length;
	
	
	average_tau_leftWall = 0;
	average_tau_rightWall = 0;
	
	for (j=0;j<numberOfNodes_yDirection;j++)
	{
		average_tau_leftWall = average_tau_leftWall + tau_leftWall[j]*dy;
		average_tau_rightWall = average_tau_rightWall + tau_rightWall[j]*dy;
	}
		
	average_tau_leftWall = average_tau_leftWall/Length;
	average_tau_rightWall = average_tau_rightWall/Length;
	
	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

wallFunction* wallFunction::New(const char *path, double **U, double **V)
{
	string wallFunctionType_Name, flowModel_Name_, usingWallFunction_Flag_;

	ifstream file(path);

	string searchForFlowModelType_Name_ = "FlowModelType";
	string searchForUsingWallFunction_Flag_ = "UsingWallFunction";
	string searchForWallFunctionType_Name = "WallFunctionName";

	string lineOfText;

	for(;;)
	{
		getline(file, lineOfText);

		if (file.eof()) break;

		if (lineOfText.find(searchForFlowModelType_Name_, 0) != string::npos)
		  file >> flowModel_Name_;
		  
		if (lineOfText.find(searchForUsingWallFunction_Flag_, 0) != string::npos)
		  file >> usingWallFunction_Flag_;
		  
		if (lineOfText.find(searchForWallFunctionType_Name, 0) != string::npos)
		  file >> wallFunctionType_Name;
    }

    file.close();
    
    wallFunction* myObject;
    
    if (usingWallFunction_Flag_ == "False" || flowModel_Name_ == "laminar")
		myObject = new notWallFunction(path, U, V);
		
	else
	{
		if (wallFunctionType_Name=="StandardWallFunction")
			myObject = new standardWallFunction(path, U, V);
			
		if (wallFunctionType_Name=="EnhancedWallFunction")
			myObject = new enhancedWallFunction(path, U, V);
	}
	    
	return myObject;
}
