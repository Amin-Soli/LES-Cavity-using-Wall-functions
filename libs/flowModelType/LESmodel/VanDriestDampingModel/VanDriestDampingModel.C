#include "VanDriestDampingModel.H"

//////////////// constructor ///////////////////////

VanDriestDampingModel::VanDriestDampingModel (const char *path, double **U, double **V , double **mu_SGS): LESmodel(path, U, V, mu_SGS)
{
	
	readInputData(path);
	
	l_Smagorinsky = Smagorinsky_coeff*sqrt(dx*dy);
	
	Yplus_bottomWall = new double [numberOfNodes_xDirection];
	Yplus_leftWall = new double [numberOfNodes_yDirection];
	Yplus_rightWall = new double [numberOfNodes_yDirection];
	
	mu = new double *[numberOfNodes_yDirection];
	mixingLength = new double *[numberOfNodes_yDirection];
	yFromBottomWall = new double *[numberOfNodes_yDirection];
	yFromLeftWall = new double *[numberOfNodes_yDirection];
	yFromRightWall = new double *[numberOfNodes_yDirection];
	min_y = new double *[numberOfNodes_yDirection];
	Yplus = new double *[numberOfNodes_yDirection];
		
    for(int j = 0; j <numberOfNodes_yDirection; j++)
    {
		mu[j] = new double[numberOfNodes_xDirection];
		mixingLength[j] = new double[numberOfNodes_xDirection];
		yFromBottomWall[j] = new double[numberOfNodes_xDirection];
		yFromLeftWall[j] = new double[numberOfNodes_xDirection];
		yFromRightWall[j] = new double[numberOfNodes_xDirection];
		min_y[j] = new double[numberOfNodes_xDirection];
		Yplus[j] = new double[numberOfNodes_xDirection];
	}
	
	for(int j=0;j<numberOfNodes_yDirection;j++)
        for(int i=0;i<numberOfNodes_xDirection;i++)
            mu[j][i] = mu_;
	
	
	printInputData();
	
}

//////////////// destructor  ///////////////////////

VanDriestDampingModel::~VanDriestDampingModel ()
{
	
	for(int j = 0; j < numberOfNodes_yDirection; ++j)
	{
		 delete[] mu[j];
         delete[] mixingLength[j];
         delete[] yFromBottomWall[j];
         delete[] yFromLeftWall[j];
         delete[] yFromRightWall[j];
         delete[] min_y[j];
         delete[] Yplus[j];
	}
	
	delete[] mu;
    delete[] mixingLength;
    delete[] yFromBottomWall;
    delete[] yFromLeftWall;
    delete[] yFromRightWall;
    delete[] min_y;
    delete[] Yplus;
    
    delete[] Yplus_bottomWall;   
    delete[] Yplus_leftWall;
    delete[] Yplus_rightWall;
	
}
	
//////////// public functions /////////////////////////////////////////////////////////////////////////

void VanDriestDampingModel::readInputData(const char *path)
{
	ifstream file(path);
	
	string searchForDynamicViscosity = "DynamicViscosity";
	string searchForSmagorinsky_coeff = "SmagorinskyCoeff";
	string searchForKappa = "Kappa";
	string searchForE = "E_";
	string searchForA_plus = "A_plus";

	string lineOfText;

	for(;;)
	{
		getline(file, lineOfText);

		if (file.eof()) break;

		if (lineOfText.find(searchForDynamicViscosity, 0) != string::npos)
		  file >> mu_;
		
		if (lineOfText.find(searchForSmagorinsky_coeff, 0) != string::npos)
		  file >> Smagorinsky_coeff;
		  
		if (lineOfText.find(searchForKappa, 0) != string::npos)
		  file >> kappa;
		  
		if (lineOfText.find(searchForE, 0) != string::npos)
		  file >> E;
		  
		if (lineOfText.find(searchForA_plus, 0) != string::npos)
		  file >> A_plus;

    }

    file.close();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VanDriestDampingModel::printInputData()
{
	cout << "Flow model type is: " << modelFlow_name << endl;
	cout << "LES model is: " << turbulentModel_name << endl;
	cout << "Smagorinsky Coefficient is: " << Smagorinsky_coeff << endl;
	cout << "Kappa is: " << kappa << endl;
	cout << "E is: " << E << endl;
	cout << "A_plus is: " << A_plus << endl << endl;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double VanDriestDampingModel::calculate_Yplus_forEachFace_usingNewtonRaphson(double U, double d, double rho, double mu)
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

void VanDriestDampingModel::calculate_Yplus_for_allBoundaryWalls()
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VanDriestDampingModel::calculate_mixingLengths_from_allBoundaryWalls()
{
	int i,j;
	
	// calculate all cells distance from bottom wall	    
	for(j=0;j<numberOfNodes_yDirection;j++)
		for(i=0;i<numberOfNodes_xDirection;i++)
			yFromBottomWall[j][i] = j*dy + dy/2.0 ;
	
	// calculate all cells distance from left wall	    
	for(j=0;j<numberOfNodes_yDirection;j++)
		for(i=0;i<numberOfNodes_xDirection;i++)
			yFromLeftWall[j][i] = i*dx + dx/2.0 ;
			
	// calculate all cells distance from right wall	    
	for(j=0;j<numberOfNodes_yDirection;j++)
		for(i=0;i<numberOfNodes_xDirection;i++)
			yFromRightWall[j][i] = (numberOfNodes_xDirection-1-i)*dx + dx/2.0 ;
			
	
	double min_wallDistance, y_plus;
	
	for(j=0;j<numberOfNodes_yDirection;j++)
		for(i=0;i<numberOfNodes_xDirection;i++)
		{
			min_wallDistance = yFromBottomWall[j][i];
			y_plus = Yplus_bottomWall[i];
			if (min_wallDistance > yFromLeftWall[j][i])
			{
				min_wallDistance = yFromLeftWall[j][i];
				y_plus = Yplus_leftWall[j];
			}
			if (min_wallDistance > yFromRightWall[j][i])
			{
				min_wallDistance = yFromRightWall[j][i];
				y_plus = Yplus_rightWall[j];
			}
			
			min_y[j][i] = min_wallDistance;
			Yplus[j][i] = y_plus;		
		}
	
	
	double VanDriestDamping_func;
	
	// calculate mixing length     
	for(j=0;j<numberOfNodes_yDirection;j++)
		for(i=0;i<numberOfNodes_xDirection;i++)
		{
			VanDriestDamping_func = (1.0 - exp(-Yplus[j][i]/A_plus));
			mixingLength[j][i] = kappa*min_y[j][i]*VanDriestDamping_func;				
		}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VanDriestDampingModel::calculate_delta_subgrid()
{	
	calculate_Yplus_for_allBoundaryWalls();
	
	calculate_mixingLengths_from_allBoundaryWalls();
	
	// find minimum between Smagorinsky length and mixing lengths
	for(int j=0;j<numberOfNodes_yDirection;j++)
		for(int i=0;i<numberOfNodes_xDirection;i++)
			delta_sgs[j][i] = min(l_Smagorinsky, mixingLength[j][i]);

}
