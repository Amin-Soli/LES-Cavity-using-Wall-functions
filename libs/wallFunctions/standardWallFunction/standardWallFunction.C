#include "standardWallFunction.H"

//////////////// constructor ///////////////////////

standardWallFunction::standardWallFunction (const char *path, double **U, double **V): wallFunction(path, U, V)
{print_wallFunction_information();}

	
//////////// public functions /////////////////////////////////////////////////////////////////////////

void standardWallFunction::print_wallFunction_information()
{
	cout << wallFunction_name << " is used to calculate walls shear stress." << endl << endl << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void standardWallFunction::calculate_wallsKinematicViscosity()
{
	int i,j; 
	double yPlus;
	
	calculate_Yplus_for_allBoundaryWalls();	
	
	// bottom wall
	for(i=0;i<numberOfNodes_xDirection;i++)
	{
		yPlus = Yplus_bottomWall[i];
		if ( yPlus < Yplus_laminar)
			nu_bottomWall[i] = nu[0][i];
			
		else 
			nu_bottomWall[i] = nu[0][i]*(yPlus/(1.0/kappa *log(E*yPlus)));
	}
	
	// left wall
	for(j=0;j<numberOfNodes_yDirection;j++)
	{
		yPlus = Yplus_leftWall[j];
		if ( yPlus < Yplus_laminar)
			nu_leftWall[j] = nu[j][0];
			
		else 
			nu_leftWall[j] = nu[j][0]*(yPlus/(1.0/kappa *log(E*yPlus)));
	}
	
	// right wall
	for(j=0;j<numberOfNodes_yDirection;j++)
	{
		yPlus = Yplus_rightWall[j];
		if ( yPlus < Yplus_laminar)
			nu_rightWall[j] = nu[j][numberOfNodes_xDirection-1];
			
		else 
			nu_rightWall[j] = nu[j][numberOfNodes_xDirection-1]*(yPlus/(1.0/kappa *log(E*yPlus)));
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void standardWallFunction::write_solutionField()
{
	int i,j;
	
	wallFunction::write_solutionField();
	
	stringstream ss;
	ss << time;
	string str = ss.str();
	string myName = "results/t=" + str ;	

    ofstream file1;
    string yPlus_bottomWallPath = myName + "/yPlus_bottomWall.txt";
    const char* yPlus_bottomWallPathName = yPlus_bottomWallPath.c_str();	
    file1.open(yPlus_bottomWallPathName);

    file1 << fixed << setprecision(6);
    
    for (i=0;i<numberOfNodes_xDirection;i++)
		file1 << i*dx + dx/2.0 << ", " << Yplus_bottomWall[i] << endl;

    file1.close();
    
    ofstream file2;
    string yPlus_leftWallPath = myName + "/yPlus_leftWall.txt";
    const char* yPlus_leftWallPathName = yPlus_leftWallPath.c_str();	
    file2.open(yPlus_leftWallPathName);

    file2 << fixed << setprecision(6);
    
    for (j=0;j<numberOfNodes_yDirection;j++)
		file2 << j*dy + dy/2.0 << ", " << Yplus_leftWall[j] << endl;

    file2.close();
    
    ofstream file3;
    string yPlus_rightWallPath = myName + "/yPlus_rightWall.txt";
    const char* yPlus_rightWallPathName = yPlus_rightWallPath.c_str();	
    file3.open(yPlus_rightWallPathName);

    file3 << fixed << setprecision(6);
    
    for (j=0;j<numberOfNodes_yDirection;j++)
		file3 << j*dy + dy/2.0 << ", " << Yplus_rightWall[j] << endl;

    file3.close();
}
