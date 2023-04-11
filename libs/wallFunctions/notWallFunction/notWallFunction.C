#include "notWallFunction.H"

//////////////// constructor ///////////////////////

notWallFunction::notWallFunction (const char *path, double **U, double **V): wallFunction(path, U, V)
{print_wallFunction_information();}

	
//////////// public functions /////////////////////////////////////////////////////////////////////////

void notWallFunction::print_wallFunction_information()
{
	if (flowModel == "laminar")
		cout <<"No wall function is used to calculate walls shear stress, since flow type is laminar." << endl << endl << endl;
		
	else if (usingWallFunction_flag == "False")
		cout <<"No wall function is used to calculate walls shear stress." << endl << endl << endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void notWallFunction::calculate_wallsKinematicViscosity()
{
	int i,j; 
	
	for (i=0;i<numberOfNodes_xDirection;i++)
			nu_bottomWall[i] = mu[0][i]/rho[0][i]; 
		
	for (j=0;j<numberOfNodes_yDirection;j++)
	{
		nu_leftWall[j] = mu[j][0]/rho[j][0];
		nu_rightWall[j] = mu[j][numberOfNodes_xDirection-1]/rho[j][numberOfNodes_xDirection-1];
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void notWallFunction::write_solutionField()
{
	wallFunction::write_solutionField();
}
