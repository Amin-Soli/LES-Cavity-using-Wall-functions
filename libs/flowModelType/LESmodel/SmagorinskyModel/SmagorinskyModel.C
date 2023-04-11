#include "SmagorinskyModel.H"

///////////////////////////////////////////////// constructor ///////////////////////////////////////////////////////////////////

SmagorinskyModel::SmagorinskyModel (const char *path, double **U, double **V , double **mu_SGS): LESmodel(path, U, V, mu_SGS)
{
	
	readInputData(path);
	
	l_Smagorinsky = Smagorinsky_coeff*sqrt(dx*dy);
	  
	printInputData();
	
}
	
//////////// public functions /////////////////////////////////////////////////////////////////////////

void SmagorinskyModel::readInputData(const char *path)
{
	ifstream file(path);
	
	string searchForSmagorinsky_coeff = "SmagorinskyCoeff";

	string lineOfText;

	for(;;)
	{
		getline(file, lineOfText);

		if (file.eof()) break;


		if (lineOfText.find(searchForSmagorinsky_coeff, 0) != string::npos)
		  file >> Smagorinsky_coeff;

    }

    file.close();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SmagorinskyModel::printInputData()
{
	cout << "Flow model type is: " << modelFlow_name << endl;
	cout << "LES model is: " << turbulentModel_name << endl;
	cout << "Smagorinsky Coefficient is: " << Smagorinsky_coeff << endl << endl;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SmagorinskyModel::calculate_delta_subgrid()
{
	for(int j=0;j<numberOfNodes_yDirection;j++)
		for(int i=0;i<numberOfNodes_xDirection;i++)
			delta_sgs[j][i] = l_Smagorinsky;
}
