#include "LESmodel.H"
#include "../SmagorinskyModel/SmagorinskyModel.C"
#include "../VanDriestDampingModel/VanDriestDampingModel.C"

//////////////////////////////////////////// constructor //////////////////////////////////////////////////////////////////////

LESmodel::LESmodel (const char *path, double **U, double **V , double **mu_SGS)
{
	
	readInputData(path);
	
	u = U;
	v = V;
	mu_sgs = mu_SGS;
	
	time = start_time;
	
	dx = Length/numberOfNodes_xDirection;
	dy = Length/numberOfNodes_yDirection;
	
	count_for_writingSolutionField = 0;
	
	rho = new double *[numberOfNodes_yDirection];
	delta_sgs = new double *[numberOfNodes_yDirection];
		
    for(int j = 0; j <numberOfNodes_yDirection; j++)
    {
		rho[j] = new double[numberOfNodes_xDirection];
		delta_sgs[j] = new double[numberOfNodes_xDirection];
	}
	
	u_t = new double [numberOfNodes_xDirection+1];
    u_b = new double [numberOfNodes_xDirection+1];
    
    for(int i=0;i<numberOfNodes_xDirection+1;i++)
    {
        u_b[i] = 0 ;
        u_t[i] = free_velocity ;
    }
    
    v_l = new double [numberOfNodes_yDirection+1];
    v_r = new double [numberOfNodes_yDirection+1];
    
    for(int j=0;j<numberOfNodes_yDirection+1;j++)
    {
        v_l[j] = 0 ;
        v_r[j] = 0 ;
    }

    for(int j=0;j<numberOfNodes_yDirection;j++)
        for(int i=0;i<numberOfNodes_xDirection;i++)
            rho[j][i] = rho_;
	
}

//////////////////////////////////////////////// destructor  ///////////////////////////////////////////////

LESmodel::~LESmodel ()
{

	for(int j = 0; j < numberOfNodes_yDirection; ++j)
	{
         delete[] rho[j];
         delete[] delta_sgs[j];
	}

     delete[] rho;
     delete[] delta_sgs;
     
     delete[] u_b;
     delete[] u_t;
     
     delete[] v_l;
     delete[] v_r;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void LESmodel::update_time()
{
	check_writing_solutionField();
	time = time + dt;
	count_for_writingSolutionField++;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void LESmodel::check_writing_solutionField()
{
	if(count_for_writingSolutionField % numberOfDt_for_writingSolutionField == 0)
		write_solutionField();	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void LESmodel::write_solutionField()
{
	int i,j;
	
	stringstream ss;
	ss << time;
	string str = ss.str();
	string myName = "results/t=" + str ;	

    ofstream file1;
    string mu_SGSpath = myName + "/Mu_turbulence.txt";
    const char* mu_SGSpathName = mu_SGSpath.c_str();	
    file1.open(mu_SGSpathName);

    file1 << fixed << setprecision(6);

    for (i=0;i<numberOfNodes_xDirection;i++)
        file1 << i*dx + dx/2.0 << ", " << 0*dy << ", " << mu_sgs[0][i] << endl;

    for(j=0;j<numberOfNodes_yDirection;j++)
        for (i=0;i<numberOfNodes_xDirection;i++)
            file1 << i*dx + dx/2.0 << ", " << j*dy + dy/2.0 << ", " << mu_sgs[j][i] << endl;

    for (i=0;i<numberOfNodes_xDirection;i++)
        file1 << i*dx + dx/2.0 << ", " << numberOfNodes_yDirection*dy << ", " << mu_sgs[numberOfNodes_yDirection-1][i] << endl;

    for (j=0;j<numberOfNodes_yDirection;j++)
        file1 << 0*dx << ", " << j*dy + dy/2.0 << ", " << mu_sgs[j][0] << endl;

    for (j=0;j<numberOfNodes_yDirection;j++)
        file1 << numberOfNodes_xDirection*dx << ", " << j*dy + dy/2.0 << ", " << mu_sgs[j][numberOfNodes_xDirection-1] << endl;

    file1 << 0*dx << ", " << 0*dy << ", " << mu_sgs[0][0]/rho[0][0] << endl;
    file1 << 0*dx << ", " << numberOfNodes_yDirection*dy << ", " << mu_sgs[numberOfNodes_yDirection-1][0] << endl;
    file1 << numberOfNodes_xDirection*dx << ", " << 0*dy << ", " << mu_sgs[0][numberOfNodes_xDirection-1] << endl;
    file1 << numberOfNodes_xDirection*dx << ", " << numberOfNodes_yDirection*dy << ", " << mu_sgs[numberOfNodes_yDirection-1][numberOfNodes_xDirection-1] << endl;

    file1.close();

}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void LESmodel::readInputData(const char *path)
{
	ifstream file(path);

	string searchForLength = "Length";
	string searchForStartTime = "Startime";
	string searchForNumberOfDt_for_writingSolutionFields = "NumberOfDtForWritingSolutionFields";
	string searchForDelta_t = "Delta_t";
	string searchForDensity = "Density";
	string searchForMovingWallVelocity = "MovingWallVelocity";
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
		  
		if (lineOfText.find(searchForMovingWallVelocity, 0) != string::npos)
		  file >> free_velocity;
		  
		if (lineOfText.find(searchForNumberOfDt_for_writingSolutionFields, 0) != string::npos)
		  file >> numberOfDt_for_writingSolutionField;
		  
		if (lineOfText.find(searchForDelta_t, 0) != string::npos)
		  file >> dt;

		if (lineOfText.find(searchForDensity, 0) != string::npos)
		  file >> rho_;

		if (lineOfText.find(searchForNumberOfNodesInXDirection, 0) != string::npos)
		  file >> numberOfNodes_xDirection;

		if (lineOfText.find(searchForNumberOfNodesInYDirection, 0) != string::npos)
		  file >> numberOfNodes_yDirection;

    }

    file.close();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void LESmodel::update_mu_turbulence()
{
	
	calculate_delta_subgrid();
	
    double s11, s22, s12, p_k;
	int i,j;
	    
	j=0;
	
	i=0;
	
	s11 = (u[j][i+1] - u[j][i])/dx;
	s22 = (v[j+1][i] - v[j][i])/dy;
	s12 = 0.5*( (u[j][i+1] + u[j][i] + u[j+1][i+1] + u[j+1][i])/4.0 - (u_b[i+1] + u_b[i])/2.0 )/dy
		+ 0.5*( (v[j+1][i] + v[j][i] + v[j+1][i+1] + v[j][i+1])/4.0 - (v_l[j+1] + v_l[j])/2.0 )/dx ;

	p_k = s11*s11 + s22*s22 + 2*s12*s12;

	mu_sgs[j][i] = rho[j][i]*delta_sgs[j][i]*delta_sgs[j][i]*sqrt(2*p_k);
			
	for(i=1;i<numberOfNodes_xDirection-1;i++)
	{				
		s11 = (u[j][i+1] - u[j][i])/dx;
		s22 = (v[j+1][i] - v[j][i])/dy;
		s12 = 0.5*( (u[j][i+1] + u[j][i] + u[j+1][i+1] + u[j+1][i])/4.0 - (u_b[i+1] + u_b[i])/2.0 )/dy
			+ 0.5*( (v[j+1][i] + v[j][i] + v[j+1][i+1] + v[j][i+1])/4.0 - (v[j+1][i] + v[j][i] + v[j+1][i-1] + v[j][i-1])/4.0 )/dx ;

		p_k = s11*s11 + s22*s22 + 2*s12*s12;

		mu_sgs[j][i] = rho[j][i]*delta_sgs[j][i]*delta_sgs[j][i]*sqrt(2*p_k);
	}
		
	i=numberOfNodes_xDirection-1;
	
	s11 = (u[j][i+1] - u[j][i])/dx;
	s22 = (v[j+1][i] - v[j][i])/dy;
	s12 = 0.5*( (u[j][i+1] + u[j][i] + u[j+1][i+1] + u[j+1][i])/4.0 - (u_b[i+1] + u_b[i])/2.0 )/dy
		+ 0.5*( (v_r[j+1] + v_r[j])/2.0 - (v[j+1][i] + v[j][i] + v[j+1][i-1] + v[j][i-1])/4.0 )/dx ;

	p_k = s11*s11 + s22*s22 + 2*s12*s12;

	mu_sgs[j][i] = rho[j][i]*delta_sgs[j][i]*delta_sgs[j][i]*sqrt(2*p_k);
			

	for(j=1;j<numberOfNodes_yDirection-1;j++)
	{
		i=0;
		
		s11 = (u[j][i+1] - u[j][i])/dx;
		s22 = (v[j+1][i] - v[j][i])/dy;
		s12 = 0.5*( (u[j][i+1] + u[j][i] + u[j+1][i+1] + u[j+1][i])/4.0 - (u[j][i+1] + u[j][i] + u[j-1][i+1] + u[j-1][i])/4.0 )/dy
			+ 0.5*( (v[j+1][i] + v[j][i] + v[j+1][i+1] + v[j][i+1])/4.0 - (v_l[j+1] + v_l[j])/2.0 )/dx ;

		p_k = s11*s11 + s22*s22 + 2*s12*s12;

		mu_sgs[j][i] = rho[j][i]*delta_sgs[j][i]*delta_sgs[j][i]*sqrt(2*p_k);
				
		for(i=1;i<numberOfNodes_xDirection-1;i++)
		{				
			s11 = (u[j][i+1] - u[j][i])/dx;
			s22 = (v[j+1][i] - v[j][i])/dy;
			s12 = 0.5*( (u[j][i+1] + u[j][i] + u[j+1][i+1] + u[j+1][i])/4.0 - (u[j][i+1] + u[j][i] + u[j-1][i+1] + u[j-1][i])/4.0 )/dy
				+ 0.5*( (v[j+1][i] + v[j][i] + v[j+1][i+1] + v[j][i+1])/4.0 - (v[j+1][i] + v[j][i] + v[j+1][i-1] + v[j][i-1])/4.0 )/dx ;

			p_k = s11*s11 + s22*s22 + 2*s12*s12;

			mu_sgs[j][i] = rho[j][i]*delta_sgs[j][i]*delta_sgs[j][i]*sqrt(2*p_k);
		}
			
		i=numberOfNodes_xDirection-1;
		
		s11 = (u[j][i+1] - u[j][i])/dx;
		s22 = (v[j+1][i] - v[j][i])/dy;
		s12 = 0.5*( (u[j][i+1] + u[j][i] + u[j+1][i+1] + u[j+1][i])/4.0 - (u[j][i+1] + u[j][i] + u[j-1][i+1] + u[j-1][i])/4.0 )/dy
			+ 0.5*( (v_r[j+1] + v_r[j])/2.0 - (v[j+1][i] + v[j][i] + v[j+1][i-1] + v[j][i-1])/4.0 )/dx ;

		p_k = s11*s11 + s22*s22 + 2*s12*s12;

		mu_sgs[j][i] = rho[j][i]*delta_sgs[j][i]*delta_sgs[j][i]*sqrt(2*p_k);
					
	}
	
	j=numberOfNodes_yDirection-1;
	
	i=0;
	
	s11 = (u[j][i+1] - u[j][i])/dx;
	s22 = (v[j+1][i] - v[j][i])/dy;
	s12 = 0.5*( (u_t[i+1] + u_t[i])/2.0 - (u[j][i+1] + u[j][i] + u[j-1][i+1] + u[j-1][i])/4.0 )/dy
		+ 0.5*( (v[j+1][i] + v[j][i] + v[j+1][i+1] + v[j][i+1])/4.0 - (v_l[j+1] + v_l[j])/2.0 )/dx ;

	p_k = s11*s11 + s22*s22 + 2*s12*s12;

	mu_sgs[j][i] = rho[j][i]*delta_sgs[j][i]*delta_sgs[j][i]*sqrt(2*p_k);
			
	for(i=1;i<numberOfNodes_xDirection-1;i++)
	{				
		s11 = (u[j][i+1] - u[j][i])/dx;
		s22 = (v[j+1][i] - v[j][i])/dy;
		s12 = 0.5*( (u_t[i+1] + u_t[i])/2.0 - (u[j][i+1] + u[j][i] + u[j-1][i+1] + u[j-1][i])/4.0 )/dy
			+ 0.5*( (v[j+1][i] + v[j][i] + v[j+1][i+1] + v[j][i+1])/4.0 - (v[j+1][i] + v[j][i] + v[j+1][i-1] + v[j][i-1])/4.0 )/dx ;

		p_k = s11*s11 + s22*s22 + 2*s12*s12;

		mu_sgs[j][i] = rho[j][i]*delta_sgs[j][i]*delta_sgs[j][i]*sqrt(2*p_k);
	}
		
	i=numberOfNodes_xDirection-1;
	
	s11 = (u[j][i+1] - u[j][i])/dx;
	s22 = (v[j+1][i] - v[j][i])/dy;
	s12 = 0.5*( (u_t[i+1] + u_t[i])/2.0 - (u[j][i+1] + u[j][i] + u[j-1][i+1] + u[j-1][i])/4.0 )/dy
		+ 0.5*( (v_r[j+1] + v_r[j])/2.0 - (v[j+1][i] + v[j][i] + v[j+1][i-1] + v[j][i-1])/4.0 )/dx ;

	p_k = s11*s11 + s22*s22 + 2*s12*s12;

	mu_sgs[j][i] = rho[j][i]*delta_sgs[j][i]*delta_sgs[j][i]*sqrt(2*p_k);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

LESmodel* LESmodel::New(const char *path, double **U, double **V , double **mu_SGS)
{
	string LESmodel_Name;

	ifstream file(path);

	string searchForLESmodel_Name = "LESmodelName";

	string lineOfText;

	for(;;)
	{
		getline(file, lineOfText);

		if (file.eof()) break;

		if (lineOfText.find(searchForLESmodel_Name, 0) != string::npos)
		  file >> LESmodel_Name;
    }

    file.close();
    
    LESmodel* myObject;

    if (LESmodel_Name=="Smagorinsky")
	    myObject = new SmagorinskyModel(path, U, V, mu_SGS);
	    
	if (LESmodel_Name=="SmagorinskyWithVanDriestDamping")
	    myObject = new VanDriestDampingModel(path, U, V, mu_SGS);
	    
	return myObject;
}

