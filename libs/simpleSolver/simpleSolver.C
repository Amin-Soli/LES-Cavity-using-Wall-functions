#include "simpleSolver.H"

//////////////// constructor //////////////////////////////////////////////////////////////////////

simpleSolver::simpleSolver (const char *path, double **U, double **V , double **P, double **mu_turb)
{
	
	readInputData(path);
	
	ReynoldsNumber = rho_*free_velocity*Length/mu_;
	
	u = U;
	v = V;
	p = P;
	mu_turbulence = mu_turb;
	
	time = start_time;
	
	dx = Length/numberOfNodes_xDirection;
	dy = Length/numberOfNodes_yDirection;
	
	simple_iteration = 0;
	count_for_writingSolutionFields = 0;
	
	u0 = new double *[numberOfNodes_yDirection];
	v0 = new double *[numberOfNodes_yDirection+1];
    pc = new double *[numberOfNodes_yDirection];
    rho = new double *[numberOfNodes_yDirection];
    mu = new double *[numberOfNodes_yDirection];
    
    for(int j = 0; j <numberOfNodes_yDirection+1; j++)
		v0[j] = new double[numberOfNodes_xDirection];
		
    for(int j = 0; j <numberOfNodes_yDirection; j++)
    {
		u0[j] = new double[numberOfNodes_xDirection+1];
		pc[j] = new double[numberOfNodes_xDirection];
		rho[j] = new double[numberOfNodes_xDirection];
		mu[j] = new double[numberOfNodes_xDirection];
	}
	
	Ao_X = new double *[numberOfNodes_yDirection];
    A0o_X = new double *[numberOfNodes_yDirection];
    Ae_X = new double *[numberOfNodes_yDirection];
    Aw_X = new double *[numberOfNodes_yDirection];
    An_X = new double *[numberOfNodes_yDirection];
    As_X = new double *[numberOfNodes_yDirection];
    So_X = new double *[numberOfNodes_yDirection];
    
    for(int j = 0; j <numberOfNodes_yDirection; j++)
    {
		Ao_X[j] = new double[numberOfNodes_xDirection];
		A0o_X[j] = new double[numberOfNodes_xDirection];
		Ae_X[j] = new double[numberOfNodes_xDirection];
		Aw_X[j] = new double[numberOfNodes_xDirection];
		An_X[j] = new double[numberOfNodes_xDirection];
	    As_X[j] = new double[numberOfNodes_xDirection];
		So_X[j] = new double[numberOfNodes_xDirection];
	}
	
	Ao_Y = new double *[numberOfNodes_yDirection];
    A0o_Y = new double *[numberOfNodes_yDirection];
    Ae_Y = new double *[numberOfNodes_yDirection];
    Aw_Y = new double *[numberOfNodes_yDirection];
    An_Y = new double *[numberOfNodes_yDirection];
    As_Y = new double *[numberOfNodes_yDirection];
    So_Y = new double *[numberOfNodes_yDirection];
    
    for(int j = 0; j <numberOfNodes_yDirection; j++)
    {
		Ao_Y[j] = new double[numberOfNodes_xDirection];
		A0o_Y[j] = new double[numberOfNodes_xDirection];
		Ae_Y[j] = new double[numberOfNodes_xDirection];
		Aw_Y[j] = new double[numberOfNodes_xDirection];
		An_Y[j] = new double[numberOfNodes_xDirection];
	    As_Y[j] = new double[numberOfNodes_xDirection];
		So_Y[j] = new double[numberOfNodes_xDirection];
	}
	
	Ao_P = new double *[numberOfNodes_yDirection];
    Ae_P = new double *[numberOfNodes_yDirection];
    Aw_P = new double *[numberOfNodes_yDirection];
    An_P = new double *[numberOfNodes_yDirection];
    As_P = new double *[numberOfNodes_yDirection];
    So_P = new double *[numberOfNodes_yDirection];
    
    for(int j = 0; j <numberOfNodes_yDirection; j++)
    {
		Ao_P[j] = new double[numberOfNodes_xDirection];
		Ae_P[j] = new double[numberOfNodes_xDirection];
		Aw_P[j] = new double[numberOfNodes_xDirection];
		An_P[j] = new double[numberOfNodes_xDirection];
	    As_P[j] = new double[numberOfNodes_xDirection];
		So_P[j] = new double[numberOfNodes_xDirection];
	}
	
	u_l = new double [numberOfNodes_yDirection];
    u_r = new double [numberOfNodes_yDirection];
    u_t = new double [numberOfNodes_xDirection+1];
    u_b = new double [numberOfNodes_xDirection+1];
    
    //input boundary condition for u
	
    for(int j=0;j<numberOfNodes_yDirection;j++)
    {
        u_l[j] = 0 ;
        u_r[j] = 0 ;
    }

    for(int i=0;i<numberOfNodes_xDirection+1;i++)
    {
        u_b[i] = 0 ;
        u_t[i] = free_velocity ;
    }
    
    v_l = new double [numberOfNodes_yDirection+1];
    v_r = new double [numberOfNodes_yDirection+1];
    v_t = new double [numberOfNodes_xDirection];
    v_b = new double [numberOfNodes_xDirection];
    
    //input boundary condition for v
    
    for(int i=0;i<numberOfNodes_xDirection;i++)
    {
        v_b[i] = 0 ;
        v_t[i] = 0 ;
    }

    for(int j=0;j<numberOfNodes_yDirection+1;j++)
    {
        v_l[j] = 0 ;
        v_r[j] = 0 ;
    }
    
    for(int j=0;j<numberOfNodes_yDirection;j++)
        for(int i=0;i<numberOfNodes_xDirection;i++)
        {
            rho[j][i] = rho_;
            mu[j][i] = mu_;
            mu_turbulence[j][i] = 0;
        }

	printInputData();
	
	if (start_time == 0)
	{
		// create a folder for writing results:
		
		#ifdef _WIN32

		if (opendir("results"))
			system("rmdir /Q /S results");

		if (mkdir("results") != -1)
			cout << "Directory of results was created." << endl;

		#endif
		
		#ifdef linux
		
		DIR *dir = opendir("results");
		
		if (dir)
			system("rm -rf results");

		closedir(dir);

		if (mkdir("results", 0777) != -1)
			cout << "Directory of results was created." << endl << endl;
			
		#endif
		
	}
}

//////////////////////////////////   destructor  ///////////////////////////////////////////////////

simpleSolver::~simpleSolver ()
{
	for(int j = 0; j <numberOfNodes_yDirection+1; j++)
		delete[] v0[j];
		
	for(int j = 0; j < numberOfNodes_yDirection; ++j)
	{
		 delete[] u0[j];
         delete[] pc[j];
         delete[] rho[j];
         delete[] mu[j];
	}

	 delete[] u0;
     delete[] v0;
     delete[] pc;
     delete[] rho;
     delete[] mu;
     
	for(int j = 0; j < numberOfNodes_yDirection; ++j)
	{
         delete[] Ao_X[j];
         delete[] A0o_X[j];
         delete[] Ae_X[j];
         delete[] Aw_X[j];
         delete[] An_X[j];
         delete[] As_X[j];
         delete[] So_X[j];
	}

     delete[] Ao_X;
     delete[] A0o_X;
     delete[] Ae_X;
     delete[] Aw_X;
     delete[] An_X;
     delete[] As_X;
     delete[] So_X;
     
	for(int j = 0; j < numberOfNodes_yDirection; ++j)
	{
         delete[] Ao_Y[j];
         delete[] A0o_Y[j];
         delete[] Ae_Y[j];
         delete[] Aw_Y[j];
         delete[] An_Y[j];
         delete[] As_Y[j];
         delete[] So_Y[j];
	}

     delete[] Ao_Y;
     delete[] A0o_Y;
     delete[] Ae_Y;
     delete[] Aw_Y;
     delete[] An_Y;
     delete[] As_Y;
     delete[] So_Y;
     
	for(int j = 0; j < numberOfNodes_yDirection; ++j)
	{
         delete[] Ao_P[j];
         delete[] Ae_P[j];
         delete[] Aw_P[j];
         delete[] An_P[j];
         delete[] As_P[j];
         delete[] So_P[j];
	}

     delete[] Ao_P;
     delete[] Ae_P;
     delete[] Aw_P;
     delete[] An_P;
     delete[] As_P;
     delete[] So_P;

     delete[] u_l;
     delete[] u_r;
     delete[] u_b;
     delete[] u_t;
     
     delete[] v_l;
     delete[] v_r;
     delete[] v_b;
     delete[] v_t;
}

//////////////////////////////// public functions /////////////////////////////////////

void simpleSolver::setBoundaryConditions()
{    
    //// assign boundary conditions to field variables
 
    for (int j=0;j<numberOfNodes_yDirection;j++)
    {
        u[j][0] = u_l[j];
        u[j][numberOfNodes_xDirection] = u_r[j];
    }

    for (int i=0;i<numberOfNodes_xDirection;i++)
    {
        v[0][i] = v_b[i];
        v[numberOfNodes_yDirection][i] = v_t[i];
    }
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::calculate_linkCoeffs_for_XmomentumEq()
{
   double ac_e, ac_w, mu_e, mu_w, ac_n, ac_s, mu_n, mu_s, rho_p;
   
   int i,j;

   j=0;
   for (i=1;i<numberOfNodes_xDirection;i++)
   {
	   rho_p = (rho[j][i] + rho[j][i-1])/2.0;
	   ac_e = rho[j][i]*(u[j][i] + u[j][i+1])/2.0;
	   ac_w = rho[j][i-1]*(u[j][i] + u[j][i-1])/2.0;
	   mu_e = mu[j][i] + mu_turbulence[j][i];
	   mu_w = mu[j][i-1] + mu_turbulence[j][i-1];
	   ac_n = (rho[j][i] + rho[j][i-1] + rho[j+1][i] + rho[j+1][i-1])/4.0 * (v[j+1][i] + v[j+1][i-1])/2.0;
	   mu_n = (mu[j][i] + mu[j][i-1] + mu[j+1][i] + mu[j+1][i-1] + mu_turbulence[j][i] + mu_turbulence[j][i-1] + mu_turbulence[j+1][i] + mu_turbulence[j+1][i-1])/4.0;
	   ac_s = (rho[j][i] + rho[j][i-1])/2.0 * (v[j][i] + v[j][i-1])/2.0;
	   mu_s = (mu[j][i] + mu[j][i-1] + mu_turbulence[j][i] + mu_turbulence[j][i-1])/2.0;

	   A0o_X[j][i] = rho_p*dx*dy/dt;
	   Ae_X[j][i] = -dy*(fabs(ac_e) - ac_e)/2.0 - mu_e*dy/dx ;
	   Aw_X[j][i] = -dy*(fabs(ac_w) + ac_w)/2.0 - mu_w*dy/dx ;
	   An_X[j][i] = -dx*(fabs(ac_n) - ac_n)/2.0 - mu_n*dx/dy - mu_s*dx/(3.0*dy) ;
	   As_X[j][i] = 0.0  ;
	   Ao_X[j][i] = dy*(fabs(ac_e) + ac_e)/2.0 + mu_e*dy/dx + dy*(fabs(ac_w) - ac_w)/2.0 + mu_w*dy/dx
					  + dx*(fabs(ac_n) + ac_n)/2.0 + mu_n*dx/dy + 3.0*mu_s*dx/dy + A0o_X[j][i];

	   So_X[j][i] = (p[j][i-1]-p[j][i])*dy + A0o_X[j][i]*u0[j][i] + ac_s*u_b[i]*dx + 8.0*mu_s*u_b[i]*dx/(3.0*dy);
   }

   for(j=1;j<numberOfNodes_yDirection-1;j++)
	   for (i=1;i<numberOfNodes_xDirection;i++)
	   {
		   rho_p = (rho[j][i] + rho[j][i-1])/2.0;
		   ac_e = rho[j][i]*(u[j][i] + u[j][i+1])/2.0;
		   ac_w = rho[j][i-1]*(u[j][i] + u[j][i-1])/2.0;
		   mu_e = mu[j][i] + mu_turbulence[j][i];
		   mu_w = mu[j][i-1] + mu_turbulence[j][i-1];
		   ac_n = (rho[j][i]+ rho[j][i-1] + rho[j+1][i] + rho[j+1][i-1])/4.0 * (v[j+1][i] + v[j+1][i-1])/2.0;
		   mu_n = (mu[j][i]+ mu[j][i-1] + mu[j+1][i] + mu[j+1][i-1] + mu_turbulence[j][i]+ mu_turbulence[j][i-1] + mu_turbulence[j+1][i] + mu_turbulence[j+1][i-1])/4.0;
		   ac_s = (rho[j][i] + rho[j][i-1] + rho[j-1][i] + rho[j-1][i-1])/4.0 * (v[j][i] + v[j][i-1])/2.0;
		   mu_s = (mu[j][i] + mu[j][i-1] + mu[j-1][i] + mu[j-1][i-1] + mu_turbulence[j][i] + mu_turbulence[j][i-1] + mu_turbulence[j-1][i] + mu_turbulence[j-1][i-1])/4.0;

		   A0o_X[j][i] = rho_p*dx*dy/dt;
		   Ae_X[j][i] = -dy*(fabs(ac_e) - ac_e)/2.0 - mu_e*dy/dx ;
		   Aw_X[j][i] = -dy*(fabs(ac_w) + ac_w)/2.0 - mu_w*dy/dx ;
		   An_X[j][i] = -dx*(fabs(ac_n) - ac_n)/2.0 - mu_n*dx/dy ;
		   As_X[j][i] = -dx*(fabs(ac_s) + ac_s)/2.0 - mu_s*dx/dy  ;
		   Ao_X[j][i] = dy*(fabs(ac_e) + ac_e)/2.0 + mu_e*dy/dx + dy*(fabs(ac_w) - ac_w)/2.0 + mu_w*dy/dx
						  + dx*(fabs(ac_n) + ac_n)/2.0 + mu_n*dx/dy + dx*(fabs(ac_s) - ac_s)/2.0 + mu_s*dx/dy + A0o_X[j][i];

		   So_X[j][i] = (p[j][i-1]-p[j][i])*dy + A0o_X[j][i]*u0[j][i];
	   }


   j=numberOfNodes_yDirection-1;
   for (i=1;i<numberOfNodes_xDirection;i++)
   {
	   rho_p = (rho[j][i] + rho[j][i-1])/2.0;
	   ac_e = rho[j][i]*(u[j][i] + u[j][i+1])/2.0;
	   ac_w = rho[j][i-1]*(u[j][i] + u[j][i-1])/2.0;
	   mu_e = mu[j][i] + mu_turbulence[j][i];
	   mu_w = mu[j][i-1] + mu_turbulence[j][i-1];
	   ac_n = (rho[j][i]+ rho[j][i-1])/2.0 * (v[j+1][i] + v[j+1][i-1])/2.0;
	   mu_n = (mu[j][i]+ mu[j][i-1] + mu_turbulence[j][i]+ mu_turbulence[j][i-1])/2.0;
	   ac_s = (rho[j][i] + rho[j][i-1] + rho[j-1][i] + rho[j-1][i-1])/4.0 * (v[j][i] + v[j][i-1])/2.0;
	   mu_s = (mu[j][i] + mu[j][i-1] + mu[j-1][i] + mu[j-1][i-1] + mu_turbulence[j][i] + mu_turbulence[j][i-1] + mu_turbulence[j-1][i] + mu_turbulence[j-1][i-1])/4.0;

	   A0o_X[j][i] = rho_p*dx*dy/dt;
	   Ae_X[j][i] = -dy*(fabs(ac_e) - ac_e)/2.0 - mu_e*dy/dx ;
	   Aw_X[j][i] = -dy*(fabs(ac_w) + ac_w)/2.0 - mu_w*dy/dx ;
	   An_X[j][i] = 0.0 ;
	   As_X[j][i] = -dx*(fabs(ac_s) + ac_s)/2.0 - mu_s*dx/dy - mu_n*dx/(3.0*dy) ;
	   Ao_X[j][i] = dy*(fabs(ac_e) + ac_e)/2.0 + mu_e*dy/dx + dy*(fabs(ac_w) - ac_w)/2.0 + mu_w*dy/dx
					+ dx*(fabs(ac_s) - ac_s)/2.0 + mu_s*dx/dy + 3.0*mu_n*dx/dy + A0o_X[j][i];
					
	   So_X[j][i] = (p[j][i-1]-p[j][i])*dy + A0o_X[j][i]*u0[j][i] - ac_n*u_t[i]*dx + 8.0*mu_n*u_t[i]*dx/(3.0*dy);
   }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::calculate_linkCoeffs_for_YmomentumEq()
{
   double ac_e, ac_w, mu_e, mu_w, ac_n, ac_s, mu_n, mu_s, rho_p;
   
   int i,j;

   i=0;
   for (j=1;j<numberOfNodes_yDirection;j++)
   {
	   rho_p = (rho[j][i] + rho[j-1][i])/2.0;
	   ac_e = (rho[j][i] + rho[j][i+1] + rho[j-1][i] + rho[j-1][i+1])/4.0 * (u[j-1][i+1] + u[j][i+1])/2.0;
	   ac_w = (rho[j][i] + rho[j-1][i])/2.0 * (u[j][i] + u[j-1][i])/2.0;
	   mu_e = (mu[j][i] + mu[j][i+1] + mu[j-1][i] + mu[j-1][i+1] + mu_turbulence[j][i] + mu_turbulence[j][i+1] + mu_turbulence[j-1][i] + mu_turbulence[j-1][i+1])/4.0;
	   mu_w = (mu[j][i] + mu[j-1][i] + mu_turbulence[j][i] + mu_turbulence[j-1][i])/2.0;
	   ac_n = rho[j][i] * (v[j+1][i] + v[j][i])/2.0;
	   mu_n = mu[j][i] + mu_turbulence[j][i];
	   ac_s = rho[j-1][i] * (v[j][i] + v[j-1][i])/2.0;
	   mu_s = mu[j-1][i] + mu_turbulence[j-1][i];

	   A0o_Y[j][i] = rho_p*dx*dy/dt;
	   Ae_Y[j][i] = -dy*(fabs(ac_e) - ac_e)/2.0 - mu_e*dy/dx - mu_w*dy/(3.0*dx) ;
	   Aw_Y[j][i] = 0.0 ;
	   An_Y[j][i] = -dx*(fabs(ac_n) - ac_n)/2.0 - mu_n*dx/dy ;
	   As_Y[j][i] = -dx*(fabs(ac_s) + ac_s)/2.0 - mu_s*dx/dy  ;
	   Ao_Y[j][i] = dy*(fabs(ac_e) + ac_e)/2.0 + mu_e*dy/dx + 3.0*mu_w*dy/dx
					  + dx*(fabs(ac_n) + ac_n)/2.0 + mu_n*dx/dy + dx*(fabs(ac_s) - ac_s)/2.0 + mu_s*dx/dy + A0o_Y[j][i];

	   So_Y[j][i] = (p[j-1][i]-p[j][i])*dx + A0o_Y[j][i]*v0[j][i] + ac_w*v_l[j]*dy + 8.0*mu_w*v_l[j]*dy/(3.0*dx);
   }

   for(i=1;i<numberOfNodes_xDirection-1;i++)
	   for (j=1;j<numberOfNodes_yDirection;j++)
	   {
		   rho_p = (rho[j][i] + rho[j-1][i])/2.0;
		   ac_e = (rho[j][i] + rho[j][i+1] + rho[j-1][i] + rho[j-1][i+1])/4.0 * (u[j-1][i+1] + u[j][i+1])/2.0;
		   ac_w = (rho[j][i] + rho[j][i-1] + rho[j-1][i] + rho[j-1][i-1])/4.0 * (u[j][i] + u[j-1][i])/2.0;
		   mu_e = (mu[j][i] + mu[j][i+1] + mu[j-1][i] + mu[j-1][i+1] + mu_turbulence[j][i] + mu_turbulence[j][i+1] + mu_turbulence[j-1][i] + mu_turbulence[j-1][i+1])/4.0;
		   mu_w = (mu[j][i] + mu[j][i-1] + mu[j-1][i] + mu[j-1][i-1] + mu_turbulence[j][i] + mu_turbulence[j][i-1] + mu_turbulence[j-1][i] + mu_turbulence[j-1][i-1])/4.0;
		   ac_n = rho[j][i] * (v[j+1][i] + v[j][i])/2.0;
		   mu_n = mu[j][i] + mu_turbulence[j][i];
		   ac_s = rho[j-1][i] * (v[j][i] + v[j-1][i])/2.0;
		   mu_s = mu[j-1][i] + mu_turbulence[j-1][i];

		   A0o_Y[j][i] = rho_p*dx*dy/dt;
		   Ae_Y[j][i] = -dy*(fabs(ac_e) - ac_e)/2.0 - mu_e*dy/dx ;
		   Aw_Y[j][i] = -dy*(fabs(ac_w) + ac_w)/2.0 - mu_w*dy/dx ;
		   An_Y[j][i] = -dx*(fabs(ac_n) - ac_n)/2.0 - mu_n*dx/dy ;
		   As_Y[j][i] = -dx*(fabs(ac_s) + ac_s)/2.0 - mu_s*dx/dy  ;
		   Ao_Y[j][i] = dy*(fabs(ac_e) + ac_e)/2.0 + mu_e*dy/dx + dy*(fabs(ac_w) - ac_w)/2.0 + mu_w*dy/dx
						  + dx*(fabs(ac_n) + ac_n)/2.0 + mu_n*dx/dy + dx*(fabs(ac_s) - ac_s)/2.0 + mu_s*dx/dy + A0o_Y[j][i];

		   So_Y[j][i] = (p[j-1][i]-p[j][i])*dx + A0o_Y[j][i]*v0[j][i];
	   }


   i=numberOfNodes_xDirection-1;
   for (j=1;j<numberOfNodes_yDirection;j++)
   {
	   rho_p = (rho[j][i] + rho[j-1][i])/2.0;
	   ac_e = (rho[j][i] + rho[j-1][i])/2.0 * (u[j-1][i+1] + u[j][i+1])/2.0;
	   ac_w = (rho[j][i] + rho[j][i-1] + rho[j-1][i] + rho[j-1][i-1])/4.0 * (u[j][i] + u[j-1][i])/2.0;
	   mu_e = (mu[j][i] + mu[j-1][i] + mu_turbulence[j][i] + mu_turbulence[j-1][i])/2.0;
	   mu_w = (mu[j][i] + mu[j][i-1] + mu[j-1][i] + mu[j-1][i-1] + mu_turbulence[j][i] + mu_turbulence[j][i-1] + mu_turbulence[j-1][i] + mu_turbulence[j-1][i-1])/4.0;
	   ac_n = rho[j][i] * (v[j+1][i] + v[j][i])/2.0;
	   mu_n = mu[j][i] + mu_turbulence[j][i];
	   ac_s = rho[j-1][i] * (v[j][i] + v[j-1][i])/2.0;
	   mu_s = mu[j-1][i] + mu_turbulence[j-1][i];

	   A0o_Y[j][i] = rho_p*dx*dy/dt;
	   Ae_Y[j][i] = 0.0 ;
	   Aw_Y[j][i] = -dy*(fabs(ac_w) + ac_w)/2.0 - mu_w*dy/dx - mu_e*dy/(3.0*dx);
	   An_Y[j][i] = -dx*(fabs(ac_n) - ac_n)/2.0 - mu_n*dx/dy ;
	   As_Y[j][i] = -dx*(fabs(ac_s) + ac_s)/2.0 - mu_s*dx/dy  ;
	   Ao_Y[j][i] =  3.0*mu_e*dy/dx + dy*(fabs(ac_w) - ac_w)/2.0 + mu_w*dy/dx
					  + dx*(fabs(ac_n) + ac_n)/2.0 + mu_n*dx/dy + dx*(fabs(ac_s) - ac_s)/2.0 + mu_s*dx/dy + A0o_Y[j][i];
					  
	   

	   So_Y[j][i] = (p[j-1][i]-p[j][i])*dx + A0o_Y[j][i]*v0[j][i] - ac_e*v_r[j]*dy + 8.0*mu_e*v_r[j]*dy/(3.0*dx);
   }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::solve_XmomentumEq()
{
	   int i,j;
	   double a_row[numberOfNodes_xDirection], b_row[numberOfNodes_xDirection+1], c_row[numberOfNodes_xDirection], r_row[numberOfNodes_xDirection+1], du_row[numberOfNodes_xDirection+1];
	   double a_column[numberOfNodes_yDirection-1], b_column[numberOfNodes_yDirection], c_column[numberOfNodes_yDirection-1], r_column[numberOfNodes_yDirection], du_column[numberOfNodes_yDirection];
	   double res_u = 0;
	   
	   j=0;
	   for(i=1;i<numberOfNodes_xDirection;i++)
			res_u = res_u + (Ae_X[j][i]*u[j][i+1] + Aw_X[j][i]*u[j][i-1] + An_X[j][i]*u[j+1][i]  + Ao_X[j][i]*u[j][i] - So_X[j][i])
						  * (Ae_X[j][i]*u[j][i+1] + Aw_X[j][i]*u[j][i-1] + An_X[j][i]*u[j+1][i]  + Ao_X[j][i]*u[j][i] - So_X[j][i]);  
	   
	   for(j=1;j<numberOfNodes_yDirection-1;j++)
					for(i=1;i<numberOfNodes_xDirection;i++)
						res_u = res_u + (Ae_X[j][i]*u[j][i+1] + Aw_X[j][i]*u[j][i-1] + An_X[j][i]*u[j+1][i] + As_X[j][i]*u[j-1][i] + Ao_X[j][i]*u[j][i] - So_X[j][i])
									  * (Ae_X[j][i]*u[j][i+1] + Aw_X[j][i]*u[j][i-1] + An_X[j][i]*u[j+1][i] + As_X[j][i]*u[j-1][i] + Ao_X[j][i]*u[j][i] - So_X[j][i]);
	   
	   j=numberOfNodes_yDirection-1;
	   for(i=1;i<numberOfNodes_xDirection;i++)
			res_u = res_u + (Ae_X[j][i]*u[j][i+1] + Aw_X[j][i]*u[j][i-1] + As_X[j][i]*u[j-1][i]  + Ao_X[j][i]*u[j][i] - So_X[j][i])
						  * (Ae_X[j][i]*u[j][i+1] + Aw_X[j][i]*u[j][i-1] + As_X[j][i]*u[j-1][i]  + Ao_X[j][i]*u[j][i] - So_X[j][i]);
	   
	   residual_U = sqrt(res_u);

	   for(int inner_iteration_for_u =0; inner_iteration_for_u < inner_iteration_for_momentumEq ;inner_iteration_for_u++)
	   {
		   /////////////////////////////////////////////////////////////////////////////////////////////////
		   //////////////////////////////////// row sweep //////////////////////////////////////////////////
		   /////////////////////////////////////////////////////////////////////////////////////////////////
		   
		   j=0;
		   
		   a_row[numberOfNodes_xDirection-1] = 0.0;
		   b_row[0] = 1.0;
		   b_row[numberOfNodes_xDirection] = 1.0;
		   c_row[0] = 0.0;
		   r_row[0] = 0;
		   r_row[numberOfNodes_xDirection] = 0;
		   
		   for (i=1;i<numberOfNodes_xDirection;i++)
		   {
			   a_row[i-1] = Aw_X[j][i];
			   b_row[i] = (1.0 + gama)*Ao_X[j][i];
			   c_row[i] = Ae_X[j][i];
			   r_row[i] = -An_X[j][i]*u[j+1][i] + So_X[j][i];
			   r_row[i] = r_row[i] - Ae_X[j][i]*u[j][i+1] - Aw_X[j][i]*u[j][i-1] - Ao_X[j][i]*u[j][i];
		   }

		   TDMA_algorithm(a_row,b_row,c_row,r_row,du_row,numberOfNodes_xDirection+1);

		   for (i=0;i<numberOfNodes_xDirection+1;i++)
			  u[j][i] = u[j][i] + du_row[i];
		   

		   for(j=1;j<numberOfNodes_yDirection-1;j++)
		   {
			   a_row[numberOfNodes_xDirection-1] = 0.0;
			   b_row[0] = 1.0;
			   b_row[numberOfNodes_xDirection] = 1.0;
			   c_row[0] = 0.0;
			   r_row[0] = 0;
			   r_row[numberOfNodes_xDirection] = 0;

			   for (i=1;i<numberOfNodes_xDirection;i++)
			   {
				   a_row[i-1] = Aw_X[j][i];
				   b_row[i] = (1.0 + gama)*Ao_X[j][i];
				   c_row[i] = Ae_X[j][i];
				   r_row[i] = -An_X[j][i]*u[j+1][i] - As_X[j][i]*u[j-1][i] + So_X[j][i];
				   r_row[i] = r_row[i] - Ae_X[j][i]*u[j][i+1] - Aw_X[j][i]*u[j][i-1] - Ao_X[j][i]*u[j][i];
			   }

			   TDMA_algorithm(a_row,b_row,c_row,r_row,du_row,numberOfNodes_xDirection+1);

			   for (i=0;i<numberOfNodes_xDirection+1;i++)
				  u[j][i] = u[j][i] + du_row[i];
		   }
		   
		   j=numberOfNodes_yDirection-1;
		   
		   a_row[numberOfNodes_xDirection-1] = 0.0;
		   b_row[0] = 1.0;
		   b_row[numberOfNodes_xDirection] = 1.0;
		   c_row[0] = 0.0;
		   r_row[0] = 0;
		   r_row[numberOfNodes_xDirection] = 0;
		   
		   for (i=1;i<numberOfNodes_xDirection;i++)
		   {
			   a_row[i-1] = Aw_X[j][i];
			   b_row[i] = (1.0 + gama)*Ao_X[j][i];
			   c_row[i] = Ae_X[j][i];
			   r_row[i] = - As_X[j][i]*u[j-1][i] + So_X[j][i];
			   r_row[i] = r_row[i] - Ae_X[j][i]*u[j][i+1] - Aw_X[j][i]*u[j][i-1] - Ao_X[j][i]*u[j][i];
		   }

		   TDMA_algorithm(a_row,b_row,c_row,r_row,du_row,numberOfNodes_xDirection+1);

		   for (i=0;i<numberOfNodes_xDirection+1;i++)
			  u[j][i] = u[j][i] + du_row[i];


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////// column sweep //////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////


		   for (i=1;i<numberOfNodes_xDirection;i++)
		   {
			  for(j=0;j<numberOfNodes_yDirection-1;j++)
			  {
				   a_column[j] = As_X[j+1][i];
				   b_column[j] = (1.0 + gama)*Ao_X[j][i];
				   c_column[j] = An_X[j][i];
			  }
			  b_column[numberOfNodes_yDirection-1] = (1.0 + gama)*Ao_X[numberOfNodes_yDirection-1][i];
				
			  j=0;
			  r_column[j] = -Ae_X[j][i]*u[j][i+1] - Aw_X[j][i]*u[j][i-1] + So_X[j][i];
			  r_column[j] = r_column[j] - An_X[j][i]*u[j+1][i] - Ao_X[j][i]*u[j][i];
			  	
			  for(j=1;j<numberOfNodes_yDirection-1;j++)
			  {
				   r_column[j] = -Ae_X[j][i]*u[j][i+1] - Aw_X[j][i]*u[j][i-1] + So_X[j][i];
				   r_column[j] = r_column[j] - An_X[j][i]*u[j+1][i] - As_X[j][i]*u[j-1][i] - Ao_X[j][i]*u[j][i];
			  }
			  
			  j=numberOfNodes_yDirection-1;
			  r_column[j] = -Ae_X[j][i]*u[j][i+1] - Aw_X[j][i]*u[j][i-1] + So_X[j][i];
			  r_column[j] = r_column[j] - As_X[j][i]*u[j-1][i] - Ao_X[j][i]*u[j][i];

			  TDMA_algorithm(a_column,b_column,c_column,r_column,du_column,numberOfNodes_yDirection);

			  for (j=0;j<numberOfNodes_yDirection;j++)
				u[j][i] = u[j][i] + du_column[j];

		   }

	   }    

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::solve_YmomentumEq()
{
	   int i,j;
	   double a_column[numberOfNodes_yDirection], b_column[numberOfNodes_yDirection+1], c_column[numberOfNodes_yDirection], r_column[numberOfNodes_yDirection+1], dv_column[numberOfNodes_yDirection+1];
	   double a_row[numberOfNodes_xDirection-1], b_row[numberOfNodes_xDirection], c_row[numberOfNodes_xDirection-1], r_row[numberOfNodes_xDirection], dv_row[numberOfNodes_xDirection];
	   double res_v = 0;
	   
	   i=0;
	   for(j=1;j<numberOfNodes_yDirection;j++)
			res_v = res_v + (Ae_Y[j][i]*v[j][i+1] + An_Y[j][i]*v[j+1][i] + As_Y[j][i]*v[j-1][i] + Ao_Y[j][i]*v[j][i] - So_Y[j][i])
						  * (Ae_Y[j][i]*v[j][i+1] + An_Y[j][i]*v[j+1][i] + As_Y[j][i]*v[j-1][i] + Ao_Y[j][i]*v[j][i] - So_Y[j][i]);
	  
	   for(j=1;j<numberOfNodes_yDirection;j++)
			       for(i=1;i<numberOfNodes_xDirection-1;i++)
				        res_v = res_v + (Ae_Y[j][i]*v[j][i+1] + Aw_Y[j][i]*v[j][i-1] + An_Y[j][i]*v[j+1][i] + As_Y[j][i]*v[j-1][i] + Ao_Y[j][i]*v[j][i] - So_Y[j][i])
							          * (Ae_Y[j][i]*v[j][i+1] + Aw_Y[j][i]*v[j][i-1] + An_Y[j][i]*v[j+1][i] + As_Y[j][i]*v[j-1][i] + Ao_Y[j][i]*v[j][i] - So_Y[j][i]);
		
		i=numberOfNodes_xDirection-1;
	    for(j=1;j<numberOfNodes_yDirection;j++)
			res_v = res_v + (Aw_Y[j][i]*v[j][i-1] + An_Y[j][i]*v[j+1][i] + As_Y[j][i]*v[j-1][i] + Ao_Y[j][i]*v[j][i] - So_Y[j][i])
						  * (Aw_Y[j][i]*v[j][i-1] + An_Y[j][i]*v[j+1][i] + As_Y[j][i]*v[j-1][i] + Ao_Y[j][i]*v[j][i] - So_Y[j][i]);	 
		
		residual_V = sqrt(res_v);

	   for(int inner_iteration_for_v =0; inner_iteration_for_v < inner_iteration_for_momentumEq ;inner_iteration_for_v++)
	   {
		   ///////////////////////////////////////////////////////////////////////////////////////////////////////
		   //////////////////////////////////// column sweep /////////////////////////////////////////////////////
		   ///////////////////////////////////////////////////////////////////////////////////////////////////////
		   
		   i=0;
		   
		   a_column[numberOfNodes_yDirection-1] = 0.0;
		   b_column[0] = 1.0;
		   b_column[numberOfNodes_yDirection] = 1.0;
		   c_column[0] = 0.0;
		   r_column[0] = 0;
		   r_column[numberOfNodes_yDirection] = 0;

		   for (j=1;j<numberOfNodes_yDirection;j++)
		   {
			   a_column[j-1] = As_Y[j][i];
			   b_column[j] = (1.0 + gama)*Ao_Y[j][i];
			   c_column[j] = An_Y[j][i];
			   r_column[j] = -Ae_Y[j][i]*v[j][i+1] + So_Y[j][i];
			   r_column[j] = r_column[j] - An_Y[j][i]*v[j+1][i] - As_Y[j][i]*v[j-1][i] - Ao_Y[j][i]*v[j][i];
		   }

		   TDMA_algorithm(a_column,b_column,c_column,r_column,dv_column,numberOfNodes_yDirection+1);

		   for (j=0;j<numberOfNodes_yDirection+1;j++)
			  v[j][i] = v[j][i] + dv_column[j];
		   

		   for(i=1;i<numberOfNodes_xDirection-1;i++)
		   {
			   a_column[numberOfNodes_yDirection-1] = 0.0;
			   b_column[0] = 1.0;
			   b_column[numberOfNodes_yDirection] = 1.0;
			   c_column[0] = 0.0;
			   r_column[0] = 0;
			   r_column[numberOfNodes_yDirection] = 0;

			   for (j=1;j<numberOfNodes_yDirection;j++)
			   {
				   a_column[j-1] = As_Y[j][i];
				   b_column[j] = (1.0 + gama)*Ao_Y[j][i];
				   c_column[j] = An_Y[j][i];
				   r_column[j] = -Ae_Y[j][i]*v[j][i+1] - Aw_Y[j][i]*v[j][i-1] + So_Y[j][i];
				   r_column[j] = r_column[j] - An_Y[j][i]*v[j+1][i] - As_Y[j][i]*v[j-1][i] - Ao_Y[j][i]*v[j][i];
			   }

			   TDMA_algorithm(a_column,b_column,c_column,r_column,dv_column,numberOfNodes_yDirection+1);

			   for (j=0;j<numberOfNodes_yDirection+1;j++)
				  v[j][i] = v[j][i] + dv_column[j];
		   }
		   
		   i=numberOfNodes_xDirection-1;
		   
		   a_column[numberOfNodes_yDirection-1] = 0.0;
		   b_column[0] = 1.0;
		   b_column[numberOfNodes_yDirection] = 1.0;
		   c_column[0] = 0.0;
		   r_column[0] = 0;
		   r_column[numberOfNodes_yDirection] = 0;

		   for (j=1;j<numberOfNodes_yDirection;j++)
		   {
			   a_column[j-1] = As_Y[j][i];
			   b_column[j] = (1.0 + gama)*Ao_Y[j][i];
			   c_column[j] = An_Y[j][i];
			   r_column[j] = - Aw_Y[j][i]*v[j][i-1] + So_Y[j][i];
			   r_column[j] = r_column[j] - An_Y[j][i]*v[j+1][i] - As_Y[j][i]*v[j-1][i] - Ao_Y[j][i]*v[j][i];
		   }

		   TDMA_algorithm(a_column,b_column,c_column,r_column,dv_column,numberOfNodes_yDirection+1);

		   for (j=0;j<numberOfNodes_yDirection+1;j++)
			  v[j][i] = v[j][i] + dv_column[j];

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////// row sweep ///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

		   for (j=1;j<numberOfNodes_yDirection;j++)
		   {
			  for(i=0;i<numberOfNodes_xDirection-1;i++)
			  {
				   a_row[i] = Aw_Y[j][i+1];
				   b_row[i] = (1.0 + gama)*Ao_Y[j][i];
				   c_row[i] = Ae_Y[j][i];
			  }
			  b_row[numberOfNodes_xDirection-1] = (1.0 + gama)*Ao_Y[j][numberOfNodes_xDirection-1];

			  i=0;
			  r_row[i] = -An_Y[j][i]*v[j+1][i] - As_Y[j][i]*v[j-1][i] + So_Y[j][i];
		      r_row[i] = r_row[i] - Ae_Y[j][i]*v[j][i+1] - Ao_Y[j][i]*v[j][i];
		      
			  for(i=1;i<numberOfNodes_xDirection-1;i++)
			  {
				   r_row[i] = -An_Y[j][i]*v[j+1][i] - As_Y[j][i]*v[j-1][i] + So_Y[j][i];
				   r_row[i] = r_row[i] - Ae_Y[j][i]*v[j][i+1] - Aw_Y[j][i]*v[j][i-1] - Ao_Y[j][i]*v[j][i];
			  }
			  
			  i=numberOfNodes_xDirection-1;
			  r_row[i] = -An_Y[j][i]*v[j+1][i] - As_Y[j][i]*v[j-1][i] + So_Y[j][i];
		      r_row[i] = r_row[i] - Aw_Y[j][i]*v[j][i-1] - Ao_Y[j][i]*v[j][i];

			  TDMA_algorithm(a_row,b_row,c_row,r_row,dv_row,numberOfNodes_xDirection);

			  for (i=0;i<numberOfNodes_xDirection;i++)
				v[j][i] = v[j][i] + dv_row[i];

		   }       

	   }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::calculate_linkCoeffs_for_PressureCorrectionEq() 
{
	  int i,j;
      double rho_e, rho_w, rho_n, rho_s;

	  j=0;

	  i=0;

	  rho_w = rho[j][i];
	  rho_s = rho[j][i];
	  rho_e = (rho[j][i] + rho[j][i+1])/2.0;
	  rho_n = (rho[j][i] + rho[j+1][i])/2.0;

	  Aw_P[j][i] = 0;
	  As_P[j][i] = 0;
	  Ae_P[j][i] = -rho_e*dy*dy/Ao_X[j][i+1];
	  An_P[j][i] = -rho_n*dx*dx/Ao_Y[j+1][i];
	  Ao_P[j][i] = -(Ae_P[j][i] + Aw_P[j][i] + An_P[j][i] + As_P[j][i]);
	  So_P[j][i] = -((rho_e*u[j][i+1] - rho_w*u[j][i])*dy + (rho_n*v[j+1][i] - rho_s*v[j][i])*dx);


	  for(i=1;i<numberOfNodes_xDirection-1;i++)
	  {
		  rho_w = (rho[j][i] + rho[j][i-1])/2.0;
		  rho_s = rho[j][i];
		  rho_e = (rho[j][i] + rho[j][i+1])/2.0;
		  rho_n = (rho[j][i] + rho[j+1][i])/2.0;

		  Aw_P[j][i] = -rho_w*dy*dy/Ao_X[j][i];
		  As_P[j][i] = 0;
		  Ae_P[j][i] = -rho_e*dy*dy/Ao_X[j][i+1];
		  An_P[j][i] = -rho_n*dx*dx/Ao_Y[j+1][i];
		  Ao_P[j][i] = -(Ae_P[j][i] + Aw_P[j][i] + An_P[j][i] + As_P[j][i]);
		  So_P[j][i] = -((rho_e*u[j][i+1] - rho_w*u[j][i])*dy + (rho_n*v[j+1][i] - rho_s*v[j][i])*dx);
	  }

	  i=numberOfNodes_xDirection-1;

	  rho_w = (rho[j][i] + rho[j][i-1])/2.0;
	  rho_s = rho[j][i];
	  rho_e = rho[j][i];
	  rho_n = (rho[j][i] + rho[j+1][i])/2.0;

	  Aw_P[j][i] = -rho_w*dy*dy/Ao_X[j][i];
	  As_P[j][i] = 0;
	  Ae_P[j][i] = 0;
	  An_P[j][i] = -rho_n*dx*dx/Ao_Y[j+1][i];
	  Ao_P[j][i] = -(Ae_P[j][i] + Aw_P[j][i] + An_P[j][i] + As_P[j][i]);
	  So_P[j][i] = -((rho_e*u[j][i+1] - rho_w*u[j][i])*dy + (rho_n*v[j+1][i] - rho_s*v[j][i])*dx);


	  for(j=1;j<numberOfNodes_yDirection-1;j++)
	  {
		  i=0;

		  rho_w = rho[j][i];
		  rho_s = (rho[j][i] + rho[j-1][i])/2.0;
		  rho_e = (rho[j][i] + rho[j][i+1])/2.0;
		  rho_n = (rho[j][i] + rho[j+1][i])/2.0;

		  Aw_P[j][i] = 0;
		  As_P[j][i] = -rho_s*dx*dx/Ao_Y[j][i];
		  Ae_P[j][i] = -rho_e*dy*dy/Ao_X[j][i+1];
		  An_P[j][i] = -rho_n*dx*dx/Ao_Y[j+1][i];
		  Ao_P[j][i] = -(Ae_P[j][i] + Aw_P[j][i] + An_P[j][i] + As_P[j][i]);
		  So_P[j][i] = -((rho_e*u[j][i+1] - rho_w*u[j][i])*dy + (rho_n*v[j+1][i] - rho_s*v[j][i])*dx);

		  for(i=1;i<numberOfNodes_xDirection-1;i++)
		  {
			  rho_w = (rho[j][i] + rho[j][i-1])/2.0;
			  rho_s = (rho[j][i] + rho[j-1][i])/2.0;
			  rho_e = (rho[j][i] + rho[j][i+1])/2.0;
			  rho_n = (rho[j][i] + rho[j+1][i])/2.0;

			  Aw_P[j][i] = -rho_w*dy*dy/Ao_X[j][i];
			  As_P[j][i] = -rho_s*dx*dx/Ao_Y[j][i];
			  Ae_P[j][i] = -rho_e*dy*dy/Ao_X[j][i+1];
			  An_P[j][i] = -rho_n*dx*dx/Ao_Y[j+1][i];
			  Ao_P[j][i] = -(Ae_P[j][i] + Aw_P[j][i] + An_P[j][i] + As_P[j][i]);
			  So_P[j][i] = -((rho_e*u[j][i+1] - rho_w*u[j][i])*dy + (rho_n*v[j+1][i] - rho_s*v[j][i])*dx);
		  }

		  i=numberOfNodes_xDirection-1;

		  rho_w = (rho[j][i] + rho[j][i-1])/2.0;
		  rho_s = (rho[j][i] + rho[j-1][i])/2.0;
		  rho_e = rho[j][i];
		  rho_n = (rho[j][i] + rho[j+1][i])/2.0;

		  Aw_P[j][i] = -rho_w*dy*dy/Ao_X[j][i];
		  As_P[j][i] = -rho_s*dx*dx/Ao_Y[j][i];
		  Ae_P[j][i] = 0;
		  An_P[j][i] = -rho_n*dx*dx/Ao_Y[j+1][i];
		  Ao_P[j][i] = -(Ae_P[j][i] + Aw_P[j][i] + An_P[j][i] + As_P[j][i]);
		  So_P[j][i] = -((rho_e*u[j][i+1] - rho_w*u[j][i])*dy + (rho_n*v[j+1][i] - rho_s*v[j][i])*dx);

	  }

	  j=numberOfNodes_yDirection-1;

	  i=0;

	  rho_w = rho[j][i];
	  rho_s = (rho[j][i] + rho[j-1][i])/2.0;
	  rho_e = (rho[j][i] + rho[j][i+1])/2.0;
	  rho_n = rho[j][i];

	  Aw_P[j][i] = 0;
	  As_P[j][i] = -rho_s*dx*dx/Ao_Y[j][i];
	  Ae_P[j][i] = -rho_e*dy*dy/Ao_X[j][i+1];
	  An_P[j][i] = 0;
	  Ao_P[j][i] = -(Ae_P[j][i] + Aw_P[j][i] + An_P[j][i] + As_P[j][i]);
	  So_P[j][i] = -((rho_e*u[j][i+1] - rho_w*u[j][i])*dy + (rho_n*v[j+1][i] - rho_s*v[j][i])*dx);

	  for(i=1;i<numberOfNodes_xDirection-1;i++)
	  {
		  rho_w = (rho[j][i] + rho[j][i-1])/2.0;
		  rho_s = (rho[j][i] + rho[j-1][i])/2.0;
		  rho_e = (rho[j][i] + rho[j][i+1])/2.0;
		  rho_n = rho[j][i];

		  Aw_P[j][i] = -rho_w*dy*dy/Ao_X[j][i];
		  As_P[j][i] = -rho_s*dx*dx/Ao_Y[j][i];
		  Ae_P[j][i] = -rho_e*dy*dy/Ao_X[j][i+1];
		  An_P[j][i] = 0;
		  Ao_P[j][i] = -(Ae_P[j][i] + Aw_P[j][i] + An_P[j][i] + As_P[j][i]);
		  So_P[j][i] = -((rho_e*u[j][i+1] - rho_w*u[j][i])*dy + (rho_n*v[j+1][i] - rho_s*v[j][i])*dx);
	  }

	  i=numberOfNodes_xDirection-1;

	  rho_w = (rho[j][i] + rho[j][i-1])/2.0;
	  rho_s = (rho[j][i] + rho[j-1][i])/2.0;
	  rho_e = rho[j][i];
	  rho_n = rho[j][i];

	  Aw_P[j][i] = -rho_w*dy*dy/Ao_X[j][i];
	  As_P[j][i] = -rho_s*dx*dx/Ao_Y[j][i];
	  Ae_P[j][i] = 0;
	  An_P[j][i] = 0;
	  Ao_P[j][i] = -(Ae_P[j][i] + Aw_P[j][i] + An_P[j][i] + As_P[j][i]);
	  So_P[j][i] = -((rho_e*u[j][i+1] - rho_w*u[j][i])*dy + (rho_n*v[j+1][i] - rho_s*v[j][i])*dx);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::solve_PressureCorrectionEq()
{
	  int i,j; 
	  double a_column[numberOfNodes_yDirection-1], b_column[numberOfNodes_yDirection], c_column[numberOfNodes_yDirection-1], r_column[numberOfNodes_yDirection], pc_column[numberOfNodes_yDirection];
	  double a_row[numberOfNodes_xDirection-1], b_row[numberOfNodes_xDirection], c_row[numberOfNodes_xDirection-1], r_row[numberOfNodes_xDirection], pc_row[numberOfNodes_xDirection];

	  for(j=0;j<numberOfNodes_yDirection;j++)
			for(i=0;i<numberOfNodes_xDirection;i++)
				pc[j][i] = 0;
				
	  double res_p = 0;
	  for(j=0;j<numberOfNodes_yDirection;j++)
		for(i=0;i<numberOfNodes_xDirection;i++)
			res_p = res_p + So_P[j][i]*So_P[j][i];

	  residual_P = sqrt(res_p);

	  for(int inner_iteration_for_p =0; inner_iteration_for_p < inner_iteration_for_continuityEq ;inner_iteration_for_p++)
	  {

	 //////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// row sweep ////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

		  j=0;
		  	  
		  for(i=0;i<numberOfNodes_xDirection-1;i++)
		  {
			   a_row[i] = Aw_P[j][i+1];
			   b_row[i] = Ao_P[j][i];
			   c_row[i] = Ae_P[j][i];
			   r_row[i] = So_P[j][i];
			   r_row[i] = r_row[i] - An_P[j][i]*pc[j+1][i];
		  }

		  i=numberOfNodes_xDirection-1;
		  b_row[i] = Ao_P[j][i];
		  r_row[i] = So_P[j][i];
		  r_row[i] = r_row[i] - An_P[j][i]*pc[j+1][i];

		  TDMA_algorithm(a_row,b_row,c_row,r_row,pc_row,numberOfNodes_xDirection);

		  for(i=0;i<numberOfNodes_xDirection;i++)
			 pc[j][i] = pc_row[i];
			 
		  
		  for(j=1;j<numberOfNodes_yDirection-1;j++)
		  {
			  for(i=0;i<numberOfNodes_xDirection-1;i++)
			  {
				   a_row[i] = Aw_P[j][i+1];
				   b_row[i] = Ao_P[j][i];
				   c_row[i] = Ae_P[j][i];
				   r_row[i] = So_P[j][i];
				   r_row[i] = r_row[i] - As_P[j][i]*pc[j-1][i] - An_P[j][i]*pc[j+1][i];
			  }

			  i=numberOfNodes_xDirection-1;
			  b_row[i] = Ao_P[j][i];
			  r_row[i] = So_P[j][i];
			  r_row[i] = r_row[i] - As_P[j][i]*pc[j-1][i] - An_P[j][i]*pc[j+1][i];

			  TDMA_algorithm(a_row,b_row,c_row,r_row,pc_row,numberOfNodes_xDirection);

			  for(i=0;i<numberOfNodes_xDirection;i++)
				 pc[j][i] = pc_row[i];
		  }
		  
		  
		  j=numberOfNodes_yDirection-1;
		  
		  for(i=0;i<numberOfNodes_xDirection-1;i++)
		  {
			   a_row[i] = Aw_P[j][i+1];
			   b_row[i] = Ao_P[j][i];
			   c_row[i] = Ae_P[j][i];
			   r_row[i] = So_P[j][i];
			   r_row[i] = r_row[i] - As_P[j][i]*pc[j-1][i] ;
		  }

		  i=numberOfNodes_xDirection-1;
		  b_row[i] = Ao_P[j][i];
		  r_row[i] = So_P[j][i];
		  r_row[i] = r_row[i] - As_P[j][i]*pc[j-1][i];

		  TDMA_algorithm(a_row,b_row,c_row,r_row,pc_row,numberOfNodes_xDirection);

		  for(i=0;i<numberOfNodes_xDirection;i++)
			 pc[j][i] = pc_row[i];

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////// column sweep //////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		  i=0;
		  
		  for(j=0;j<numberOfNodes_yDirection-1;j++)
		  {
			   a_column[j] = As_P[j+1][i];
			   b_column[j] = Ao_P[j][i];
			   c_column[j] = An_P[j][i];
			   r_column[j] = So_P[j][i];
			   r_column[j] = r_column[j] - Ae_P[j][i]*pc[j][i+1];
		  }

		  j=numberOfNodes_yDirection-1;
		  b_column[j] = Ao_P[j][i];
		  r_column[j] = So_P[j][i];
		  r_column[j] = r_column[j] - Ae_P[j][i]*pc[j][i+1];

		  TDMA_algorithm(a_column,b_column,c_column,r_column,pc_column,numberOfNodes_yDirection);

		  for(j=0;j<numberOfNodes_yDirection;j++)
			 pc[j][i] = pc_column[j];
		 
		 
		  for(i=1;i<numberOfNodes_xDirection-1;i++)
		  {
			  for(j=0;j<numberOfNodes_yDirection-1;j++)
			  {
				   a_column[j] = As_P[j+1][i];
				   b_column[j] = Ao_P[j][i];
				   c_column[j] = An_P[j][i];
				   r_column[j] = So_P[j][i];
				   r_column[j] = r_column[j] - Ae_P[j][i]*pc[j][i+1] - Aw_P[j][i]*pc[j][i-1];
			  }

			  j=numberOfNodes_yDirection-1;
			  b_column[j] = Ao_P[j][i];
			  r_column[j] = So_P[j][i];
			  r_column[j] = r_column[j] - Ae_P[j][i]*pc[j][i+1] - Aw_P[j][i]*pc[j][i-1];

			  TDMA_algorithm(a_column,b_column,c_column,r_column,pc_column,numberOfNodes_yDirection);

			  for(j=0;j<numberOfNodes_yDirection;j++)
				 pc[j][i] = pc_column[j];
		  }
		  
		  i=numberOfNodes_xDirection-1;
		  
		  for(j=0;j<numberOfNodes_yDirection-1;j++)
		  {
			   a_column[j] = As_P[j+1][i];
			   b_column[j] = Ao_P[j][i];
			   c_column[j] = An_P[j][i];
			   r_column[j] = So_P[j][i];
			   r_column[j] = r_column[j] - Aw_P[j][i]*pc[j][i-1];
		  }

		  j=numberOfNodes_yDirection-1;
		  b_column[j] = Ao_P[j][i];
		  r_column[j] = So_P[j][i];
		  r_column[j] = r_column[j] - Aw_P[j][i]*pc[j][i-1];

		  TDMA_algorithm(a_column,b_column,c_column,r_column,pc_column,numberOfNodes_yDirection);

		  for(j=0;j<numberOfNodes_yDirection;j++)
			 pc[j][i] = pc_column[j];      

	   }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::update_solutionFields_inNewTimeStep()
{
   int i,j;
   
   for(j=0;j<numberOfNodes_yDirection;j++)
  	  for(i=0;i<numberOfNodes_xDirection;i++)
			p[j][i] = p[j][i] + omega_p * pc[j][i];

   for(j=0;j<numberOfNodes_yDirection;j++)
		for(i=1;i<numberOfNodes_xDirection;i++)
				u[j][i] = u[j][i] + omega_uv * (pc[j][i-1] - pc[j][i])*dy/Ao_X[j][i];

   for(j=1;j<numberOfNodes_yDirection;j++)
		for(i=0;i<numberOfNodes_xDirection;i++)
				v[j][i] = v[j][i] + omega_uv * (pc[j-1][i] - pc[j][i])*dy/Ao_Y[j][i];
				
	simple_iteration++;
				
	cout << "Time = " << time << '\t' << '\t' << "iteration = " << simple_iteration << endl << endl;
    cout << "residual for X-momentum Eq is = " << residual_U << endl;
    cout << "residual for Y-momentum Eq is = " << residual_V << endl;
    cout << "residual for Continuity Eq is = " << residual_P << endl;

    cout << endl << "................................................." << endl << endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::update_solutionFields_inPreviousTimeStep()
{
   for(int j=0;j<numberOfNodes_yDirection;j++)
			for(int i=0;i<numberOfNodes_xDirection+1;i++)
				u0[j][i] = u[j][i];
				
		for(int j=0;j<numberOfNodes_yDirection+1;j++)
			for(int i=0;i<numberOfNodes_xDirection;i++)
				v0[j][i] = v[j][i];

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::update_time()
{
	check_writing_solutionFields();
	time = time + dt;
	simple_iteration = 0;
	count_for_writingSolutionFields++;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool simpleSolver::runtime()
{
	if(time > end_time)
		return false;
		
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::check_writing_solutionFields()
{
	if(count_for_writingSolutionFields % numberOfDt_for_writingSolutionFields == 0)
		write_solutionFields();	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
bool simpleSolver::check_convergency()
{
	if(residual_U < error_max && residual_V < error_max && residual_P < error_max)
		return false;

	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void simpleSolver::write_solutionFields()
{
	int i,j;
	
	stringstream ss;
	ss << time;
	string str = ss.str();
	string myName = "results/t=" + str ;
	const char* fileName = myName.c_str();	
		
	// create a folder for writing results in specified time:
	
	#ifdef _WIN32
	
	string winPath = "rmdir /Q /S " + myName;
	const char* winNamePath = winPath.c_str();	
	
	if (opendir(fileName))
			system(winNamePath);

	if (mkdir(fileName) != -1)
        {}
        
	#endif
	
	#ifdef linux
	
	string linPath = "rm -rf " + myName;
	const char* linNamePath = linPath.c_str();
	
	DIR *dir = opendir(fileName);
	
	if (dir)
		system(linNamePath);

	closedir(dir);

	if (mkdir(fileName, 0777) != -1)
        {}
        
    #endif
       
    ofstream file1;
    string Upath = myName + "/U.txt";
    const char* UpathName = Upath.c_str();	
    file1.open(UpathName);

    file1 << fixed << setprecision(6);

    for(j=0;j<numberOfNodes_yDirection;j++)
        for (i=0;i<numberOfNodes_xDirection+1;i++)
            file1 << i*dx << ", " << j*dy + dy/2.0 << ", " << u[j][i] << endl;
            
    for (i=0;i<numberOfNodes_xDirection+1;i++)
        file1 << i*dx << ", " << 0*dy << ", " << u_b[i] << endl;

    for (i=0;i<numberOfNodes_xDirection+1;i++)
        file1 << i*dx << ", " << numberOfNodes_yDirection*dy << ", " << u_t[i] << endl;

    file1.close();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ofstream file2;
    string Vpath = myName + "/V.txt";
    const char* VpathName = Vpath.c_str();	
    file2.open(VpathName);

    file2 << fixed << setprecision(6);

    for(i=0;i<numberOfNodes_xDirection;i++)
        for (j=0;j<numberOfNodes_yDirection+1;j++)
            file2 << i*dx + dx/2.0 << ", " << j*dy << ", " << v[j][i] << endl;
            
    for (j=0;j<numberOfNodes_yDirection+1;j++)
        file2 << 0*dx << ", " << j*dy << ", " << v_l[j] << endl;

    for (j=0;j<numberOfNodes_yDirection+1;j++)
        file2 << numberOfNodes_xDirection*dx << ", " << j*dy << ", " << v_r[j] << endl;

    file2.close();
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    
    ofstream file3;
    string Ppath = myName + "/P.txt";
    const char* PpathName = Ppath.c_str();	
    file3.open(PpathName);

    file3 << fixed << setprecision(6);

    for(j=0;j<numberOfNodes_yDirection;j++)
        for (i=0;i<numberOfNodes_xDirection;i++)
            file3 << i*dx + dx/2.0 << ", " << j*dy + dy/2.0 << ", " << p[j][i] << endl;
            
    for (i=0;i<numberOfNodes_xDirection;i++)
        file3 << i*dx + dx/2.0 << ", " << 0*dy << ", " << p[0][i] << endl;

    for (i=0;i<numberOfNodes_xDirection;i++)
        file3 << i*dx + dx/2.0 << ", " << numberOfNodes_yDirection*dy << ", " << p[numberOfNodes_yDirection-1][i] << endl;

    for (j=0;j<numberOfNodes_yDirection;j++)
        file3 << 0*dx << ", " << j*dy + dy/2.0 << ", " << p[j][0] << endl;

    for (j=0;j<numberOfNodes_yDirection;j++)
        file3 << numberOfNodes_xDirection*dx << ", " << j*dy + dy/2.0 << ", " << p[j][numberOfNodes_xDirection-1] << endl;

    file3 << 0*dx << ", " << 0*dy << ", " << p[0][0] << endl;
    file3 << 0*dx << ", " << numberOfNodes_yDirection*dy << ", " << p[numberOfNodes_yDirection-1][0] << endl;
    file3 << numberOfNodes_xDirection*dx << ", " << 0*dy << ", " << p[0][numberOfNodes_xDirection-1] << endl;
    file3 << numberOfNodes_xDirection*dx << ", " << numberOfNodes_yDirection*dy << ", " << p[numberOfNodes_yDirection-1][numberOfNodes_xDirection-1] << endl;

    file3.close();

     //////////////////////////////////////////////////////////////////////////////////////////////////////
     
    ofstream file4;
    string Streamlinespath = myName + "/Streamlines.txt";
    const char* StreamlinespathName = Streamlinespath.c_str();
	file4.open(StreamlinespathName);

	file4 << fixed << setprecision(6);
	
	// interior nodes:
	for(j=0;j<numberOfNodes_yDirection;j++)
		for (i=0;i<numberOfNodes_xDirection;i++)
			file4 << i*dx + dx/2.0 << ", " << j*dy + dy/2.0 << ", " << (u[j][i] + u[j][i+1])/2.0 << ", " << (v[j][i] + v[j+1][i])/2.0 << endl;
	
	// bottom wall:		
	for (i=0;i<numberOfNodes_xDirection;i++)
		file4 << i*dx + dx/2.0 << ", " << 0*dy << ", " << (u_b[i] + u_b[i+1])/2.0 << ", " << v_b[i] << endl;

	// top wall:
	for (i=0;i<numberOfNodes_xDirection;i++)
		file4 << i*dx + dx/2.0 << ", " << numberOfNodes_yDirection*dy << ", " << (u_t[i] + u_t[i+1])/2.0  << ", " << v_t[i] << endl;
		
	// left wall:		
	for (j=0;j<numberOfNodes_yDirection;j++)
		file4 << 0*dx << ", " << j*dy + dy/2.0 << ", " << u_l[j] << ", " << (v_l[j] + v_l[j+1])/2.0 << endl;
		
	// right wall:		
	for (j=0;j<numberOfNodes_yDirection;j++)
		file4 << numberOfNodes_xDirection*dx << ", " << j*dy + dy/2.0 << ", " << u_r[j] << ", " << (v_r[j] + v_r[j+1])/2.0 << endl;
		
	// four points at edges of wall:
	file4 << 0*dx << ", " << 0*dy << ", " << (u_b[0] + u_l[0])/2.0 << ", " << (v_b[0] + v_l[0])/2.0 << endl;
	file4 << 0*dx << ", " << numberOfNodes_yDirection*dy << ", " << (u_t[0] + u_l[numberOfNodes_yDirection-1])/2.0 << ", " << (v_t[0] + v_l[numberOfNodes_yDirection])/2.0 << endl;
	file4 << numberOfNodes_xDirection*dx << ", " << 0*dy << ", " << (u_b[numberOfNodes_xDirection] + u_r[0])/2.0 << ", " << (v_b[numberOfNodes_xDirection-1] + v_r[0])/2.0 << endl;
	file4 << numberOfNodes_xDirection*dx << ", " << numberOfNodes_yDirection*dy << ", " << (u_t[numberOfNodes_xDirection] + u_r[numberOfNodes_yDirection-1])/2.0 << ", " << (v_t[numberOfNodes_xDirection-1] + v_r[numberOfNodes_yDirection])/2.0 << endl;

	file4.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::readInputData(const char *path)
{
	ifstream file(path);

	string searchForLength = "Length";
	string searchForEndTime = "EndTime";
	string searchForStartTime = "Startime";
	string searchForNumberOfDt_for_writingSolutionFields = "NumberOfDtForWritingSolutionFields";
	string searchForDelta_t = "Delta_t";
	string searchForMovingWallVelocity = "MovingWallVelocity";
	string searchForDensity = "Density";
	string searchForDynamicViscosity = "DynamicViscosity";
	string searchForNumberOfIterationForMomentumEq = "NumberOfIterationForMomentumEq";
	string searchForNumberOfIterationForContinuityEq = "NumberOfIterationForPressureCorrectionEq";
    string searchForGamma = "Gamma";
    string searchForOmega_uv = "Omega_uv";
    string searchForOmega_p = "Omega_p";
	string searchForNumberOfNodesInXDirection = "NumberOfNodesInXDirection";
	string searchForNumberOfNodesInYDirection = "NumberOfNodesInYDirection";
	string searchForMaxError = "MaximumErrorForSolving";


	string lineOfText;

	for(;;)
	{
		getline(file, lineOfText);

		if (file.eof()) break;

		if (lineOfText.find(searchForLength, 0) != string::npos)
		  file >> Length;
		  
		if (lineOfText.find(searchForStartTime, 0) != string::npos)
		  file >> start_time;

		if (lineOfText.find(searchForEndTime, 0) != string::npos)
		  file >> end_time;
		  
		if (lineOfText.find(searchForNumberOfDt_for_writingSolutionFields, 0) != string::npos)
		  file >> numberOfDt_for_writingSolutionFields;
		  
		if (lineOfText.find(searchForDelta_t, 0) != string::npos)
		  file >> dt;

		if (lineOfText.find(searchForMovingWallVelocity, 0) != string::npos)
		  file >> free_velocity;

		if (lineOfText.find(searchForDensity, 0) != string::npos)
		  file >> rho_;

		if (lineOfText.find(searchForDynamicViscosity, 0) != string::npos)
		  file >> mu_;

		if (lineOfText.find(searchForNumberOfIterationForMomentumEq, 0) != string::npos)
		  file >> inner_iteration_for_momentumEq;

		if (lineOfText.find(searchForNumberOfIterationForContinuityEq, 0) != string::npos)
		  file >> inner_iteration_for_continuityEq;

		if (lineOfText.find(searchForGamma, 0) != string::npos)
		  file >> gama;

		if (lineOfText.find(searchForOmega_uv, 0) != string::npos)
		  file >> omega_uv;

		if (lineOfText.find(searchForOmega_p, 0) != string::npos)
		  file >> omega_p;

		if (lineOfText.find(searchForNumberOfNodesInXDirection, 0) != string::npos)
		  file >> numberOfNodes_xDirection;

		if (lineOfText.find(searchForNumberOfNodesInYDirection, 0) != string::npos)
		  file >> numberOfNodes_yDirection;

		if (lineOfText.find(searchForMaxError, 0) != string::npos)
		  file >> error_max;
    }

    file.close();

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::printInputData()
{
	cout << "Input data for the simulation: " << endl << endl;
	cout << "Reynolds Number is: " << ReynoldsNumber << endl;
	cout << "Length is: " << Length << endl;
	cout << "Moving Wall Velocity is: " << free_velocity << endl;
	cout << "Density is: " << rho_ << endl;
	cout << "DynamicViscosity is: " << mu_ << endl;
	cout << "Number of iteration for solving Momentum Eqs is: " << inner_iteration_for_momentumEq << endl;
	cout << "Number of iteration for solving Pressure Correction Eq is: " << inner_iteration_for_continuityEq << endl;
	cout << "Gamma for stability of Momentum Eqs is: " << gama << endl;
	cout << "Omega_uv for relaxation factor of Momentum Eqs is: " << omega_uv << endl;
	cout << "Omega_p for relaxation factor of Pressure Correction Eq is: " << omega_p << endl;
	cout << "Number of nodes in X direction is: " << numberOfNodes_xDirection << endl;
	cout << "Number of nodes in Y direction is: " << numberOfNodes_yDirection << endl;
	cout << "Start time is: " << start_time << endl;
	cout << "End time is: " << end_time << endl;
	cout << "Delta_t is: " << dt << endl;
	cout << "Number of dt for writing solution fields is: " << numberOfDt_for_writingSolutionFields << endl;
	cout << "Maximum error for solving Momentum Pressure Correction Eqs is: " << error_max << endl << endl;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void simpleSolver::TDMA_algorithm(double* a, double* b, double* c, double* r, double* x, int n)
{
    for (int i = 0; i < n-1; i++)
    {
        b[i+1]  = b[i+1] - (a[i]/b[i])*c[i];
        r[i+1] = r[i+1] - (a[i]/b[i])*r[i];
    }

    x[n-1] = 1.0*r[n-1]/b[n-1];

    for (int i = n-2; i >= 0; i--)
        x[i] = (r[i] - c[i]*x[i+1])/b[i];
}
