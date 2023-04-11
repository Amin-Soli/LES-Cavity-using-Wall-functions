
class simpleSolver {
	
	private:
	
		int numberOfNodes_xDirection, numberOfNodes_yDirection;
		
		int simple_iteration, count_for_writingSolutionFields, numberOfDt_for_writingSolutionFields;
		
		double Length, dx, dy;
		
		double start_time, end_time, time, dt;
		
		double free_velocity, rho_, mu_, error_max;
		
		int inner_iteration_for_momentumEq, inner_iteration_for_continuityEq;  // number of iteration for solving Momentum and Pressure Correction Eqs
		
		double gama, omega_uv, omega_p;      // relaxtation factors
		
		double residual_U, residual_V, residual_P;   
		
		double **p, **pc, **u, **v, **u0, **v0, **rho, **mu, **mu_turbulence;       //field variables
		
		double **Ao_X, **A0o_X, **Ae_X, **Aw_X, **An_X, **As_X, **So_X;      // link Coeffs for X-momentum Eq
			   
		double **Ao_Y, **A0o_Y, **Ae_Y, **Aw_Y, **An_Y, **As_Y, **So_Y;    // link Coeffs for Y-momentum Eq
		
		double **Ao_P, **Ae_P, **Aw_P, **An_P, **As_P, **So_P;       //link Coeffs for Pressure correction Eq;
		
		double *u_l, *u_r, *u_t, *u_b;     // boundary variable for U
		
		double *v_b, *v_t, *v_r, *v_l;     // boundary variable for V 
		
		double ReynoldsNumber;
    
	       
   public:
   
	  simpleSolver (const char* , double** , double** , double** ,  double** );
	  
	  ~simpleSolver ();
	  
      void calculate_linkCoeffs_for_XmomentumEq(); 
      
      void solve_XmomentumEq(); 
      
      void calculate_linkCoeffs_for_YmomentumEq(); 
      
      void solve_YmomentumEq(); 
      
      void calculate_linkCoeffs_for_PressureCorrectionEq(); 
      
      void solve_PressureCorrectionEq(); 
      
      void update_solutionFields_inPreviousTimeStep();
      
      void update_solutionFields_inNewTimeStep();
      
      void update_time();
      
      void setBoundaryConditions(); 
      
      void readInputData(const char* );
      
      void printInputData();
      
      void write_solutionFields();
      
      void TDMA_algorithm(double* , double* , double* , double* , double* , int );
      
      bool check_convergency();
      
      bool runtime();
      
      void check_writing_solutionFields();

};
