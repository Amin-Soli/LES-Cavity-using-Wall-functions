/////////////////////// abstract base class ///////////////////////////

class wallFunction{
	
   protected:
   
       int numberOfNodes_xDirection, numberOfNodes_yDirection;
       
       int count_for_writingSolutionField, numberOfDt_for_writingSolutionField;
       
       double start_time, time, dt;
		
	   double Length, dx, dy, rho_, mu_, free_velocity;
   
	   double *Yplus_bottomWall, *Yplus_leftWall, *Yplus_rightWall;
	   
	   double **u, **v, **mu, **rho, **nu;
	   
	   double *tau_bottomWall, *tau_leftWall, *tau_rightWall;     // walls shear stress
	   
	   double *nu_bottomWall, *nu_leftWall, *nu_rightWall;     // walls kinematic viscosity
		
	   double average_tau_bottomWall, average_tau_leftWall, average_tau_rightWall;     // average walls shear stress
	   
	   double kappa, E;
	   
	   string flowModel, usingWallFunction_flag;
	   
	   
   public:
   
	   wallFunction (const char* , double** , double** );
	   
	   virtual ~wallFunction();
	   
	   void update_time();
	   
	   void check_writing_solutionField();
	   
	   virtual void write_solutionField();
	      
       virtual void calculate_wallsKinematicViscosity() = 0;
       
       void calculate_Yplus_for_allBoundaryWalls();
       
       void readInputData(const char* );
       
       double calculate_Yplus_forEachFace_usingNewtonRaphson(double , double , double , double );
       
       void calculate_shearStress_for_allBoundaryWalls();
       
       static wallFunction* New(const char* , double** , double** );
		
};

