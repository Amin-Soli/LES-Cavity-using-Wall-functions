
/////////////////////// abstract derived class ///////////////////////////

class LESmodel: public flowModelType {
	
	protected:
	
		int numberOfNodes_xDirection, numberOfNodes_yDirection;
		
		int count_for_writingSolutionField, numberOfDt_for_writingSolutionField;
		
		double Length, dx, dy, rho_, free_velocity;
		
		double start_time, time, dt;
		
		double **u, **v, **mu_sgs, **rho, **delta_sgs;
		
		double *v_r, *v_l, *u_t, *u_b; 
  
	       
   public:
   
      string modelFlow_name = "turbulent";
   
	  LESmodel (const char* , double** , double** , double** );
	  
	  ~LESmodel ();
	  
	  void update_mu_turbulence();
	  
	  virtual void calculate_delta_subgrid() = 0;
      
      void readInputData(const char* );
      
      virtual void printInputData() = 0;
      
      void write_solutionField();
      
      void update_time();
      
      void check_writing_solutionField();
      
      static LESmodel* New(const char* , double** , double** , double** );

};
