
/////////////////////// derived class ///////////////////////////

class LaminarModel: public flowModelType {
		       
   public:
   
      string modelFlow_name = "laminar";
      
      LaminarModel (const char *path, double **U, double **V , double **mu_SGS)
      { printInputData(); }
     
	  void update_mu_turbulence() 
	  {}
      
      void write_solutionField()
      {}
      
      void update_time()
      {}
      
      void printInputData()
	  { cout << "Flow model type is: " << modelFlow_name << endl << endl; }

};
