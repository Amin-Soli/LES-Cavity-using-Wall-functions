
/////////////////////// abstract class ///////////////////////////

class flowModelType {
		       
   public:
   
	  virtual ~flowModelType()
	  {}
     
	  virtual void update_mu_turbulence() = 0;
      
      virtual void write_solutionField() = 0;
      
      virtual void update_time() = 0;
      
      static flowModelType* New(const char* , double** , double** , double** );

};
