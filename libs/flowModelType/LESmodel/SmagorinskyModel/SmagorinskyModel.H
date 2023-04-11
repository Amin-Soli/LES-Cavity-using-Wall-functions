
/////////////////////// derived class ///////////////////////////

class SmagorinskyModel: public LESmodel {
	
   private:
   
	   double Smagorinsky_coeff, l_Smagorinsky;
	           
   public:
   
       string turbulentModel_name = "Smagorinsky";
   
	   SmagorinskyModel (const char* , double** , double** , double** );
	      
       void calculate_delta_subgrid();
       
       void readInputData(const char* );
       
       void printInputData();
		
};

