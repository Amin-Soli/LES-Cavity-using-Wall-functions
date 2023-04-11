/////////////////////// derived class ///////////////////////////

class enhancedWallFunction: public wallFunction{
		   
   public:
   
	   string wallFunction_name = "Enhanced Wall Function";
   
	   enhancedWallFunction (const char* , double** , double** );
	      
       void calculate_wallsKinematicViscosity();
       
       void print_wallFunction_information();
       
       void write_solutionField();
		
};

