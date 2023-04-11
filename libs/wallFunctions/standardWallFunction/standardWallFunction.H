/////////////////////// derived class ///////////////////////////

class standardWallFunction: public wallFunction{
		   
   public:
   
	   string wallFunction_name = "Standard Wall Function";
	   
	   double Yplus_laminar = 11.25;
   
	   standardWallFunction (const char* , double** , double** );
	      
       void calculate_wallsKinematicViscosity();
       
       void print_wallFunction_information();
       
       void write_solutionField();
		
};

