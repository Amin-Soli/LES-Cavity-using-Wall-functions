/////////////////////// derived class ///////////////////////////

class notWallFunction: public wallFunction{
		   
   public:
   
	   notWallFunction (const char* , double** , double** );
	      
       void calculate_wallsKinematicViscosity();
       
       void print_wallFunction_information();
       
       void write_solutionField();
		
};

