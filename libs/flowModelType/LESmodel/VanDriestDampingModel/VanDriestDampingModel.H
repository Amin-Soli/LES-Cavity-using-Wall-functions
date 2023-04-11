
/////////////////////// derived class ///////////////////////////

class VanDriestDampingModel: public LESmodel {
	
   private:
   
       double Smagorinsky_coeff, l_Smagorinsky, kappa, E, A_plus, mu_;
   
	   double *Yplus_bottomWall, *Yplus_leftWall, *Yplus_rightWall;
	   
	   double **mu, **mixingLength;
	   
	   double **yFromBottomWall, **yFromLeftWall, **yFromRightWall, **min_y, **Yplus;
	   
	   
   public:
   
       string turbulentModel_name = "Smagorinsky with Van Driest damping function";
   
	   VanDriestDampingModel (const char* , double** , double** , double** );
	   
	   ~VanDriestDampingModel();
	      
       void calculate_delta_subgrid();
       
       void calculate_Yplus_for_allBoundaryWalls();
       
       void calculate_mixingLengths_from_allBoundaryWalls();
       
       void readInputData(const char* );
       
       void printInputData();
       
       double calculate_Yplus_forEachFace_usingNewtonRaphson(double , double , double , double );
		
};

