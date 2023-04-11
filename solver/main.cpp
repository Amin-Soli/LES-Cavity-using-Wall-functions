#include "allHeaders.H"

int main() 
{	
	#include "createFields.H"
    
    SIMPLE_Object.setBoundaryConditions(); 
   
    while (SIMPLE_Object.runtime())
    {	
		    wallFunction_Object -> calculate_shearStress_for_allBoundaryWalls();
		    
		    FlowModel_Object -> update_mu_turbulence();
				
			SIMPLE_Object.update_time();
			
			FlowModel_Object -> update_time();	
			
			wallFunction_Object -> update_time();
			
			SIMPLE_Object.update_solutionFields_inPreviousTimeStep();
			
			do
			{
				   SIMPLE_Object.calculate_linkCoeffs_for_XmomentumEq();
				   
				   SIMPLE_Object.calculate_linkCoeffs_for_YmomentumEq();
				   
				   SIMPLE_Object.solve_XmomentumEq();
				   
				   SIMPLE_Object.solve_YmomentumEq();
				   
				   SIMPLE_Object.calculate_linkCoeffs_for_PressureCorrectionEq();
				   
				   SIMPLE_Object.solve_PressureCorrectionEq();
				   
				   SIMPLE_Object.update_solutionFields_inNewTimeStep();;
			   
			} while(SIMPLE_Object.check_convergency());			
	}
	
    wallFunction_Object -> calculate_shearStress_for_allBoundaryWalls();
		    
	FlowModel_Object -> update_mu_turbulence();
	
	SIMPLE_Object.write_solutionFields();
	
	FlowModel_Object -> write_solutionField();
	
	wallFunction_Object -> write_solutionField();
	
	#include "destroyFields.H"
    
    cout << "End program." << endl << endl;
	
	return 0;
}
