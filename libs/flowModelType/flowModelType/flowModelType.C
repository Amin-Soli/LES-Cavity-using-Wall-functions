#include "flowModelType.H"
#include "../LESmodel/LESmodel/LESmodel.C"
#include "../LaminarModel/LaminarModel.H"

flowModelType* flowModelType::New(const char *path, double **U, double **V , double **mu_SGS)
{
	string flowModelType_Name;

	ifstream file(path);

	string searchForFlowModelType_Name = "FlowModelType";

	string lineOfText;

	for(;;)
	{
		getline(file, lineOfText);

		if (file.eof()) break;

		if (lineOfText.find(searchForFlowModelType_Name, 0) != string::npos)
		  file >> flowModelType_Name;
    }

    file.close();
    
    flowModelType* myObject;

    if (flowModelType_Name=="turbulent")
	    myObject = LESmodel::New(path, U, V, mu_SGS);
	    
	if (flowModelType_Name=="laminar")
	    myObject = new LaminarModel(path, U, V, mu_SGS);
	    
	return myObject;
}

