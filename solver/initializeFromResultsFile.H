
stringstream ss;
ss << start_time;
string str = ss.str();

int z;

// initialize u
string U_Name = "results/t=" + str + "/U.txt";
const char* U_fileName = U_Name.c_str();

ifstream U_file (U_fileName); 

string U_line;
z = numberOfNodes_yDirection*(numberOfNodes_xDirection+1); 
string U[z]; 
string U_word; 

j=0;
if(U_file.is_open())
	while (getline(U_file, U_line) && j<z) 
	{
		stringstream U_str(U_line);
		i=0;
		while (getline(U_str, U_word, ','))
		{
			if(i==2)
				U[j] = U_word; 
			i++;
		}
		j++;
	}
	
z=0;
for(j=0;j<numberOfNodes_yDirection;j++)
	for(i=0;i<numberOfNodes_xDirection+1;i++)
	{
		u[j][i] = stod(U[z]); 
		z++;
	}
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// initialize v
string V_Name = "results/t=" + str + "/V.txt";
const char* V_fileName = V_Name.c_str();

ifstream V_file (V_fileName); 

string V_line;
z = numberOfNodes_xDirection*(numberOfNodes_yDirection+1); 
string V[z]; 
string V_word; 

j=0;
if(V_file.is_open())
	while (getline(V_file, V_line) && j<z) 
	{
		stringstream V_str(V_line);
		i=0;
		while (getline(V_str, V_word, ','))
		{
			if(i==2)
				V[j] = V_word; 
			i++;
		}
		j++;
	}
	
z=0;
for(i=0;i<numberOfNodes_xDirection;i++)
    for (j=0;j<numberOfNodes_yDirection+1;j++)
	{
		v[j][i] = stod(V[z]); 
		z++;
	}
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// initialize p
string P_Name = "results/t=" + str + "/P.txt";
const char* P_fileName = P_Name.c_str();

ifstream P_file (P_fileName); 

string P_line;
z = numberOfNodes_yDirection*numberOfNodes_xDirection; 
string P[z]; 
string P_word; 

j=0;
if(P_file.is_open())
	while (getline(P_file, P_line) && j<z) 
	{
		stringstream P_str(P_line);
		i=0;
		while (getline(P_str, P_word, ','))
		{
			if(i==2)
				P[j] = P_word; 
			i++;
		}
		j++;
	}
	
z=0;
for(j=0;j<numberOfNodes_yDirection;j++)
	for(i=0;i<numberOfNodes_xDirection;i++)
	{
		p[j][i] = stod(P[z]); 
		z++;
	}
