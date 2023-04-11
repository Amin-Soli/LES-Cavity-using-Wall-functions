
delete FlowModel_Object;

delete wallFunction_Object;

for(int j = 0; j <numberOfNodes_yDirection+1; j++)
	delete[] v[j];

for(int j = 0; j < numberOfNodes_yDirection; ++j)
{
	 delete[] u[j];
	 delete[] p[j];
	 delete[] mu_turbulence[j];
}

delete[] u;
delete[] v;
delete[] p;
delete[] mu_turbulence;
