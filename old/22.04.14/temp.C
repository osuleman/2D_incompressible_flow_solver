














void IC(double **A,M,N)

	for (int k = 0 ; k <= nx - 1 ; k++)
	{
		for (int l = 0 ; l <= ny - 1 ; l++)
		{
			x[k][l]   = x_0 + k*dx;
			y[k][l]   = y_0 + l*dy;
			phi[k][l] = sin(x[k][l]) + sin(y[k][l]);
			RHS[k][l] = 1; //rand();  // this is just a temporary RHS for testing.
		}
	}
 
