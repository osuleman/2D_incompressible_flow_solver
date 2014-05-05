class GridRefinementStudy
{

	//static const int nGridMax = 3;
	int nGrids, nMin, nMax, dn;	
	Grid grid;	 


public:
	GridRefinementStudy(int, int, int);
	void test_poisson(double omega, double tol);
	void test_derivative1();
        void test_derivative2();


};






