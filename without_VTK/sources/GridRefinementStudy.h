class GridRefinementStudy
{

	//static const int nGridMax = 3;
	int nGrids, nMin, nMax, dn;	
	Grid grid;	 


public:
	GridRefinementStudy(int, int, int);
	void test_poisson_sor(double,double);
	void test_poisson_adi(double);
	void test_derivative1();
        void test_derivative2();


};






