class Grid
{

friend class TestCase;
friend class TimeStepper;
protected:
	int nx , ny, nz, dims, iterations;
	double dx, dy, dz, l, w, t, lRef, uRef;
	double residualPsi[2], errorPsi[2], uBC[4];
	static const int nMax = 200;	
	FluidElement nodes[nMax][nMax];

public:
	Grid();
	Grid(int, int, double, double, double, double);
	void set_grid(int, int, double, double, double, double);
	FluidElement *get_node(int,int);
	void set(std::string,int,int,double);
	double get(std::string,int,int);
	void set(std::string,double);
	double get(std::string);	
	void set_psi_bc(std::string,double, double, double, double);
	void set_zeta_bc(std::string, double, double, double, double);
	void set_boundary_velocities(double,double,double,double);
	void set_psi_bc();	
	void set_zeta_bc();	
	double diff(int, int, std::string, std::string, int);	
	void calc_psi_error();
	int test_poisson(double,double);	
	int solve_poisson(double,double);
	double solve_dzeta_dt(int, int, double);
	void solve_zeta_tp1(double, double);
	void rk4_vorticity(double, double);	
	void solve(double,double,double,double,double,double,double,std::string);
	void create_grid_params_header(std::string);
	void save_grid_params(std::string);
	void save_vtk_field(std::string);
	void save_txt_field(std::string filenameSpec);

};

