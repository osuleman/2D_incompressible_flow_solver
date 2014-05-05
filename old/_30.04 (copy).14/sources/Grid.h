class Grid
{


protected:
	int nx , ny, nz, dims, iterations;
	double dx, dy, dz, l, w;
	double residualPsi[2], errorPsi[2];	
	FluidElement nodes[200][200];

public:

	Grid(int,int,double,double);
	FluidElement *get_node(int,int);
	void set(std::string,int,int,double);
	double get(std::string,int,int);
	double get(std::string);	
	void calc_psi_error();
	int solve_poisson(double,double);
	void create_grid_params_header(std::string);
	void save_grid_params(std::string);
	void save_vtk_scaler_field(std::string);
	void set_psi_bc(double, double, double, double);

	


};

