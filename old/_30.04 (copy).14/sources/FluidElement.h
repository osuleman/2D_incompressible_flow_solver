class FluidElement
{
friend class Grid;
protected:
	double pos[3];
	double psi, zeta, rhs, psiManu;

public:
	// set value function	
	void set_position(double, double, double);
	void set_psi(double);	
	void set_zeta(double); 
	void set_rhs(double);
	void set_psiManu(double);	

	// get value function
	double get_psi();	
	double get_zeta();
	double get_rhs();
	double get_psiManu();
	double get_x();
	double get_y();
	double get_z();
	// clear1 FluidElement Variables (psi, zeta, rhs, psiManu)
	void clear();

};
