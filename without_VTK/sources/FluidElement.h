class FluidElement
{
friend class Grid;
protected:
	double pos[3];
	double psi, zeta, zetaTp1, u, v, w, rhs, psiManu;

public:
	//constructor	
	FluidElement();
	//FluidElement *operator= (const FluidElement&);	
	// set value function	
	void set_position(double, double, double);
	void set_psi(double);	
	void set_zeta(double); 
	void set_rhs(double);
	void set_psiManu(double);
	void set_u(double);	
	void set_v(double);
	void set_w(double);

	// get value function
	double get_psi();	
	double get_zeta();
	double get_u();
	double get_v();
	double get_w();
	double get_rhs();
	double get_psiManu();
	double get_x();
	double get_y();
	double get_z();

	
	// clear1 FluidElement Variables (psi, zeta, rhs, psiManu)
	void clear();

};
