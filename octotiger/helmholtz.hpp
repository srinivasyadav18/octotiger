#ifndef ___HELMHOLTZ___C
#define ___HELMHOLTZ___C

typedef struct {
	double p, e, cs, cv, abar, zbar, rho, T, s;
} eos_t;

void helmholtz_reset_counters();
void helmholtz_eos( eos_t* );
double helmholtz_max_iters();
void helmholtz_compute_T( eos_t* );
void helmholtz_ztwd( double*, double*, double, double );
void helmholtz_set_cgs_units( double cm, double g, double s, double T );
double helmholtz_iters_per_call();
double helmholtz_pressure_from_energy(double rho, double E, double abar, double zbar);
double helmholtz_T_from_energy(double rho, double E, double abar, double zbar);
std::pair<double, double> helmholtz_pressure_and_soundspeed(double rho, double E, double abar, double zbar);
void helmholtz_initialize();

#endif
