  
double *fd3d_diffusion_3D ( int x_num, double x[], int y_num, double y[],  int z_num, double z[], double t, double dt, 
  double cflx, double cfly, double cflz, double *porosity , double *humidity, double value_Bond, double *coeff ( double *p , double *h, int nx, int ny, int nz ), 
  void bc (int x_num, int y_num, int z_num, double *c, double value ), double *h );
  
double fd3d_diffusion_cfl ( double k, int t_num, double t_min, double t_max, 
  int x_num, double x_min, double x_max, int y_num, double y_min, double y_max, int z_num, double z_min, double z_max );
 
 
double get_dtdx2 ( int t_num, double t_min, double t_max, 
  int x_num, double x_min, double x_max );
  
  
double *r8vec_linspace_new ( int n, double a_first, double a_last );

void timestamp ( );


double *porosity_co2(double *por_ini , double *por_current, double *por_last, double *co2, int m, int n);

double *porosity_cl(double *por_ini , double *por_current, double *por_last, double *co2, int m, int n);

double *ic_3D ( int x_num, int y_num, int z_num );

void write_vtk3D ( string output_filename, int n, double x[], int m, double y[], int h, double z[], double *c, string variable  );


double Phi0(double Pa, double Pc, double Pw, double a_c, double agr_c) ;
double Phit(double Phi0, double t ) ;



double calculate_rpo(double a_c);
double calculate_rpc(double a_c);

double Ccsh0(double Pcsh, double Ccao0);
double deltaPhic(double Pcsh, double Ccao0);



struct posCoord {
    int posX;
    int posY;
};
