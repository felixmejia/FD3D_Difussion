#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream> 
#include <string> 
#include <stdlib.h>
#include <stdio.h>
#include <ctime>

using namespace std;

#include "fd3d_diffusion.hpp"
#include "tic_toc.h"

int main ();
void fd3d_diffusion();

void bc_3D (  int x_num, int y_num, int z_num, double *c, double value );

double *ic_3D ( int x_num, int y_num, int z_num, double value );

double *coeff_diff_co2( double *p , double *h, int nx, int ny, int nz);

double *coeff_dif_cl( double *p , double *h, int nx, int ny, int nz);


double *set_porosity_random( int x_num, int y_num, int z_num, double min, double max );
double *set_value_3D( int x_num, int y_num, int z_num, double value );

double *rp( int nx, int ny, int nz, double rp_ref, double value_rpo, double value_rpc, double a_c, double * carbonation);

 double *set_porosity_ini( double a_c, double ag_c, double da, double dag, double dc, int x_num, int y_num, int z_num);
 
 double tortuosity(double phi);
 double constrictivity(double rp);


double h_2(double Uc, double R, double Tref, double T);
double h_3(double tref, double te, double m);

double Dfc_ref(double a_c);
double Dfc(double a_c, double h, double hc, double Uc, double R, double Tref, double T, double tref, double te, double m );

 double * porosity_time(double * porosity_ini, double time, int nx, int ny, int nz );
 

double *phi_we(double *h, double *T, double te, double a_c);


double *multiply(int nx, int ny, int nz, double *A, double *B);

 
 
string IntToString ( int number )
{
  std::ostringstream oss;

  // Works just like cout
  oss<< number;

  // Return the underlying string
  return oss.str();
  

}




//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FD2D_DIFFUSION_PRB.
//
//  Discussion:
//
//    FD2D_DIFFUSION_TEST tests the FD2D_DIFFUSION library.
//
//  
//  Modified:
//
//    26 October 2017
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
{
  cout << "\n";
  cout << "FD3D_DIFFUSION_TEST:\n";
  cout << "  C++ version.\n";
  cout << "  Test the FD3D_DIFFUSION library.\n";

  fd3d_diffusion( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FD3D_DIFFUSION_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void fd3d_diffusion( )
//****************************************************************************80
//
//  Purpose:
//
//    FD2D_DIFFUSION_TEST.
////
//  
//  Modified:
//
//    26 October 2017
//
//  Author

//
//    Félix Armando Mejía Cajicá
//


{
  double cflx;
  double cfly;
  double cflz;
  
  double dt;
  double *chloride, *carbonation;
  double *chloride_new;
  double *porosity;
  double *humidity, *humidity_new;
  double *temperature, *temperature_new;

double *tortuosity_ini;
double *constructivity_ini;
double *rel_tortini_constrini;
double *Dh, *Dh_car;
double *Dfc, *Dfc_car;

  int i, h;
  int j;
  int vt;
  double k;
  double *t;
  double t_max;
  double t_min;
  int t_num;
  double *x;
  double x_max;
  double x_min;
  int x_num;


  double *y;
  double y_max;
  double y_min;
  int y_num;
  
  double *z;
  double z_max;
  double z_min;
  int z_num;
  
  
 
  
  cout << "\n";
  cout << "FD3D_DIFFUSION:\n";
  cout << "  Compute an approximate solution to the diffusion-dependent\n";
  cout << "  three dimensionals diffusion equation:\n";
  cout << "\n";
  cout << "    dU/dt - K * d2U/dx2 = f(x,t)\n";
  cout << "\n";
  cout << "  Run a simple test case.\n";

//
//  X_NUM is the number of equally spaced nodes to use between Xmin and Xmax.
//
  x_num = 256;
  x_min = 0.0;
  x_max = 200.0;
  x = r8vec_linspace_new ( x_num, x_min, x_max );
  
  
//  Y_NUM is the number of equally spaced nodes to use between Ymin and Ymax.
//
  y_num = 256;
  y_min = 0.0;
  y_max = 200.0;
  y = r8vec_linspace_new ( y_num, y_min, y_max );
//

//  Z_NUM is the number of equally spaced nodes to use between Ymin and Ymax.
//
  z_num = 400;
  z_min = 10.0;
  z_max = 440.0;
  z = r8vec_linspace_new ( z_num, z_min, z_max );
//

// Num Beam
	int numBeam = 4;
	int radius = 8; // millimeters
    posCoord *posBeam;
    posBeam = ( posCoord * ) malloc ( numBeam * sizeof ( posCoord ) );
    
    posBeam[0].posX = 50; //mm
    posBeam[0].posY = 50; //mm
    
    posBeam[1].posX = 50; //mm
    posBeam[1].posY = y_max - 50; //mm
    
    posBeam[2].posX = x_max - 50; //mm
    posBeam[2].posY = y_max - 50; //mm
    
    posBeam[3].posX = x_max - 50; //mm
    posBeam[3].posY = 50; //mm
    
    posCoord posVector;
    posVector = get_Coordenadas (numBeam, posBeam, x_num, x_min, x_max, y_num, y_min, y_max);
     
//
//
//  T_NUM is the number of equally spaced time points between tmin and tmax.
//
  t_num = 15501;
  t_min = 0.0;
  t_max = 80.0;
  dt = ( t_max - t_min ) / ( double ) ( t_num - 1 );
  t = r8vec_linspace_new ( t_num, t_min, t_max );
//
//  Get the CFL coefficient.
//
	cflx = get_dtdx2 ( t_num, t_min, t_max, x_num, x_min, x_max);
	cfly = get_dtdx2 ( t_num, t_min, t_max, y_num, y_min, y_max);
	cflz = get_dtdx2 ( t_num, t_min, t_max, z_num, z_min, z_max);
	
	k=1.0;
// Verify the stability criterion value
  double verify = fd3d_diffusion_cfl( k, t_num, t_min, t_max, x_num, x_min, x_max, y_num, y_min, y_max, z_num, z_min, z_max);


  //*******************************************************************************
  // Set initial value of the concentration
  //*******************************************************************************
  double CInitial = 0.0;
  double CBound_chloride= 0.045;
  double CBound_carbonation= 4.0;

  chloride = ic_3D ( x_num, y_num, z_num );
  bc_3D ( x_num, y_num, z_num, c , CBound_chloride);

  carbonation = ic_3D ( x_num, y_num, z_num );
  bc_3D ( x_num, y_num, z_num, c , CBound_carbonation);
  
  
  //*************************************************************************
  // Valores Iniciales 	
  //*************************************************************************
  double a_c=0.4;
  double ag_c=6;
  double da = 0.975 ;
  double dag = 111;
  double dc = 111;
  double value_hum = 0.05;
  double humedad_rel = 0.75;
	double dens_concrete=3450.0;
	double calor_especifico = 5.25;
    	double conductivity=1.5;
	double temp_bound=295;

  // Parámetros para la humedad
  double hum_bound=0.75;
  double Dh_ref = 3.0*pow(10.0,-10.0);  //Coeficiente de difusion de la humedad de referencia m2/s
  double Dfc_ref = Dfc_reference(a_c);
  double alfa0 = 0.05;
  double hc = 0.75; /// Humedad relativa del poro
  double n = 11.0;  // Exponente de extension de la caida Dh
  double Um = 20.3;  /// Energia de activacion por hidratación
  double R = 8.314*pow(10.0,-3.0);  /// Constante de los gases  kJ/K. mol
  double Tref = 296.0 ; // Temperatura de referencia [K]

  double Uc = 41.8;  /// Uc Energia de Activación   41.8, 44.6 o 32 kJ/mol para a/c 0.4, 0.5 y 0.6
  double tref=1.0;  // tiempo en lo que se ha medido la difusividad del concreto Dfc_ref
  double m = 0.2;  /// m factor de reduccion de la edad m=0.2

  porosity = set_porosity_ini( a_c, ag_c, da, dag, dc, x_num, y_num, z_num);
  humidity = set_value_3D( x_num, y_num, z_num, value_hum);

  double rpo = calculate_rpo(a_c);
  double rpc = calculate_rpc(a_c);
  double *grad_carbonatacion=ic_3D ( x_num, y_num, z_num );

  double *rp_ini = rp(nx, ny, nz, rpo, rpc, grad_carbonatacion);
  double *rp;
  double *constrictivity = constrictivity(nx, ny, nz, rp_ini);
  double *tortuosity = tortuosity(nx, ny, nz, porosity);
  double *rel_const_tortu_ini = divide(nx, ny, nz, tortuosity, constrictivity);
  double *rel_const_tortu;
  double *influ_constort;
	
  free(constrictivity);
  free(tortuosity);
  free(rp_ini);

  bc_3D ( x_num, y_num, z_num, humidity , humedad_rel);

  vt=0;
  string strNumber = IntToString(vt);
  string fileName = "Simulacion_" + strNumber + ".vtk";
  
 // write_vtk3D ("Porosity.vtk", x_num,  x, y_num, y, z_num, z, porosity, "porosity" );
  write_vtk3D(fileName, x_num,  x, y_num, y, z_num, z, c, "chloride");
  for ( vt = 1; vt < t_num; vt++ )
  {
	tic();
	cout << " time: " <<  vt <<"/" << t_num << "  \n";

	// Calculation of the Temperature and Humidity in the time
	//
	//
	temperature_new = fd3d_difussion_Temp_3D(x_num, x, y_num, y, z_num, z, t[vt], dt, 
  			 cflx, cfly, cflz, dens_concrete, calor_especifico, conductivity, temperature, temp_bound, bc_3D);


	Dh = Dh_3D(nx, ny, nz, Dh_ref, alfa0, humidity, hc, Um, R, Tref, temperature,  t[vt]);

	Dfc = Dfc_3D(nx, ny, nz, Dfc_ref, humidity, hc, Uc, R, Tref, temperature, tref, t[vt], m);


	rp = rp(nx, ny, nz, rpo, rpc, grad_carbonatacion);

  	constrictivity = constrictivity(nx, ny, nz, rp);
	tortuosity = tortuosity(nx, ny, nz, porosity);
	rel_const_tortu = divide(nx, ny, nz, constrictivity, tortuosity);
	influ_constort = multiply(nx, ny, nz, rel_const_tortu_ini, rel_const_tortu);

	Dh_car = multiply(nx, ny, nz, influ_constort, Dh);
 	Dfc_car = multiply(nx, ny, nz, influ_constort, Dfc);
	

	humidity_new = fd3d_difussion_Humidity_3D ( x_num, x, y_num, y, z_num, z, t[vt], dt, 
  			 cflx, cfly, cflz, Dh_car, humidity, hum_bound, bc_3D);
	
	


	c_new = fd3d_diffusion_3D ( x_num, x, y_num, y, z_num, z, t[vt-1], dt, cflx, cfly, cflz, porosity, humidity, CBoundary, coeff_diff_co2, bc_3D, c );
 
    for ( h = 0; h < z_num; h++ )
    	{
		for ( i = 0; i < x_num; i++ )
		{
			for ( j = 0; j < y_num; j++ )
			{
				c[i*y_num+j+h*(x_num*y_num)] = c_new[i*y_num+j+h*(x_num*y_num)];
			}
		}
	}
	toc();
	//cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl; 
    strNumber = IntToString(vt);
    fileName = "Simulacion_" + strNumber + ".vtk";
    //
	//  Write the data to VTK files.
	//
  //  write_vtk3D(fileName, x_num,  x, y_num, y, z_num, z, c, "chloride");
    free(c_new);    
    toc();
	
  }
  free (c);
  delete [] t;
  delete [] x;
  delete [] y;
  return;
}



//****************************************************************************80

void bc_3D ( int x_num, int y_num, int z_num, double *V, double value )

//****************************************************************************80
//
//  Purpose:
//
//    BC_3D evaluates the boundary conditions for problem.
//
//  Modified:
//
//    26 December 2017
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//
//  Parameters:
//
//    Input, int X_NUM, the number of nodes in X.
//
//    Input, int Y_NUM, the number of nodes in Y.
//
//    Input, int ZNUM, the number of nodes in Z.
//
//    Input, double C[X_NUM,Y_NUM,Z_NUM], the current concentration values, after boundary
//    conditions have been imposed.
//
//    Input, double value, the current time.
//
{
	
	int i,j, k;
	
 
	for ( k = 0; k < z_num; k++ )
	{
		for ( i = 0; i < x_num; i++ )
		{
			c[i*y_num+k*(x_num*y_num)] = value;
			c[i*y_num+y_num-1+k*(x_num*y_num)] = value;
			//cout << "C[" <<  i << "," << x_num-1 << "," << k << "]= " << h[i*x_num+x_num-1+k*(x_num*y_num)] << "\n"  ;
			//cout << "C[" <<  i << "," << 0 << "," << k << "]= " << h[i*x_num+k*(x_num*y_num)] << "\n"  ;
		}
		for ( j = 1; j < y_num-1; j++ )
		{
			c[j+k*(x_num*y_num)] = value;
			c[(x_num-1)*y_num+j+k*(x_num*y_num)] = value;
			//cout << "C[" <<  0 << "," << j << "," << k << "]= " << h[j+k*(x_num*y_num)] << "\n"  ;
			//cout << "C[" <<  y_num-1 << "," << j << "," << k << "]= " << h[(y_num-1)*x_num+j+k*(x_num*y_num)] << "\n"  ;
		}
		
	}
  return;
}



//****************************************************************************80

double *ic_3D ( int x_num, int y_num, int z_num )

//****************************************************************************80
//
//  Purpose:
//
//    IC_3D evaluates the initial condition for problem.
//
//  Modified:
//
//    15 December 2017
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Parameters:
//
//    Input, int X_NUM, the number of nodes in X.
//
//    Input, int Y_NUM, the number of nodes in Y.
//
//    Input, int ZNUM, the number of nodes in Z.
//
//    Input, double T, the initial time.
//
//    Output, double C[X_NUM,Y_NUM,Z_NUM], the heat values at the initial time.
//
{
  int i,j, k;
  double *c;

  c = ( double * ) malloc ( x_num * y_num * (z_num) *sizeof ( double ) );
 
	for ( k = 0; k <z_num; k++ )
	{
		for ( i = 0; i < x_num; i++ )
		{
			for  ( j = 0; j < y_num; j++ )
			{
				c[i*y_num+j+k*(x_num*y_num)] = 0.0;
				//cout << "C[ " << i << "," << j << "," << k << "]= " << c[i*x_num+j+k*(x_num*y_num)] << "\n"  ;
				//getchar();
			}
		}
	}
  return c;
}



//****************************************************************************80

 double *set_porosity_random( int x_num, int y_num, int z_num, double min, double max )
 
//****************************************************************************80
//
//  Purpose:
//
//    SET_POROSITY evaluates the right hand side for problem 1.
//
//  Modified:
//
//    15 December 2017
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Parameters:
//
//    Input, int X_NUM, the number of nodes in X.
//
//    Input, int Y_NUM, the number of nodes in Y.
//
//    Input, int ZNUM, the number of nodes in Z.
//
//    Input, double min, the min Value .
//
//    Input, double max, the max Value .
//
//    Output, double p[X_NUM,Y_NUM,Z_NUM], the porosity.
//
{
  int i,j, k;
  double *c;
  double r;
  double deltaM=(max - min)/RAND_MAX ;

	

  c = ( double * ) malloc ( x_num * y_num * (z_num) *sizeof ( double ) );
	for ( i = 0; i < x_num; i++ )
	{
		for  ( j = 0; j < y_num; j++ )
		{
			for ( k = 0; k < z_num; k++ )
			{
				
				r = std::rand()*1.0/RAND_MAX  ;
				
				if(r<0.3)
				{
					c[i*y_num+j+k*(x_num*y_num)] = min + (r * deltaM*RAND_MAX);
					//cout << "r = " << c[i*y_num+j+k*(x_num*y_num)] << "\n" ;
				}
				else
				{
					//cout << "piedra " << p << "\n" ;
					c[i*y_num+j+k*(x_num*y_num)] = 0.0;	
				}
				//c[i*x_num+j+k*(x_num*y_num)]=0.0;
				//cout << "Porosity[ " << i << "," << j << "," << k << "]= " << c[i*x_num+j+k*(x_num*y_num)] << "\n"  ;
			}
		}
	}
  return c;
}


//****************************************************************************80

 double *set_porosity_ini( double a_c, double ag_c, double da, double dag, double dc, int x_num, int y_num, int z_num)
 
//****************************************************************************80
//
//  Purpose:
//
//    SET_POROSITY ini evaluates the right hand side for problem 1.
//
//  Modified:
//
//    15 December 2017
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Parameters:
//
//    Input, double a_c, relation water/cement.
//
//    Input, double ag_c, relation aggregates/cement.
//
//    Input, double da, density of water.
//
//    Input, double dag, density of aggregates.
//
//    Input, int X_NUM, the number of nodes in X.
//
//    Input, int Y_NUM, the number of nodes in Y.
//
//    Input, int ZNUM, the number of nodes in Z.
//
//    Input, double min, the min Value .
//
//    Input, double max, the max Value .
//
//    Output, double p[X_NUM,Y_NUM,Z_NUM], the porosity.
//
{
  int i,j, k;
  double *c;
  double r;

	c = ( double * ) malloc ( x_num * y_num * (z_num) *sizeof ( double ) );
	r=(a_c*dc/da)/(a_c*dc/da+ag_c*dc/dag+1.0);
	for ( i = 0; i < x_num; i++ )
	{
		for  ( j = 0; j < y_num; j++ )
		{
			for ( k = 0; k < z_num; k++ )
			{
				random = std::rand()*1.0/RAND_MAX  ;
				c[i*y_num+j+k*(x_num*y_num)] = r*random;
			}
		}
	}
  return c;
}


//****************************************************************************80

 double *set_value_3D( int x_num, int y_num, int z_num, double value )
 
//****************************************************************************80
//
//  Purpose:
//
//    set_value_3D stablish the value to all cube 3D.
//
//  Modified:
//
//    15 January 2018
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Parameters:
//
//    Input, int X_NUM, the number of nodes in X.
//
//    Input, int Y_NUM, the number of nodes in Y.
//
//    Input, int ZNUM, the number of nodes in Z.
//
//    Input, double min, the min Value .
//
//    Input, double max, the max Value .
//
//    Output, double pointer to Vector of  Cube 3D.
//
{
	int i,j, k;
	double *c;
 
  	c = ( double * ) malloc ( x_num * y_num * z_num *sizeof ( double ) );
	for ( i = 0; i < x_num; i++ )
	{
		for  ( j = 0; j < y_num; j++ )
		{
			for ( k = 0; k < z_num; k++ )
			{
				c[i*y_num+j+k*(x_num*y_num)] = value;
			}
		}
	}
  	return c;
}

//****************************************************************************80

double *coeff_diff_co2( double *p , double *h, int nx, int ny, int nz)

//****************************************************************************80
//
//  Purpose:
//
//    coefficient_diffusion_co2 calculate the coefficient diffusion of carbon dioxide.
//
//
//  Modified:
//
//    26 October 2017
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Parameters:
//
//    Input, double p, the porosity of hardered binder.
//
//    Input, double h, the humidity.
//
//    Output, double coefficient_diffusion_co2, Coefficient Diffusion.
//
{
	int i,j,k;
	double *Dco2;
	cout << "Entrando coefficentes" << endl;
	Dco2 = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) );
	cout << "Apartada Memoria coefficentes" << endl;
	#pragma acc kernels copyin(p[0:nx * ny * nz], h[0:nx * ny * nz]) copyout(Dco2[0:nx * ny * nz])
	
	#pragma acc loop collapse(3) independent	
	for(i = 0; i < nx ; i++ )
	{
	  for ( j = 0; j < ny ; j++ )
	  {
		for ( k = 0; k < nz ; k++ )
		{
	      	
		  Dco2[i*ny+j+k*(nx*ny)] = 1.64 * pow(10.0, -6.0)*pow(p[i*ny+j+k*(nx*ny)], 1.8)*pow(1.0-h[i*ny+j+k*(nx*ny)], 2.2); 
		  //cout << "Dco2[ " << i << "," << j << "," << k << "]= " << Dco2[i*ny+j+k*(nx*ny)] << "\n"  ;
		
		}
		  
	  }
	}
	cout << "Calculados coefficentes" << endl;
	return Dco2;
}





//****************************************************************************80

double *coeff_dif_cl( double *p , double *h, int nx, int ny, int nz)

//****************************************************************************80
//
//  Purpose:
//
//    coefficient_diffusion_co2 calculate the coefficient diffusion of carbon dioxide.
//
//
//  Modified:
//
//    26 October 2017
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Parameters:
//
//    Input, double p, the porosity of hardered binder.
//
//    Input, double h, the humidity.
//
//    Output, double coefficient_diffusion_co2, Coefficient Diffusion.
//
{
	int i,j,k;
	double *Dco2;
	Dco2 = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) );
	
	for(i = 0; i < nx ; i++ )
	{
	  for ( j = 0; j < ny ; j++ )
	  {
		for ( k = 0; k < nz ; k++ )
		{
		  Dco2[i*ny+j+k*(nx*ny)] = 1.64 * pow(10.0, -6.0)*pow(p[i*ny+j+k*(nx*ny)], 1.8)*pow(1.0-h[i*ny+j+k*(nx*ny)], 2.2); 
		}
		  
	  }
	}
	return Dco2;
}



//****************************************************************************80

//La porosidad del concreto (P) o total de huecos en el material
//compuesto, se ha modelado como una función de: (a/c) la relación
//agua/cemento, el grado de hidratación del cemento (h), el volumen
//de aire atrapado (A), las cantidades de agregados fino (arena, Af) y
//grueso (grava, Ag), y del cemento (c); y las gravedades específicas
//de los agregados (ρf y ρg), en el modelo matemático de la Ecuación
//(2) de (Solís & Moreno, 2006) para concretos con agregados
//comunes de densidades normales, cuando la relación
//agua/cemento es alta y el grado de hidratación es bajo, la pasta de
//cemento tendrá una alta porosidad capilar; contendrá un número
//relativamente grande de poros amplios y bien conectados, y por lo
//tanto, su coeficiente de permeabilidad será alto. En cuanto avanza
//la hidratación, la mayoría de los poros serán reducidos a un
//tamaño pequeño (100 nm o menos) y también perderán sus
//interconexiones, de manera que la permeabilidad se abatirá.
double porosidad_ini(double a_c, double hidratacion, double air_c, double agFin, double denagFin, 
					double grava, double dengrava)
//****************************************************************************80
//
//  Purpose:
//
//    porosidad_ini calculate the initial porosity.
//
//
//  Modified:
//
//    16 December 2017
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Parameters:
//
//    Input, double a_c, the relation water/cement.
//
//    Input, double hidratacion, the hidratation.
//
//    Input, double agFin, the hidratation.
//
//    Input, double denagFin, the hidratation.
//
//    Input, double grava, the hidratation.
//
//    Input, double dengrava, the hidratation.
//
//    Output, double porosidad_ini, Coefficient Diffusion.
//					
{
	
	double P;
	
	P = (a_c - 0.36 * hidratacion + air_c)/(0.317 + (1/denagFin)*agFin+(1/dengrava)*grava+air_c);
	
	return P;


}



//****************************************************************************80
/// Funcion Calcula el radio pico de los poros segun el nivel de carbonatacion
/////////////////////////////////////
/// nx, ny, nz  
/// value_rpo depende de la a/c 
/// value_rpc depende de a/c
/// grad_carb Grado de la carbonatacion
//////////

double *rp(int nx, int ny, int nz, double value_rpo, double value_rpc, double * grad_carb)
{
	double *value_rp;
	int i,j,k;
	
	value_rp = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) );
	for(i = 0; i < nx ; i++ )
	{
		for ( j = 0; j < ny ; j++ )
		{
			for ( k = 0; k < nz ; j++ )
			{
				value_rp[i*ny+j+k*(nx*ny)]=5.0*pow(10.0,-8.0)*((value_rpc-value_rpo)*grad_carb[i*ny+j+k*(nx*ny)] + value_rpo);
			}
		}
	}
	return value_rp;
}

double *tortuosity(int nx, int ny, int nz, double *porosity)
{
	double *result;
	double b1, b2, b3, b4;
	b1=1.5; 
	b2= 8.0;
	b3=0.25;
	b4=2.5;
	
	result = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) );
	
	for(i = 0; i < nx ; i++ )
	{
		for ( j = 0; j < ny ; j++ )
		{
			for ( k = 0; k < nz ; j++ )
			{
				result[i*ny+j+k*(nx*ny)] = -b1 * tanh(b2*(porosity[i*ny+j+k*(nx*ny)]-b3)) + b4;
			}
		}
	}
	return result;
	
}


double *constrictivity(int nx, int ny, int nz, double *rp)
{
	double *result;
	double c1, c2, c3, c4;
	c1=0.395;
	c2=4.0;
	c3=5.95;
	c4=0.405;

	result = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) );
	
	for(i = 0; i < nx ; i++ )
	{
		for ( j = 0; j < ny ; j++ )
		{
			for ( k = 0; k < nz ; j++ )
			{
				result[i*ny+j+k*(nx*ny)] = c1*tanh(c2*(log(rp[i*ny+j+k*(nx*ny)])+c3))+c4;
			}
		}
	}
	return result;
} 

double *multiply(int nx, int ny, int nz, double *A, double *B)
{
	double *result;
	
	result = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) );
	
	for(i = 0; i < nx ; i++ )
	{
		for ( j = 0; j < ny ; j++ )
		{
			for ( k = 0; k < nz ; j++ )
			{
				result[i*ny+j+k*(nx*ny)] = A[i*ny+j+k*(nx*ny)]*B[i*ny+j+k*(nx*ny)];
			}
		}
	}
	return result;
} 

double *divide(int nx, int ny, int nz, double *tortuosity, double *constrictivity)
{
	double *result;
	
	result = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) );
	
	for(i = 0; i < nx ; i++ )
	{
		for ( j = 0; j < ny ; j++ )
		{
			for ( k = 0; k < nz ; j++ )
			{
				result[i*ny+j+k*(nx*ny)] = tortuosity[i*ny+j+k*(nx*ny)]/constrictivity[i*ny+j+k*(nx*ny)];
			}
		}
	}
	return result;
} 

double h_1(double h, double hc)
{
	double res;
	res = 1.0+ (pow(1-h, 4.0)/pow(1-hc, 4.0));
	return res;

}

double h_2(double Uc, double R, double Tref, double T)
{
	double res;
	res = exp(Uc/(R*Tref)-Uc/(R*T));
	return res;

}

double h_3(double tref, double te, double m)
{
	double res;
	res = pow(tref/te, m);
	return res;

}


// Coeficiente de Difusion de cloruro de referencia
// a_c relacion agua cemento

double Dfc_reference(double a_c)
{
	double res;
	double exponente = 1.776 + 1.364*a_c;
	double denominador = 3.1536*pow(10.0,13.0);
	
	res = pow(10.0, exponente)/denominador;
	return res;
	
}


// Coeficiente de Difusion de cloruro 


double Dfc(double a_c, double h, double hc, double Uc, double R, double Tref, double T, double tref, double te, double m )
{
/// Uc Energia de Activación   41.8, 44.6 32 kJ/mol paqra a/c 0.4, 0.5 y 0.6
/// R constante de los gases 8.314 * 10^-3 kJ/Kmol
/// Tc,ref y tref tempersatura y tiempo en lo que se ha medido la difusividad del concreto Dfc,ref
/// m factor de reduccion de la edad m=0.2
/// hc es el nivel de humedad se toma 0.75
/// T es la temperatura actual en el concreto
/// h es la humedad actual en el poro de concreto


	double res;
	
	res = Dfc_ref(a_c)*h_1(h,hc)*h_2(Uc,R,Tref,T)*h_3(tref,te,m);
	return res;
	
}

double Dfc_3D(int nx, int ny, int nz, double Dfc_ref, double *h, double hc, double Uc, double R, double Tref, double *T, double tref, double te, double m )
{
/// Uc Energia de Activación   41.8, 44.6 32 kJ/mol paqra a/c 0.4, 0.5 y 0.6
/// R constante de los gases 8.314 * 10^-3 kJ/Kmol
/// Tc,ref y tref tempersatura y tiempo en lo que se ha medido la difusividad del concreto Dfc,ref
/// m factor de reduccion de la edad m=0.2
/// hc es el nivel de humedad se toma 0.75
/// T es la temperatura actual en el concreto
/// h es la humedad actual en el poro de concreto
	
	double *result;
	
	result = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) );
	double h3=h_3(tref,te,m);
	for(i = 0; i < nx ; i++ )
	{
		for ( j = 0; j < ny ; j++ )
		{
			for ( k = 0; k < nz ; j++ )
			{
				result[i*ny+j+k*(nx*ny)] = Dfc_ref*h_1(h[i*ny+j+k*(nx*ny)],hc)*h_2(Uc,R,Tref,T[i*ny+j+k*(nx*ny)])*h3;
			}
		}
	}
	return result;
}





//CALCULO DE COEFICIENTE DE DIFUSION DE LA HUMEDAD


/// Funcion de la influencia de la humedad relativa del Poro
/////
///  alfa0 valor entre 0.025 y 0.1 se tomo 0.05 para este modelo
/// h es la humedad
/// hc la humedad relativa se toma 0.75
/// n es un valor entre 6 y 16 que caracteriza la extension de la caida del Dh


double f_1(double alfa0, double h, double hc, double n=11.0)
{
	double res;
	
	
	res = alfa0 + (1.0-alfa0)/(1.0+(pow(1.0-h, n)/pow(1.0-hc, n))
	return res;

}


///  Funcion de la influencia de la Temoperatura 
//////
/// Um es la energia de activacion por hidratacion 20.3 kJ/mol
/// R constante de los gases 8.314*.10^-3 kJ/Kmol 
/// Tref es la temperatura de referencia 296K
double f_2(double Um, double R, double Tref, double T)
{
	double res;
	res = exp(Um/(R*Tref)-Um/(R*T));
	return res;

}


/// Funcion de influencia de la edad
/////
/// te edad del concreto en dias
////
double f_3(double te)
{
	double res;
	res = 0.3+sqrt(13.0/te);
	return res;

}



/// Dh_ref coeficiente de difusion de la humedad de refrencia se toma 3*10^-10 m2/s
/// 
double Dh(double h, double hc, double Um, double R, double Tref, double T,  double te)
{
	double res;
	double n=11.0;
	
	res = Dh_ref*f_1(alfa0, h, hc , n)*f_2(Um,R,Tref,T)*f_3(te);
	return res;
	
}

/// Dh_ref coeficiente de difusion de la humedad de refrencia se toma 3*10^-10 m2/s
/// 
double *Dh_3D(int nx, int ny, int nz, double Dh_ref, double alfa0, double n, double *h, double hc, double Um, double R, double Tref, double *T,  double te)
{
	double *result;
	
	
	result = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) );
	double f3=f_3(te);
	for(i = 0; i < nx ; i++ )
	{
		for ( j = 0; j < ny ; j++ )
		{
			for ( k = 0; k < nz ; j++ )
			{

				result[i*ny+j+k*(nx*ny)] = Dh_ref*f_1(alfa0, h[i*ny+j+k*(nx*ny)], hc , n)*f_2(Um,R,Tref,T[i*ny+j+k*(nx*ny)])*f3;
			}
		}
	}
	return result;
}



//*******************************************************************************
//// Calculo de Porosity in the time
//*******************************************************************************
//
//  Purpose:
//
//    porosity_time calculate the Porosity of concrete in the time
//
//
//  Modified:
//
//    15 March 2018
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Parameters:
//

double * porosity_time(double * porosity_ini, double time, int nx, int ny, int nz )
{
	double bheta, n_min;
	double *res;
	
	n_min=17.8  ;  // mol/m3dia;  1.78*pow(10.0,-5)  mol/cm3dia;
	bheta= 4.19 * pow(10.0, -6.0);   // m3/mol
	res = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) );
	for(i = 0; i < nx ; i++ )
	{
		for ( j = 0; j < ny ; j++ )
		{
			for ( k = 0; k < nz ; j++ )
			{
				res[i*ny+j+k*(nx*ny)]=porosity_ini[i*ny+j+k*(nx*ny)]*exp(-bheta*n_min*time);
			}
		}
	}
	return res;
}
//*******************************************************************************
//// Calculo de Phi_We, es la fraccion de volumen de agua en poro
//*******************************************************************************
//
//  Purpose:
//
//    Phi_We calculate the fraction of water volume in poro.
//
//
//  Modified:
//
//    26 February 2018
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Parameters:
//
//    Input, double p, the porosity of hardered binder.
//
//    Input, double h, the humidity.
//
//    Input, double T, the Temperature.
//
//    Input, double te, time Hidratation.
//
//    Input, double nx, relation water cement.
//
//    Input, int nx, the number of nodes in X.
//
//    Input, int ny, the number of nodes in Y.
//
//    Input, int nz, the number of nodes in Z.
//
//    Output, double phi_we, fraction of water volume in poro
//
double *phi_we(double *h, double *T, double te, double a_c, int nx, int ny, int nz)
{
	double *result;
	result = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) );
	
	double nw=(2.5+15.0/te)*(0.33+2.2*a_c);
	double Vm=(0.068-0.22/te)*(0.85+0.45*a_c);
	double C, Km;

	for(i = 0; i < nx ; i++ )
	{
	  for ( j = 0; j < ny ; j++ )
	  {
		for ( k = 0; k < nz ; j++ )
		{
			C=exp(855.0/T[i*ny+j+k*(nx*ny)]);
			Km=(C*(1.0-1.0/nw)-1.0)/(C-1.0);
			result[i*ny+j+k*(nx*ny)]=(C*Km*Vm*h[i*ny+j+k*(nx*ny)])/((1.0-Km*h[i*ny+j+k*(nx*ny)])*(1.0+(C-1.0))*Km*h[i*ny+j+k*(nx*ny)]));
		}
	  }
	}
	
	return result;

}

