# include <cstdlib>
#include <iostream>
#include <string>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>


using namespace std;

# include "fd3d_diffusion.hpp"


//**************************************************************************************************************************
//******  Calculate diffussion of Temperature   
//**************************************************************************************************************************

double *fd3d_difussion_Temp_3D ( int x_num, double x[], int y_num, double y[], int z_num, double z[], double t, double dt, 
  double cflx, double cfly,  double cflz, double dens, double calesp, double conductivity, double *Temperature, double valueBond,
  void bc ( int x_num, int y_num, int z_num, double *V, double value  ) )

//////
//////

{
	double *T_new;
	int j, i, k;
	double mult;

	mult = conductivity/(dens*calesp);

	T_new = ( double * ) malloc ( x_num * y_num * (z_num) *sizeof ( double ) );
	
	k=0;
	#pragma acc kernels copyin(Temperature[0:x_num * y_num * z_num], h[0:x_num * y_num * z_num]) copyout(T_new[0:x_num * y_num * z_num])
	#pragma acc loop collapse(2) independent 
	for(i = 1; i < x_num - 1; i++ )
	{
		for ( j = 1; j < y_num - 1; j++ )
		{
			T_new[i*y_num+j+k*(x_num*y_num)] = Temperature[i*y_num+j] 
		  			+ mult * ( cfly * ( Temperature[i*y_num+j-1]
					- 2.0 * Temperature[i*y_num+j]
					+       Temperature[i*y_num+j+1] ) + cflx * (Temperature[(i-1)*y_num+j]
					- 2.0 * Temperature[i*y_num+j]
					+       Temperature[(i+1)*y_num+j] )+ cflz * (1.0*Temperature[i*y_num+j+(k+2)*(x_num*y_num)]
					- 2.0 * Temperature[i*y_num+j+(k)*(x_num*y_num)]
					+       Temperature[i*y_num+j+(k+1)*(x_num*y_num)] ));
		}
	}
	
	
		//cout << " K0 \n";
	
	#pragma acc loop collapse(3) independent 
	//#pragma acc kernels loop present(h_new, coeff_diff, h)
	for(k = 1; k < z_num -1; k++ )
	{
		for(i = 1; i < x_num - 1; i++ )
		{
			for ( j = 1; j < y_num - 1; j++ )
			{
				T_new[i*y_num+j+k*(x_num*y_num)] = Temperature[i*y_num+j+k*(x_num*y_num)]  
			  + mult * (  cfly* (Temperature[i*y_num+j-1+k*(x_num*y_num)]
						- 2.0 * Temperature[i*y_num+j+k*(x_num*y_num)]
						+       Temperature[i*y_num+j+1+k*(x_num*y_num)] ) + cflx * (Temperature[(i-1)*y_num+j+k*(x_num*y_num)]
						- 2.0 * Temperature[i*y_num+j+k*(x_num*y_num)]
						+       Temperature[(i+1)*y_num+j+k*(x_num*y_num)] ) + cflz * (Temperature[i*y_num+j+(k-1)*(x_num*y_num)]
						- 2.0 * Temperature[i*y_num+j+k*(x_num*y_num)]
						+       Temperature[i*y_num+j+(k+1)*(x_num*y_num)] ));	
			}
		}
	}
	
	k=z_num-1;
	
	#pragma acc loop collapse(2) independent 
	//#pragma acc kernels loop present(h_new, coeff_diff, h)
	for(i = 1; i < x_num - 1; i++ )
	{
		for ( j = 1; j < y_num - 1; j++ )
		{
			T_new[i*y_num+j+k*(x_num*y_num)] = Temperature[i*y_num+j+k*(x_num*y_num)] +
			  + mult * (cfly* (Temperature[i*y_num+j-1+k*(x_num*y_num)]
						- 2.0 * Temperature[i*y_num+j+k*(x_num*y_num)]+  Temperature[i*y_num+j+1+k*(x_num*y_num)] ) 
			  + cflx * ( Temperature[(i-1)*y_num+j+k*(x_num*y_num)]
						- 2.0 * Temperature[i*y_num+j+k*(x_num*y_num)]+  Temperature[(i+1)*y_num+j+k*(x_num*y_num)] ) 
			  + cflz * ( Temperature[i*y_num+(k-1)*(x_num*y_num)]
						- 2.0 * Temperature[i*y_num+j+(k)*(x_num*y_num)]	+   1.0*  Temperature[i*y_num+j+(k-2)*(x_num*y_num)] ));
		}
	}
			
	bc ( x_num, y_num, z_num, T_new, valueBond );

  return T_new;
}



//**************************************************************************************************************************
//******  Calculate diffussion of Humidity   
//**************************************************************************************************************************

double *fd3d_difussion_Humidity_3D ( int x_num, double x[], int y_num, double y[], int z_num, double z[], double t, double dt, 
  double cflx, double cfly,  double cflz, double * Coeff_Diff, double *Humidity, double valueBond,
  void bc ( int x_num, int y_num, int z_num, double *V, double value  ) )

//////
//////

{
	double *H_new;
	int j, i, k;
	
	H_new = ( double * ) malloc ( x_num * y_num * (z_num) *sizeof ( double ) );
	
	k=0;
	#pragma acc kernels copyin(Humidity[0:x_num * y_num * z_num], h[0:x_num * y_num * z_num]) copyout(H_new[0:x_num * y_num * z_num])
	#pragma acc loop collapse(2) independent 
	for(i = 1; i < x_num - 1; i++ )
	{
		for ( j = 1; j < y_num - 1; j++ )
		{
			H_new[i*y_num+j+k*(x_num*y_num)] = Humidity[i*y_num+j] 
		  			+ Coeff_Diff[i*y_num+j] * ( cfly * ( Humidity[i*y_num+j-1]
					- 2.0 * Humidity[i*y_num+j]
					+       Humidity[i*y_num+j+1] ) + cflx * (Humidity[(i-1)*y_num+j]
					- 2.0 * Humidity[i*y_num+j]
					+       Humidity[(i+1)*y_num+j] )+ cflz * (1.0*Humidity[i*y_num+j+(k+2)*(x_num*y_num)]
					- 2.0 * Humidity[i*y_num+j+(k)*(x_num*y_num)]
					+       Humidity[i*y_num+j+(k+1)*(x_num*y_num)] ));
		}
	}
	
	
		//cout << " K0 \n";
	
	#pragma acc loop collapse(3) independent 
	//#pragma acc kernels loop present(h_new, coeff_diff, h)
	for(k = 1; k < z_num -1; k++ )
	{
		for(i = 1; i < x_num - 1; i++ )
		{
			for ( j = 1; j < y_num - 1; j++ )
			{
				H_new[i*y_num+j+k*(x_num*y_num)] = Humidity[i*y_num+j+k*(x_num*y_num)]  
			  + mult * (  cfly* (Humidity[i*y_num+j-1+k*(x_num*y_num)]
						- 2.0 * Humidity[i*y_num+j+k*(x_num*y_num)]
						+       Humidity[i*y_num+j+1+k*(x_num*y_num)] ) + cflx * (Humidity[(i-1)*y_num+j+k*(x_num*y_num)]
						- 2.0 * Humidity[i*y_num+j+k*(x_num*y_num)]
						+       Humidity[(i+1)*y_num+j+k*(x_num*y_num)] ) + cflz * (Humidity[i*y_num+j+(k-1)*(x_num*y_num)]
						- 2.0 * Humidity[i*y_num+j+k*(x_num*y_num)]
						+       Humidity[i*y_num+j+(k+1)*(x_num*y_num)] ));	
			}
		}
	}
	
	k=z_num-1;
	
	#pragma acc loop collapse(2) independent 
	//#pragma acc kernels loop present(h_new, coeff_diff, h)
	for(i = 1; i < x_num - 1; i++ )
	{
		for ( j = 1; j < y_num - 1; j++ )
		{
			H_new[i*y_num+j+k*(x_num*y_num)] = Humidity[i*y_num+j+k*(x_num*y_num)] +
			  + mult * (cfly* (Humidity[i*y_num+j-1+k*(x_num*y_num)]
						- 2.0 * Humidity[i*y_num+j+k*(x_num*y_num)]+  Humidity[i*y_num+j+1+k*(x_num*y_num)] ) 
			  + cflx * ( Humidity[(i-1)*y_num+j+k*(x_num*y_num)]
						- 2.0 * Humidity[i*y_num+j+k*(x_num*y_num)]+  Humidity[(i+1)*y_num+j+k*(x_num*y_num)] ) 
			  + cflz * ( Humidity[i*y_num+(k-1)*(x_num*y_num)]
						- 2.0 * Humidity[i*y_num+j+(k)*(x_num*y_num)]	+   1.0*  Humidity[i*y_num+j+(k-2)*(x_num*y_num)] ));
		}
	}
			
	bc ( x_num, y_num, z_num, H_new, valueBond );

  return T_new;
}



double *fd3d_diffusion_3D ( int x_num, double x[], int y_num, double y[], int z_num, double z[], double t, double dt, 
  double cflx, double cfly,  double cflz, double *porosity, double *humidity, double value_Bond, double *coeff ( double *p , double *h, int nx, int ny, int nz ), 
  void bc ( int x_num, int y_num, int z_num, double *c, double value  ), double *h )

//****************************************************************************80
//
//  Purpose:
//
//    FD3D_DIFFUSION: Finite difference solution of 3D Diffusion equation.
//
//  Discussion:
//
//    This program takes one time step to solve the 3D Diffusion equation 
//    with an explicit method.
//
//    This program solves
//
//      dUdT - k * d2UdX2 = F(X,T)
//
//    over the interval [A,B] with boundary conditions
//
//      U(A,T) = UA(T),
//      U(B,T) = UB(T),
//
//    over the time interval [T0,T1] with initial conditions
//
//      U(X,T0) = U0(X)
//
//    The code uses the finite difference method to approximate the
//    second derivative in space, and an explicit forward Euler approximation
//    to the first derivative in time.
//
//    The finite difference form can be written as
//
//      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
//      ------------------  = F(X,T) + k *  ------------------------------------
//               dt                                   dx * dx
//
//    or, assuming we have solved for all values of U at time T, we have
//
//      U(X,T+dt) = U(X,T) 
//        + cfl * ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) ) + dt * F(X,T) 
//
//    Here "cfl" is the Courant-Friedrichs-Loewy coefficient:
//
//      cfl = k * dt / dx / dx
//
//    In order for accurate results to be computed by this explicit method,
//    the CFL coefficient must be less than 0.5.
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
//    Input, int X_NUM, the number of points to use in the 
//    spatial dimension.
//
//    Input, double X(X_NUM), the coordinates of the nodes.
//
//    Input, double T, the current time.
//
//    Input, double DT, the size of the time step.
//
//    Input, double CFL, the Courant-Friedrichs-Loewy coefficient,
//    computed by FD1D_HEAT_EXPLICIT_CFL.
//
//    Input, double H[X_NUM], the solution at the current time.
//
//    Input, double *coeff ( int x_num, double x[], double t ), the function 
//    which evaluates the right hand side.
//
//    Input, void BC ( int x_num, double x[], double t, double h[] ), 
//    the function which evaluates the boundary conditions.
//
//    Output, double FD1D_HEAT_EXPLICIT[X_NUM)], the solution at time T+DT.
//
{
	double *coeff_diff;
	double *h_new;
	int j, i, k;

	coeff_diff = coeff ( porosity, humidity, x_num, y_num, z_num);
	cout <<  x_num << " - " << y_num << " - " << z_num << " \n";
	
	h_new = ( double * ) malloc ( x_num * y_num * (z_num) *sizeof ( double ) );
	cout << " New MAtriz \n";
	
	//#pragma acc kernels copyin(coeff_diff[0:x_num * y_num * z_num], h[0:x_num * y_num * z_num]) copyout(h_new[0:x_num * y_num * z_num])
	
	k=0;
	//#pragma acc loop collapse(2) independent 
	//#pragma acc data copyin(coeff_diff[0:x_num * y_num * z_num], h[0:x_num * y_num * z_num]) copyout(h_new[0:x_num * y_num * z_num])
	
	//#pragma acc kernels loop present(h_new, coeff_diff, h)
	#pragma acc kernels copyin(coeff_diff[0:x_num * y_num * z_num], h[0:x_num * y_num * z_num]) copyout(h_new[0:x_num * y_num * z_num])
	
	#pragma acc loop collapse(2) independent 
	for(i = 1; i < x_num - 1; i++ )
		{
		  for ( j = 1; j < y_num - 1; j++ )
		  {
			h_new[i*y_num+j+k*(x_num*y_num)] = h[i*y_num+j] 
			  + coeff_diff[i*y_num+j+k*(x_num*y_num)] * ( cfly * ( h[i*y_num+j-1]
						- 2.0 * h[i*y_num+j]
						+       h[i*y_num+j+1] ) + cflx * ( h[(i-1)*y_num+j]
						- 2.0 * h[i*y_num+j]
						+       h[(i+1)*y_num+j] )+ cflz * (1.0*h[i*y_num+j+(k+2)*(x_num*y_num)]
						- 2.0 * h[i*y_num+j+(k)*(x_num*y_num)]
						+       h[i*y_num+j+(k+1)*(x_num*y_num)] ));
		  }
		}
	
	
		//cout << " K0 \n";
	
	#pragma acc loop collapse(3) independent 
	//#pragma acc kernels loop present(h_new, coeff_diff, h)
	for(k = 1; k < z_num -1; k++ )
	{
		for(i = 1; i < x_num - 1; i++ )
		{
		  for ( j = 1; j < y_num - 1; j++ )
		  {
			h_new[i*y_num+j+k*(x_num*y_num)] = h[i*y_num+j+k*(x_num*y_num)]  
			  + coeff_diff[i*y_num+j+k*(x_num*y_num)] * (  cfly* (   h[i*y_num+j-1+k*(x_num*y_num)]
						- 2.0 * h[i*y_num+j+k*(x_num*y_num)]
						+       h[i*y_num+j+1+k*(x_num*y_num)] ) + cflx * (h[(i-1)*y_num+j+k*(x_num*y_num)]
						- 2.0 * h[i*y_num+j+k*(x_num*y_num)]
						+       h[(i+1)*y_num+j+k*(x_num*y_num)] ) + cflz * (h[i*y_num+j+(k-1)*(x_num*y_num)]
						- 2.0 * h[i*y_num+j+k*(x_num*y_num)]
						+       h[i*y_num+j+(k+1)*(x_num*y_num)] ));
			//cout << i*x_num+j+(k+1)*(x_num*y_num) << " \n";
			//cout << "C[" <<  i << "," << j << "," << k << "]= " << h[i*x_num+j+(k+1)*(x_num*y_num)]<< "__"  ;
		  }
		}
		//cout << " K = " << k << " \n";
	}
	
	k=z_num-1;
	//cout << " \n K = " << k << " \n";
	
	#pragma acc loop collapse(2) independent 
	//#pragma acc kernels loop present(h_new, coeff_diff, h)
	for(i = 1; i < x_num - 1; i++ )
		{
		  for ( j = 1; j < y_num - 1; j++ )
		  {
			h_new[i*y_num+j+k*(x_num*y_num)] = h[i*y_num+j+k*(x_num*y_num)] +
			  + coeff_diff[i*y_num+j+k*(x_num*y_num)] * (cfly* (h[i*y_num+j-1+k*(x_num*y_num)]
						- 2.0 * h[i*y_num+j+k*(x_num*y_num)]+  h[i*y_num+j+1+k*(x_num*y_num)] ) 
			  + cflx * ( h[(i-1)*y_num+j+k*(x_num*y_num)]
						- 2.0 * h[i*y_num+j+k*(x_num*y_num)]+  h[(i+1)*y_num+j+k*(x_num*y_num)] ) 
			  + cflz * ( h[i*y_num+(k-1)*(x_num*y_num)]
						- 2.0 * h[i*y_num+j+(k)*(x_num*y_num)]	+   1.0*  h[i*y_num+j+(k-2)*(x_num*y_num)] ));
		  }
		}
			
	bc ( x_num, y_num, z_num, h_new, value_Bond );

  delete [] coeff_diff;
  return h_new;
}


double *fd3d_diffusion_pri ( int x_num, double x[], int y_num, double y[], int z_num, double z[], double t, double dt, 
  double cflx, double cfly,  double cflz, double *porosity, double *humidity, double value_Bond, double *coeff ( double *p , double *h, int nx, int ny, int nz ), 
  void bc ( int x_num, int y_num, int z_num, double *c, double value  ), double *h )

//****************************************************************************80
//
//  Purpose:
//
//    FD3D_DIFFUSION: Finite difference solution of 3D Diffusion equation.
//
//  Discussion:
//
//    This program takes one time step to solve the 3D Diffusion equation 
//    with an explicit method.
//
//    This program solves
//
//      dUdT - k * d2UdX2 = F(X,T)
//
//    over the interval [A,B] with boundary conditions
//
//      U(A,T) = UA(T),
//      U(B,T) = UB(T),
//
//    over the time interval [T0,T1] with initial conditions
//
//      U(X,T0) = U0(X)
//
//    The code uses the finite difference method to approximate the
//    second derivative in space, and an explicit forward Euler approximation
//    to the first derivative in time.
//
//    The finite difference form can be written as
//
//      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
//      ------------------  = F(X,T) + k *  ------------------------------------
//               dt                                   dx * dx
//
//    or, assuming we have solved for all values of U at time T, we have
//
//      U(X,T+dt) = U(X,T) 
//        + cfl * ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) ) + dt * F(X,T) 
//
//    Here "cfl" is the Courant-Friedrichs-Loewy coefficient:
//
//      cfl = k * dt / dx / dx
//
//    In order for accurate results to be computed by this explicit method,
//    the CFL coefficient must be less than 0.5.
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
//    Input, int X_NUM, the number of points to use in the 
//    spatial dimension.
//
//    Input, double X(X_NUM), the coordinates of the nodes.
//
//    Input, double T, the current time.
//
//    Input, double DT, the size of the time step.
//
//    Input, double CFL, the Courant-Friedrichs-Loewy coefficient,
//    computed by FD1D_HEAT_EXPLICIT_CFL.
//
//    Input, double H[X_NUM], the solution at the current time.
//
//    Input, double *coeff ( int x_num, double x[], double t ), the function 
//    which evaluates the right hand side.
//
//    Input, void BC ( int x_num, double x[], double t, double h[] ), 
//    the function which evaluates the boundary conditions.
//
//    Output, double FD1D_HEAT_EXPLICIT[X_NUM)], the solution at time T+DT.
//
{
	double *coeff_diff;
	double *h_new;
	int j, i, k;

	coeff_diff = coeff ( porosity, humidity, x_num, y_num, z_num);
	cout <<  x_num << " - " << y_num << " - " << z_num << " \n";
	
	h_new = ( double * ) malloc ( x_num * y_num * (z_num) *sizeof ( double ) );
	cout << " New MAtriz \n";
	
	//#pragma acc kernels copyin(coeff_diff[0:x_num * y_num * z_num], h[0:x_num * y_num * z_num]) copyout(h_new[0:x_num * y_num * z_num])
	
	k=0;
	//#pragma acc loop collapse(2) independent 
	//#pragma acc data copyin(coeff_diff[0:x_num * y_num * z_num], h[0:x_num * y_num * z_num]) copyout(h_new[0:x_num * y_num * z_num])
	
	//#pragma acc kernels loop present(h_new, coeff_diff, h)
	#pragma acc kernels copyin(coeff_diff[0:x_num * y_num * z_num], h[0:x_num * y_num * z_num]) copyout(h_new[0:x_num * y_num * z_num])
	
	#pragma acc loop collapse(2) independent 
	for(i = 1; i < x_num - 1; i++ )
		{
		  for ( j = 1; j < y_num - 1; j++ )
		  {
			h_new[i*y_num+j+k*(x_num*y_num)] = h[i*y_num+j] + coeff_diff[i*y_num+j+k*(x_num*y_num)] *
			   (cfly* ( h[i*y_num+j+1]-h[i*y_num+j-1] ) + cflx * ( h[(i+1)*y_num+j]
						-  h[(i-1)*y_num+j] )+ cflz* (h[i*y_num+j+(k+2)*(x_num*y_num)]
						- h[i*y_num+j+(k)*(x_num*y_num)]));
		  }
		}
	
	
		//cout << " K0 \n";
	
	#pragma acc loop collapse(3) independent 
	//#pragma acc kernels loop present(h_new, coeff_diff, h)
	for(k = 1; k < z_num -1; k++ )
	{
		for(i = 1; i < x_num - 1; i++ )
		{
		  for ( j = 1; j < y_num - 1; j++ )
		  {
			h_new[i*y_num+j+k*(x_num*y_num)] = h[i*y_num+j+k*(x_num*y_num)]  
			  + coeff_diff[i*y_num+j+k*(x_num*y_num)] * (  cfly *( -h[i*y_num+j-1+k*(x_num*y_num)]
						+       h[i*y_num+j+1+k*(x_num*y_num)] ) + cflx* (-h[(i-1)*y_num+j+k*(x_num*y_num)]
						+       h[(i+1)*y_num+j+k*(x_num*y_num)] ) + cflz * (-h[i*y_num+j+(k-1)*(x_num*y_num)]
						+       h[i*y_num+j+(k+1)*(x_num*y_num)] ));
			//cout << i*x_num+j+(k+1)*(x_num*y_num) << " \n";
			//cout << "C[" <<  i << "," << j << "," << k << "]= " << h[i*x_num+j+(k+1)*(x_num*y_num)]<< "__"  ;
		  }
		}
		//cout << " K = " << k << " \n";
	}
	
	k=z_num-1;
	//cout << " \n K = " << k << " \n";
	
	#pragma acc loop collapse(2) independent 
	//#pragma acc kernels loop present(h_new, coeff_diff, h)
	for(i = 1; i < x_num - 1; i++ )
		{
		  for ( j = 1; j < y_num - 1; j++ )
		  {
			h_new[i*y_num+j+k*(x_num*y_num)] = h[i*y_num+j+k*(x_num*y_num)] +
			  + coeff_diff[i*y_num+j+k*(x_num*y_num)] * (cfly*( - h[i*y_num+j-1+k*(x_num*y_num)]
						+  h[i*y_num+j+1+k*(x_num*y_num)] ) 
			  + cflx* ( -h[(i-1)*y_num+j+k*(x_num*y_num)]
						+  h[(i+1)*y_num+j+k*(x_num*y_num)] ) 
			  + cflz * ( -h[i*y_num+(k-1)*(x_num*y_num)]
						+   h[i*y_num+j+(k)*(x_num*y_num)] ));  // (k-2)
		  }
		}
			
	bc ( x_num, y_num, z_num, h_new, value_Bond );

  delete [] coeff_diff;
  return h_new;
}




//****************************************************************************80

double *ic_3D ( int x_num, int y_num, int z_num, double value )

//****************************************************************************80
//
//  Purpose:
//
//    IC_TEST01 evaluates the initial condition for problem 1.
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
//    Input, int X_NUM, the number of nodes.
//
//    Input, double X[X_NUM], the node coordinates.
//
//    Input, double T, the initial time.
//
//    Output, double C[X_NUM, Y_NUM, Z_NUM], the heat values at the initial time.
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
				c[i*y_num+j+k*(x_num*y_num)] = value;
				//cout << "C[ " << i << "," << j << "," << k << "]= " << c[i*x_num+j+k*(x_num*y_num)] << "\n"  ;
				//getchar();
			}
		}
	}
  return c;
}



//****************************************************************************80

double fd3d_diffusion_cfl ( double k, int t_num, double t_min, double t_max, 
  int x_num, double x_min, double x_max, int y_num, double y_min, double y_max, int z_num, double z_min, double z_max )

//****************************************************************************80
//
//  Purpose:
//
//    FD2D_DIFFUSION_CFL: compute the Courant-Friedrichs-Loewy coefficient.
//
//  Discussion:
//
//    The equation to be solved has the form:
//
//      dUdT - k * d2UdX2 = F(X,T)
//
//    over the interval [X_MIN,X_MAX][Y_MIN,Y_MAX] with boundary conditions
//
//      U(X_MIN,T) = U_X_MIN(T),
//      U(X_MIN,T) = U_X_MAX(T),
//
//    over the time interval [T_MIN,T_MAX] with initial conditions
//
//      U(X,T_MIN) = U_T_MIN(X)
//
//    The code uses the finite difference method to approximate the
//    second derivative in space, and an explicit forward Euler approximation
//    to the first derivative in time.
//
//    The finite difference form can be written as
//
//      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
//      ------------------  = F(X,T) + k *  ------------------------------------
//               dt                                   dx * dx
//
//    or, assuming we have solved for all values of U at time T, we have
//
//      U(X,T+dt) = U(X,T) 
//        + cfl * ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) ) + dt * F(X,T) 
//
//    Here "cfl" is the Courant-Friedrichs-Loewy coefficient:
//
//      cfl = k * dt / dx / dx
//
//    In order for accurate results to be computed by this explicit method,
//    the CFL coefficient must be less than 0.5!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 December 2017
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Reference:
//
//    George Lindfield, John Penny,
//    Numerical Methods Using MATLAB,
//    Second Edition,
//    Prentice Hall, 1999,
//    ISBN: 0-13-012641-1,
//    LC: QA297.P45.
//
//  Parameters:
//
//    Input, double K, the heat conductivity coefficient.
//
//    Input, int T_NUM, the number of time values, including 
//    the initial value.
//
//    Input, double T_MIN, T_MAX, the minimum and maximum times.
//
//    Input, int X_NUM, the number of nodes.
//
//    Input, double X_MIN, X_MAX, the minimum and maximum spatial 
//    coordinates.
//
//    Output, double FD1D_HEAT_EXPLICIT_CFL, the Courant-Friedrichs-Loewy coefficient.
//
{
  double cfl;
  double dx, dy, dz;
  double dt;

  dx = ( x_max - x_min ) / ( double ) ( x_num - 1 );
  dy = ( y_max - y_min ) / ( double ) ( y_num - 1 );
  dz = ( z_max - z_min ) / ( double ) ( z_num - 1 );
  dt = ( t_max - t_min ) / ( double ) ( t_num - 1 );
//
//  Check the CFL condition, print out its value, and quit if it is too large.
//
  cfl = k * dt / dx / dx + k * dt / dy / dy + k * dt / dz / dz;

  cout << "\n";
  cout << " CFL stability criterion value = " << cfl << "\n";

  if ( 0.5 <= cfl )
  {
    cerr << "\n";
    cerr << "FD3D_DIFFUSION_CFL - Fatal error!\n";
    cerr << "  CFL condition failed.\n";
    cerr << "  0.5 <= K * dT / dX / dX + K * dT / dY / dY  = CFL.\n";
    exit ( 1 );
  }

  return cfl;
}

double get_dtdx2 ( int t_num, double t_min, double t_max, 
  int x_num, double x_min, double x_max)
{
  double cfl;
  double dx;
  double dt;

  dx = ( x_max - x_min ) / ( double ) ( x_num - 1 );
  dt = ( t_max - t_min ) / ( double ) ( t_num - 1 );
  cfl =  dt / dx/dx;
  return cfl;
}


posCoord *get_Coordenadas (int numBeam, posCoord *posBeam, int x_num, double x_min, double x_max, int y_num, double y_min, double y_max)
{
  double cfl;
  double dx;
  double dy;
	posCoord *posVerify;
		
	
  dx = ( x_max - x_min ) / ( double ) ( x_num - 1 );
  dy = ( y_max - y_min ) / ( double ) ( y_num - 1 );
  
  posVerify = ( posCoord * ) malloc ( numBeam * sizeof ( posCoord ) );
  for(int i;i<numBeam; i++)
  {
  		posverify[i].posX= posBeam[i].posX/dx;	
  		posverify[i].posY= posBeam[i].posY/dy;	
  		cout << "posX=" << posverify[i].posX;
  		cout << "  posY=" << posverify[i].posY << "\n";
  }
  return posVerify;
}

//****************************************************************************80

double *r8vec_linspace_new ( int n, double a_first, double a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 November 2017
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
{
  double *a;
  int i;

  a = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return a;
}

//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 December 2017
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}




//****************************************************************************80

double *porosity_co2( double *por_ini , double *por_current, double *por_last, double *co2, int n, int m)

//****************************************************************************80
//
//  Purpose:
//
//    POROSITY calculate the current porosity of concrete.
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
//    Input, double *por_ini, the porosity initial.
//
//    Input, double h, the humidity.
//
//    Output, double coefficient_diffusion_co2, Coefficient Diffusion.
//
{
	double *por;
	int i,j;
	
	por = ( double * ) malloc ( n * m * sizeof ( double ) );
	
	for(i = 0; i < n ; i++ )
	{
	  for ( j = 0; j < m ; j++ )
	  {
		por[i*m+j] = por_ini[i*m+j] - (por_current[i*m+j]-por_last[i*m+j])*co2[i*m+j];  
	  }
	}
	return por;
}


//********************************************************************************************

void write_vtk3D ( string output_filename, int n, double x[], int m, double y[], int h, double z[], double *c, string variable )

//*******************************************************************************************
//
//  Purpose:
//
//    WRITE_VTK writes an VTK Unstructured file.
//
//  Modified:
//
//    26 October 2017
//
//  Author:
//
//    Félix Armando Mejía
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int N, the number of points in X.
//
//    Input, double X[N], pointS in X.
//
//    Input, int M, the number of points in Y.
//
//    Input, double Y[M], pointS in Y.
//    
//    Input, double *c, Matrix of Chloride Concentration .
//
{
  int j, i, k;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "WRITE_VTK - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//

	output << "# vtk DataFile Version 2.0\n";
	output << "Simula \n";
	output << "ASCII\n";
	output << "DATASET UNSTRUCTURED_GRID\n";
	output << "POINTS "<< n*m*h << " double\n" ;
	

  for ( k = 0; k < h; k++ )
  {
	  for ( i = 0; i < n; i++ )
	  {
		  for ( j = 0; j < m; j++ )
		  {
			/*SwapEnd(x[i]);
			output.write((char*)&x[i], sizeof(double));
			SwapEnd(y[j]);
			output.write((char*)&y[j], sizeof(double));
			SwapEnd(z[k]);
			output.write((char*)&z[k], sizeof(double));
			output << x[i] <<  " " << y[j] <<  " " << z[k]  <<" \n";*/
			output << x[i] <<  " " << y[j] <<  " " << z[k]  <<" \n";

		  }
	  }
  }
  output << "CELLS "<< (n-1)*(m-1)*(h-1) << "   " << (n-1)*(m-1)*(h-1)*9  << " \n";
  int numPoints = 8; 
  for ( k = 0; k < h-1; k++ )
  {
	  for ( i = 0; i < n-1; i++ )
	  {
		  for ( j = 0; j < m -1 ; j++ )
		  {
			//output <<  numPoints <<  " " << i*m+j+k*(n*m) <<  " " << i*m+j+1+k*(n*m) << " " << (i+1)*m+j+k*(n*m) << " " << (i+1)*m+j+1+k*(n*m);
			//output <<  " " << i*m+j+(k+1)*(n*m) <<  " " << i*m+j+1+(k+1)*(n*m) << " " << (i+1)*m+j+(k+1)*(n*m) << " " << (i+1)*m+j+1+(k+1)*(n*m) << "\n";

			output <<  numPoints <<  " " << i*m+j+k*(n*m) <<  " " << (i+1)*m+j+k*(n*m)  << " " <<  (i+1)*m+j+1+k*(n*m) << " " << i*m+j+1+k*(n*m);
			output <<  " " << i*m+j+(k+1)*(n*m) <<  " " <<  (i+1)*m+j+(k+1)*(n*m)  << " " << (i+1)*m+j+1+(k+1)*(n*m) << " " << i*m+j+1+(k+1)*(n*m) << "\n";
		  }
	  }
  }
  
  output << "CELL_TYPES "<< (n-1)*(m-1)*(h-1) << " \n";
  int cellType = 12; 
  for ( i = 0; i < (n-1)*(m-1)*(h-1); i++ )
  {
	  output <<  cellType << " \n";
  }
  output << "POINT_DATA "<< n*m*h << "\n" ;
  output << "SCALARS " << variable << " double" << " \n";
  output << "LOOKUP_TABLE default\n";
  for ( k = 0; k < h; k++ )
  {
	  for ( i = 0; i < n; i++ )
	  {
		  for ( j = 0; j < m; j++ )
		  {
			output << c[i*m+j+k*(n*m)] << "\n";
		  }
	  }
  }
//
//  Close the file.
//
  output.close ( );

  return;
}



//****************************************************************************80

double *porosity_cl( double *por_ini , double *por_current, double *por_last, double *co2, int n, int m)

//****************************************************************************80
//
//  Purpose:
//
//    POROSITY calculate the current porosity of concrete.
//
//
//  Modified:
//
//    26 November 2017
//
//  Author:
//
//    Félix Armando Mejía Cajicá
//
//  Parameters:
//
//    Input, double *por_ini, the porosity initial.
//
//    Input, double h, the humidity.
//
//    Output, double coefficient_diffusion_co2, Coefficient Diffusion.
//
{
	double *por;
	int i,j;
	
	por = ( double * ) malloc ( n * m * sizeof ( double ) );
	
	for(i = 0; i < n ; i++ )
	{
	  for ( j = 0; j < m ; j++ )
	  {
		por[i*m+j] = por_ini[i*m+j] - (por_current[i*m+j]-por_last[i*m+j])*co2[i*m+j];  
	  }
	}
	return por;
}



//****************************************************************************80
double calculate_rpo(double a_c)
{
	double result;
	result = -4.66*pow(a_c,2.0)+8.72*a_c-1.78;
	
	return result;
	
}


//****************************************************************************80
double calculate_rpc(double a_c)
{
	double result;
	result = calculate_rpo(a_c)*(10.59*pow(a_c,2.0)-11.36*a_c+4.29);
	
	return result;
		
}


//****************************************************************************80
//   Porosity of the non-carbonated concrete
//****************************************************************************80
double Phi0(double Pa, double Pc, double Pw, double a_c, double agr_c) 
{
	double result;
	double Pc_Pa = Pc/Pw; 
	double Pc_Pagr = Pc/Pa;
	result = (a_c*Pc_Pa)/(a_c*Pc_Pa+agr_c*Pc_Pagr+1.0);
	
	return result;
		
}


//****************************************************************************80
//   Porosity in the time carbonated concrete
//****************************************************************************80
// t in days
//****************************************************************************80

double Phit(double Phi0, double t ) 
{
	double result;
	double nmin = 17.8  ;  // mol/m3dia;  1.78*pow(10.0,-5)  mol/cm3dia;
	double beta = 4.19 * pow(10.0, -6.0);   // m3/mol
	double exponente = -beta * nmin  * t;
	result = Phi0*exp(exponente);
	
	return result;
		
}


//****************************************************************************80
//   Ultimate Reduction of the porosity when the concrete is fully
//   carbonated.
//****************************************************************************80
double Ccsh0(double Pcsh, double Ccao0)
{
	double result;
	result = Pcsh*Ccao0/3.0;
	return result;
}


//****************************************************************************80
//   Ultimate Reduction of the porosity when the concrete is fully
//   carbonated.
//****************************************************************************80
double deltaPhic(double Pcsh, double Ccao0) 
{
	double result;
	double deltaVch=3.85*pow(10.0, -6.0); // m3/mol
	double deltaVcsh=15.39*pow(10.0, -6.0); // m3/mol
	double vCcsh0=Ccsh0(Pcsh, Ccao0);
	result = (Ccao0-3.0*vCcsh0)*deltaVch+vCcsh0*deltaVcsh;
	
	return result;
		
}



//****************************************************************************80
//   Grade of carbonation 
//****************************************************************************80

double *Alpha_c(double Phi0, double deltaPhic, double * porosity_time, int nx, int ny, int nz) 
{
	double *result;
	result = ( double * ) malloc ( nx * ny * nz * sizeof ( double ) );
	

	for(i = 0; i < nx ; i++ )
	{
	  for ( j = 0; j < ny ; j++ )
	  {
		for ( k = 0; k < nz ; j++ )
		{
			result[i*ny+j+k*(nx*ny)]=(Phi0 - porosity_time[i*ny+j+k*(nx*ny)])/(deltaPhic);
		}
	  }
	}
	
	return result;

		
}





