#include "cmdstuff.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>

using namespace std;

static const double hatocm=219474.63067;

EXTERN void FORTRAN(vinit)();
EXTERNC void gauleg(double x1,double x2,double *x,double *w,int n);
EXTERNC double plgndr(int j,int m,double x);
vector thetagrid(int nsize,vector &weights);
double PjNormal(int j,int m, double x);
double Pj0(int j,double x);


int main(int argc,char **argv) {
  int i,j,k,j1,j1p,j2,j2p,jp,n1,n2,n1p,n2p;

  int sizej=atoi(argv[1]);
  int jmax=2*(sizej-1);
  
  // ******************************
  // Set up basis for theta 
  // ******************************
  int size_theta=2.*jmax+10;

  // define grid points for theta 
  // initialize vectors of weights
  vector weights_theta(size_theta);
  vector grid_theta=thetagrid(size_theta,weights_theta);

  ofstream logout("log");
  // test orthonormality of theta and phi basis
  logout<<"  // test orthonormality of theta "<<endl;
  for (j=0;j<=jmax;j++){
	if (j%2) continue; // skip iteration if j odd
  	for (jp=0;jp<=jmax;jp++){
		if (jp%2) continue; 
		double sum=0.;
		for (i=0;i<size_theta;i++){
	  		double theta=grid_theta(i);
	  		sum+=weights_theta(i)*Pj0(j,cos(theta))*Pj0(jp,cos(theta));
		}
		logout<<j<<" "<<jp<<" "<<sum<<endl;
      	}
  }

// read potential data
// 10 angles
  double r,theta,pot;
  int size_r=134;
  int size_theta_pot=19;
  std::ifstream potfile("RbH2.dat");
  std::ofstream potfileout("pot.dat");
  matrix V2d(size_r,size_theta_pot);
  vector theta_pot_grid(size_theta_pot);
  vector r_pot_grid(size_r);
  for (i=0;i<size_r;i++) {
	for (j=0;j<size_theta_pot;j++) {
		potfile>>r>>theta>>pot;
		if (i==0) theta_pot_grid(j)=theta*M_PI/180.;
		if (j==0) r_pot_grid(i)=r;
		V2d(i,j)=pot*0.695;
	}
  }
  for (i=0;i<size_r;i++) {
	potfileout<<r_pot_grid(i)<<" ";
	for (j=0;j<size_theta_pot;j++) {
		potfileout<<V2d(i,j)<<" ";
	}
	potfileout<<endl;
  }

  // **************************************
  // Define kinetic energy info
  // **************************************

  // the rovibrational energy levels for the monomers
  vector E0j(3);
  vector E1j(3);

  E0j(0)=-36118.074;
  E0j(1)=-35763.701;
  E0j(2)=-34949.276;

  E1j(0)=-31956.927;
  E1j(1)=-31620.254;
  E1j(2)=-30846.711;

  double B0=(E0j(1)-E0j(0))/6.;
  double B1=(E1j(1)-E1j(0))/6.;

  std::ofstream potfinefileout("potfine.dat");
  for (i=0;i<size_r;i++) {
	vector pot_at_r(size_theta_pot);
	for (j=0;j<size_theta_pot;j++) pot_at_r(j)=V2d(i,j);
  	Interp V0_func(size_theta_pot,theta_pot_grid,pot_at_r);
   	for (k=0;k<size_theta;k++){
		double theta=grid_theta(k);
		double potvalue=0.;
		if (theta <= M_PI/2.)
			potvalue=V0_func.interp(theta);
		else
			potvalue=V0_func.interp((M_PI-theta));
		potfinefileout<<theta<<" "<<potvalue<<endl;;
	}
	potfinefileout<<endl;
  }
  // setup H
  for (i=0;i<size_r;i++) {
	vector pot_at_r(size_theta_pot);
	for (j=0;j<size_theta_pot;j++) pot_at_r(j)=V2d(i,j);
  	Interp V0_func(size_theta_pot,theta_pot_grid,pot_at_r);
  	matrix H(sizej,sizej);
// potential matrix elements
	for (j=0;j<sizej;j++) {
		int jvalue=2*j;
		for (jp=0;jp<sizej;jp++) {
			int jpvalue=2*jp;

			double sum=0.;
       			for (k=0;k<size_theta;k++){
	                        double theta=grid_theta(k);
				double potvalue=0.;
				if (theta <=M_PI/2.)
					potvalue=V0_func.interp(theta);
				else
					potvalue=V0_func.interp((M_PI-theta));
       		                sum+=weights_theta(k)*Pj0(jvalue,cos(theta))*Pj0(jpvalue,cos(theta))*potvalue;
                	}
			H(j,jp)=sum;
		}
		H(j,j)+=B0*(double)(jvalue*(jvalue+1));
	}
  	vector ev=diag(H);
	cout<<r_pot_grid(i)<<" "<<ev(0)<<endl;
  }
}
vector thetagrid(int nsize,vector &weights)
{
  int i;
  vector grid(nsize);
  double *x=new double[nsize];
  double *w=new double[nsize];
  double x1=-1.;
  double x2=1.;
  gauleg(x1,x2,x,w,nsize);
  for (i=0;i<nsize;i++) {
    grid(i)=acos(x[i]);
    weights(i)=w[i]; // in this case weights_theta=weights since the argument uses a reference operator
  }
  return grid;
}

double Pj0(int j,double x)
{
  return sqrt((double)j+.5)*plgndr(j,0,x);
}
