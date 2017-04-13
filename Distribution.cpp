#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <memory>

#include "Distribution.h"
#include "Point.h"
#include "Random.h"

double uniform_distribution::sample() { return a + Urand() * ( b - a ); }

double linear_distribution::sample() {
  double r1 = Urand(), r2 = Urand();
  double p  = 2.0 * std::fmin( fa, fb ) / ( fa + fb );
  if ( r1 < p ) { return a + r2 * ( b - a ); }
  else {
    if ( fb > fa ) { return a + ( b - a ) * std::sqrt( r2 ); }
    else           { return a + ( b - a ) * ( 1.0 - std::sqrt( r2 )); }
  }
}

double exponential_distribution::sample() { return -std::log( Urand() ) / lambda; }

double normal_distribution::sample() 
{ return mu + sigma * std::sqrt( -2.0 * std::log( Urand() ) ) * std::cos( twopi * Urand() ); }

double HenyeyGreenstein_distribution::sample() {
  if ( a != 0.0 ) {
    return ( 1.0 + a*a - pow( ( 1.0 - a*a )/( 1.0 + a*(2.0*Urand() - 1.0) ), 2 ) ) / ( 2.0 * a ); 
  }
  else {
    return 2.0 * Urand() - 1.0;
  }
}

int meanMultiplicity_distribution::sample() {
  return (int) std::floor( nu + Urand() );
}

TerrellFission_distribution::TerrellFission_distribution( std::string label, double p1, double p2, double p3 ) 
    : distribution(label), nubar(p1), sigma(p2), b(p3) {
  double c  = 0.0;
  int nu = 0;
  while ( c < 1.0 - 1.0e-12 ) {
    double a  = ( nu - nubar + 0.5 + b ) / sigma;
    c = 0.5 * ( 1 + erf( a / sqrt(2.0) ) ) ;

    cdf.push_back(c);
    nu += 1;
  }
  cdf.push_back(1.0);
}

int TerrellFission_distribution::sample() {
  double r  = Urand();
  int nu;
  for ( int i = 0 ; i < cdf.size() ; i++ ) {
    if ( r < cdf[i] ) {
      nu = i;
      break;
    }
  }
 return nu; 
}

point isotropicDirection_distribution::sample() {
  // sample polar cosine and azimuthal angle uniformly
  double mu  = 2.0 * Urand() - 1.0;
  double azi = twopi * Urand();

  // convert to Cartesian coordinates
  double c = std::sqrt( 1.0 - mu * mu );
  point p;
  p.x = std::cos( azi ) * c;
  p.y = std::sin( azi ) * c;
  p.z = mu;

  return p;
}

point anisotropicDirection_distribution::sample() {
  double mu  = dist_mu->sample(); 
  double azi = twopi * Urand();
  double cos_azi = std::cos(azi);
  double sin_azi = std::sin(azi);

  // rotate the local particle coordinate system aligned along the incident direction
  // to the global problem (x,y,z) coordinate system 
  double sin_t0 = std::sqrt( 1.0 - mu * mu );
  double c = sin_t0 / sin_t;

  point p;
  p.x = axis.x * mu + ( axis.x * axis.z * cos_azi - axis.y * sin_azi ) * c;
  p.y = axis.y * mu + ( axis.y * axis.z * cos_azi + axis.x * sin_azi ) * c;
  p.z = axis.z * mu - cos_azi * sin_t0 * sin_t;
  return p;
}

point independentXYZ_distribution::sample() {
  return point( dist_x->sample(), dist_y->sample(), dist_z->sample() );
}

point uniform_spherical_dist::sample() {
	double r0_cb = r0 * r0 * r0;
	double R_cb = R * R * R;
	double diff = R_cb - r0_cb;

	double r = std::cbrt( r0_cb + ( diff * Urand() ) );
	double eta = ( 2 * Urand() ) - 1;
	double phi = twopi * Urand();

	double eta_sin = std::sqrt( 1 - ( eta * eta ) );

	return point( center.x + (r * eta_sin * std::cos(phi) ), center.y + (r * eta_sin * std::sin(phi) ), center.z + ( r * eta ) );
}

point uniform_disk_dist::sample() {
	double r0_sq = r0 * r0;
	double R_sq = R * R;
	double diff = R_sq - r0_sq;

	double r = std::sqrt( r0_sq + ( diff * Urand() ) );
	double phi = twopi * Urand();

	if( axis == "x" )
		return point( center.x, center.y + ( r * std::cos(phi) ), center.z + ( r * std::sin(phi) ) );
	else if( axis == "y" )
		return point( center.x + ( r * std::cos(phi) ), center.y, center.z + ( r * std::sin(phi) ) );
	else if( axis == "z" )
		return point( center.x + ( r * std::cos(phi) ), center.y + ( r * std::sin(phi) ), center.z );
	else
		return point( std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max() );
}
double angulardirection_distribution::sample() {
	
	// sample polar cosine angle 
	double u = 2.0*Urand()-1.0 ;
	//sample from 0 to 1 
	if (u > 0){
		 double mu = Urand();
	     if (Urand() < 0.79){
		 return mu;
	 }
	 else {
		   return std::pow(6.0*Urand(),1.0/3.0)-1.0 * copysign(1.0,mu);
	    }
	}
	//sample from 0 to -1 
	 if (u < 0){
		 double mu = Urand()-1.0;
		 if (Urand() < 0.29){
		 return mu;
		}
	 else {
		   return std::pow(6.0*Urand(),1.0/3.0)-1.0 * copysign(1.0,mu);
		}
	}
 }

 