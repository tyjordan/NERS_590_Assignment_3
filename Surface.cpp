#include <cmath>
#include <limits>
#include <cassert>

#include "Point.h"
#include "QuadSolver.h"
#include "Surface.h"

// evaluates the surface equation w.r.t. to point p
double plane::eval( point p ) {
  return a * p.x  +  b * p.y  +  c * p.z  - d;
}

// determines the mininum positive distance to intersection for a ray r
// (returns a very large number if no intersection along ray for ease of calculation down the line)
double plane::distance( ray r ) {
  point p = r.pos;
  point u = r.dir;

  double denom = a * u.x  +  b * u.y  +  c * u.z;
  if ( std::fabs( denom ) > 100.0 * std::numeric_limits<double>::epsilon() ) {
    double dist = ( d - a * p.x - b * p.y - c * p.z ) / denom;
    if ( dist > 0.0 ) { return dist; }
    else { return std::numeric_limits<double>::max(); }
  }
  else {
    // moving in a direction that is (or is very close to) parallel to the surface
    return std::numeric_limits<double>::max();
  }
}

// get new reflected direction
point plane::reflect( ray r ) {
   assert( std::fabs( eval( r.pos ) ) < std::numeric_limits<float>::epsilon() );

   point  u = r.dir;
   double t = 2.0 * ( a * u.x  +  b * u.y  +  c * u.z ) / ( a*a + b*b + c*c );
   return point( u.x - a*t, u.y - b*t, u.z - c*t );
}


double sphere::eval( point p ) {
  return std::pow( p.x - x0, 2 ) + std::pow( p.y - y0, 2 ) + std::pow( p.z - z0, 2 )  - rad*rad;
}

double sphere::distance( ray r ) {
  point p = r.pos;
  point u = r.dir;

  // difference between each coordinate and current point
  point q( p.x - x0, p.y - y0, p.z - z0 );

  // put into quadratic equation form: a*s^2 + b*s + c = 0, where a = 1
  double b = 2.0 * ( q.x * u.x  +  q.y * u.y  +  q.z * u.z );
  double c = eval( p );

  return quad_solve( 1.0, b, c );
}

point sphere::reflect( ray r ) {
   assert( std::fabs( eval( r.pos ) ) < std::numeric_limits<float>::epsilon() );

   point p = r.pos;
   point u = r.dir;

   point q( p.x - x0, p.y - y0, p.z - z0 );

   double t = 2.0 * ( q.x * u.x  +  q.y * u.y  +  q.z * u.z ) / ( rad*rad );
   return point( u.x - q.x * t,  u.y - q.y * t,  u.z - q.z * t );
}

double x_cylinder::eval( point p ) {
	double dy = p.y - pos.y;
	double dz = p.z - pos.z;
	return ( dy * dy ) + ( dz * dz ) - ( rad * rad );
}

double x_cylinder::distance( ray r ) {
	double dy = r.pos.y - pos.y;
	double dz = r.pos.z - pos.z;

	double a = ( r.dir.y * r.dir.y ) + ( r.dir.z * r.dir.z );
	double b = 2 * ( ( r.dir.y * dy ) + ( r.dir.z * dz ) );
	double c = eval(r.pos);

	return quad_solve( a, b, c );
}

point x_cylinder::reflect( ray r ) {
	double dy = r.pos.y - pos.y;
	double dz = r.pos.z - pos.z;

	double fy = ( 2 * ( ( dy * r.dir.y ) + ( dz * r.dir.z ) ) * dy ) / ( rad * rad );
	double fz = ( 2 * ( ( dy * r.dir.y ) + ( dz * r.dir.z ) ) * dz ) / ( rad * rad );

	return point( r.dir.x, ( r.dir.y - fy ), ( r.dir.z - fz ) );
}

double y_cylinder::eval( point p ) {
	double dx = p.x - pos.x;
	double dz = p.z - pos.z;
	return ( dx * dx ) + ( dz * dz ) - ( rad * rad );
}

double y_cylinder::distance( ray r ) {
	double dx = r.pos.x - pos.x;
	double dz = r.pos.z - pos.z;

	double a = ( r.dir.x * r.dir.x ) + ( r.dir.z * r.dir.z );
	double b = 2 * ( ( r.dir.x * dx ) + ( r.dir.z * dz ) );
	double c = eval(r.pos);

	return quad_solve( a, b, c );
}

point y_cylinder::reflect( ray r ) {
	double dx = r.pos.x - pos.x;
	double dz = r.pos.z - pos.z;

	double fx = ( 2 * ( ( dx * r.dir.x ) + ( dz * r.dir.z ) ) * dx ) / ( rad * rad );
	double fz = ( 2 * ( ( dx * r.dir.x ) + ( dz * r.dir.z ) ) * dz ) / ( rad * rad );

	return point( ( r.dir.x - fx ), r.dir.y, ( r.dir.z - fz ) );
}

double z_cylinder::eval( point p ) {
	double dx = p.x - pos.x;
	double dy = p.y - pos.y;
	return ( dx * dx ) + ( dy * dy ) - ( rad * rad );
}

double z_cylinder::distance( ray r ) {
	double dx = r.pos.x - pos.x;
	double dy = r.pos.y - pos.y;

	double a = ( r.dir.x * r.dir.x ) + ( r.dir.y * r.dir.y );
	double b = 2 * ( ( r.dir.x * dx ) + ( r.dir.y * dy ) );
	double c = eval(r.pos);

	return quad_solve( a, b, c );
}

point z_cylinder::reflect( ray r ) {
	double dx = r.pos.x - pos.x;
	double dy = r.pos.y - pos.y;

	double fx = ( 2 * ( ( dx * r.dir.x ) + ( dy * r.dir.y ) ) * dx ) / ( rad * rad );
	double fy = ( 2 * ( ( dx * r.dir.x ) + ( dy * r.dir.y ) ) * dy ) / ( rad * rad );

	return point( ( r.dir.x - fx ), ( r.dir.y - fy ), r.dir.z );
}

double x_cone::eval( point p ) {
 return std::pow( (p.y - pos.y), 2 ) + std::pow( (p.z - pos.z), 2 )  - m * m * std::pow( (p.x - pos.x), 2 );
}

double x_cone::distance( ray r ) {
 point a = r.pos;
 point c = r.dir;
 point t( a.x - pos.x, a.y - pos.y, a.z - pos.z );
 
	double a1 = ( c.y * c.y + c.z * c.z - m * m*c.x * c.x );
	double b = 2.0 * ( t.y * c.y  +  t.z * c.z - m * m*t.x * c.x);
	double c1 = eval( a );

	return quad_solve(a1, b, c1);
}
