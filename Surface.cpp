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
double x_cone::eval(point p){
	double dx = p.x - pos.x;
	double dy = p.y - pos.y;
	double dz = p.z - pos.z; 
	return (-(rad * rad)*(dx * dx)) +(dy * dy) + (dz * dz);
}
double x_cone::distance( ray r) {
	double dx = r.pos.x - pos.x; 
	double dy = r.pos.y - pos.y;
	double dz = r.pos.z - pos.z;
	
	double a = -(rad * rad)*(r.dir.x * r.dir.x) + (r.dir.y * r.dir.y) + (r.dir.z * r.dir.z);
	double b = 2 * ( (rad*rad)*(r.pos.x*r.dir.x - pos.x*r.dir.x) - (r.pos.y*r.dir.z) - (r.pos.z*r.dir.z) );
	double c = eval(r.pos);
	
	return quad_solve(a, b, c);
}
point x_cone::reflect(ray r){
	double dx = r.pos.x - pos.x;
	double dy = r.pos.y - pos.y; 
	double dz = r.pos.z - pos.z; 
	
	double fx = (2 * (dx*r.dir.x + dy*r.dir.y - std::pow(rad,2.0)*dz*r.dir.z))/((rad * rad)*(1 + rad * rad)*(dz * dz));
    double fy = (2 * (dx*r.dir.x + dy*r.dir.y - std::pow(rad,2.0)*dz*r.dir.z)* dy)/((rad * rad)*(1 + rad * rad)*(dz * dz));
    double fz = (2 * (dx*r.dir.x + dy*r.dir.y - std::pow(rad,2.0)*dz*r.dir.z)* dz)/((rad * rad)*(1 + rad * rad)*(dz * dz));
	
	return point(r.dir.x - fx,r.dir.y - fy, r.dir.z - fz);
}
double y_cone::eval(point p){
	double dx = p.x - pos.x;
	double dy = p.y - pos.y;
	double dz = p.z - pos.z; 
	return (dx * dx) - ((rad * rad)*(dy * dy)) + (dz * dz);
}
double y_cone::distance( ray r) {
	double dx = r.pos.x - pos.x; 
	double dy = r.pos.y - pos.y;
	double dz = r.pos.z - pos.z;
	
	double a = -(rad * rad)*(r.dir.x * r.dir.x) + (r.dir.y * r.dir.y) + (r.dir.z * r.dir.z);
	double b = 2 * ((dx * dx) - ((rad * rad)*(dy * dy)) + (dz * dz));
	double c = eval(r.pos);
	
	return quad_solve(a, b, c);
}
point y_cone::reflect(ray r){
	double dx = r.pos.x - pos.x;
	double dy = r.pos.y - pos.y; 
	double dz = r.pos.z - pos.z; 
	
	double fx = (2 * (dx*r.dir.x + dy*r.dir.y - std::pow(rad,2.0)*dz*r.dir.z)*dx )/((rad * rad)*(1 + rad * rad)*(dz * dz));
    double fy = (2 * (dx*r.dir.x + dy*r.dir.y - std::pow(rad,2.0)*dz*r.dir.z))/((rad * rad)*(1 + rad * rad)*(dz * dz));
    double fz = (2 * (dx*r.dir.x + dy*r.dir.y - std::pow(rad,2.0)*dz*r.dir.z)* dz)/((rad * rad)*(1 + rad * rad)*(dz * dz));
	
	return point(r.dir.x - fx,r.dir.y - fy, r.dir.z - fz);
}
double z_cone::eval(point p){
	double dx = p.x - pos.x;
	double dy = p.y - pos.y;
	double dz = p.z - pos.z; 
	return (dx * dx)+(dy * dy) -((rad * rad)*(dz * dz));
}
double z_cone::distance( ray r) {
	double dx = r.pos.x - pos.x; 
	double dy = r.pos.y - pos.y;
	double dz = r.pos.z - pos.z;
	
	double a = (r.dir.x * r.dir.x) + (r.dir.y * r.dir.y) - ((rad * rad)*(r.dir.z * r.dir.z));
	double b = 2 * ((dx * dx) + (dy * dy) - ((rad * rad)*(dz * dz)));
	double c = eval(r.pos);
	
	return quad_solve(a, b, c);
}
point z_cone::reflect(ray r){
	double dx = r.pos.x - pos.x;
	double dy = r.pos.y - pos.y; 
	double dz = r.pos.z - pos.z; 
	
	double fx = (2 * (dx*r.dir.x + dy*r.dir.y - std::pow(rad,2.0)*dz*r.dir.z)*dx )/((rad * rad)*(1 + rad * rad)*(dz * dz));
    double fy = (2 * (dx*r.dir.x + dy*r.dir.y - std::pow(rad,2.0)*dz*r.dir.z)*dy)/((rad * rad)*(1 + rad * rad)*(dz * dz));
    double fz = (2 * (dx*r.dir.x + dy*r.dir.y - std::pow(rad,2.0)*dz*r.dir.z))/((rad * rad)*(1 + rad * rad)*(dz * dz));
	
	return point(r.dir.x - fx,r.dir.y - fy, r.dir.z - fz);
}
