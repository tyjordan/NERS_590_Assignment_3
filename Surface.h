#ifndef _SURFACE_HEADER_
#define _SURFACE_HEADER_

#include <iostream>
#include <string>
#include <vector>
#include <limits>

#include "Point.h"
#include "Particle.h"
#include "Estimator.h"

class surface {
  private:
    std::string surface_name;
    bool        reflect_bc;
	std::vector< std::shared_ptr< estimator > > surface_estimators;
  public:

     surface( std::string label ) : surface_name(label) { reflect_bc = false; };
    ~surface() {};

    virtual std::string name()    final { return surface_name; };
    virtual void makeReflecting() final { reflect_bc = true; };

    virtual void attachEstimator( std::shared_ptr< estimator > E ) final {
      surface_estimators.push_back( E );
    };

    void crossSurface( particle* p, double d ) { 
      // score estimators
      for ( auto e : surface_estimators ) { e->score( p, d ); }
      // reflect if needed
      if ( reflect_bc ) { 
        point d = reflect( p->getRay() );
        p->setDirection( d );
      }

      // advance particle off the surface
      p->move( std::numeric_limits<float>::epsilon() );
    };

    virtual double eval( point p )   = 0;
    virtual double distance( ray r ) = 0;
    virtual point  reflect( ray r )  = 0;
};

class plane : public surface {
  private:
    double a, b, c, d;
  public:
     plane( std::string label, double p1, double p2, double p3, double p4 ) : 
       surface(label), a(p1), b(p2), c(p3), d(p4) {};
    ~plane() {};

     double eval( point p );
     double distance( ray r );
     point  reflect( ray r );
};

class sphere : public surface {
  private:
    double x0, y0, z0, rad;
  public:
     sphere( std::string label, double p1, double p2, double p3, double p4 ) : 
       surface(label), x0(p1), y0(p2), z0(p3), rad(p4) {};
    ~sphere() {};

     double eval( point p );
     double distance( ray r );
     point  reflect( ray r );
};

//cylinders centered on axes:

class x_cylinder : public surface {
	private:
		point pos;
		double rad;
	public:
		x_cylinder( std::string label, point P, double R ) : surface(label), pos(P), rad(R) {};
		~x_cylinder() {};

		double eval( point p );
		double distance( ray r );
		point reflect( ray r );
};

class y_cylinder : public surface {
	private:
		point pos;
		double rad;
	public:
		y_cylinder( std::string label, point P, double R ) : surface(label), pos(P), rad(R) {};
		~y_cylinder() {};

		double eval( point p );
		double distance( ray r );
		point reflect( ray r );
};

class z_cylinder : public surface {
	private:
		point pos;
		double rad;
	public:
		z_cylinder( std::string label, point P, double R ) : surface(label), pos(P), rad(R) {};
		~z_cylinder() {};

		double eval( point p );
		double distance( ray r );
		point reflect( ray r );
};

class x_cone : public surface {
	private:
		point pos;
		double m;
	public:
		x_cone( std::string label, point P, double M ) : surface(label), pos(P), m(M) {};
		~x_cone() {};

		double eval( point p );
		double distance( ray r );
		point reflect( ray r ) { return r.pos; };
};

#endif
