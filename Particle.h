#ifndef _PARTICLE_HEADER_
#define _PARTICLE_HEADER_

#include <memory>

#include "Point.h"

class cell;

class particle {
  private:
    point  p_pos, p_dir;
    double p_wgt, p_eng, p_mass, p_time;
    bool   exist;
    std::shared_ptr< cell > p_cell;
  public:
     particle( point p, point d, double e );
    ~particle() {};

    point pos() { return p_pos; };      // return particle position
    point dir() { return p_dir; };      // return particle direction 
    double wgt() { return p_wgt; };     // return particle weight
    bool alive() { return exist; };     // return particle state flag
    double energy() { return p_eng; };  // return particle energy
    double mass() { return p_mass };    // return particle mass
    double time() { return p_time };    // return time of particle life

    ray getRay() { return ray( p_pos, p_dir ); }

    std::shared_ptr< cell > cellPointer() { return p_cell; }

    void move( double s );
    void scatter( double mu0 );
    void kill();
    void setDirection( point p );
    void adjustWeight( double f );
    void recordCell( std::shared_ptr< cell > cel );
    void setEnergy( double new_energy );
    void setMass( double particle_mass);
    void updateTime( double s );
};

#endif
