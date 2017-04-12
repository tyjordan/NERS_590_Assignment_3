#ifndef _SOURCE_HEADER_
#define _SOURCE_HEADER_

#include <stack>
#include <memory>

#include "Point.h"
#include "Distribution.h"
#include "Particle.h"

class source {
  private:
    std::shared_ptr< distribution<point> > dist_pos;
    std::shared_ptr< distribution<point> > dist_dir;
    std::shared_ptr< distribution<double>> energy;
  public:
     source( std::shared_ptr< distribution<point> > pos, 
     	       std::shared_ptr< distribution<point> > dir, 
     	       std::shared_ptr< distribution<double>> eng )
           : dist_pos(pos), dist_dir(dir), energy(eng) {};
    ~source() {};

    std::stack<particle> sample();
};

#endif
