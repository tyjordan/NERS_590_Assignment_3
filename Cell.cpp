#include <iostream>
#include <vector>
#include <utility>
#include <memory>
#include <cassert>
#include <cmath>
#include <limits>

#include "Cell.h"
#include "Particle.h"

// this is a test
// take in a pair of surface pointer and integer describing sense (must not be zero!)
// and append to vector of surfaces
void cell::addSurface( std::shared_ptr< surface > S, int sense ) {

  assert( sense != 0 );
  int sgn = std::copysign( 1, sense );

  surfaces.push_back( std::make_pair( S, sgn ) );
}

void cell::attachEstimator( std::shared_ptr <estimator> E) {
	cell_estimators.push_back( E );
	if( ( E->estimator_type() == "pathLengthFlux" ) || ( E->estimator_type() == "time_binned_pLF" )) {
		std::vector< std::pair< std::shared_ptr< nuclide >, double > > nuc = getMaterial()->getNuclides();
		for(auto n: nuc) {
			std::vector< std::shared_ptr< reaction > > rxn = n.first->getReactions();
			for(auto r: rxn) {
				if( ( r->name() == E->reaction_type() ) || ( E->reaction_type() == "total" ) ) {
					E->add_to_list( r, getMaterial()->atom_density() * n.second );
				}
			}
		}
	}
	return;
}

// test if point p inside the current cell
bool cell::testPoint( point p ) {

  // loop over surfaces in cell, if not on correct side return false
  // if on correct side of all surfaces, particle is in the cell and return true
  for ( auto s : surfaces ) {
    // first = surface pointer, second = +/- 1 indicating sense
    if ( s.first->eval( p ) * s.second < 0 ) { return false; }  
  }
  return true;
}

// find first intersecting surface of ray r and distance to intersection
std::pair< std::shared_ptr< surface >, double > cell::surfaceIntersect( ray r ) {

  double dist = std::numeric_limits<double>::max();
  std::shared_ptr< surface > S = nullptr;
  for ( auto s : surfaces ) {
    // first is a surface pointer; distance is always positive or huge if invalid
    double d = s.first->distance( r );
    if ( d < dist ) {
      // current surface intersection is closer
      dist = d;
      S    = s.first;
    }
  }
  return std::make_pair( S, dist );
}

void cell::moveParticle( particle* p, double s ) {
  p->move( s );
}

void cell::sampleCollision( particle* p, std::stack<particle>* bank ) {
  cell_material->sample_collision( p, bank );   
}

