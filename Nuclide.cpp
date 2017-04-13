#include <vector>
#include <memory>
#include <cassert>

#include "Random.h"
#include "Nuclide.h"

// add a new reaction to the current nuclide
void nuclide::addReaction( std::shared_ptr< reaction > R ) { rxn.push_back( R ); }

// return the total microscopic cross section
double nuclide::total_xs(double E) {
  double xs = 0.0;
  for ( auto r : rxn ) { xs += r->xs(E); }
  return xs;
}

// randomly sample a reaction type from this nuclide
std::shared_ptr< reaction > nuclide::sample_reaction(double E) {
  double u = total_xs(E) * Urand();
  double s = 0.0;
  for ( auto r : rxn ) {
    s += r->xs(E);
    if ( s > u ) { return r; }
  }
  assert( false ); // should never reach here
  return nullptr;
}
