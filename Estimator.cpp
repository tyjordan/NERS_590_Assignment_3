#include <cmath>
#include <iostream>

#include "Estimator.h"
#include "Material.h"
#include "Particle.h"

void surface_current_estimator::score( particle* p, double null ) { tally_hist += p->wgt() ; }

double cell_pathLengthFlux_estimator::total_rxn_xs(double E) {
	double xs = 0.0;
	for (auto r: rxn_list) { xs += r.first->xs(E) * r.second; }
	return xs;
}

void cell_pathLengthFlux_estimator::score( particle* p, double path_length ) {
	tally_hist += p->wgt() * path_length * total_rxn_xs( p->energy() );
	return;
}

void counting_estimator::score( particle* p, double null ) { count_hist++; }

void counting_estimator::endHistory() {
  if ( tally.size() < count_hist + 1 ) { tally.resize( count_hist + 1, 0.0 ); }
  tally[ count_hist ] += 1.0;
  nhist++;
  count_hist = 0;
}

void counting_estimator::report( int T ) {
  std::cout << name() << std::endl;
  double s1 = 0.0, s2 = 0.0;
  for ( int i = 0 ; i < tally.size() ; i++ ) {
    double p = tally[i] / nhist;
    std::cout << i << " " << p << "   " << std::sqrt( p * ( 1.0 - p ) / nhist ) / p <<  std::endl;
    s1 += p * i;
    s2 += p * i * i;
  }
  std::cout << "   mean = " << s1 << std::endl;
  std::cout << "   var  = " << s2 - s1*s1 << std::endl;
}
