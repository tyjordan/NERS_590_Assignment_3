#include <cmath>
#include <iostream>
#include <string>

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

time_binned_pLF_est::time_binned_pLF_est( std::string label, std::string rt, double min, double max, int binnum ) :
	estimator(label), rxn_type(rt) {

	est_type = "time_binned_pLF";

	double binwidth = (max - min) / binnum;
	double t = min + binwidth;
	upper_bin_edges.push_back(t);
	while(t < max) {
		t += binwidth;
		upper_bin_edges.push_back(t);
	}

	for(int i = 0; i < upper_bin_edges.size(); i++) {
		cell_pathLengthFlux_estimator b ( std::to_string( upper_bin_edges.at(i) ), rt );
		bins.push_back(b);
	}
}

void time_binned_pLF_est::add_to_list( std::shared_ptr <reaction> R, double N ) {
	for(auto b: bins) { b.add_to_list(R, N); }
	return;
}

void time_binned_pLF_est::score( particle* p, double path_length) {

	double start = p->time();
	double v = std::sqrt( 2 * p->energy() * 1.602e-13 / p->mass() ) * 100;
	double stop = start + (path_length / v);

	for(int i = 0; i < upper_bin_edges.size(); i++) {
		if( start < upper_bin_edges.at(i) ) {
			if( stop < upper_bin_edges.at(i) ) {
				bins.at(i).score(p, path_length);
			}
			else {
				bins.at(i).score( p, v * (upper_bin_edges.at(i) - start) );
				for(int j = (i + 1); j < upper_bin_edges.size(); j++) {
					if( stop < upper_bin_edges.at(j) ) {
						bins.at(j).score( p, v * ( stop - upper_bin_edges.at(j - 1) ) );
						break;
					}
					else {
						bins.at(j).score( p, v * ( upper_bin_edges.at(j) - upper_bin_edges.at(j - 1) ) );
					}
				}
			}
			break;
		}
	}
}
						
void time_binned_pLF_est::endHistory() {
	for(auto b: bins) { b.endHistory(); }
	return;
}

void time_binned_pLF_est::report( int T ) {
	std::cout << name() << std::endl;
	for(auto b: bins) { b.report(T); }
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
