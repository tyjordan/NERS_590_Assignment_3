#ifndef _ESTIMATOR_HEADER_
#define _ESTIMATOR_HEADER_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <typeinfo>
#include <utility>
#include <memory>

#include "Particle.h"
#include "Material.h"
#include "Reaction.h"

class estimator {
  private:
    std::string estimator_name;
  protected:
    unsigned long long nhist;
	std::string est_type;
  public:
     estimator( std::string label ) : estimator_name(label) {};
    ~estimator() {};

    virtual std::string name() final { return estimator_name; };
	virtual std::string estimator_type() final { return est_type; };

    virtual void score( particle*, double ) = 0;

    virtual void endHistory() = 0;
    virtual void report( int T ) = 0;

	virtual std::string reaction_type() { return "null"; };
	virtual void add_to_list( std::shared_ptr <reaction> R, double N ) { return; };

};

class single_valued_estimator : public estimator {
  private:

  protected:
    double tally_hist, tally_sum, tally_squared;
  public:

     single_valued_estimator(std::string label ) : estimator(label) { 
       nhist         = 0;
       tally_hist    = 0.0;   
       tally_sum     = 0.0; 
       tally_squared = 0.0;
     };
    ~single_valued_estimator() {};

     virtual void endHistory()    final { 
       nhist++;
       tally_sum     += tally_hist;
       tally_squared += tally_hist * tally_hist;
       tally_hist = 0.0; }

     virtual void score( particle*, double ) = 0;

     virtual void report( int T ) final { //T is the number of particle tracks for FOM calculation
       double mean = tally_sum / nhist;
       double var  = ( tally_squared / nhist - mean*mean ) / nhist;
		if( T == 0 ) {
       		std::cout << name() << "   " << mean << "   " << std::sqrt( var ) / mean << std::endl;
		}
		else {
			std::cout << name() << "   " << mean << "   " << std::sqrt( var ) / mean << std::endl;
			std::cout << "		FOM: " << 1 / ( ( var / ( mean * mean ) ) * T ) << std:: endl;
		}
     };
};

class surface_current_estimator : public single_valued_estimator {
  private:

  public:
     surface_current_estimator( std::string label ) : single_valued_estimator(label)
		{ est_type = "surface_current"; };
    ~surface_current_estimator() {};

    void score( particle*, double );
};

class cell_pathLengthFlux_estimator : public single_valued_estimator {
  private:
	std::string rxn_type;
	std::vector < std::pair < std::shared_ptr <reaction>, double > > rxn_list;
	double total_rxn_xs(double E);
  public:
     cell_pathLengthFlux_estimator( std::string label, std::string rt ) : 
       single_valued_estimator(label), rxn_type(rt) { est_type = "pathLengthFlux"; };
    ~cell_pathLengthFlux_estimator() {};

	std::string reaction_type() { return rxn_type; };
	void add_to_list( std::shared_ptr <reaction> R, double N )
		{ rxn_list.push_back( std::make_pair(R, N) ); return; };

    void score( particle*, double );
};

class time_binned_pLF_est : public estimator {
	private:
		std::string rxn_type;
		std::vector <double> upper_bin_edges;
		std::vector < std::shared_ptr <cell_pathLengthFlux_estimator> > bins;
	public:
		time_binned_pLF_est( std::string label, std::string rt, double min, double max, int binnum );
		~time_binned_pLF_est() {};

		std::string reaction_type() { return rxn_type; };
		void add_to_list( std::shared_ptr <reaction> R, double N );

		void score( particle*, double );

    	void endHistory();
    	void report( int T );
};

class counting_estimator : public estimator {
  private:
    int count_hist;
    std::vector< double > tally;
  public:
     counting_estimator( std::string label ) : estimator(label)
		{ count_hist = 0; est_type = "counting"; };
    ~counting_estimator() {};

    void score( particle*, double);
    void endHistory();
    void report( int T );
};

#endif
