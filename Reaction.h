#ifndef _REACTION_HEADER_
#define _REACTION_HEADER_

#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include <stack>
#include <utility>

#include "Particle.h"
#include "Distribution.h"
#include "Caffeine.h"

class reaction {
  protected:
    std::string rxn_name;
  public:
     reaction() {};
    ~reaction() {};

	virtual double xs(double E) = 0;
    virtual std::string name() final { return rxn_name; };
    virtual void sample( particle* p, std::stack<particle>* bank ) = 0;
};

class capture_reaction : public reaction {
  private:
    
  public:
     capture_reaction() { rxn_name = "capture"; };
    ~capture_reaction() {};

	virtual double xs(double E) = 0;
    virtual void sample( particle* p, std::stack<particle>* bank ) = 0;
};

class SE_capture_reaction : public capture_reaction {
  private:
    double rxn_xs;
  public:
    SE_capture_reaction( double x ) : rxn_xs( x ) {};
   ~SE_capture_reaction() {}; 

	double xs(double E) { return rxn_xs; };
    void sample(particle* p, std::stack<particle>* bank );
};

class CE_capture_reaction : public capture_reaction {
  private:
	std::shared_ptr <caffeine> energy_dependence;
  public:
    CE_capture_reaction( std::shared_ptr< caffeine > eng_dep )
    : energy_dependence(eng_dep) {}; 
   ~CE_capture_reaction() {};

	double xs(double E) { return energy_dependence->sample(E); };
    void sample( particle* p, std::stack<particle>* bank );
};

class scatter_reaction : public reaction {
  protected:
    std::shared_ptr< distribution<double> > scatter_dist; 
  public:
     scatter_reaction( std::shared_ptr< distribution<double> > D ) :
       scatter_dist(D) { rxn_name = "scatter"; };
    ~scatter_reaction() {};

	 virtual double xs(double E) = 0;
     virtual void sample( particle* p, std::stack<particle>* bank ) = 0;
};

class SE_scatter_reaction : public scatter_reaction {
  private:
	double rxn_xs;
  public:
    SE_scatter_reaction( double x, std::shared_ptr< distribution<double> > D ) 
    : rxn_xs(x), scatter_reaction(D) {};
   ~SE_scatter_reaction() {};

	double xs(double E) { return rxn_xs; };
    void sample( particle* p, std::stack<particle>* bank );
};

class CE_scatter_reaction : public scatter_reaction {
  private:
	double A;
	std::shared_ptr <caffeine> energy_dependence;
  public:
    CE_scatter_reaction( std::shared_ptr< caffeine > eng_dep, std::shared_ptr< distribution<double> > D, double a ) 
    :  scatter_reaction(D), energy_dependence(eng_dep), A(a) {};
   ~CE_scatter_reaction() {};

	double xs(double E) { return energy_dependence->sample(E); };
    void sample( particle*p, std::stack<particle>* bank );
};


class fission_reaction : public reaction {
  protected:
    std::shared_ptr< distribution<int> >   multiplicity_dist; 
    std::shared_ptr< distribution<point> > isotropic;
  public:
     fission_reaction( std::shared_ptr< distribution<int> > D ) :
       multiplicity_dist(D) { 
         rxn_name = "fission";
         isotropic = std::make_shared< isotropicDirection_distribution > ( "isotropic" ); 
       };
    ~fission_reaction() {};

	virtual double xs(double E) = 0;
    virtual void sample( particle* p, std::stack<particle>* bank ) = 0;
};

class SE_fission_reaction : public fission_reaction {
  private:
    double rxn_xs;
  public:
     SE_fission_reaction( double x, std::shared_ptr< distribution<int> > D ) :
       rxn_xs(x), fission_reaction( D ) {}; 
    ~SE_fission_reaction() {};

	double xs(double E) { return rxn_xs; };
    void sample( particle* p, std::stack<particle>* bank );
};

class CE_fission_reaction : public fission_reaction {
  private:
    std::shared_ptr <caffeine> energy_dependence;
	std::shared_ptr < distribution <double> > fis_eng_dist;
  public:
     CE_fission_reaction( std::shared_ptr< caffeine > eng_dep, std::shared_ptr< distribution<int> > D, std::shared_ptr< distribution<double> > DE ) :
       energy_dependence(eng_dep), fission_reaction( D ), fis_eng_dist(DE) {}; 
    ~CE_fission_reaction() {};

	double xs(double E) { return energy_dependence->sample(E); };
    void sample( particle* p, std::stack<particle>* bank );
};


#endif
