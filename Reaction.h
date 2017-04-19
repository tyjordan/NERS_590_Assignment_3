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
	double rxn_xs;
    std::string rxn_name;
  public:
     reaction( double x ) : rxn_xs(x) {};
    ~reaction() {};

    virtual std::string name() final { return rxn_name; };
    virtual double xs(double E) = 0;
    virtual void sample( particle* p, std::stack<particle>* bank ) = 0;
};

class capture_reaction : public reaction {
  private:
    
  public:
     capture_reaction( double x ) : reaction(x) { rxn_name = "capture"; };
    ~capture_reaction() {};

	virtual double xs(double E) = 0;
    virtual void sample( particle* p, std::stack<particle>* bank ) = 0;
};

class SE_capture_reaction : public capture_reaction {
  private:

  public:
    SE_capture_reaction( double x ) : capture_reaction( x ) {};
   ~SE_capture_reaction() {}; 

	double xs(double E) { return rxn_xs; };
    void sample(particle* p, std::stack<particle>* bank );
};

class CE_capture_reaction : public capture_reaction {
  private:
	std::shared_ptr <caffeine> energy_dependence;
  public:
    CE_capture_reaction( double x, std::shared_ptr< caffeine > eng_dep )
    : capture_reaction( x ), energy_dependence(eng_dep) {}; 
   ~CE_capture_reaction() {};

	double xs(double E) { return rxn_xs * energy_dependence->sample(E); };
    void sample( particle* p, std::stack<particle>* bank );
};

class scatter_reaction : public reaction {
  protected:
    std::shared_ptr< distribution<double> > scatter_dist; 
  public:
     scatter_reaction( double x, std::shared_ptr< distribution<double> > D ) :
       reaction(x), scatter_dist(D) { rxn_name = "scatter"; };
    ~scatter_reaction() {};

	 virtual double xs(double E) = 0;
     virtual void sample( particle* p, std::stack<particle>* bank ) = 0;
};

class SE_scatter_reaction : public scatter_reaction {
  private:

  public:
    SE_scatter_reaction( double x, std::shared_ptr< distribution<double> > D ) 
    :  scatter_reaction( x, D ) {};
   ~SE_scatter_reaction() {};

	double xs(double E) { return rxn_xs; };
    void sample( particle* p, std::stack<particle>* bank );
};

class CE_scatter_reaction : public scatter_reaction {
  private:
	double A;
	std::shared_ptr <caffeine> energy_dependence;
  public:
    CE_scatter_reaction( double x, std::shared_ptr< caffeine > eng_dep, std::shared_ptr< distribution<double> > D, double a ) 
    :  scatter_reaction( x, D ), energy_dependence(eng_dep), A(a) {};
   ~CE_scatter_reaction() {};

	double xs(double E) { return rxn_xs * energy_dependence->sample(E); };
    void sample( particle*p, std::stack<particle>* bank );
};



class fission_reaction : public reaction {
  protected:
    std::shared_ptr< distribution<int> >   multiplicity_dist; 
    std::shared_ptr< distribution<point> > isotropic;
  public:
     fission_reaction( double x, std::shared_ptr< distribution<int> > D ) :
       reaction(x), multiplicity_dist(D) { 
         rxn_name = "fission";
         isotropic = std::make_shared< isotropicDirection_distribution > ( "isotropic" ); 
       };
    ~fission_reaction() {};

	virtual double xs(double E) = 0;
    virtual void sample( particle* p, std::stack<particle>* bank ) = 0;
};

class SE_fission_reaction : public fission_reaction {
  private:
    
  public:
     SE_fission_reaction( double x, std::shared_ptr< distribution<int> > D ) :
       fission_reaction( x, D ) {}; 
    ~SE_fission_reaction() {};

	double xs(double E) { return rxn_xs; };
    void sample( particle* p, std::stack<particle>* bank );
};

class CE_fission_reaction : public fission_reaction {
  private:
    std::shared_ptr <caffeine> energy_dependence;
	std::shared_ptr < distribution <double> > fis_eng_dist;
  public:
     CE_fission_reaction( double x, std::shared_ptr< caffeine > eng_dep, std::shared_ptr< distribution<int> > D, std::shared_ptr< distribution<double> > DE ) :
       fission_reaction( x, D ), energy_dependence(eng_dep), fis_eng_dist(DE) {}; 
    ~CE_fission_reaction() {};

	double xs(double E) { return rxn_xs * energy_dependence->sample(E); };
    void sample( particle* p, std::stack<particle>* bank );
};


#endif
