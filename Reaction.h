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

class reaction {
  private:
    double rxn_xs;
  protected:
    std::string rxn_name;
  public:
     reaction( double x ) : rxn_xs(x) {};
    ~reaction() {};

    virtual std::string name() final { return rxn_name; };
    virtual double xs() final { return rxn_xs; };
    virtual void sample( particle* p, std::stack<particle>* bank ) = 0;
};

class capture_reaction : public reaction {
  private:
    
  public:
     capture_reaction( double x ) : reaction(x) { rxn_name = "capture"; };
    ~capture_reaction() {};

    virtual void sample( particle* p, std::stack<particle>* bank );
};

class SE_capture_reaction : public capture_reaction {
  private:

  public:
    SE_capture_reaction( double x ) : capture_reaction( x );
   ~SE_capture_reaction() {}; 

    void sample(particle* p, std::stack<particle>* bank );
};

class CE_capture_reaction : public capture_reaction {
  private:

  public:
    CE_capture_reaction( double x, std::shared_ptr< Caffeine > eng_dep )
    : capture_reaction( x ), energy_dependence(eng-dep) ; 
   ~CE_capture_reaction() {};

   void sample( particle* p, std::stack<particle>* bank );
};


class scatter_reaction : public reaction {
  private:
    std::shared_ptr< distribution<double> > scatter_dist; 
  public:
     scatter_reaction( double x, std::shared_ptr< distribution<double> > D ) :
       reaction(x), scatter_dist(D) { rxn_name = "scatter"; };
    ~scatter_reaction() {};

     virtual void sample( particle* p, std::stack<particle>* bank );
};

class SE_scatter_reaction : scatter_reaction {
  private:

  public:
    SE_scatter_reaction( double x, std::shared_ptr< distribution<double> > D ) 
    :  scatter_reaction( double x, std::shared_ptr< distribution<double> > D );
   ~SE_scatter_reaction() {};

    void sample( particle*p, std::stack<particle>* p );
};

class CE_scatter_reaction : scatter_reaction {
  private:

  public:
    CE_scatter_reaction( double x, std::shared_ptr< distribution<double> > D, std::shared_ptr< Caffeine > eng_dep ) 
    :  scatter_reaction( double x, std::shared_ptr< distribution<double> > D ), energy_dependence(eng_dep);
   ~CE_scatter_reaction() {};

    void sample( particle*p, std::stack<particle>* p );
};



class fission_reaction : public reaction {
  private:
    std::shared_ptr< distribution<int> >   multiplicity_dist; 
    std::shared_ptr< distribution<point> > isotropic;
  public:
     fission_reaction( double x, std::shared_ptr< distribution<int> > D ) :
       reaction(x), multiplicity_dist(D) { 
         rxn_name = "fission";
         isotropic = std::make_shared< isotropicDirection_distribution > ( "isotropic" ); 
       };
    ~fission_reaction() {};

    void sample( particle* p, std::stack<particle>* bank );
};

class SE_fission_reaction : public reaction {
  private:
    
  public:
     SE_fission_reaction( double x, std::shared_ptr< distribution<int> > D ) :
       fission_reaction( double x, std::shared_ptr< distribution<int> > D ); 
    ~SE_fission_reaction() {};

    void sample( particle* p, std::stack<particle>* bank );
};

class CE_fission_reaction : public reaction {
  private:
    
  public:
     CE_fission_reaction( double x, std::shared_ptr< distribution<int> > D, std::shared_ptr< Caffeine > eng_dep ) :
       fission_reaction( double x, std::shared_ptr< distribution<int> > D ), energy_dependence(eng_dep); 
    ~CE_fission_reaction() {};

    void sample( particle* p, std::stack<particle>* bank );
};


#endif
