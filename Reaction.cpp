#include "Reaction.h"
#include "Particle.h"
#include "Distribution.h"

void  SE_capture_reaction::sample( particle* p, std::stack< particle >* bank ) {
  // kill the particle and leave the bank unmodified
  p->kill();
}

void  CE_capture_reaction::sample( particle* p, std::stack< particle >* bank ) {
  // kill the particle and leave the bank unmodified
  p->kill();
}


void  SE_scatter_reaction::sample( particle* p, std::stack< particle >* bank ) {
  // scatter the particle and leave the bank unmodified
  double mu0 = scatter_dist->sample();
  p->scatter( mu0 );
}

void  CE_scatter_reaction::sample( particle* p, std::stack< particle >* bank ) {

  double vin = std::sqrt(2.0*p->energy()*1.0e-6*1.6022e-19/p->mass());

  double c = vin*(1/(1+A));
  point vcm ( (c * p->dir().x), (c * p->dir().y), (c * p->dir().z) );

  // sample scattering angle in CoM frame
  double mu0com = scatter_dist->sample();
  p->scatter( mu0com );

  double voutx = vcm.x + vin*(A/(A+1))*p->dir().x;
  double vouty = vcm.y + vin*(A/(A+1))*p->dir().y;
  double voutz = vcm.z + vin*(A/(A+1))*p->dir().z;

  p->setEnergy(0.5*p->mass()*( voutx*voutx + vouty*vouty + voutz*voutz ));
  
}


void  SE_fission_reaction::sample( particle* p, std::stack< particle >* bank ) {
  // create random number of secondaries from multiplicity distributon and
  // push all but one of them into the bank, and set working particle to the last one
  // if no secondaries, kill the particle

  int n = multiplicity_dist->sample();
  if ( n <= 0 ) {
    p->kill();
  }
  else {
    // bank all but last particle (skips if n = 1)
    for ( int i = 0 ; i < (n - 1) ; i++ ) {
      particle q( p->pos(), isotropic->sample(), p->energy() );
      q.recordCell( p->cellPointer() );
      bank->push( q );
    }
    // set working particle to last one
    particle q( p->pos(), isotropic->sample(), p->energy() );
    q.recordCell( p->cellPointer() );
    *p = q;
  }
}

void  CE_fission_reaction::sample( particle* p, std::stack< particle >* bank ) {
  // create random number of secondaries from multiplicity distributon and
  // push all but one of them into the bank, and set working particle to the last one
  // if no secondaries, kill the particle

  int n = multiplicity_dist->sample();
  if ( n <= 0 ) {
    p->kill();
  }
  else {
    // bank all but last particle (skips if n = 1)
    for ( int i = 0 ; i < (n - 1) ; i++ ) {
      particle q( p->pos(), isotropic->sample(), p->energy() );
      q.recordCell( p->cellPointer() );
      bank->push( q );
    }
    // set working particle to last one
    particle q( p->pos(), isotropic->sample(), p->energy() );
    q.recordCell( p->cellPointer() );
    *p = q;
  }
}

