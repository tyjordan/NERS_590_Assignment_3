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

  double mu0 = scatter_dist->sample();
  p->scatter( mu0 );
  double energy_Cm = p->energy()*std::pow((A/(A+1.0)),2.0);
  double energy_labframe = energy_Cm + ((2.0*mu0*(A+1.0)*std::sqrt(p->energy()*energy_Cm))+p->energy())/std::pow(A+1.0,2.0);
  // sample scattering angle in CoM frame
  double mu0com = mu0*std::sqrt(energy_Cm/energy_labframe)+std::sqrt(p->energy()/energy_labframe)/(A+1);;
  p->scatter( mu0com );
  p->setEnergy(energy_labframe);
  
  return;
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

