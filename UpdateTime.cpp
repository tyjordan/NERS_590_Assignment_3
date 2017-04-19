#include "UpdateTime.h"

void UpdateTime(bool time_tracking, double time_cutoff, double* transDist, particle* p)
{
	if(time_tracking) {
		p->updateTime(*transDist);

		if( p->time() > time_cutoff ) {
			p->setTime(time_cutoff);
			*transDist -= ( p->time() - time_cutoff ) * std::sqrt( 2.0 * p->energy() * 1.6022e-13 / p->mass() ) * 100.0;
			p->kill();
		}
	}
	return;
}