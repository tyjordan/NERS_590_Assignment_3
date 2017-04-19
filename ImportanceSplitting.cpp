#include "ImportanceSplitting.h"

void ImportanceSplitting
(bool split_roulette, double prevImportance, std::shared_ptr <cell> currentCell,
 particle* p, std::stack<particle>* bank)
{
	if( split_roulette ) {

		double r = ( currentCell->getImportance() ) / prevImportance;

		if( r < 1 ) { //roulette
			if( Urand() < r ) {
				p->adjustWeight( 1 / r );
			}
			else {
				p->kill();
			}
		}
		else if(r > 1) { //split
			int split_num = std::floor( r + Urand() );
			for(int i = 0; i < split_num; i++) {
				bank->emplace( p->pos(), p->dir(), p->energy() );
				bank->top().adjustWeight( p->wgt() / split_num );
			}
			p->kill();
		}
	}
	else if( currentCell->getImportance() == 0 ) { //allows for problem truncation when splitting and rouletting is disabled
		p->kill();
	}
	return;
}