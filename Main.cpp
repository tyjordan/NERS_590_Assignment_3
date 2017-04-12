#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <cmath>
#include <vector>
#include <iomanip>
#include <cfloat>
#include <limits>
#include <memory>
#include <cassert>
#include <time.h>
#include <stack>

#include "pugixml.hpp"
#include "Distribution.h"
#include "Particle.h"
#include "Surface.h"
#include "Cell.h"
#include "Estimator.h"
#include "Source.h"
#include "Reaction.h"
#include "Material.h"
#include "Nuclide.h"
#include "Point.h"
#include "Input_XML.h"

// Find what cell a particle is in
void FindCurrentCell ( particle* p , std::vector < std::shared_ptr < cell > >* vec ) {
    std::shared_ptr< cell > cel = nullptr;
    for ( int v = 0; v < vec->size(); v++ ) { // std::shared_ptr<cell> v : vec
        if ( vec->at(v)->testPoint( p->pos() ) ) { // v->testPoint(p->pos())
            cel = vec->at(v);
            break;
        } 
    }
    assert( cel ); // terminates program if cell is null
    p->recordCell( cel );
} // end find cell function




int main() 
{
    //Initialize problem variables
    double err = 1000.0 * std::numeric_limits < double > ::epsilon();

    unsigned long long NSamples;
    bool continuous_eng;
    bool split_roulette;

    // distributuions
    std::vector< std::shared_ptr< distribution<double> > >  double_distributions;
    std::vector< std::shared_ptr< distribution<int>    > >  int_distributions;
    std::vector< std::shared_ptr< distribution<point>  > >  point_distributions;
	std::vector< std::shared_ptr< caffeine > > eng_dependences;

    std::vector< std::shared_ptr<nuclide> > nuclides;
    std::vector< std::shared_ptr<material> > materials;
    std::vector< std::shared_ptr< surface > > surfaces;
    std::vector< std::shared_ptr< cell > > cells;
    std::vector< std::shared_ptr< estimator > > estimators;

    std::shared_ptr< source > src;

	//Read XML input
    Input_Problem_Data( &NSamples , &continuous_eng , &split_roulette , &double_distributions ,
                        &int_distributions , &point_distributions , &eng_dependences , &nuclides ,
                        &materials , &surfaces , &cells , &estimators , &src);

  
	//Transport loop:

    clock_t init, final;
    init=clock();
    int track_count = 0;	//Count number of tracks for Figure of Merit

    // Begin loop over histories
    for ( int n = 1; n <= NSamples; n++ ) {

        std::stack < particle > bank = src->sample();

        while ( ! bank.empty() ) {
            particle p = bank.top() ; bank.pop();
            // find cell particle is in
            FindCurrentCell ( &p , &cells );
            std::shared_ptr < cell > currentCell = p.cellPointer();

            while ( p.alive() ) {
                // find distance to nearest boundary
                std::pair < std::shared_ptr < surface > , double > rayIntersect = 
                    currentCell->surfaceIntersect( p.getRay() );
                double distToBound  = rayIntersect.second;

                // find distance to collision
                double distToCollision = -std::log( Urand() )/(currentCell->macro_xs());
                double transDist = std::fmin( distToCollision , distToBound );
                // move particle to new location
                currentCell->moveParticle( &p , transDist );
				track_count++;
				currentCell->scoreEstimators( &p, transDist );

                // determine if boundary interaction
                if ( transDist == distToBound ) { // boundary interaction
                    // advance to surface and apply crossSurface method from Surface class
					double prevImportance = currentCell->getImportance();
                    rayIntersect.first->crossSurface( &p, transDist );
                    FindCurrentCell ( &p , &cells );
                    currentCell = p.cellPointer();

					//Particle splitting and rouletting based on cell importances
					if( split_roulette ) {

						double r = ( currentCell->getImportance() ) / prevImportance;

						if( r < 1 ) { //roulette
							if( Urand() < r ) {
							 	p.adjustWeight( 1 / r );
							}
							else {
								p.kill();
							}
						}
						else if(r > 1) { //split
							int split_num = std::floor( r + Urand() );
							for(int i = 0; i < split_num; i++)
							{
								bank.emplace( p.pos(), p.dir() );
								bank.top().adjustWeight( p.wgt() / split_num );
							}
							p.kill();
						}
					}
					else if( currentCell->getImportance() == 0 ) { //allows for problem truncation when
							//splitting and rouletting is disabled
						p.kill();
					}
                }
                else { // collision
                    // sample collision
                    currentCell->sampleCollision( &p , &bank );
                }

            } // end while particle alive loop

        } // end while bank !empty loop
        for (auto e: estimators) { e->endHistory(); } //resolve estimators and the end of each history
 
		//simulation progress output
		if( ( n % (int) std::floor( NSamples * 0.1 ) ) == 0 ) {
			std::cout << "	" << n << " histories completed" << std::endl;
		}


    } // end for loop over number of histories
	std::cout << std::endl;
	
    for (auto e: estimators) { e->report( track_count ); }

    final=clock()-init;
    std::cout << "Run time: " << (double)final / ((double)CLOCKS_PER_SEC) << " seconds" << std::endl;

    return 0;
}