#ifndef _IMPORTANCE_SPLITTING_HEADER_
#define _IMPORTANCE_SPLITTING_HEADER_

#include <memory>
#include <stack>

#include "Particle.h"
#include "Cell.h"
#include "Random.h"

void ImportanceSplitting
(bool split_roulette, double prevImportance, std::shared_ptr <cell> currentCell,
 particle* p, std::stack<particle>* bank);

#endif