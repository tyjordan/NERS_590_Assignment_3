#include <iostream>
#include <string>
#include <cmath>

#include "Caffeine.h"


double constant_dependence::sample(double E) { return 1.0; }

double inverse_sqrt_dependence::sample(double E) { return 1.0/std::sqrt(E); }

