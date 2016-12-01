#ifndef SQUFOF_H
#define SQUFOF_H


#include <gmpxx.h>
#include <list>


bool squfof(mpz_class &N, mpz_class &p, float time_allocated,
        clock_t t);

bool Shanks_factoring(mpz_class N, std::list<mpz_class> &factors,
        int factoring_power, float time_allocated, clock_t t);


#endif
