#ifndef DIXON_H
#define DIXON_H

#include <gmpxx.h> 
#include <list>
#include <vector>

// base factor
//const int pi_B = 4; // 8 (prime numbers < 20) + 1 (special prime -1)
//const mpz_class factor_base[9] = {-1,2,3,5,7,11,13,17,19};


void seive_erastothen(int B, std::vector<int> & factor_base);


bool  dixon_factoring(mpz_class N, std::list<mpz_class> & factors, 
        int B, const std::vector<int> & factor_base, 
        float time_allocated, clock_t t);

bool dixon_recursive(mpz_class N, std::list<mpz_class> & factors, 
        int b, const std::vector<int> & factor_base, 
        float time_allocated, clock_t t);

bool ce_factoring(mpz_class & N, mpz_class & p, int b, 
        const std::vector<int> &  factor_base, 
        float time_allocated, clock_t t);

bool isBsmooth(mpz_class N, mpz_class Q, mpz_class &q, 
        std::list< std::vector<int> > & Q_factors, int b, 
        std::vector<int> & factor_base);

bool gauss_elimination( 
        const std::list< std::vector<int> > & Q_factors_
        , int b, std::list<int> &dependent_rows);
 


#endif
