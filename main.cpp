#include <iostream>  // cin, cout 
#include <ctime>     // clock_t, clock(), CLOCKS_PER_SEC 
#include <gmpxx.h>   // GMP
#include <list>      // list
#include <cstring>   // stdcmp()
#include <cassert>   // assert()
#include <algorithm>
#include <stdlib.h>
#include <cmath> 

/* IMPLEMENTED ALGORITHM  */

#include "squfof.h"
#include "dixon.h"


#define checkTime(t,time_allocated) (((float)clock()-(t))/CLOCKS_PER_SEC > (time_allocated)) 
#define myTime(t0) ((float)clock()-(t0))/CLOCKS_PER_SEC
#define max(a,b) (a>b? a:b)
#define min(a,b) (a<b? a:b)

/*
 *  Authors: Yassir Jedra (jedra@kth.se) and Lucas Rodes Guirao (hello@lucasrodesguirao.com)
 *  
 *  There are several comments throughout the code. They might reffer to different issues,
 *  hence, in order to differentiate its nature we use different symbology:
 */  
    /* This explains the code below and should be used as a guideline to understand the code */
    // This are print outs that are usefull when debugging.
    //// This contain parts of the code that might be commented or uncommented



/* Choose which algorithm to run: 
      1 (Trial divisions), 2 (Fermat), 3 (Rho), 4 (SHank's SQUFOF), 5(DIXON) */
int ALGORITHM = 3; 
/* Randomized versions*/
int RANDOMIZED = 0;

/* Variable needed to generate random numbers */
gmp_randstate_t STATE;

/* Declaration of the algorithm functions */
bool naive_factoring(mpz_class &number, std::list<mpz_class> & factors, 
        float time_allocated, clock_t t);
bool fermat_method(mpz_class & N, mpz_class & p, mpz_class & q, 
        float time_allocated, clock_t t);
bool fermat_factoring(mpz_class N, std::list<mpz_class> & factors, 
        float time_allocated, clock_t t);
bool rho_factoring(mpz_class number, std::list<mpz_class> &factors,
        float time_allocated, clock_t t1, int depth);
/* Algorithms involving randomization */
bool rand_naive_factoring(mpz_class & number, std::list<mpz_class> & factors,
        float time_allocated, clock_t t);

bool rand_fermat_factoring(mpz_class N, std::list<mpz_class> & factors, 
        float time_allocated, clock_t t);
bool fermat_fact(mpz_class & N, mpz_class & p, mpz_class & q, float time_allocated, clock_t t);
 

/* Main */
int main(int argc, char* argv[]){
    srand(time(NULL));
    /* The starting time */
    clock_t t0 = clock();
    /* Read all the input numbers and store them in a list*/
    std::list<mpz_class> numbers;
    mpz_class number;
    while(std::cin >> number)
        numbers.push_back(number); 
    
    /* Initilizations required to generate a random number in 
       POLLARD'S  RHO ALGORITHM */
    long seed;
    gmp_randinit_mt(STATE);
    mpz_class x, y, z, g, limit; g = 0;
    time(&seed);
    gmp_randseed_ui (STATE, seed);
    
    int B_MAX = 100000;
    
    /* construct an array of primes < B_MAX */
    std::vector<int> factor_base;
    if(ALGORITHM == 5)
        seive_erastothen(B_MAX, factor_base);

    /* While time limit not exceeded (i.e time < 15 ) keep working on */
    while (myTime(t0) < 14.9){
        std::list<mpz_class>::iterator it_n;
        uint visited_numbers = 0;

        /* Loop through all the numbers in the input */ 
        for(it_n = numbers.begin() ; it_n != numbers.end(); ++it_n){
            number = *it_n;
            
            /* Update the time remaining, and the time allocated to number */                    
            clock_t t1 = clock();
            float time_left = 14.95 - ((float)t1 - t0)/CLOCKS_PER_SEC; 
            /* Generate a normally distributed variable using Box–Muller transform */
            float std_dev = 0.1;
            float mean = 0.2 ;
            float a = ((double) rand() / (RAND_MAX));
            float b = ((double) rand() / (RAND_MAX));

            float extra_time = 0;
            ////float extra_time = std_dev*(sqrt(-2*log(a))*cos(2*M_PI*b))+mean; //DECOMENT THIS FOR RANDOM TIME ALLOCATION
            // std::cout << extra_time << std::endl;
            // std::cout << std::endl;
            /* Allocate the time according to the remaining time and 
               some randomness given by the normally distributed variable 
               extra_time */
            float time_allocated = time_left/(numbers.size()-visited_numbers) + extra_time;

            visited_numbers++;
            
            std::list<mpz_class> factors;            
            bool fail = false;
            
            // std::cout << "Number: " << number << std::endl;


            ////naive_factoring(number, factors, time_allocated, t1); 
 

            /* Make the number odd */
            while (number%2==0){
                factors.push_back(2);
                number = number/2; 
                //std::cout << "loop" << std::endl;
                //std::cout << "2" << ", number = " << number << std::endl;
            }

            /* Prevent the number to be divisible by 25, it caused some issues in Pollard's*/
            while (number%25==0){
                factors.push_back(5);
                factors.push_back(5);
                number = number/25; 
                // std::cout << "loop" << std::endl;
                // std::cout << "2" << ", number = " << number << std::endl;
            }
            mpz_class number_sum = 1;

            /* If the input number is 1, output 1. If the remaining number (after dividing
               by 2 and 25 is 1, then print the stored factors */
            if (number <= 1){	
                if (*it_n <= 1)
                    factors.push_back(number);
            }

           /* Create a list containing the factors of number */
            else{
                /* Factorize according to the chosen ALGORITHM */
                switch(ALGORITHM){
                    /* NAIVE ALGORITHM */
                    case 1 :    if(!RANDOMIZED){
                                    fail = naive_factoring(number, factors, 
                                            time_allocated, t1); 
                                    break;
                                } else{
                                    fail = rand_naive_factoring(number, factors,
                                            time_allocated, t1); 
                                    break;
                                }

                    /* FERMAT ALGORITHM */
                    case 2 :    if(!RANDOMIZED){
                                    fail = fermat_factoring(number, factors,
                                            time_allocated, t1); 
                                    break;
                                } else{
                                    fail = rand_fermat_factoring(number, 
                                            factors, time_allocated, t1); 
                                    break;
                                }

                    case 3 :    fail = rho_factoring(number, factors, 
                                        time_allocated, t1, 0); 
                                break;
                    /* SHANKS' ALGORITHM */
                    case 4 :    fail = Shanks_factoring(number, factors, 1, 
                                        time_allocated, t1); break;

                    /* DIXON'S ALGORITHM */            
                    case 5 :    fail = dixon_factoring(number, factors,
                                        B_MAX, factor_base, 
                                        time_allocated, t1);
                                break;
                }
            }
            // std::cout << " " << std::endl;

            /* Output result for current number: (i) either "fail" or (ii) the list of found factors*/
            /* (i) */
            if (fail)
                std::cout << "fail \n";
            // (ii) */
            else {
                std::list<mpz_class>::iterator it;
                number_sum = 1;
                for(it = factors.begin(); it != factors.end(); ++it)
                    number_sum *= *it;
                if(number_sum == *it_n){
                    for(it = factors.begin(); it != factors.end(); ++it)
                        std::cout << *it << "\n";
                }else
                    /* Found factors are not the correct ones! We checked in Kattis and the algorithm
                     never enters this condition. However we leave it in case of future improvements
                     */
                    std::cout << "fail \n";
            }
            std::cout << std::endl;
        }
        
        /* Write something after the execution of the file, i.e. "./file X" to display the runing time */
        if(argc == 2 && strcmp(argv[0],"1"))
            std::cout << "time :" << myTime(t0) << "\n";
        /* if visited numbers != numbers.size() raise an error */
        assert(visited_numbers == numbers.size());
        return 0;
    }            
}


/* THHE FOLLOWING SECTION CONTAINS THE FUNCTIONS IMPLEMENTING THE ALGORITHMS USED
   If the reader is not confident with the theory, we encourage her/him to have a
   look at our report
*/


/* NAIVE FACTORIZATION ALGORITHM (brute force) */
bool naive_factoring(mpz_class & number, std::list<mpz_class> & factors,
        float time_allocated, clock_t t){

    if(1 <= mpz_probab_prime_p(number.get_mpz_t(),50)){
        factors.push_back(number);
    } else {
        mpz_class d(3); /* denominator*/  
        mpz_class q(number); /* quotient */
        mpz_class r;         /* remainder */
        mpz_class root;      /* square root */  

        /* mpz_root <=> root = square_root(number)*/
        mpz_root(root.get_mpz_t(),number.get_mpz_t(),2);            

        /* Find factors */

        ////while (q != 1 && d<102){
        while (q != 1){
            /* Check if we have exceeded the time allocated */
            if (myTime(t)>time_allocated) 
                return true;
            
            /* case1: 
                   d > square root of the current quotient q */
            if(d > root){
                factors.push_back(q);
                break;
            }
            
            /* mpz_mod <=> r = q mod d */
            mpz_mod(r.get_mpz_t(),q.get_mpz_t(),d.get_mpz_t());


            /* case2: 
                   d <= square root of the current quotient q*/
            
            /* check if d is divisor*/ 
            if (r == 0){
                factors.push_back(d);
                q = q/d;
            }else{

                ////mpz_nextprime(d.get_mpz_t(), d.get_mpz_t());
                d = d+1; 
                mpz_root(root.get_mpz_t(),q.get_mpz_t(),2);            
            }
        }

        ////if (1<=mpz_probab_prime_p(q.get_mpz_t(),50))
        ////    number = 1;
        ////else
        ////    number = q;
    } 
    /* False if no time limit exeeded */
    return false;
}

/* RANDOMIZED NAIVE FACTORIZATION ALGORITHM */
bool rand_naive_factoring(mpz_class & number, std::list<mpz_class> & factors,
        float time_allocated, clock_t t){

    if(1 <= mpz_probab_prime_p(number.get_mpz_t(),50)){
        factors.push_back(number);
    } else {
        mpz_class root;      /* square root */  
        /* mpz_root <=> root = square_root(number)*/
        mpz_root(root.get_mpz_t(),number.get_mpz_t(),2);            

        /* Randomly obtain a prime number within {3,...,sqrt(N)} */
        mpz_class d;
        mpz_urandomm(d.get_mpz_t(), STATE, root.get_mpz_t()); 
        d = d+2;
        // std::cout << "selected: " << d << std::endl;
        while(mpz_probab_prime_p(d.get_mpz_t(), 50)<1){
            d= (d+1)%root+2;
            // std::cout << d << std::endl;/* Check if we have exceeded the time allocated */
            if (myTime(t)>time_allocated) 
                return true;
        }
        /* Alternatively, one might prefer to initialize it without randomness */
        // std::cout << "exited while loop, d = " << d << std::endl;
        mpz_class q(number); /* quotient */
        mpz_class r;         /* remainder */


        /* Find factors */
        while (q != 1){

            // std::cout << d << std::endl;

            /* Check if we have exceeded the time allocated */
            if (myTime(t)>time_allocated) 
                return true;
            
            /* case1: 
                   d > square root of the current quotient q */
            if(d > root){
                factors.push_back(q);
                break;
            }
            
            /* mpz_mod <=> r = q mod d */
            mpz_mod(r.get_mpz_t(),q.get_mpz_t(),d.get_mpz_t());


            /* case2: 
                   d <= square root of the current quotient q*/
            
            /* check if d is divisor*/ 
            if (r == 0){
                factors.push_back(d);
                q = q/d;
            }else{
                mpz_nextprime(d.get_mpz_t(), d.get_mpz_t()); 

                if (d > root)
                    d = 3;
                mpz_root(root.get_mpz_t(),q.get_mpz_t(),2);            
            }
        }
    }
    /* false if no time limit exeeded */
    return false;
}


/* FERMAT's FACTORIZATION ALGORITHM */
bool fermat_method(mpz_class & N, mpz_class & p, mpz_class & q, 
        float time_allocated, clock_t t){

    /* N must be odd */
    assert(0 != mpz_odd_p(N.get_mpz_t()));
    
    /* Compute the root of N */
    mpz_class root_N,rem;
    mpz_sqrtrem(root_N.get_mpz_t(),rem.get_mpz_t(),N.get_mpz_t());
    
    /* x = ceiling(root of N) */
    mpz_class x = (rem == 0) ? root_N : root_N+1;
    /* x*x = Y mod N */
    mpz_class Y = x*x - N;
    
    /* y = sqrt(Y) */
    mpz_class y;
    mpz_sqrt(y.get_mpz_t(), Y.get_mpz_t());
       
    /* While Y is not a perfect square */
    while (0 == mpz_perfect_square_p(Y.get_mpz_t()) /*&& x < N-1 && (x-y)*(x+y) != N*/) {
        
        /* Check if we have exceeded the time allocated */
        if(myTime(t)>time_allocated)
            return true; 
        x = x+1;
        Y = x*x - N;
        
        /* y = sqrt(Y) */
        mpz_sqrt(y.get_mpz_t(), Y.get_mpz_t());
    }
     
    /* Compute the factors */
    p = x+y;
    q = x-y; 
    
    /* N must = p*q */
    assert (N == p*q);
    
    /* p must not be prime */
    assert (N != p);
    
    /* Update N */
    N = N/(p*q);
    
    /* no time limit exeeded */
    return false;
}

bool fermat_factoring(mpz_class N, std::list<mpz_class> & factors, 
        float time_allocated, clock_t t){
    
    /* Check if we have exceeded the time allocated */
    if (myTime(t) > time_allocated)
        return true;

    mpz_class p, q;

  
    /* case 1: N is prime */
    if (1 <= mpz_probab_prime_p(N.get_mpz_t(),50)){
        factors.push_back(N);
        N = 1;
        return false;
    }


    /* case 2: N is odd
            => use fermat's method to find p and q */
    if (fermat_method(N,p,q,time_allocated,t))
        return true;

    /* N must be 1 otherwise p,q is wrong */ 
    assert(N == 1);


    /* Case 2.1: p is not prime */
    if ( fermat_factoring(p,factors,time_allocated,t))
        return true;

    /* Case 2.2: q is not prime and q!=p
             => use fermat's method on q */
    if (fermat_factoring(q,factors,time_allocated,t))
        return true;
    
    /* No time limit exeeded */
    return false;
}


/* POLLARD's RHO ALGORIHTM FACTORIZATION */
bool rho_factoring(mpz_class number, std::list<mpz_class> &factors,
        float time_allocated, clock_t t1, int depth){
    /* Initialize main parameters: x and y = numbers in the generated sequence; z = x - y ; limit = limit  of 
       the sequence size; g = GCD computed in the algorithm */
    mpz_class limit, x, y, z, g;
    
    bool factorized = false;
    int a ; 
    // std::cout << "entered level " << depth << std::endl;
    
    while(!factorized){
        // std::cout << number << std::endl;        
        /* Primality test */
        int is_prime = mpz_probab_prime_p(number.get_mpz_t(), 50);
               
        if(is_prime==0){
            // std::cout << "not prime" << std::endl;
            /* Obtain the sequence size n^0.25 */
            mpz_root(limit.get_mpz_t(), number.get_mpz_t(), 4);
            limit+=3;
            /* Define the seeds */
            mpz_urandomm(x.get_mpz_t(), STATE, number.get_mpz_t()); 
            a = 1;
            x = x + a;
            y = x;
            /* Factor search */ 
            // std::cout << "looking for factor" << std::endl;
            for(int i = 0; i < limit; i++){
                if (myTime(t1) > time_allocated){ 
                    return true;
                }
                /* Update x_i */
                x = (x*x + a) % number;
                /* Update y_i */
                y = (((y*y + a)%number)*((y*y +a)%number) + a) % number;
                /* Compute the difference */
                z = x - y;
                /* Obtain GCD(z, number) */ 
                mpz_gcd(g.get_mpz_t(), z.get_mpz_t(), number.get_mpz_t()); 
      
                // std::cout << x << "," << y << "factor = " << g << std::endl;
                if (g != 1 and g != number){
                    /* Rho factoring on g */
                    rho_factoring(g, factors, time_allocated, t1, depth+1);
                    number = number/g;
                    break;
                }else if (g == number)
                    a = a+1;
            }
            //// count++;
            //// if (count > lim_iter){
            ////    return true;
            ////} 

        }
        /* If number is prime, it is a factor itself and hence we are done. 
         In addition, if number is 1, we break the loop */
        else if (is_prime != 0 or number == 1){
            // std::cout << "prime" << std::endl;
            factors.push_back(number);
            break;
        }

    }
   //  std::cout << "skipping level " << depth <<  std::endl;
   return false;
}    

 

/* RANDOMIZED FERMAT FACTORIZATION ALGORITHM*/
bool fermat_fact(mpz_class & N, mpz_class & p, mpz_class & q, float time_allocated, clock_t t){
 
    mpz_class root, rem;
    mpz_class x,y,a,b,X,Y;
    
    mpz_sqrtrem(root.get_mpz_t(), rem.get_mpz_t(), N.get_mpz_t());

   if(rem == 0){        
        p = root;
        N = N/p;
        return false;
    }
    
    x = root;
    int r = 1;
    int k = 0;
    //mpz_class a,b,ra,rb;
    bool not_yet_found = true;
    while (not_yet_found){        
        
        if( myTime(t)> time_allocated){
            return true;
        }       
        mpz_class R; 
        if(k < 100){
            x = root + k;
            k++;
        }else {
            r = rand() % 1000;
            R = r*N; 
            k = -100;
        }
        
         //x = x+1;
        
        mpz_urandomm(x.get_mpz_t(), STATE, N.get_mpz_t());
        if (x<= 1)
            continue;
        
        mpz_gcd(p.get_mpz_t(), x.get_mpz_t(), N.get_mpz_t());
        if(p>1){
            N = N/p;
            return false;
        }
        
        // X = x²
        X = x*x;        
        
        // x² = y² mod N
        mpz_mod(Y.get_mpz_t(),X.get_mpz_t(),N.get_mpz_t());
        
        // if Y = y²
        if (0 != mpz_perfect_square_p(Y.get_mpz_t())){

            mpz_sqrt(y.get_mpz_t(), Y.get_mpz_t());
            // check if they won't lead to non trivial factors
            //mpz_mod(ra.get_mpz_t(),a.get_mpz_t(),N.get_mpz_t());
            //mpz_mod(rb.get_mpz_t(),b.get_mpz_t(),N.get_mpz_t());
            
            if ((x+y) % N != 0 and (x-y) % N != 0) {
                a = x+y;
                b = x-y;
                mpz_gcd(p.get_mpz_t(),a.get_mpz_t(),N.get_mpz_t());
                N = N/p;
                return false;
                not_yet_found = false;
            }
        } 
    }
    return false;
}

bool rand_fermat_factoring(mpz_class N, std::list<mpz_class> &factors,
        float time_allocated, clock_t t){

    if(N == 1)
        return false;

    if(1 <= mpz_probab_prime_p(N.get_mpz_t(),20)){
        factors.push_back(N);
        return false;
    }

    mpz_class p,q;
    if(fermat_fact(N,p,q,time_allocated,t))
        return true;

    if(1 <= mpz_probab_prime_p(p.get_mpz_t(),20)){
        factors.push_back(p);
        while(N%p == 0){
            N = N/p;
            factors.push_back(p);
        }
    } else if (rand_fermat_factoring(p,factors,time_allocated,t))
        return true;

    if(rand_fermat_factoring(N,factors,time_allocated,t))
        return true;

    return false;
}


