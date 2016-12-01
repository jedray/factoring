#include <list>     // lists
#include <vector>   // vector
#include <iostream> //
#include <gmpxx.h>  // mpz_class and large numbers
#include <assert.h> 

#define myTime(t0) ((float)clock()-(t0))/CLOCKS_PER_SEC



void seive_erastothen(int B, std::vector<int> & factor_base){
    bool isPrime[B];
    
    for(int i=0;i<B;i++)
        isPrime[i] = true;

    isPrime[0] = false;
    isPrime[1] = false;
    for(int p=2; p*p < B; p++){
        if(isPrime[p]){
            int k = 2;
            while (k*p < B){
                isPrime[k*p] = false;
                k++;
            }
        }
    }


    int b = 0;
    for(int i=0; i<B; i++){
        if(isPrime[i])
            b++;
    }

    factor_base.resize(b+1);
    factor_base[0] = -1;
    int j = 1;
    for(int i=2; i<B; i++){
        if(isPrime[i]){
            factor_base[j] = i;
            j++;
        }
    }
}

bool gauss_elimination( const std::list< std::vector<int> > & Q_factors_
        , int b, std::list<int> &dependent_rows){
    
    int m = b + 1; // m columns
    int n = Q_factors_.size(); // n rows
    

    bool mark[n];     // row marker
    for(int i=0; i< n; i++)
        mark[i] = false;
    std::list< std::vector<int> > Q_factors = Q_factors_;
    std::list< std::vector<int> >::iterator i,l;
    int d;
    // triangulize the matrix 
    for(int j=0; j<m; j++){
        i = Q_factors.begin();
        d = 0;
        while (i != Q_factors.end()){
            if((*i)[j] == 1)
                break;
            ++i;
            ++d;
        }
        if (i != Q_factors.end()){
            mark[d] = true;
            for(int k=0; k<m; k++){
                if(k!=j and (*i)[k] == 1){
                    for(l=Q_factors.begin(); l!=Q_factors.end(); ++l){
                        (*l)[k] = ((*l)[k] + (*l)[j])%2;
                    }
                }

            }
        }
    }

//    std::cout << " after elimination \n";
//    for(int i = 0; i<n; i++){
//        for(int j = 0; j<m; j++){
//            std::cout << matrix[i][j];
//            if (j != m-1)
//                std::cout << " , ";
//            else 
//                std::cout << std::endl;
//        }
//    }




    // find dependent rows
    
    // find unmarked row
    i = Q_factors.begin();
    d = 0;
    while (d<n && mark[d]){
        ++i;
        ++d;
    }
    if(d>=n)
        return false;
    dependent_rows.push_back(d);

    // find linear combination for unmarked row
    for (int j = 0; j<m; j++){
        if((*i)[j] == 1){
            d = 0;
            for(l = Q_factors.begin(); l != Q_factors.end(); ++l){
                if(l != i and (*l)[j] == 1){
                    // save row k
                    dependent_rows.push_back(d);
                    break;
                }
                d++;
            }
        }
    }
    return true;
}




// checks whether Q² mod N is factorizable in the factor base
bool isBsmooth(mpz_class N, mpz_class Q, mpz_class &q, 
        std::list< std::vector<int> > & Q_factors, int b, 
        const std::vector<int> & factor_base){
    std::cout << " b=" << b;
    std::vector<int> f(b+1); // b+1 is the size of the base factor
    std::cout << "why ? \n";    
    // compute Q² mod N 
    // take smallest absolute value
    q = Q*Q % N;
    mpz_class q_ = N - (Q*Q % N);
    if (q_< q){
        q = q_;
        f[0] = 1;
        
    }
    int sum = f[0];

    // find factors in the factor base
    for(int i=1; i<b+1; i++){
        // check if p_i in base factor divides 
        //std::cout << "here \n";
        while(q % (factor_base[i]) == 0){
           q /= (factor_base[i]);
           f[i] = (f[i] + 1) % 2;
           std::cout << factor_base[i] <<", ";
        }  
        sum += f[i];
        //std::cout << "or here \n";
    } 

    //std::cout << "q=" << q << "\n";
    if (q != 1 or sum == 0)
        return false;
    
    std::cout << " not even once \n";
    // add vector to list
    Q_factors.push_back(f);
    return true;
}


// main factoring function
bool ce_factoring(mpz_class & N, mpz_class & p,  
        int b, const std::vector<int> & factor_base, float time_allocated, clock_t t){
    
    // init search 
    mpz_class Q,q;
    mpz_class prod_Q, prod_q;
    mpz_sqrt(Q.get_mpz_t(), N.get_mpz_t());
    
    std::list< std::vector<int> > Q_factors; // list of vectors<b+1>
    std::list< std::vector<int> >::iterator l;
    std::list< mpz_class > Qs,qs;               // list of saved Qis
    std::list< mpz_class >::iterator i,j;
    
    gmp_randstate_t STATE_;
    gmp_randinit_mt(STATE_);
    
    
    mpz_class a_(1), b_(0), root_a_N, a_N(N);
    mpz_sqrt (root_a_N.get_mpz_t(),a_N.get_mpz_t());
    std::cout << "[3.1.1] Entering main loop ... \n";
    while (true){
         
        if(myTime(t)>time_allocated)
            return true;
        
        // generate random Q in the range 0..N-1
        // mpz_urandomm(Q.get_mpz_t(), STATE_, N.get_mpz_t());
        
        
        if(b_ < 15){
            Q = root_a_N + b_;
            b_++;
        }else{
            a_++;
            a_N = a_*N;
            mpz_sqrt (root_a_N.get_mpz_t(),a_N.get_mpz_t());
            b = -15;
        }
        
        if(Q == 0)
            continue;

        // if Q is B-smooth than add it to our lists
                    
        std::cout << "[3.1.1.1] is it B-smooth ? \n";
        if( isBsmooth(N, Q, q, Q_factors, b, factor_base) ){
            Qs.push_back(Q);
            Qs.push_back(q);
        }

        // if picked B-smooth numbers are enough then proceed 

        if(Q_factors.size()> 2){
            std::list<int> dependent_rows;
            std::list<int>::iterator k;
            int d;
            std::cout << "starting gauss elimination .... \n";
            // if dependent rows exist then proceed
            if (gauss_elimination(Q_factors,b,dependent_rows)){
                std::cout << "gauss elimination done \n ";
                prod_Q = 1;
                prod_q = 1;
                k = dependent_rows.begin();
                i = Qs.begin();
                j = qs.begin();
                d = 0;
                // generate the perfect squares x² and y²
                while(k != dependent_rows.end()){
                    if(d == (*k)){
                        prod_Q *= (*i)*(*i);
                        prod_q *= (*j);
                        ++k;
                    }
                    ++d;
                    ++i;
                    ++j;
                }
                mpz_class x,y;
                mpz_sqrt(x.get_mpz_t(), prod_Q.get_mpz_t());
                mpz_sqrt(y.get_mpz_t(), prod_q.get_mpz_t());
                
                // compute the gcd(x+y,N) or gcd(x-y)
                mpz_class candidate;
                candidate = x+y;
                mpz_gcd(p.get_mpz_t(),candidate.get_mpz_t(), N.get_mpz_t());
                // if success return factor gcd(x+y,N)
                if(p>1 and p<N){
                    std::cout << "successefuly factorized \n";
                    N = N/p;
                    return false;
                }
                candidate = x-y;
                mpz_gcd(p.get_mpz_t(),candidate.get_mpz_t(), N.get_mpz_t());
                // if success return factor gcd(x-y,N)
                if(p>1 and q<N){
                    N = N/p;
                    return false;
                }

                // if failure remove first Q in dependent row
                i = Qs.begin();
                j = qs.begin();
                k = dependent_rows.begin();
                l = Q_factors.begin();
                d = 0;
                while(d<(*k)){
                    ++i;
                    ++j;
                    ++l;
                    ++d;
                }
                Qs.erase(i);
                qs.erase(j);
                Q_factors.erase(l);
            } 
        }
    }
    return true;
}



bool dixon_recursive(mpz_class N, std::list<mpz_class> & factors, 
        int b, const std::vector<int> & factor_base, float time_allocated, clock_t t){
    
    if(myTime(t)> time_allocated)
        return true;

    if(N==1)
        return false;

    if(1 <= mpz_probab_prime_p(N.get_mpz_t(),20)){
        factors.push_back(N);
        return false;
    }
    std::cout << "[3.1] Entering ... \n";
    mpz_class p;
    if(ce_factoring(N,p,b,factor_base,time_allocated,t))
        return true;
    std::cout << "[3.1] done. \n";
    mpz_class b_;    
    mpz_root(b_.get_mpz_t(),p.get_mpz_t(),12);
    b = mpz_get_si(b_.get_mpz_t());
    std::cout << "[3.2] Entering recursion over p ... \n";
    if(dixon_recursive(p,factors,b,factor_base,time_allocated,t))
        return true;
    std::cout << "[3.2] Done. \n";
    mpz_root(b_.get_mpz_t(),N.get_mpz_t(),12);
    b = mpz_get_si(b_.get_mpz_t()); 
    std::cout << "[3.3] Entering recursion ever q ... \n";
    if(dixon_recursive(N,factors,b,factor_base,time_allocated,t))
        return true;
    std::cout << "[3.3] Done. \n";

    return false;



}

bool  dixon_factoring(mpz_class N, std::list<mpz_class> & factors, 
        int B, const std::vector<int> & factor_base, float time_allocated, clock_t t){

    

    std::cout << "[DIXON]   Running  .... \n"; 
    // we choose pi(B) to be the 12th rooth of n such that the 
    // gauss elimnation wo'nt take more than n^4
    mpz_class b;
    mpz_root(b.get_mpz_t(),N.get_mpz_t(),4);

    //std::cout << "Dixon Factoring START ..... \n";
    // step 1:
    // first check if the primes in the factor base are factors of N
    // and if find size of the factor base 
    int i = 1;

    std::cout << "[1]   trial devisions .... \n";
    while( factor_base[i] < B){
        while( N % factor_base[i] == 0){
            N = N/factor_base[i];
            factors.push_back(factor_base[i]);
        }
        i++;
        if(i>b)
            break;
    } 
    
    //return false;
    // step 2:
    // apply dixon's factoring
    

    std::cout << "[2]   recursivity starts .... \n";
    return(dixon_recursive(N,factors,mpz_get_si(b.get_mpz_t()), factor_base, 
                time_allocated,t));
}















