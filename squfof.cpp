#include <gmpxx.h>
#include <list>
#include <iostream>
#include <cassert>
#include <utility>
#define myTime(t0) ((float)clock()-(t0))/CLOCKS_PER_SEC

//    these are the nmultipliers considered to increase N 
// in case of failure
const int mult[] = {1, 3, 5, 7, 11, 3*5, 17, 19, 3*7, 3*11, 5*7, 5*11, 7*11};



/*          SHANK'S SQUFOF FACTORIZATION whithout decomposition   */

bool squfof(mpz_class &N, mpz_class &p, float time_allocated,
        clock_t t){
        
    mpz_class P0, P1; 
    mpz_class Q0, Q1, Q2, sqrt_Q1;
    mpz_class b0, b1;
    mpz_class R, sqrt_R;
    int k = 1;
    
    while (true){
        R = mult[k]*N;            
        mpz_sqrt(sqrt_R.get_mpz_t(), R.get_mpz_t());
        if(sqrt_R*sqrt_R == R){
            p = sqrt_R;
            N = N/p;
            return false;
        }
        // forward cycling untill Q is found to be square

        P0 = sqrt_R;
        Q0 = 1;
        Q1 = R - (P0*P0);
        //int i = 1;
        while (0 == mpz_perfect_square_p(Q1.get_mpz_t()) /*or 
               i % 2 == 1*/){
            // check time
            if( myTime(t)> time_allocated){
                return true;
            } 
            //compute
            b1 =  (sqrt_R + P0)/Q1;
            P1 = (b1*Q1) - P0;
            Q2 = Q0 + (b1*(P0-P1)); 
            
            //update 
            P0 = P1;
            Q0 = Q1;
            Q1 = Q2;
            //i++;
            //if (Q1 == 0)
            //    return true;
        } 
        // backward steps
        mpz_sqrt(sqrt_Q1.get_mpz_t(), Q1.get_mpz_t());
        b1 =  (sqrt_R - P0)/sqrt_Q1;
        P0 = b1*sqrt_Q1 + P0;
        Q0 = sqrt_Q1;
        Q1 = (R - (P0*P0))/Q0;
        while (true){
            // check time  
            if( myTime(t)> time_allocated){
                return true;
            } 
            
            // compute
            b1 = (sqrt_R + P0)/Q1;
            P1 = b1*Q1 - P0;
            if (P1 == P0)
                break;
            Q2 = Q0 + (b1*(P0-P1));
            
            // stop when P1 = P0
            
            // update
            P0 = P1;
            Q0 = Q1;
            Q1 = Q2;

            //if (Q1 == 0)
            //    return false;
        }
      
        mpz_gcd(p.get_mpz_t(), P1.get_mpz_t(), N.get_mpz_t()); 
        if(p>1 and p<N){
            N=N/p;
            // std::cout << "not cool  \n";
            return false;
        }
        else{
            k++;
            continue;
        }
    }
    return true;
}





bool Shanks_factoring(mpz_class N, std::list<mpz_class> &factors,
        int factoring_power, float time_allocated, clock_t t){

    if(N == 1)
        return false;

    if(1 <= mpz_probab_prime_p(N.get_mpz_t(),20)){
        for(int i=0; i<factoring_power; i++)
            factors.push_back(N);
        return false;
    }

    // find factoring power (i.e N = n^p)
    // we know that N < 2^100
    // so p<100    

    mpz_class temp, n = 0;
    int i = 2, power = 1;
    while(i<100){ 
        if(mpz_root(temp.get_mpz_t(),N.get_mpz_t(),i)){
            power = i;
            n = temp;
        }
        if(temp < 2){
            break;
        }
        i++;
    }
    if(n>1){
        factoring_power *= power;
        N = n;
    }

    // check again if n is prime
    if(1 <= mpz_probab_prime_p(N.get_mpz_t(),20)){
        for(int i=0; i<factoring_power; i++)
            factors.push_back(N);
        return false;
    }


    mpz_class p,q;
    if(squfof(N,p,time_allocated,t))
        return true;

    /*if(1 <= mpz_probab_prime_p(p.get_mpz_t(),20)){ 
        for(int i=0; i<factoring_power; i++)
            factors.push_back(p);
    } else*/
    if (Shanks_factoring(p,factors,factoring_power,time_allocated,t))
        return true;

    if(Shanks_factoring(N,factors,factoring_power,time_allocated,t))
        return true;

    return false;
}


/*    Shanks' SQUare FOrms Factorization with decompostion  */
/*                                                          */
/*              NOT FULLY IMPLEMENTED THOUGH                */
/*                                                          */
/*                                                          */

// Shanks' SQUare FOrms Factorization
//bool squfof(mpz_class &N, mpz_class &p, float time_allocated,
//        clock_t t){
//        
//    mpz_class P0, P1; 
//    mpz_class Q0, Q1, Q2, sqrt_Q1;
//    mpz_class b0, b1;
//    mpz_class R, sqrt_R;
//    mpz_class L, L1, B;
//    mpz_class r;
//    // i-th
//    int i;
//    int k = 1;
//    while (true){
//        
//        // N should not be even
//        assert(N%2 != 0);
//        
//         
//        // forward steps
//        if(N%4 == 1)
//            R = 2*N;            
//        else 
//            R = N;
//        //std::cout << " R=" << R <<"\n";
//        mpz_sqrt(sqrt_R.get_mpz_t(), R.get_mpz_t());
//        L1 = 2*sqrt_R;
//        mpz_sqrt(L.get_mpz_t(),L1.get_mpz_t());
//        L = 2*L;
//        B = 2*L;
//      
//        P0 = sqrt_R;
//        Q0 = 1;
//        Q1 = R - (P0*P0);
//        i = 2;
//        
//         std::cout << " L=" << L << "\n"; 
//         std::cout << " B=" << B << "\n"; 
//         std::cout << " P=S=sqrt_R=" << sqrt_R  << "\n"; 
//         std::cout << " Q=" << Q1 << "\n"; 
//
//
//        
//        std::list< std::pair<mpz_class,mpz_class> > queue;
//        while (i<B){
//                
//            // check time
//            if( myTime(t)> time_allocated){
//                return true;
//            } 
//            //compute
//            b1 =  (sqrt_R + P0)/Q1;
//            P1 = (b1*Q1) - P0;
//            if (Q1 <= L){
//                if (Q1 % 2 == 0){
//                    mpz_class l = Q1/2;
//                    mpz_class t = P0 % l;
//                    std::pair< mpz_class, mpz_class > temp(l,t);
//                    queue.push_back(temp);
//                }else{
//                    mpz_class l = Q1;
//                    mpz_class t = P0 % l;
//                    std::pair<mpz_class, mpz_class> temp(l,t);
//                    queue.push_back(temp);
//                }
//
//            }
//            Q2 = Q0 + (b1*(P0-P1)); 
//            
//            //update 
//            P0 = P1;
//            Q0 = Q1;
//            Q1 = Q2;
//                
//            if(i%2 == 1){
//                i = i+1;
//                continue;
//            }
//            
//            if(0 == mpz_perfect_square_p(Q1.get_mpz_t())){
//                //std::cout << " Q=" << Q1 << "\n"; 
//                i = i+1;
//                continue;
//            }else{
//                mpz_sqrt(r.get_mpz_t(), Q1.get_mpz_t());
//                //std::cout << " r=" << r << "\n"; 
//                std::list< std::pair<mpz_class,mpz_class> >::iterator it;
//                std::pair<mpz_class, mpz_class> saved;
//                for(it=queue.begin(); it!=queue.end(); ++it){
//                    //std::cout << "time \n";
//                    if( r == (*it).first /*and t == (*it).second)*/ 
//                            and ((P1-Q2) % ((*it).first) == 0)){
//                        if(r>1){
//                            //removal = true;
//                            //it_r = queue.begin();
//                            it++;
//                            queue.erase(queue.begin(), it);
//                            i++;
//                            continue;
//                        }
//                        if(r==1){
//                            k = k+1;
//                            return true;
//                        }
//                    }
//                }
//                break;
//            }
//        }
//        if(i>B)
//            return true;
//        // backward steps
//        //mpz_sqrt(sqrt_Q1.get_mpz_t(), Q1.get_mpz_t());
//        Q0 = r;
//        b1 =  (sqrt_R - P0)/r;
//        P0 = r*b1 + P0;
//        Q1 = (R - (P0*P0))/Q0;
//        while (true){
//            // check time  
//            if( myTime(t)> time_allocated){
//                return true;
//            } 
//            
//            // compute
//            b1 = (sqrt_R + P0)/Q1;
//            P1 = b1*Q1 - P0;
//            if(P1 == P0){
//                break;
//            }
//            Q2 = Q0 + (b1*(P0-P1));
//            // update
//            P0 = P1;
//            Q0 = Q1;
//            Q1 = Q2;
//
//            //if (Q1 == 0)
//            //    return false;
//        }
//      
//        if(Q1%2 == 0){
//            Q1 = Q1/2;
//        }
//        p = Q1;
//        N = N/p;
//        return false;
//    }
//    return true;
//}
//



