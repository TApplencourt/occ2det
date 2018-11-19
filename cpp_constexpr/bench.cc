#include <benchmark/benchmark.h>

#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>


unsigned long n_combinaison(unsigned n, unsigned k) {
    //Try to avoid overflow by interleaving * and /
    //https://blog.plover.com/math/choose.html
    unsigned long r = 1;
    unsigned d;
    if (k > n) return 0;
    for (d=1; d <= k; d++, n--) {
            r = (r*n) / d;
     }
     return r;
}

void gen_permutation0(unsigned const n_orbital, unsigned const n_int, unsigned const size_orbital_bucket, 
                    unsigned const * const occ, unsigned const n_alpha, 
                    unsigned long * const n_det, unsigned **l_det_alpha, unsigned **l_det_beta)
{
    // Will store the value 
    unsigned * const det_patern =  (unsigned*) malloc ( n_int * sizeof(unsigned));
    unsigned * const occ_if = (unsigned*) malloc ( n_orbital * sizeof(unsigned));
 
    unsigned i_int;
    for (i_int=0; i_int < n_int; i_int++){
        det_patern[i_int] = 0;
    }

    unsigned n_alpha_f = n_alpha;
    // Number of single occupied alpha orbital
    unsigned n_single_orbital = 0;

    unsigned idx;
    for (idx = 0 ; idx < n_orbital; idx++){
        switch(occ[idx]) {
            case 1 :
            {
                occ_if[n_single_orbital] = idx;
                n_single_orbital += 1;
                break;
            }
            case 2 :
            {
                unsigned const det_i = idx / size_orbital_bucket;
                unsigned const det_p = idx % size_orbital_bucket;

                det_patern[det_i] |= 1 << det_p; // set bit to one
                n_alpha_f -= 1;
                break;
            } 
        } 
   }
    assert (n_single_orbital <= (sizeof(unsigned)*8));

    *n_det = n_combinaison(n_single_orbital, n_alpha_f);
  
    *l_det_beta =  (unsigned*) malloc ( (*n_det) * n_int * sizeof(unsigned));
    *l_det_alpha =  (unsigned*) malloc ( (*n_det) * n_int * sizeof(unsigned));

    // p represent the permutation of n_alpha_f orbital in the single occupated orbitals
    // Initialize it with the first valid occupations
    // who correspond to the bitmask with `n_alpha_f` right bits set
    // aka p = \sum_{n=0}^m 2^n =  2^{m+1} - 1 == (1 << m) - 1. 
    unsigned p = (1 << n_alpha_f ) - 1;

    unsigned long idx_d;
    unsigned i;

    // Set the value of the determinant
    for (idx_d= 0 ; idx_d < (*n_det); idx_d++){

        unsigned *alpha = *l_det_alpha+idx_d*n_int;
        unsigned *beta =  *l_det_beta+idx_d*n_int; 

        for (i_int=0; i_int < n_int; i_int++){
            beta[i_int] = det_patern[i_int];
            alpha[i_int] = det_patern[i_int];
        }
        
        for ( i = 0 ; i < n_single_orbital; i++){
            bool const o = ( p >> i ) & 1; //index i of permutation p
            unsigned * const det = (o==1) ? alpha: beta; //Determinant are zero-ed .Modify only if need to set the bit to 1

            idx = occ_if[i];
            unsigned const det_i = idx / size_orbital_bucket;
            unsigned const det_p = idx % size_orbital_bucket;

            det[det_i]  |= 1 << det_p ; //Set the bit to one
        }

        // Compute the permutation
        // https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        unsigned const t = p | (p - 1); 
        p = (t + 1) | (((~t & (t+1)) - 1) >> (__builtin_ctz(p) + 1));
    }
}

static void no_constexp_divmod(benchmark::State& state) {

    const unsigned size_orbital_bucket = sizeof(unsigned)*8;

    const unsigned n_orbital = 36;
    const unsigned n_alpha = 15;

    const unsigned n_int = n_orbital  / size_orbital_bucket  + 1 ;

    //unsigned *occ = (unsigned*) malloc ( n_orbital * sizeof(unsigned));
    unsigned occ[36] = {1,1,2,2,1,1,1,1,1,1,1,0,2,1,0,1,2,1,1,1,0,2,1,1,2,1,1,2,1,1,1,1,2,1,0,1};
 
    unsigned long n_det;
    unsigned *l_det_alpha;
    unsigned *l_det_beta;


  // Code inside this loop is measured repeatedly
  for (auto _ : state) {
    //unsigned long n =  n_combinaison(2000,10);
    // Make sure the variable is not optimized away by compiler
    gen_permutation0(n_orbital, n_int, size_orbital_bucket, occ, n_alpha, &n_det, &l_det_alpha, &l_det_beta);
    benchmark::DoNotOptimize(n_det);
  }
}
// Register the function as a benchmark
BENCHMARK(no_constexp_divmod);

void gen_permutation1(unsigned const n_orbital, unsigned const n_int,
                    unsigned const * const occ, unsigned const n_alpha, 
                    unsigned long * const n_det, unsigned **l_det_alpha, unsigned **l_det_beta)
{
    constexpr unsigned size_orbital_bucket = sizeof(unsigned)*8;

    // Will store the value 
    unsigned * const det_patern =  (unsigned*) malloc ( n_int * sizeof(unsigned));
    unsigned * const occ_if = (unsigned*) malloc ( n_orbital * sizeof(unsigned));
 
    unsigned i_int;
    for (i_int=0; i_int < n_int; i_int++){
        det_patern[i_int] = 0;
    }

    unsigned n_alpha_f = n_alpha;
    // Number of single occupied alpha orbital
    unsigned n_single_orbital = 0;

    unsigned idx;
    for (idx = 0 ; idx < n_orbital; idx++){
        switch(occ[idx]) {
            case 1 :
            {
                occ_if[n_single_orbital] = idx;
                n_single_orbital += 1;
                break;
            }
            case 2 :
            {
                unsigned const det_i = idx / size_orbital_bucket;
                unsigned const det_p = idx % size_orbital_bucket;

                det_patern[det_i] |= 1 << det_p; // set bit to one
                n_alpha_f -= 1;
                break;
            } 
        } 
   }
    assert (n_single_orbital <= (sizeof(unsigned)*8));

    *n_det = n_combinaison(n_single_orbital, n_alpha_f);
  
    *l_det_beta =  (unsigned*) malloc ( (*n_det) * n_int * sizeof(unsigned));
    *l_det_alpha =  (unsigned*) malloc ( (*n_det) * n_int * sizeof(unsigned));

    // p represent the permutation of n_alpha_f orbital in the single occupated orbitals
    // Initialize it with the first valid occupations
    // who correspond to the bitmask with `n_alpha_f` right bits set
    // aka p = \sum_{n=0}^m 2^n =  2^{m+1} - 1 == (1 << m) - 1. 
    unsigned p = (1 << n_alpha_f ) - 1;

    unsigned long idx_d;
    unsigned i;

    // Set the value of the determinant
    for (idx_d= 0 ; idx_d < (*n_det); idx_d++){

        unsigned *alpha = *l_det_alpha+idx_d*n_int;
        unsigned *beta =  *l_det_beta+idx_d*n_int; 

        for (i_int=0; i_int < n_int; i_int++){
            beta[i_int] = det_patern[i_int];
            alpha[i_int] = det_patern[i_int];
        }
        
        for ( i = 0 ; i < n_single_orbital; i++){
            bool const o = ( p >> i ) & 1; //index i of permutation p
            unsigned * const det = (o==1) ? alpha: beta; //Determinant are zero-ed .Modify only if need to set the bit to 1

            idx = occ_if[i];
            unsigned const det_i = idx / size_orbital_bucket;
            unsigned const det_p = idx % size_orbital_bucket;

            det[det_i]  |= 1 << det_p ; //Set the bit to one
        }

        // Compute the permutation
        // https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        unsigned const t = p | (p - 1); 
        p = (t + 1) | (((~t & (t+1)) - 1) >> (__builtin_ctz(p) + 1));
    }
}


static void constexp_divmod(benchmark::State& state) {

    constexpr unsigned size_orbital_bucket = sizeof(unsigned)*8;

    const unsigned n_orbital = 36;
    const unsigned n_alpha = 15;

    const unsigned n_int = n_orbital  / size_orbital_bucket  + 1 ;

    //unsigned *occ = (unsigned*) malloc ( n_orbital * sizeof(unsigned));
    unsigned occ[36] = {1,1,2,2,1,1,1,1,1,1,1,0,2,1,0,1,2,1,1,1,0,2,1,1,2,1,1,2,1,1,1,1,2,1,0,1};
 
    unsigned long n_det;
    unsigned *l_det_alpha;
    unsigned *l_det_beta;


  // Code inside this loop is measured repeatedly
  for (auto _ : state) {
    //unsigned long n =  n_combinaison(2000,10);
    // Make sure the variable is not optimized away by compiler
    gen_permutation1(n_orbital, n_int, occ, n_alpha, &n_det, &l_det_alpha, &l_det_beta);
    benchmark::DoNotOptimize(n_det);
  }
}
// Register the function as a benchmark
BENCHMARK(constexp_divmod);


void gen_permutation2(unsigned const n_orbital, unsigned const n_int,
                    unsigned const * const occ, unsigned const n_alpha, 
                    unsigned long * const n_det, unsigned **l_det_alpha, unsigned **l_det_beta)
{

    constexpr unsigned size_orbital_bucket = sizeof(unsigned)*8;
    constexpr unsigned log_size_orbital_bucket = __builtin_ctz(size_orbital_bucket);

    // Will store the value 
    unsigned * const det_patern =  (unsigned*) malloc ( n_int * sizeof(unsigned));
    unsigned * const occ_if = (unsigned*) malloc ( n_orbital * sizeof(unsigned));

 
    unsigned i_int;
    for (i_int=0; i_int < n_int; i_int++){
        det_patern[i_int] = 0;
    }

    unsigned n_alpha_f = n_alpha;
    // Number of single occupied alpha orbital
    unsigned n_single_orbital = 0;

    unsigned idx;
    for (idx = 0 ; idx < n_orbital; idx++){
        switch(occ[idx]) {
            case 1 :
            {
                occ_if[n_single_orbital] = idx;
                n_single_orbital += 1;
                break;
            }
            case 2 :
            {
                unsigned const det_i = idx >> log_size_orbital_bucket;
                unsigned const det_p = idx - (det_i  -1 << log_size_orbital_bucket);

                det_patern[det_i] |= 1 << det_p; // set bit to one
                n_alpha_f -= 1;
                break;
            } 
        } 
   }
    assert (n_single_orbital <= (sizeof(unsigned)*8));

    *n_det = n_combinaison(n_single_orbital, n_alpha_f);
  
    *l_det_beta =  (unsigned*) malloc ( (*n_det) * n_int * sizeof(unsigned));
    *l_det_alpha =  (unsigned*) malloc ( (*n_det) * n_int * sizeof(unsigned));

    // p represent the permutation of n_alpha_f orbital in the single occupated orbitals
    // Initialize it with the first valid occupations
    // who correspond to the bitmask with `n_alpha_f` right bits set
    // aka p = \sum_{n=0}^m 2^n =  2^{m+1} - 1 == (1 << m) - 1. 
    unsigned p = (1 << n_alpha_f ) - 1;

    unsigned long idx_d;
    unsigned i;

    // Set the value of the determinant
    for (idx_d= 0 ; idx_d < (*n_det); idx_d++){

        unsigned *alpha = *l_det_alpha+idx_d*n_int;
        unsigned *beta =  *l_det_beta+idx_d*n_int; 

        for (i_int=0; i_int < n_int; i_int++){
            beta[i_int] = det_patern[i_int];
            alpha[i_int] = det_patern[i_int];
        }
        
        for ( i = 0 ; i < n_single_orbital; i++){
            bool const o = ( p >> i ) & 1; //index i of permutation p
            unsigned * const det = o ? alpha: beta; //Determinant are zero-ed .Modify only if need to set the bit to 1

            idx = occ_if[i];
            unsigned const det_i = idx >> log_size_orbital_bucket;
            unsigned const det_p = idx - (det_i  -1 << log_size_orbital_bucket);

            det[det_i]  |= 1 << det_p ; //Set the bit to one
        }

        // Compute the permutation
        // https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        unsigned const t = p | (p - 1); 
        p = (t + 1) | (((~t & (t+1)) - 1) >> (__builtin_ctz(p) + 1));
    }
}


static void constexp_no_divmod(benchmark::State& state) {

    constexpr unsigned size_orbital_bucket = sizeof(unsigned)*8;

    const unsigned n_orbital = 36;
    const unsigned n_alpha = 15;

    const unsigned n_int = n_orbital  / size_orbital_bucket  + 1 ;

    //unsigned *occ = (unsigned*) malloc ( n_orbital * sizeof(unsigned));
    unsigned occ[36] = {1,1,2,2,1,1,1,1,1,1,1,0,2,1,0,1,2,1,1,1,0,2,1,1,2,1,1,2,1,1,1,1,2,1,0,1};
 
    unsigned long n_det;
    unsigned *l_det_alpha;
    unsigned *l_det_beta;


  // Code inside this loop is measured repeatedly
  for (auto _ : state) {
    //unsigned long n =  n_combinaison(2000,10);
    // Make sure the variable is not optimized away by compiler
    gen_permutation2(n_orbital, n_int, occ, n_alpha, &n_det, &l_det_alpha, &l_det_beta);
    benchmark::DoNotOptimize(n_det);
  }
}
// Register the function as a benchmark
BENCHMARK(constexp_no_divmod);


void gen_permutation3(unsigned const n_orbital, unsigned const n_int,
                    unsigned const * const occ, unsigned const n_alpha, 
                    unsigned long * const n_det, unsigned **l_det_alpha, unsigned **l_det_beta)
{

    const unsigned size_orbital_bucket = sizeof(unsigned)*8;
    const unsigned log_size_orbital_bucket = __builtin_ctz(size_orbital_bucket);

    // Will store the value 
    unsigned * const det_patern =  (unsigned*) malloc ( n_int * sizeof(unsigned));
    unsigned * const occ_if = (unsigned*) malloc ( n_orbital * sizeof(unsigned));

 
    unsigned i_int;
    for (i_int=0; i_int < n_int; i_int++){
        det_patern[i_int] = 0;
    }

    unsigned n_alpha_f = n_alpha;
    // Number of single occupied alpha orbital
    unsigned n_single_orbital = 0;

    unsigned idx;
    for (idx = 0 ; idx < n_orbital; idx++){
        switch(occ[idx]) {
            case 1 :
            {
                occ_if[n_single_orbital] = idx;
                n_single_orbital += 1;
                break;
            }
            case 2 :
            {
                unsigned const det_i = idx >> log_size_orbital_bucket;
                unsigned const det_p = idx - ( (det_i  -1) << log_size_orbital_bucket);

                det_patern[det_i] |= 1 << det_p; // set bit to one
                n_alpha_f -= 1;
                break;
            } 
        } 
   }
    assert (n_single_orbital <= (sizeof(unsigned)*8));

    *n_det = n_combinaison(n_single_orbital, n_alpha_f);
  
    *l_det_beta =  (unsigned*) malloc ( (*n_det) * n_int * sizeof(unsigned));
    *l_det_alpha =  (unsigned*) malloc ( (*n_det) * n_int * sizeof(unsigned));

    unsigned i;

    unsigned *l_det_i =  (unsigned*) malloc (n_single_orbital * sizeof(unsigned));
    unsigned *l_det_p =  (unsigned*) malloc (n_single_orbital * sizeof(unsigned));

    // p represent the permutation of n_alpha_f orbital in the single occupated orbitals
    // Initialize it with the first valid occupations
    // who correspond to the bitmask with `n_alpha_f` right bits set
    // aka p = \sum_{n=0}^m 2^n =  2^{m+1} - 1 == (1 << m) - 1. 
    unsigned p = (1 << n_alpha_f ) - 1;

    unsigned long idx_d;

    // Set the value of the determinant
    for (idx_d= 0 ; idx_d < (*n_det); idx_d++){

        unsigned *alpha = *l_det_alpha+idx_d*n_int;
        unsigned *beta =  *l_det_beta+idx_d*n_int; 

        for (i_int=0; i_int < n_int; i_int++){
            beta[i_int] = det_patern[i_int];
            alpha[i_int] = det_patern[i_int];
        }
        
        for ( i = 0 ; i < n_single_orbital; i++){
            bool const o = 1; //( p >> i ) & 1; //index i of permutation p
            unsigned * const det = o ? alpha: beta; //Determinant are zero-ed .Modify only if need to set the bit to 1

            idx = occ_if[i];
            unsigned const det_i = idx >> log_size_orbital_bucket;
            unsigned const det_p = idx - ( (det_i  -1) << log_size_orbital_bucket);

            det[det_i]  |= 1 << det_p ; //Set the bit to one
        }

        // Compute the permutation
        // https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        unsigned const t = p | (p - 1); 
        p = (t + 1) | (((~t & (t+1)) - 1) >> (__builtin_ctz(p) + 1));
    }
}

static void no_constexp_no_divmod(benchmark::State& state) {

    constexpr unsigned size_orbital_bucket = sizeof(unsigned)*8;

    const unsigned n_orbital = 36;
    const unsigned n_alpha = 15;

    const unsigned n_int = n_orbital  / size_orbital_bucket  + 1 ;

    //unsigned *occ = (unsigned*) malloc ( n_orbital * sizeof(unsigned));
    unsigned occ[36] = {1,1,2,2,1,1,1,1,1,1,1,0,2,1,0,1,2,1,1,1,0,2,1,1,2,1,1,2,1,1,1,1,2,1,0,1};
 
    unsigned long n_det;
    unsigned *l_det_alpha;
    unsigned *l_det_beta;


  // Code inside this loop is measured repeatedly
  for (auto _ : state) {
    //unsigned long n =  n_combinaison(2000,10);
    // Make sure the variable is not optimized away by compiler
    gen_permutation3(n_orbital, n_int, occ, n_alpha, &n_det, &l_det_alpha, &l_det_beta);
    benchmark::DoNotOptimize(n_det);
  }
}
// Register the function as a benchmark
BENCHMARK(no_constexp_no_divmod);

BENCHMARK_MAIN();