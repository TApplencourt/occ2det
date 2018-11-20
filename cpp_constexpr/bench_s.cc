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

    const unsigned n_orbital = 20;
    const unsigned n_alpha = 10;

    const unsigned n_int = n_orbital  / size_orbital_bucket  + 1 ;

    unsigned occ[20] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

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

    const unsigned n_orbital = 20;
    const unsigned n_alpha = 10;

    const unsigned n_int = n_orbital  / size_orbital_bucket  + 1 ;

    unsigned occ[20] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

    unsigned long n_det;
    unsigned *l_det_alpha;
    unsigned *l_det_beta;


  // Code inside this loop is measured repeatedly
  for (auto _ : state) {
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

    const unsigned n_orbital = 20;
    const unsigned n_alpha = 10;

    const unsigned n_int = n_orbital  / size_orbital_bucket  + 1 ;

    unsigned occ[20] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

    unsigned long n_det;
    unsigned *l_det_alpha;
    unsigned *l_det_beta;


  // Code inside this loop is measured repeatedly
  for (auto _ : state) {
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
            bool const o = ( p >> i ) & 1; //index i of permutation p
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

    const unsigned n_orbital = 20;
    const unsigned n_alpha = 10;

    const unsigned n_int = n_orbital  / size_orbital_bucket  + 1 ;

    unsigned occ[20] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

    unsigned long n_det;
    unsigned *l_det_alpha;
    unsigned *l_det_beta;


  // Code inside this loop is measured repeatedly
  for (auto _ : state) {
    gen_permutation3(n_orbital, n_int, occ, n_alpha, &n_det, &l_det_alpha, &l_det_beta);
    benchmark::DoNotOptimize(n_det);
  }
}
// Register the function as a benchmark
BENCHMARK(no_constexp_no_divmod);


void gen_permutation4(unsigned const n_orbital, unsigned const n_int,
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

    // p represent the permutation of n_alpha_f orbital in the single occupated orbitals
    // Initialize it with the first valid occupations
    // who correspond to the bitmask with `n_alpha_f` right bits set
    // aka p = \sum_{n=0}^m 2^n =  2^{m+1} - 1 == (1 << m) - 1.
    unsigned p = (1 << n_alpha_f ) - 1;

    unsigned long idx_d;
    const unsigned n_det_2 = (*n_det)/2; 

    // Set the value of the determinant
    for (idx_d= 0 ; idx_d < n_det_2; idx_d++){

        unsigned *alpha = *l_det_alpha+(idx_d+1)*n_int;
        unsigned *beta =  *l_det_beta+(idx_d+1)*n_int;

        for (i_int=0; i_int < n_int; i_int++){
            beta[i_int] = det_patern[i_int];
            alpha[i_int] = det_patern[i_int];
        }

        for ( i = 0 ; i < n_single_orbital; i++){
            bool const o = ( p >> i ) & 1; //index i of permutation p

            idx = occ_if[i];
            unsigned const det_i = idx >> log_size_orbital_bucket;
            unsigned const det_p = idx - ( (det_i  -1) << log_size_orbital_bucket);

            if (o) {
                alpha[det_i]  |= 1 << det_p ; //Set the bit to one
                beta[det_i+1] |= 1 << det_p ; //Set the bit to one
            }
            else {
                beta[det_i]    |= 1 << det_p ; //Set the bit to one
                alpha[det_i+1] |= 1 << det_p ; //Set the bit to one
            }
        }

        // Compute the permutation
        // https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        unsigned const t = p | (p - 1);
        p = (t + 1) | (((~t & (t+1)) - 1) >> (__builtin_ctz(p) + 1));
    }
}

static void no_constexp_no_divmod_sym(benchmark::State& state) {

    constexpr unsigned size_orbital_bucket = sizeof(unsigned)*8;

    const unsigned n_orbital = 20;
    const unsigned n_alpha = 10;

    const unsigned n_int = n_orbital  / size_orbital_bucket  + 1 ;

    unsigned occ[20] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

    unsigned long n_det;
    unsigned *l_det_alpha;
    unsigned *l_det_beta;


  // Code inside this loop is measured repeatedly
  for (auto _ : state) {
    gen_permutation4(n_orbital, n_int, occ, n_alpha, &n_det, &l_det_alpha, &l_det_beta);
    benchmark::DoNotOptimize(n_det);
  }
}

BENCHMARK(no_constexp_no_divmod_sym);


void gen_permutation5(unsigned const n_orbital, unsigned const n_int,
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
                unsigned const det_i =  (idx >> log_size_orbital_bucket);
                unsigned const det_p = idx - ( (det_i  -1) << log_size_orbital_bucket);

                det_patern[det_i] |= 1 << det_p; // set bit to one
                n_alpha_f -= 1;
                break;
            }
        }
   }
    assert (n_single_orbital <= (sizeof(unsigned)*8));

    unsigned n_beta_f = n_single_orbital - n_alpha_f;
    //printf("n_single_orbital  %u\n",n_single_orbital);
    //printf("n_alpha_f %u\n",n_alpha_f);
    //printf("n_beta_f %u\n",n_beta_f);

    unsigned  n_pair = n_beta_f;
    //printf("n_pair %u\n", n_pair);
    unsigned n_peel = n_alpha_f - n_beta_f;
    //printf("n_peel %u\n", n_peel);

    const unsigned n_det_pair = n_combinaison(n_single_orbital, n_alpha_f);

    *n_det = n_det_pair;

    *l_det_beta =  (unsigned*) malloc ( (*n_det) * n_int * sizeof(unsigned));
    *l_det_alpha =  (unsigned*) malloc ( (*n_det) * n_int * sizeof(unsigned));

    // p represent the permutation of n_alpha_f orbital in the single occupated orbitals
    // Initialize it with the first valid occupations
    // who correspond to the bitmask with `n_alpha_f` right bits set
    // aka p = \sum_{n=0}^m 2^n =  2^{m+1} - 1 == (1 << m) - 1.
    unsigned p = (1 << n_alpha_f ) - 1;

    unsigned long idx_d;
    unsigned i;

    unsigned delta_d = 0;

    // Set the value of the determinant
     for (idx_d= 0 ; idx_d < n_det_pair/2; idx_d++){

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
            unsigned const det_i =  (idx >> log_size_orbital_bucket);
            unsigned const det_p = idx - ( (det_i  -1) << log_size_orbital_bucket);

            if (o){
                alpha[det_i]  |= 1 << det_p ; //Set the bit to one
                //beta[det_i+(n_det_pair/2)] |= 1 << det_p ; //Set the bit to one
            }
            else{
                beta[det_i]  |= 1 << det_p ; //Set the bit to one
                //alpha[det_i+(n_det_pair/2)] |= 1 << det_p ; //Set the bit to one
            }
        }

        // Compute the permutation
        // https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        unsigned const t = p | (p - 1);
        p = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(p) + 1));
    }

    
    for (idx_d= 0 ; idx_d < n_det_pair/2; idx_d++){
        unsigned *alpha = *l_det_alpha;
        unsigned *beta =  *l_det_beta;

        for (i_int=0; i_int < n_int; i_int++){
            alpha[n_det_pair - idx_d-1 + i_int*n_int] = beta[idx_d  + i_int*n_int];
            beta[n_det_pair - idx_d -1+ i_int*n_int] = alpha[idx_d + i_int*n_int];
        }
    }
    
}

static void no_constexp_no_divmod_sym2(benchmark::State& state) {

    constexpr unsigned size_orbital_bucket = sizeof(unsigned)*8;

    const unsigned n_orbital = 20;
    const unsigned n_alpha = 10;

    const unsigned n_int = n_orbital  / size_orbital_bucket  + 1 ;

    unsigned occ[20] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

    unsigned long n_det;
    unsigned *l_det_alpha;
    unsigned *l_det_beta;


  // Code inside this loop is measured repeatedly
  for (auto _ : state) {
    gen_permutation5(n_orbital, n_int, occ, n_alpha, &n_det, &l_det_alpha, &l_det_beta);
    benchmark::DoNotOptimize(n_det);
  }
}

BENCHMARK(no_constexp_no_divmod_sym2);



BENCHMARK_MAIN();
