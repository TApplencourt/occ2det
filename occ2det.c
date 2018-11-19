#define PRINTF_BINARY_PATTERN_INT8 "%c%c%c%c%c%c%c%c"
#define PRINTF_BYTE_TO_BINARY_INT8(i)    \
    (((i) & 0x80ll) ? '1' : '0'), \
    (((i) & 0x40ll) ? '1' : '0'), \
    (((i) & 0x20ll) ? '1' : '0'), \
    (((i) & 0x10ll) ? '1' : '0'), \
    (((i) & 0x08ll) ? '1' : '0'), \
    (((i) & 0x04ll) ? '1' : '0'), \
    (((i) & 0x02ll) ? '1' : '0'), \
    (((i) & 0x01ll) ? '1' : '0')

#define PRINTF_BINARY_PATTERN_INT16 \
    PRINTF_BINARY_PATTERN_INT8              PRINTF_BINARY_PATTERN_INT8
#define PRINTF_BYTE_TO_BINARY_INT16(i) \
    PRINTF_BYTE_TO_BINARY_INT8((i) >> 8),   PRINTF_BYTE_TO_BINARY_INT8(i)
#define PRINTF_BINARY_PATTERN_INT32 \
    PRINTF_BINARY_PATTERN_INT16             PRINTF_BINARY_PATTERN_INT16
#define PRINTF_BYTE_TO_BINARY_INT32(i) \
    PRINTF_BYTE_TO_BINARY_INT16((i) >> 16), PRINTF_BYTE_TO_BINARY_INT16(i)
#define PRINTF_BINARY_PATTERN_INT64    \
    PRINTF_BINARY_PATTERN_INT32             PRINTF_BINARY_PATTERN_INT32
#define PRINTF_BYTE_TO_BINARY_INT64(i) \
    PRINTF_BYTE_TO_BINARY_INT32((i) >> 32), PRINTF_BYTE_TO_BINARY_INT32(i)

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

void gen_permutation(unsigned const n_orbital, unsigned const n_int, unsigned const log_size_orbital_bucket, 
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
                unsigned const det_i =  (idx >> log_size_orbital_bucket);
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
            unsigned const det_i =  (idx >> log_size_orbital_bucket);
            unsigned const det_p = idx - ( (det_i  -1) << log_size_orbital_bucket); 

            det[det_i]  |= 1 << det_p ; //Set the bit to one
        }

        // Compute the permutation
        // https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        unsigned const t = p | (p - 1); 
        //p = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(p) + 1));
        p = (t + 1) | (((~t & (t+1)) - 1) >> (__builtin_ctz(p) + 1));
    }
}

int main(int argc, char **argv)
{
    const unsigned size_orbital_bucket = sizeof(unsigned)*8;
    const unsigned log_size_orbital_bucket = __builtin_ctz(size_orbital_bucket);
    const unsigned mode = atoi(argv[1]); 
    const unsigned n_orbital = atoi(argv[2]);;
    const unsigned n_alpha =  atoi(argv[3]);

    const unsigned n_int = n_orbital  / size_orbital_bucket  + 1 ;
    printf("n_int: %u\n", n_int);

    unsigned *occ = (unsigned*) malloc ( n_orbital * sizeof(unsigned));

    unsigned d, i = 0;
    for (i= 0 ; i < n_orbital; i++) {
        occ[i] = atoi(argv[i+4]);
    }
 
    unsigned long n_det;
    unsigned *l_det_alpha;
    unsigned *l_det_beta;

    gen_permutation(n_orbital, n_int, log_size_orbital_bucket, occ, n_alpha, &n_det, &l_det_alpha,  &l_det_beta);

    printf("n_det: %lu\n", n_det);
  
    // 0 == Production  mode ; 1 == Test mode 
    d = (mode) ? 0 : n_det; 

    for (d  ; d < n_det; d++){
        for (i =0 ; i < n_int; i++){
            printf(PRINTF_BINARY_PATTERN_INT32"\n", PRINTF_BYTE_TO_BINARY_INT32(l_det_alpha[d*n_int + i]));
        }
        for (i =0 ; i < n_int; i++){
             printf(PRINTF_BINARY_PATTERN_INT32"\n", PRINTF_BYTE_TO_BINARY_INT32(l_det_beta[d*n_int + i]));
        }
    }

    return 0;    
}
