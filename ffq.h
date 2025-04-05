/* ffq.h, Chris Monico.
     Code for synthesis of and basic arithmetic over small finite fields.
   There is nothing at all fancy here - this is minimal functionality
   for reasonably fast arithmetic over small finite fields, say with
   orders upto a few thousand or so.

     This probably should all work for fields or order up to 2^{31} or so;
   it seems to be fine at nextprime(2^{31}) = 2147483659. But it
   definitely fails for prime fields whose order is a bit bigger than that,
   because of multiplicative overflow. It was anyway intended for small
   finite fields.

   NOT thread-safe.
*/
#ifndef _FFQ_H
#define _FFQ_H

#ifndef MAX
#define MAX(_a, _b) ((_a)>(_b)?(_a):(_b))
#define MIN(_a, _b) ((_a)<(_b)?(_a):(_b))
#endif


#include <stdint.h>

typedef struct {
  int      deg;
  int64_t *coefs;
  int      alloc_deg;
} poly_t;



typedef struct {
  long p;
  int  k;
  poly_t poly;
  poly_t prim_elt;
} ffq_t;

int   poly_alloc(poly_t *f, int deg);
void  poly_dealloc(poly_t *f);
void  poly_cp(poly_t *res, poly_t *op);
char *poly_str(char *str, poly_t *f);

/*******************************************************/
int eratosthenes(int *p, int p_size, int max_p);
  /* Find all primes in [2, max_p], and store them in the given array.
     Return numprimes on success, negative on failure (e.g., p was not large enough). */
int distinct_prime_divisors(long *p, int p_size, long n);


/*******************************************************/
void  ffq_init(ffq_t *F);
void  ffq_clear(ffq_t *F);
long  ffq_size(ffq_t *F);
char *ffq_str(char *str, poly_t *op, ffq_t *F);
int   ffq_is_zero(poly_t *op, ffq_t *F);
int   ffq_is_one(poly_t *op, ffq_t *F);
void  ffq_set_zero(poly_t *res, ffq_t *F);
void  ffq_set_one(poly_t *res, ffq_t *F);
void  ffq_set_x(poly_t *res, ffq_t *F);
void  ffq_rand(poly_t *res, ffq_t *F);

void ffq_add(poly_t *res, poly_t *op1, poly_t *op2, ffq_t *F);

void ffq_sub(poly_t *res, poly_t *op1, poly_t *op2, ffq_t *F);

void ffq_mul_c(poly_t *res, poly_t *op1, int64_t c, ffq_t *F);

void ffq_mul(poly_t *res, poly_t *op1, poly_t *op2, ffq_t *F);

int  ffq_inv(poly_t *res, poly_t *op1, ffq_t *F);

void ffq_pow(poly_t *res, poly_t *op1, int64_t pow, ffq_t *F);

int64_t ffq_to_enum(poly_t *op, ffq_t *F);
  /* Fix an enumeration of 'F', and return the value of 'op' under that enumeration. */
void ffq_from_enum(poly_t *op, int64_t n, ffq_t *F);
  /* Inverse of the above. */


int ffq_elt_is_primitive(poly_t *op, ffq_t *F);
int ffq_create(ffq_t *F, int64_t p, int k);

#endif
