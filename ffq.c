/* ffq.c, Chris Monico.
   Code for synthesis of and basic arithmetic over small finite fields.
   There is nothing at all fancy here - this is minimal functionality
   for reasonably fast arithmetic over small finite fields, say with
   orders upto a few thousand or so.
*/
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ffq.h"
#include "common.h"



/*******************************************************/
long gcd_ext(long *a, long *b, long op1, long op2)
/* Set (a,b) such that g = gcd(op1, op2) = a*op1 + b*op2, and return g. */
{ long a1, b1, a2, b2, tmp, q;

  a1 = 1; b1 = 0;
  a2 = 0; b2 = 1;
  if (op1 < 0) { op1 = -op1; a1 = -1; }
  if (op2 < 0) { op2 = -op2; b2 = -1; }
  while (op2) {
    /* op1 <-- op1 % op2, and corresponding adjustments to ai,bi. */
    q = op1/op2;
    op1 -= q*op2;
    a1 -= q*a2;
    b1 -= q*b2;
    /* Swap (op1, a1, b1) <--> (op2, a2, b2). */
    tmp = op1; op1 = op2; op2 = tmp;
    tmp = a1; a1 = a2; a2 = tmp;
    tmp = b1; b1 = b2; b2 = tmp;
  }
  *a = a1;
  *b = b1;
  return op1;
}
/*******************************************************/
long inv_mod_n(long x, long n)
{ long g, a, b;

  g = gcd_ext(&a, &b, x, n);
  if (g==1)
    return (a+n)%n;
  return 0;
}



/**********************************************************/
int poly_alloc(poly_t *f, int deg)
{
  if (!(f->coefs = (int64_t *)malloc((deg+1)*sizeof(int64_t)))) {
    f->deg = f->alloc_deg = -1;
    return -1;
  }
  f->alloc_deg = deg;
  for (int i=0; i<=deg; i++)
    f->coefs[i] = 0;
  f->deg = 0;
  return 0;
}
  
/**********************************************************/
void poly_dealloc(poly_t *f)
{ 
  if (f->coefs) 
    free(f->coefs);
  f->coefs = NULL;
  f->deg = f->alloc_deg = -1;
}

/**********************************************************/
int poly_realloc(poly_t *f, int deg)
{
  if (deg <= f->alloc_deg) return 0;
  poly_dealloc(f);
  return poly_alloc(f, deg);
}

/**********************************************************/
void poly_fixdeg(poly_t *f)
{ 
  while ((f->deg > 0) && (f->coefs[f->deg]==0))
    f->deg -= 1;
}

/**********************************************************/
void poly_cp(poly_t *res, poly_t *op)
{
  poly_realloc(res, op->deg);
  memcpy(res->coefs, op->coefs, (op->deg+1)*sizeof(int64_t));
  res->deg = op->deg;
}

/**********************************************************/
char *poly_str(char *str, poly_t *f)
{ char s[64];
  int  i;

  sprintf(str, "(");
  for (i=0; i<=f->deg; i++) {
    sprintf(s, "%ld%c", f->coefs[i], (i<f->deg)?',':')');
    strcat(str, s);
  }
  return str;
}


/**********************************************************/
void poly_add(poly_t *res, poly_t *op1, poly_t *op2, int64_t modulus)
{ int deg, i, m;
  poly_t tmp;

  deg = MAX(op1->deg, op2->deg);
  poly_alloc(&tmp, deg);

  i=0;
  m = MIN(op1->deg, op2->deg);
  while (i<=m) {
    tmp.coefs[i] = (op1->coefs[i] + op2->coefs[i])%modulus;
    i++;
  }
  while (i<=op1->deg) {
    tmp.coefs[i] = op1->coefs[i]%modulus;
    i++;
  }
  while (i<=op2->deg) {
    tmp.coefs[i] = op2->coefs[i]%modulus;
    i++;
  }
  tmp.deg = deg;
  poly_fixdeg(&tmp);
  poly_cp(res, &tmp);
  poly_dealloc(&tmp);
}


/**********************************************************/
void poly_sub(poly_t *res, poly_t *op1, poly_t *op2, int64_t modulus)
{ int deg, i, m;
  poly_t tmp;

  deg = MAX(op1->deg, op2->deg);
  poly_alloc(&tmp, deg);

  i=0;
  m = MIN(op1->deg, op2->deg);
  while (i<=m) {
    tmp.coefs[i] = (op1->coefs[i] - op2->coefs[i] + modulus)%modulus;
    i++;
  }
  while (i<=op1->deg) {
    tmp.coefs[i] = op1->coefs[i]%modulus;
    i++;
  }
  while (i<=op2->deg) {
    tmp.coefs[i] = (modulus - op2->coefs[i])%modulus;
    i++;
  }
  tmp.deg = deg;
  poly_fixdeg(&tmp);
  poly_cp(res, &tmp);
  poly_dealloc(&tmp);
}

/**********************************************************/
void poly_mul(poly_t *res, poly_t *op1, poly_t *op2, int64_t modulus)
{ int  deg;
  int  i, j;
  poly_t R;

  deg = op1->deg + op2->deg;
  poly_alloc(&R, deg);
 
  for (i=0; i<=deg; i++)
    R.coefs[i] = 0;

  for (i=0; i<=op1->deg; i++) 
    for (j=0; j<=op2->deg; j++)  
      R.coefs[i+j] = (R.coefs[i+j] + op1->coefs[i]*op2->coefs[j])%modulus;
  R.deg = deg;
  poly_fixdeg(&R);

  poly_cp(res, &R);
  poly_dealloc(&R);
}


/**********************************************************/
int poly_div(poly_t *q, poly_t *r, poly_t *a, poly_t *b, int64_t modulus)
/* Get q,r such that a = q*b + r  (mod p). 
   Arguments {q,r} must be disjoint from {a,b}! */
{ int     da, db, i, j, dq;
  int64_t c_b, c;

  poly_fixdeg(a);
  poly_fixdeg(b);

  if ((b->deg == 0) && (b->coefs[0]==0)) 
    return -1; /* Division by zero. */

  dq = MAX(a->deg-b->deg, 0);
  poly_realloc(q, dq);
  poly_realloc(r, a->deg); /* It will be smaller, but we start it here. */

  q->deg = dq;
  for (i=0; i<=q->deg; i++)
    q->coefs[i] = 0;
  poly_cp(r, a);


  c_b = inv_mod_n(b->coefs[b->deg], modulus);

  for (i=r->deg; i>=b->deg; i--) {
    c = (r->coefs[i] * c_b) % modulus;
    /* Do r <-- r - c*b*x^{i - deg(b)}. */
    for (j=b->deg; j>=0; j--) 
      r->coefs[j+i-b->deg] = (modulus + r->coefs[j+i-b->deg] - (c*b->coefs[j])%modulus)%modulus;
    /* And q <-- q + c*x^{i-deg(b)}. */
    q->coefs[i-b->deg] = c;
  }

  poly_fixdeg(q);
  //r->deg = MAX(b->deg-1, 0);
  poly_fixdeg(r);
  return 0;
}


/**********************************************************/
int poly_mod(poly_t *res, poly_t *op1, poly_t *op2, int64_t modulus)
{ poly_t tmpQ, tmpR;
  int    retval;

  poly_alloc(&tmpQ, op1->deg);
  poly_alloc(&tmpR, op1->deg);
  retval = poly_div(&tmpQ, &tmpR, op1, op2, modulus);
  poly_cp(res, &tmpR);
  poly_dealloc(&tmpQ);
  poly_dealloc(&tmpR);
  return retval;
}

/**********************************************************/
int  poly_gcd_ext(poly_t *g, poly_t *a, poly_t *b, poly_t *P1, poly_t *P2, int64_t modulus)
  /* Get  (g,a,b) such that gcd(op1, op2) = g = a*P1 + b*P2 (mod p). */
{ poly_t *a1, *b1, *a2, *b2, *op1, *op2, *ptr;
  poly_t _a1, _b1, _a2, _b2, _op1, _op2, tmp, q, r;
  int    d = MAX(P1->deg, P2->deg);  
  int64_t c;

  poly_alloc(&_a1, d); poly_alloc(&_b1, d);
  poly_alloc(&_a2, d); poly_alloc(&_b2, d);
  poly_alloc(&_op1, P1->deg);
  poly_alloc(&_op2, P1->deg);
  poly_cp(&_op1, P1);
  poly_cp(&_op2, P2);

  a1 = &_a1; b1 = &_b1;
  a2 = &_a2; b2 = &_b2;
  op1 = &_op1; op2 = &_op2;

  a1->deg = b1->deg = 0;
  a1->coefs[0] = 1; b1->coefs[0] = 0;
  a2->deg = b2->deg = 0;
  a2->coefs[0] = 0; b2->coefs[0] = 1;

  poly_alloc(&q, d);
  poly_alloc(&r, d);
  poly_alloc(&tmp, d);

  while ((op2->deg > 0) || (op2->coefs[0])) {
    /* op1 <-- op1 - q*op2 = (op1 % op2). */
    poly_div(&q, &r, op1, op2, modulus);
    /* op1 <-- op1 - q*op2. */
    poly_cp(op1, &r);
    /* a1 <-- a1 - q*a2. */
    poly_mul(&tmp, &q, a2, modulus);
    poly_sub(a1, a1, &tmp, modulus);
    /* b1 <-- b1 - q*b2. */
    poly_mul(&tmp, &q, b2, modulus);
    poly_sub(b1, b1, &tmp, modulus);

    /* Swap (op1, a1, b1) <--> (op2, a2, b2). */
    ptr = op1; op1 = op2; op2 = ptr;
    ptr = a1; a1 = a2; a2 = ptr;
    ptr = b1; b1 = b2; b2 = ptr;
  }

  poly_cp(g, op1);
  poly_cp(a, a1);
  poly_cp(b, b1);

  poly_dealloc(&q); poly_dealloc(&r); poly_dealloc(&tmp);
  poly_dealloc(&_op1); poly_dealloc(&_op2);
  poly_dealloc(&_a1); poly_dealloc(&_a2);    
  poly_dealloc(&_b1); poly_dealloc(&_b2);    

  /* If g is constant, we should scale it to 1. */
  if (g->deg == 0) {
    c = inv_mod_n(g->coefs[0], modulus);
    for (int i=0; i<=a->deg; i++)
      a->coefs[i] = (a->coefs[i]*c)%modulus;
    for (int i=0; i<=b->deg; i++)
      b->coefs[i] = (b->coefs[i]*c)%modulus;
    g->coefs[0] = 1;
  }
  return 0;
}  

/**********************************************************/
int poly_gcd(poly_t *g, poly_t *op1, poly_t *op2, int64_t modulus)
{ poly_t a, b;
  int    d = MAX(op1->deg, op2->deg), retval;

  poly_alloc(&a, d);
  poly_alloc(&b, d);
  retval = poly_gcd_ext(g, &a, &b, op1, op2, modulus);
  poly_dealloc(&a);
  poly_dealloc(&b);
  return retval;
}

/**********************************************************/
void poly_powmod(poly_t *res, poly_t *base, int64_t e, poly_t *M, int64_t modulus)
{ poly_t R, pow, tmp;

  poly_alloc(&R, M->deg); 
  poly_alloc(&pow, M->deg);
  poly_alloc(&tmp, M->deg);
  R.deg = 0;
  R.coefs[0] = 1;
  poly_cp(&pow, base);
  while (e>0) {
    if (e%2) {
      poly_mul(&tmp, &R, &pow, modulus);
      poly_mod(&R, &tmp, M, modulus);
    }
    e = e/2;
    poly_mul(&tmp, &pow, &pow, modulus);
    poly_mod(&pow, &tmp, M, modulus);
  }
  poly_cp(res, &R);
  poly_dealloc(&R);
  poly_dealloc(&pow);
  poly_dealloc(&tmp);
}

/**********************************************************/
int poly_is_irreducible(poly_t *f, int64_t p)
/* This is not particularly efficient, but should be good enough
   for the problem at hand. We use the fact that f is irreducible over Z/pZ
   iff 
        x^{p^{deg}} == x (mod f), and
   and for all primes q | deg,
        gcd(x^{p^{deg/q}} - x, f) = 1.
  For the problem at hand, we are generating finite fields whose
  elements can be easily enumerated. We assume, in particular,
  that p^{deg} < 2^{32}.
*/
{ long   power;
  int    i, j, deg = f->deg, retval;
  int    smallprimes[] = {2,3,5,7,11,13,17,19,23};
  int    num_sp = sizeof(smallprimes)/sizeof(int);
  poly_t res, x, g;

  retval = 1;
  poly_alloc(&res, deg);
  poly_alloc(&x, deg);
  poly_alloc(&g, deg);

  power = p;
  for (i=2; i<= deg; i++)
    power *= p;
  x.deg = 1;
  x.coefs[0] = 0;
  x.coefs[1] = 1;

  poly_powmod(&res, &x, power, f, p);

  /* If res != x, then f is not irreducible. */
  if ((res.deg != 1) || (res.coefs[0]) || (res.coefs[1] != 1)) {
    retval = 0;
    goto IRR_DONE;
  }

 for (j=0; j<num_sp; j++) {
    if ((deg % smallprimes[j])==0) {
      power = p;
      for (i=2; i<= deg/smallprimes[j]; i++)
        power *= p;
      poly_powmod(&res, &x, power, f, p);
      /* res <-- res - x : */
      res.coefs[1] = (res.coefs[1] - 1 + p)%p;
      if (res.deg == 0) res.deg = 1;
      poly_fixdeg(&res);
      poly_gcd(&g, &res, f, p);
      if ((g.deg > 0) || (g.coefs[0] != 1)) {
        retval = 0;
        goto IRR_DONE;
      }
    }
  }

IRR_DONE:
  poly_alloc(&res, deg);
  poly_alloc(&x, deg);
  poly_alloc(&g, deg);

  return retval;
}


int *_small_primes = NULL;
int  _num_small_primes=0;

/************************************************************/
int eratosthenes(int *Primes, int max_primes, int max_p)
/* Find all primes in [2, max_p], and store them in the given array.
   Return numprimes on success, negative on failure (e.g., p was not large enough).
   This version will not be used for large values of max_p, so 
   we hash using bytes instead of bits to minimize the arithmetic
   at the expense of memory. This is not meant to be blazing fast at all,
   because we will only every use it for max_p upto a few hundred, maybe,
   and I didn't want to just hardcode the list. 
*/
{ uint8_t  *h;
  int       h_size = max_p+1;
  int       p, p_index, h_index, done=0, p_size;

  if (max_p < 2) {
    return 0;
  }

  h = (uint8_t *)malloc(h_size*sizeof(uint8_t));
  memset(h, 0x00, h_size*sizeof(uint8_t));
  /* h[j] will correspond to the positive integer 2*j+1. */

  p_index=1;
  while (done==0) {
    while ((p_index < h_size) && h[p_index])
      p_index++;
    p = 2*p_index + 1;
    if ((p_index >= h_size)||(p*p>max_p)) {
      done = 1;
    } else {
      h_index = p_index+p;
      while (h_index < h_size) {
        h[h_index] = 1;
        h_index += p;
      }
    }
    p_index++;
  }

  Primes[0] = 2;
  p_size = 1;
  h_index = 1;
  while (2*h_index+1 <= max_p) {
    if (h[h_index]==0) {
      if (p_size < max_primes-1) 
        Primes[p_size++] = 2*h_index+1;
      else {
        free(h);
        return -1;
      }
     }
     h_index++;
  }
  free(h);
  return p_size;
}

/*******************************************************/
int _init_globals()
{ int lim = 16777216; /* 2^{24}. */
  int size = 1077872; /* pi(2^{24}) + 1. */
  static int initialized = 0;

  if (initialized) return 0;
  
  if (_small_primes == NULL) {
    _small_primes = (int *)malloc(size*sizeof(int));
    _num_small_primes = eratosthenes(_small_primes,  size, lim);
  }
  initialized=1;
}

/**************************************************************/
int distinct_prime_divisors(long *p, int p_size, long n)
  /* Return a list of the distinct primes dividing n.
     number of such divisors on success, negative on failure. */
{ long max_p, k, index, d;
  long remain;

  _init_globals();

  max_p = (int)(1 + sqrt((double)n));
  k = 0;
  remain = n;
  index = 0; 
  while ((remain>1) && (index < _num_small_primes)) {
    d = 0;
    while (remain % _small_primes[index] == 0) {
      if (k >= p_size) {
        fprintf(stderr, "distinct_prime_divisors() error: n=%ld, p_size=%d is too small!\n", n, p_size);
        return -1;
      }
      p[k] = _small_primes[index];
      remain /= _small_primes[index];
      d = 1;
    }
    k += d;
    index++;
  }
  if (remain>1) {
    p[k] = remain;
    k++;
  }
  return k;
}


/******************************/
/* Start of ffq_t stuff.      */
/******************************/

/*******************************************************/
void ffq_init(ffq_t *F)
{
  poly_alloc(&F->poly, 16);
  poly_alloc(&F->prim_elt, 16);
  F->p = 2;
  F->poly.coefs[0] = 0;
  F->poly.coefs[1] = 1;
}

/*******************************************************/
void ffq_clear(ffq_t *F)
{
  poly_dealloc(&F->poly);
  poly_dealloc(&F->prim_elt);
  F->p = 0;
}

/*******************************************************/
long  ffq_size(ffq_t *F)
{ int i;
  long q;

  q = 1;
  for (i=1; i<=F->k; i++)
    q *= F->p;
  return q;
}

/*******************************************************/
void ffq_add(poly_t *res, poly_t *op1, poly_t *op2, ffq_t *F)
/* Safe for (output arg) = (input arg). */
{ 
  poly_add(res, op1, op2, F->p);
}

/*******************************************************/
void ffq_sub(poly_t *res, poly_t *op1, poly_t *op2, ffq_t *F)
/* Safe for (output arg) = (input arg). */
{ int i;

  poly_sub(res, op1, op2, F->p);
}

/*******************************************************/
void ffq_mul_c(poly_t *res, poly_t *op1, int64_t c, ffq_t *F)
/* Safe for (output arg) = (input arg). */
{ int i;
  long cc;

  poly_realloc(res, op1->deg);
  cc = ((c%F->p) + F->p) % F->p;
  for (i=0; i<=op1->deg; i++) 
    res->coefs[i] = (op1->coefs[i]*cc)%F->p;
}

/*******************************************************/
void ffq_mul(poly_t *res, poly_t *op1, poly_t *op2, ffq_t *F)
/* Safe for (output arg) = (input arg). */
{ int  d=F->k-1;
  poly_t prod;

  poly_alloc(&prod, op1->deg+op2->deg);

  poly_mul(&prod, op1, op2, F->p);
  poly_mod(res, &prod, &F->poly, F->p);
  poly_dealloc(&prod);
}


/*******************************************************/
int ffq_inv(poly_t *res, poly_t *op1, ffq_t *F)
/* Safe for (output arg) = (input arg). */
{ int  d=F->k, i;
  poly_t g, a, b;

  poly_alloc(&g, d);
  poly_alloc(&a, d);
  poly_alloc(&b, d);

  poly_fixdeg(op1);
  if ((op1->deg==0) && (op1->coefs[0]==0)) 
    return -1;

  poly_gcd_ext(&g, &a, &b, op1, &F->poly, F->p);
  poly_cp(res, &a);
  poly_dealloc(&g);
  poly_dealloc(&a);
  poly_dealloc(&b);
  return 0;
}


/*******************************************************/
void ffq_pow(poly_t *_res, poly_t *op1, int64_t pow, ffq_t *F)
{ int   d=F->k, i;
  poly_t res, p2;

  poly_alloc(&res, d);
  poly_alloc(&p2, d);

  res.coefs[0] = 1;
  res.deg = 0;
  poly_cp(&p2, op1);

  if (pow < 0) {
    ffq_inv(&p2, op1, F);
    pow = -pow;
  } 

  while (pow > 0) {
    if (pow%2) 
      ffq_mul(&res, &res, &p2, F);
    ffq_mul(&p2, &p2, &p2, F);
    pow = pow/2;
  }
  poly_cp(_res, &res);
  poly_dealloc(&res);
  poly_dealloc(&p2);
}

/*******************************************************/
int ffq_is_zero(poly_t *op, ffq_t *F)
{ int i;

  poly_fixdeg(op);
  return ((op->deg==0) && ((op->coefs[0]%F->p)==0));
}

/*******************************************************/
int ffq_is_one(poly_t *op, ffq_t *F)
{ int i;

  poly_fixdeg(op);
  return ((op->deg==0) && ((op->coefs[0]%F->p)==1));
}

/*******************************************************/
char *ffq_str(char *str, poly_t *op, ffq_t *F)
{ 
  return poly_str(str, op);
}

/*******************************************************/
void ffq_rand(poly_t *res, ffq_t *F)
{ int i;

  poly_realloc(res, F->k);
  for (i=0; i<F->k; i++)
    res->coefs[i] = lrand48()%F->p;
  res->deg = F->k-1;
  poly_fixdeg(res);
}

/*******************************************************/
void ffq_set_zero(poly_t *res, ffq_t *F)
{ int i;

  res->deg = 0;
  res->coefs[0] = 0;
}

/*******************************************************/
void ffq_set_one(poly_t *res, ffq_t *F)
{ 
  res->deg = 0;
  res->coefs[0] = 1;
}

/*******************************************************/
void ffq_set_x(poly_t *res, ffq_t *F)
{ int i;

  poly_realloc(res, 1);
  res->deg = 1;
  res->coefs[0] = 0;
  res->coefs[1] = 1;
}
  


/********************************************************/
long ffq_to_enum(poly_t *op, ffq_t *F)
/* Fix an enumeration of 'F', and return the value of 'op' under that enumeration. 
   It will have the properties that 0_F |--> 0 and 1_F |--> 1.
*/
{ long retval, pow;
  int  i;

  if (op->deg >= F->k) 
    poly_mod(op, op, &F->poly, F->p);

  retval = 0;
  pow = 1;
  for (i=0; i<=op->deg; i++) {
    retval += (op->coefs[i]%F->p)*pow;
    pow *= F->p;
  }
  return retval;
}

/********************************************************/
void ffq_from_enum(poly_t *op, long n, ffq_t *F)
/* Inverse of the above. */
{ int i;

  poly_realloc(op, F->k);
  for (i=0; i<F->k; i++) {
    op->coefs[i] = n%F->p;
    n = (n-n%F->p)/F->p;
  }
  op->deg = F->k-1;
  poly_fixdeg(op);
}

/********************************************************/
int ffq_elt_is_primitive(poly_t *op, ffq_t *F)
{ long   q1, cofactor;
  int    i, num_pdivs;
  long   pdivs[128];
  poly_t op_pow;

  poly_alloc(&op_pow, F->k);

  _init_globals();
  q1 = 1;
  for (i=0; i<F->k; i++)
    q1 *= F->p;
  q1 -= 1;

  num_pdivs = distinct_prime_divisors(pdivs, 128, q1);

  ffq_pow(&op_pow, op, q1, F);
  if (ffq_is_one(&op_pow, F) == 0) {
    /* In this case, either F->poly is not irreducible or op is zero. */
    poly_dealloc(&op_pow);
    return 0;
  }

  for (i=0; i<num_pdivs; i++) {
    cofactor = q1/pdivs[i];
    ffq_pow(&op_pow, op, cofactor, F);
    if (ffq_is_one(&op_pow, F)) {
      poly_dealloc(&op_pow);
      return 0;
    }
  }
  poly_dealloc(&op_pow);
  return 1;
}

/********************************************************/
int ffq_create(ffq_t *F, long p, int k)
{ int    i;
  poly_t f;

  poly_alloc(&f, k);
  F->k = k;
  F->p = p;
  f.deg = k;
  if (k==1) {
    f.coefs[0] = 0;
    f.coefs[1] = 1;
  } else {
    do {
      f.coefs[k] = 1; /* Monic. */
      f.coefs[0] = 1 + (lrand48()%(p-1)); /* Nonzero random. */
      for (i=1; i<k; i++)
        f.coefs[i] = lrand48()%p; /* Random. */
    } while (poly_is_irreducible(&f, p)==0);
  }
  poly_cp(&F->poly, &f);
  poly_dealloc(&f);

  /* Now find a primitive element of the form x+c. */
  poly_realloc(&F->prim_elt, k);
  for (i=0; i<=k-1; i++)
    F->prim_elt.coefs[i] = 0;
  
  F->prim_elt.deg = 1;
  if (k>1)
    F->prim_elt.coefs[1] = 1;
  for (i=0; i<p; i++) {
    F->prim_elt.coefs[0] = i;
    if (ffq_elt_is_primitive(&F->prim_elt, F))
      return 0;
  }
  return 0;
}

/********************************************************/


#if 0
/********************************************************/
int main(int argC, char *args[])
{ long  p = 7, f[16];
  int   deg=4, k, res;
  char  str[256];
  ffq_t F;
  double t0, t1;

  ffq_init(&F);

  srand48(time(NULL));


  if (argC == 3) {
    p = atol(args[1]);
    k = atoi(args[2]);
    t0 = sTime();
    res = ffq_create(&F, p, k);
    t1 = sTime();
    if (res==0) {
      printf("F = GF(%ld^%d) created.\n", p, k);
      printf("poly: %s\n", poly_str(str, &F.poly));
      printf("prim: %s\n", ffq_str(str, &F.prim_elt, &F));
    }
    printf("Elapsed time: %1.8lf seconds.\n", t1-t0);
    return 0;
  }

}
#endif
