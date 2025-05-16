/* payley.c, 7/31/2024.

   Code for doing some calculations accompanying,
     ``A Mathon-type construction for digraphs and improved lower
       bounds for Ramsey numbers''. D. McCarthy, C. Monico.
       https://arxiv.org/pdf/2408.04067

   ---------------------------------------------------------------------
   Brief summary:
     Let k>=2 be an even integer. Let q be a prime power such that
               q == k+1 (mod 2k).
     Let S_k be the subgroup of F_q^* of order (q-1)/k consisting of 
     the k-th power residues. 
     Let G_k(q) be the directed graph with 
         * vertices: elements of F_q,
         * an edge a-->b iff b-a\in S_k.
     (Note: the conditions on q ensure that -1\not\in S_k).
     Let K_m(G) denote the number of transitive subtournaments of order m
     contained in a digraph G.
     
       The goal is to find the largest q for which K_m(G_k(q)) = 0,
     for various values of m and k. The reason these quantities are of 
     interest is that there is a relationship with Ramsey numbers. See the paper
     above for details. 

     What this code actually does is, for a given k and q, find the largest 
     transitive subtournament of G_k(F_q). If that largest trans. sub-tour. has
     size t, we record the triple q,k,t to file. From that, we can extract what
     we want.

   ---------------------------------------------------------------------
   About the code:
     * This code is written with speed in mind, not safeness or re-usability.
       Nearly all of the functions here are wildly unsafe, with no argument verification
       or safety checks or anything that would slow things down.

     * It's also not thread-safe. I had considered making this multi-threaded,
       but since it was one-time use code, we took the simpler approach of launching
       multiple separate jobs.

     * This version which uses a precomputed 0/1 edge matrix in the k=2 case.
       It is faster for k=2, but could run into memory issues for larger primes 
       which can only be done for larger k values. So the larger k values just determine
       those edges on-the-fly, as needed.

     * This version also has an early-abort condition added to extend_chain which speeds 
       things up pretty good.

     * All of this is using machine precision, for speed. So things will fail for larger
       integers. I'm not sure exactly where it will fail but it would be somewhere around
       2^{31}, which is well above the range of interest (signed mult. overflow). 

     * We often used enumerations of the elements of F_q, to make various tasks
       a bit faster (i.e., using lookup tables and such).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ffq.h"
#include "common.h"

#ifndef BIT
#define BIT(_i) ((0x00000001)<<(_i))
#endif

#define USAGE "[OPTIONS]\n"\
"OPTIONS:\n"\
"---------------------------------------------------------------\n"\
"  -q0 <int> : first prime power to check (default: 2).\n"\
"  -q1 <int> : last prime power to check (default: 512).\n"\
"  -k  <int> : k value to use (default: 512).\n"\
"  -numqonly : just report the number of q values that would be \n"\
"              checked without actually doing the calculations.\n\n"\
" Raw results will be appended to the file 'results.k<val>'. \n"\
" Each line in that file is a triple q,k,t, where:\n"\
" q is a prime power, q==k+1 (mod 2k), and t is the size of the largest\n"\
" transitive subtournament of G_k(q).\n"


/* 0-1 matrices: */
typedef struct {
  uint32_t *data;
  int       rows, cols;
  int       stride;
} bitmat_t;

/* KISS: digraphs will be simply represented as a list of edges. */
typedef int64_t edge_t[2];


/* A `housekeeping' structure: these things are needed
   by multiple functions, and we don't want to be passing
   10 arguments trillions of times. */
typedef struct {
  ffq_t   F;  /* A finite field F_q. */
  int64_t q;  /* Order of the field F_q. */
  int     k;  /* q == k+1 (mod 2k). */
  poly_t *Sk; /* S_k, k-th power residues in F_q. */
  int     Sk_size;
  int64_t *Ske; /* Enumerations (non-neg. integer representations) of elements Sk, sorted. */
  int8_t  *H;   /* Characteristic function (hash) of Ske; H[j]=1 iff j is in Ske. */
  int64_t *out_nbrs_1; /* Out neighbors of 1 (enumerated values). */
  int      on1_size;   /* Size of the above. */
  edge_t  *H1edges;    /* Edges in the subgraph H_k^1, as described in the paper. */
  int      num_H1edges;
  int      *on1_successors; 
           /* on1_successors[i] = n means that the out edges of out_nbrs_1[i] in H_k^1
              are precisely H1edges[n], H1edges[n+1],..., H1edges[m-1], where m=on1_successors[i+1].
           */
  bitmat_t LM; /* Entry (i,j)=1 iff poly_from_enum(i) --> poly_from_enum(j), in H_k^1. */
} gen_payley_prob_t;


/****** Globals. Sorry-not-sorry. ******/
int  *small_primes; /* Will be initialized by init_globals(), based on the next two values. */
int   max_primes=1077871, max_p=1<<24;


/******        Prototypes.        ******/
int enum_is_less(int64_t a, int64_t b, gen_payley_prob_t *P);


/*************************************************/
int bitmat_init(bitmat_t *M, int rows, int cols)
  /* Initialize the bitmat_t M, allocating and setting 
     entries to zero. */
{ int stride;

  stride = cols/32 + ((cols%32)!=0);
  if (!(M->data = (uint32_t *)malloc(rows*stride*sizeof(int32_t))))
    return -1;
  memset(M->data, 0x00, rows*stride*sizeof(int32_t));
  M->rows = rows;
  M->cols = cols;
  M->stride = stride;
  return 0;
}

/*************************************************/
void bitmat_clear(bitmat_t *M)
  /* Deallocate. */
{
  if (M->data)
    free(M->data);
  M->data = NULL;
}

/*************************************************/
int init_globals()
{ static int initialized=0;

  if (initialized) return 0;

  small_primes = (int *)malloc(max_primes * sizeof(int));
  eratosthenes(small_primes, max_primes, max_p);
  initialized = 1;
  return 0;
}

/************************************************************/
int get_kpow_residues(poly_t **res, ffq_t *F, int k)
  /* Compute S_k, as described above. */
{ int     size, i;
  int64_t t;
  poly_t  w, wi;

  /* How many are there? */
  t = ffq_size(F) - 1;
  size = (int)(t/k);
  
  /* Allocate for them: */
  *res = (poly_t *)malloc(size*sizeof(poly_t));

  poly_alloc(&w, F->k);
  poly_alloc(&wi, F->k);
  ffq_pow(&w, &F->prim_elt, k, F);
  poly_cp(&wi, &w);
  for (i=0; i<size; i++) {
    poly_alloc(*res+i, F->k);
    poly_cp(*res+i, &wi);
    ffq_mul(&wi, &wi, &w, F);
  }
  poly_dealloc(&w);
  poly_dealloc(&wi);
  return size;
}

/************************************************************/
int cmp_int64(const void *a, const void *b)
  /* Comparison function, used for sorting. */
{ int64_t *A = (int64_t *)a, *B = (int64_t *)b;

  if (*A < *B) return -1;
  if (*A > *B) return 1;
  return 0;
}

/************************************************************/
int get_kpow_residues_enum(int64_t **res, ffq_t *F, int k)
  /* Get enumerations of elements of S_k. */
{ poly_t *Sk;
  int     size, i;

  size = get_kpow_residues(&Sk, F, 2);
  *res = (int64_t *)malloc(size*sizeof(int64_t));
  if (*res == NULL) return -1;
  for (i=0; i<size; i++)
    (*res)[i] = ffq_to_enum(&Sk[i], F);
  free(Sk);
  /* Sort them. */
  qsort(*res, size, sizeof(int64_t), cmp_int64);
  return size;
}

/************************************************************/
int get_kpow_residues_hash(int8_t **H, ffq_t *F, int k)
{ poly_t *Sk;
  int     size, i, q, index;
  
  size = get_kpow_residues(&Sk, F, 2);
  q = (int)ffq_size(F);
  *H = (int8_t *)malloc(q*sizeof(int8_t));
  if (*H == NULL) return -1;
  memset(*H, 0x00, q*sizeof(int8_t));
  for (i=0; i<size; i++) {
    index = (int)ffq_to_enum(&Sk[i], F);
    (*H)[index] = 1;
  }
  free(Sk);
  return 0;
}

/************************************************************/
int out_neighbors_of_1(int64_t **res, int8_t *Sk_hash, ffq_t *F)
{ int    size, i, j, q, t, res_size;
  poly_t one, this_residue, sum;

  q = (int)ffq_size(F);
  size = (q-1)/2; /* An upper bound. */
  *res = (int64_t *)malloc(size*sizeof(int64_t));
  if (*res == NULL) return -1;
  
  poly_alloc(&one, F->k);
  poly_alloc(&this_residue, F->k);
  poly_alloc(&sum, F->k);

  ffq_set_one(&one, F);
  res_size = 0;
  for (i=1; i<q; i++) {
    if (Sk_hash[i]) {
      /* ffq_to_enum(i) is a residue. Check if adding 1 also yields a residue. */
      ffq_from_enum(&this_residue, (int64_t)i, F);
      ffq_add(&sum, &this_residue, &one, F);     
      t = (int)ffq_to_enum(&sum, F);
      if (Sk_hash[t]) {
        (*res)[res_size] = t;
        res_size++;
      }
    }
  }
  poly_dealloc(&one);
  poly_dealloc(&this_residue);
  poly_dealloc(&sum);

  qsort(*res, res_size, sizeof(int64_t), cmp_int64);
  return res_size;
}

/**********************************************************************/
int generate_H1k(gen_payley_prob_t *P)
/* Given a gen_payley_prob_t for which all of the previous fields
   are completed, determine the edges in the induced subgraph H_k^1.
   Let H_k be the subgraph induced by S_k.
   Then H_k^1 is the subgraph of H_k induced by the out-neighbors of 1.
   Thus, (1) the vertex set is out_nbrs_1 = { v : v in S_k and v-1 in S_k}.
         (2) There is an edge from v1 to v2 iff v2-v1 is in Sk. 
*/
{ int     bound, i, j, n;
  poly_t  v1, v2, diff;
  int     de;

  bound = P->on1_size * (P->on1_size - 1) / 2;
  P->H1edges = (edge_t *)malloc(bound * sizeof(edge_t));
  P->on1_successors = (int *)malloc((1+P->on1_size)*sizeof(int));
  P->num_H1edges = 0;
  
  poly_alloc(&v1, P->F.k);
  poly_alloc(&v2, P->F.k);
  poly_alloc(&diff, P->F.k);

  n = 0;
  for (i=0; i<P->on1_size; i++) {
    P->on1_successors[i] = n;
    ffq_from_enum(&v1, P->out_nbrs_1[i], &(P->F));
    for (j=0; j<P->on1_size; j++) {
      ffq_from_enum(&v2, P->out_nbrs_1[j], &(P->F));
      ffq_sub(&diff, &v2, &v1, &(P->F));
      de = (int)ffq_to_enum(&diff, &(P->F));
      if (P->H[de]) {
        if (n >= bound) {
          /* This should be impossible. */
          printf("generate_H1k() error! n=%d, bound=%d. Bailing...\n", n, bound);
          return -1;
        }
        /* There is an edge from v1 to v2. */
        P->H1edges[n][0] = P->out_nbrs_1[i];
        P->H1edges[n][1] = P->out_nbrs_1[j];
        n++;
      }
    }
  }
  P->num_H1edges = n;
  P->on1_successors[i] = n;

  poly_dealloc(&v1);
  poly_dealloc(&v2);
  poly_dealloc(&diff);
}

/**********************************************************************/
void dealloc_gpp(gen_payley_prob_t *P)
{ int i;

  for (i=0; i<P->Sk_size; i++)
    poly_dealloc(P->Sk+i);
  free(P->Sk); P->Sk = NULL;
  free(P->Ske); P->Ske = NULL;
  free(P->H); P->H = NULL;
  free(P->out_nbrs_1); P->out_nbrs_1 = NULL;
  free(P->on1_successors);  P->on1_successors = NULL;

  if (P->k==2)
    bitmat_clear(&P->LM);
}

/************************************************************************/
int init_instance(gen_payley_prob_t *P, int field_p, int field_e, int k)
/* All fields of P are assumed to be uninitialized. So clear them between re-uses! */
{ int i, index;

  ffq_init(&(P->F));
  ffq_create(&(P->F), field_p, field_e);
  P->q = ffq_size(&(P->F));
  if ((P->q%(2*k)) != k+1) {
    printf("init_instance() error : q=%ld == %ld (mod %ld). Bailing...\n",
            P->q, (P->q)%((long)2*k), (long)2*k);
    return -1;
  }
  P->k = k;
  P->Sk_size = get_kpow_residues(&(P->Sk), &(P->F), P->k);

  /* Generate the enumerations of Sk elements, and sort them: */
  P->Ske = (int64_t *)malloc(P->Sk_size*sizeof(int64_t));
  for (i=0; i<P->Sk_size; i++)
    P->Ske[i] = ffq_to_enum(P->Sk+i, &(P->F));
  qsort(P->Ske, P->Sk_size, sizeof(int64_t), cmp_int64);

  /* The hash table/characteristic function of Ske: */
  P->H = (int8_t *)malloc(P->q*sizeof(int8_t));
  memset(P->H, 0x00, P->q*sizeof(int8_t));
  for (i=0; i<P->Sk_size; i++) {
    index = P->Ske[i];
    P->H[index] = 1;
  }
  
  P->on1_size = out_neighbors_of_1(&P->out_nbrs_1, P->H, &(P->F));

  generate_H1k(P);

  if (k==2) {
    bitmat_init(&(P->LM), P->q, P->q);
    for (int i=1; i<P->q; i++)
      for (int j=1; j<P->q; j++) 
        if (enum_is_less(i, j, P))
          P->LM.data[P->LM.stride*i + (j/32)] ^= BIT(j%32);
  } 

  return 0;
}

/*********************************************************/
int enum_is_less(int64_t a, int64_t b, gen_payley_prob_t *P)
{ int de=0, pow;
  int64_t p = P->F.p;

  pow = 1;
  for (int i=0; i<P->F.k; i++) {
    de += ((p+(b-a)%p)%p)*pow;
    a = a/p;
    b = b/p;
    pow *= p;
  } 
  return P->H[de];
}

/*********************************************************/
int extend_chain(int64_t *chain, int chain_size, int max_chain_size, 
                 int64_t *possible_successors, int num_possible_successors,
                 int      current_max,
                 gen_payley_prob_t *P)
/* Consider the relation < defined by a < b iff there is an edge from a --> b.
   Given a totally ordered chain:
      chain[0] < chain[1] < ... < chain[chain_size-1], with chain_size>=1,
   attempt to extend it to such a totally ordered chain of larger size.
   possible_successors: the set {t : c < t for all c in chain}.
   Return value: the maximum length to which this chain is extendible.
*/
{ int64_t *next_successors, t, poss;
  int      i, j, n, max, len, r, M;

  if ((num_possible_successors == 0)||(chain_size >= max_chain_size))
    return chain_size; /* It cannot be extended further. */

  /* An early-abort condition. In this case, we can't possibly beat the
     longest one which has been so far found. */
  if (chain_size + num_possible_successors < current_max)
    return chain_size; 

  next_successors = (int64_t *)malloc((num_possible_successors-1)*sizeof(int64_t));
  max  = chain_size;
  for (i=0; i<num_possible_successors; i++) {
    t = possible_successors[i];
    /* Try this one. */
    chain[chain_size] = t;
    /* Set possible_successors <-- (possible_successors-{t}) \cap successors(t). */
    for (n=j=0; j<num_possible_successors; j++) {
      poss = possible_successors[j];
      if ((poss != t)) {
        if ( ((P->k==2) && (P->LM.data[P->LM.stride*t + (poss/32)] & BIT(poss%32)))  ||
             ((P->k>2) && (enum_is_less(t, poss, P)))) {
          next_successors[n++] = poss;
        }
      }
    }
    M = MAX(max, current_max);
    len = extend_chain(chain, chain_size+1, max_chain_size, next_successors, n, M, P);
    max = MAX(max, len);
  }
  free(next_successors);
  return max;
}

/*********************************************************/
int max_trans_subtourney(gen_payley_prob_t *P)
{ int      max_chain_size = 32;
  int64_t  chain[max_chain_size], *possible_successors;
  int      chain_size, num_possible_successors, i, j, n0, n1;
  int      len, max=0;

  /* This is using Lemma 4.2(c) in the McCarthy, Springfield paper 
     to find the maximum length transitive subtournament of G_k(q);
     let H_k(q) be the subgraph of G_k(q) induced by S_k. Then let H_k^1(q)
     be the subgraph of H_k(q) induced by the set of out-neighbors of 1 in
     H_k(q). By that lemma, G_k(q) has a transitive subtournament of order m
     iff H_k^1(q) has a transitive subtournament of order m-2.

     Beyond that lemma, there is nothing clever going on. It's essentially
     a brute force searching for a maximal 'totally ordered' subset
     under the 'order' a<b iff a-->b. Here, we consider each possible starting
     point for such a subset (chain). We then call a function which will recursively
     consider all possible extensions of it.
  */
     
  possible_successors = (int64_t *)malloc(P->on1_size*sizeof(int64_t));

  for (i=0; i<P->on1_size; i++) {
    printTmp("q=%d ... %2.1lf%% done (m=%d)", P->q, 100.0*i/(double)P->on1_size, max);
    chain[0] = P->out_nbrs_1[i];
    num_possible_successors = 0;

    n0 = P->on1_successors[i];    
    n1 = P->on1_successors[i+1];
    for (j=n0; j<n1; j++) 
      possible_successors[num_possible_successors++] = P->H1edges[j][1];

    len = extend_chain(chain, 1, max_chain_size, 
                       possible_successors, num_possible_successors, max, P);
    max = MAX(max, len);
  }
  printf("\n");
  free(possible_successors);
  return max+2;
}

/*********************************************************/
void print_payley_info(FILE *fp, gen_payley_prob_t *P)
{ int i;
  char str[256];

  fprintf(fp, "F_%ld = F_%ld[x] / <%s>\n", P->q, P->F.p, poly_str(str, &(P->F.poly)));
  for (i=0; i<P->Sk_size; i++)
    fprintf(fp, "residue %d: %s, enum:%ld\n", i, poly_str(str, P->Sk+i), ffq_to_enum(P->Sk+i, &(P->F))); 

  fprintf(fp, "Just enumerations:\n");
  for (i=0; i<P->Sk_size; i++) 
    fprintf(fp, "%d, ", P->Ske[i]);
  fprintf(fp, "\n");

  fprintf(fp, "Out neighbors of 1:\n");
  for (i=0; i<P->on1_size; i++) 
    fprintf(fp, "%d, ", P->out_nbrs_1[i]);
  fprintf(fp, "\n");

  fprintf(fp, "Edges of H_k^1  (num_H1edges=%d):\n", P->num_H1edges);
  for (i=0; i<P->num_H1edges; i++) {
    fprintf(fp, "(%ld, %ld)  ", P->H1edges[i][0], P->H1edges[i][1]);
    if ((i%8==7) || (i == P->num_H1edges-1)) 
      fprintf(fp, "\n");
  }
}

/*********************************************************/
int factor_prime_power(long *p, int *e, long q)
{ int i;

  *p = 0;
  *e = 0;
  i = 0;
  while ((small_primes[i]*small_primes[i] < q) && (q%small_primes[i]))
    i++;
  if (q%small_primes[i]==0) {
    *p = small_primes[i];
    while ((q%small_primes[i])==0) {
      *e += 1;
      q = q/small_primes[i];
    }
    if (q != 1) return -1; /* Not a power of small_primes[i] ! */
    return 0;
  }
  *p = q;
  *e = 1;
  return 0;
}

/*********************************************************/
long next_prime_power(long n)
/* Return the least prime power q >= n.
   This is not particularly efficient, but it's not a bottleneck. 
*/
{ long p[32], q;

  q = MAX(2, n);
  
  while (distinct_prime_divisors(p, 32, (long)q)>1) 
    q++;
  return q;
}

/*********************************************************/
long next_admissible_q(long n, int k)
/* Return the least q>=n such that:
   (i) q is a prime power,
   (ii) q == k+1 (mod 2k).
   Again, this is not particularly efficient, but it's also
   not a bottleneck.
*/
{ long q, p[32];

  q = n + (2*k - (n%(2*k)) + k+1)%(2*k);
  while (distinct_prime_divisors(p, 32, (long)q)>1) 
    q += 2*k;
  return q;
}

/*********************************************************/
int main(int argC, char *args[])
{ int               i, e, k=2, m;
  long              q0=2, q1=512, q, p;
  char              str[256];
  gen_payley_prob_t P;
  double            t0, t1, t0_job, t1_job;
  char              filename[256];
  FILE             *fp;
  int               NumCasesOnly=0, numcases=0, count=0;

  init_globals();
  /* Parse command-line args. */
  for (i=1; i<argC; i++) {
    if (strcmp(args[i], "-q0")==0) {
      if ((++i) < argC)
        q0 = atol(args[i]);
    } else if (strcmp(args[i], "-q1")==0) {
      if ((++i) < argC)
        q1 = atol(args[i]);
    } else if (strcmp(args[i], "-k")==0) {
      if ((++i) < argC)
        k = atoi(args[i]);
    } else if (strcmp(args[i], "-numqonly")==0) {
     NumCasesOnly=1;
    } else if (strcmp(args[i], "-help")==0) {
      printf("%s %s\n", args[0], USAGE);
      return 0;
    }
  }

  sprintf(filename, "results.k%d", k);
  
  /* Determine the number of q values to be tested. */
  q = next_admissible_q(q0, k);
  printf("Counting cases...\n");
  while (q <= q1) {
    numcases++;
    q = next_admissible_q(q+2*k, k);
  }
  printf("%d values of q in [%ld, %ld] with q==%d (mod %d).\n", 
         numcases, q0, q1, k+1, 2*k);
  if (NumCasesOnly) {
    return 0;
  }

  /* Do the actual work. */
  q = next_admissible_q(q0, k);
  t0_job = sTime();
  while (q <= q1) {
    if (factor_prime_power(&p, &e, q)) {
      printf("Error: q=%d is not a prime power!\n", q);
      return -1;
    }
    if (init_instance(&P, p, e, k)) {
      printf("init_instance() failed. Bailing...\n");
      return -1;
    }

    count++;
    printf("q=%ld (%d of %d), k=%d ...\n", q, count, numcases, k); 
    t0 = sTime();
    m = max_trans_subtourney(&P);
    t1 = sTime();
    printf("q=%ld, k=%d, max. transitive subtournament length: %d\n", q, k, m);
    printf("elapsed time: %1.4lf\n", t1-t0);
    if ((fp = fopen(filename, "a"))) {
      fprintf(fp, "%d, %d, %d\n", q, k, m);
      fclose(fp);
    }
    dealloc_gpp(&P);
    q = next_admissible_q(q+2*k, k);
  }
  t1_job = sTime();
  printf("Total job time: %1.6lf seconds.\n", t1_job - t0_job);
}
