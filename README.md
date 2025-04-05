# GenPayley

This is the code used to perform the calculations in the paper,
[A Mathon-type construction for digraphs and improved lower bounds for Ramsey numbers]
(https://arxiv.org/pdf/2408.04067).

It has no external dependencies and doesn't use anything fancy in the code,
so it should build easily on a reasonable platform just using make.


##    Brief summary:

	```
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
	```

##    About the code:

     [!CAUTION]
	```
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

     * In many places, the code uses an enumeration of the elements of F_q, to make various tasks
       a bit faster (i.e., using lookup tables and such).
	```
