# GenPayley

This is the code used to perform the calculations in the paper,
[A Mathon-type construction for digraphs and improved lower bounds for Ramsey numbers]
(https://arxiv.org/pdf/2408.04067).

It has no external dependencies and doesn't use anything fancy in the code,
so it should build easily on a reasonable platform just using make.


##    Brief summary:


Let $k\ge2$ be an even integer. Let $q$ be a prime power such that
       $$q \equiv k+1 (\mathrm{mod}\\,\\, 2k).$$
Let $S_k$ be the subgroup of $\mathbb{F}_q^*$ of order $(q-1)/k$ consisting of 
the $k$-th power residues. 
Let $G_k(q)$ be the directed graph with 
 * vertices: elements of $\mathbb{F}_q$,
 * an edge $a\to b$ iff $b-a\in S_k$.
(Note: the conditions on q ensure that $-1\not\in S_k$).
Let $K_m(G)$ denote the number of transitive subtournaments of order $m$
contained in a digraph $G$.

The goal is to find the largest $q$ for which $K_m(G_k(q)) = 0$,
for various values of $m$ and $k$. The reason these quantities are of 
interest is that there is a relationship with Ramsey numbers. See the paper
above for details. 

What this code actually does is, for a given $k$ and $q$, find the largest 
transitive subtournament of $G_k(\mathbb{F}_q)$. If that largest trans. sub-tour. has
size $t$, we record the triple $q,k,t$ to file. From that, we can extract what
we want.

The way we find the size of that largest transitive subtournament is
by using Lemma 4.2(d) in the McCarthy, Springfield paper: 
Let $H_k(q)$ be the subgraph of $G_k(q)$ induced by $S_k$. Then let $H_k^1(q)$
be the subgraph of $H_k(q)$ induced by the set of out-neighbors of 1 in
$H_k(q)$. By that lemma, $G_k(q)$ has a transitive subtournament of order $m$
iff $H_k^1(q)$ has a transitive subtournament of order $m-2$. So it suffices
to find the largest transitive subtournament of $H_k^1(q)$. For that,
use a straightforward recursive approach : 

extend_chain( (a1,...,ar), S):
  * Input: a totally ordered chain $a1 < a2 < ... < ar$, and the set
    $S = \\{ x\in H_k^1(q) \\,:\\, a_j < x \\,\\, \mbox{for all} \\,\\, 1\le j\le r\\}$.
  * Set $M \gets r$.
  * For each $x\in S$, set $S' \gets S\cap \mathrm{Successors}(x)$, set 
    $M \gets \mathrm{max}\\{M, \mathrm{extend\\_chain}( (a1,...,ar,x), S')\\}$.
  * return $M$.

We then do extend_chain( (a1), Successors(a1)) for each $a_1\in H_k^1(q)$
to find the maximum length of such a chain.

	

##    About the code:

> [!CAUTION]
> This code is written with speed in mind, not safeness or re-usability.
> Nearly all of the functions here are wildly unsafe, with no argument verification
> or safety checks or anything that would slow things down.

* It's also not thread-safe. I had considered making this multi-threaded,
but since it was one-time use code, we took the simpler approach of launching
multiple separate jobs.

* This version uses a precomputed 0/1 edge matrix in the $k=2$ case.
It is faster for $k=2$, but could run into memory issues for larger primes 
which can only be done for larger $k$ values. So the larger $k$ values just determine
those edges on-the-fly, as needed.

* This version also has an early-abort condition added to extend_chain which speeds 
things up pretty good.

* All of this is using machine precision, for speed. So things will fail for larger
integers. I'm not sure exactly where it will fail but it would be somewhere around
$2^{31}$, which is well above the range of interest (signed mult. overflow). 

* In many places, the code uses an enumeration of the elements of $\mathbb{F}_q$, to make various tasks
a bit faster (i.e., using lookup tables and such).
