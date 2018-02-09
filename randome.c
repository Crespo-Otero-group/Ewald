/*  This function computes pseudo-random numbers in the
range from 0 to 1 using one of the methods proposed by
Pierre L'Ecuyer in "Efficient and Portable Combined Random
Number Generators", Communications of the ACM, Vol 31,
pp 742-749, 774, 1988.

PREVIOUS METHODS:
The most commonly employed method is the Lehmer linear
congruential generator (LCG): s(i+1) = [A * s(i) + C] MOD M.
This can be thought of as a roulette wheel with M consecutive
numbers which is spun an amount [A * s(i) + C].  When C = 0, this
becomes the multiplicative linear congruential generator (MLCG),
which has the maximum possible period P = M-1 if M is prime and
A is a primitive element modulo M (i. e. A^N MOD M = 1 is true
only for N = M - 1). By using a technique described by Schrage,
it is possible in 31-bit arithmetic to perform the MLCG when M is
as large as 2^31 - 1 and A is as large as 1.4 x 2^15. However, the
MLCG has the well-known shortcoming that if successive numbers
are plotted in a multi-dimensional space, the points will lie only
on a set of parallel hyperplanes. Such correlations are unacceptable
in Monte Carlo computations where several random numbers are
needed to define each event. One solution (Bayes and Durham) is
the use of a large (typically 97 element) "state table" which is
randomly accessed to scramble and break up these correlations.
An added advantage is that the period is increased to a very large
number. These techniques are described in the book "Numerical
Recipes" and are the basis for many commonly used random number
generators. Of special note is the ANSI C standard generator rand()
which has only 32k values, a period of 32k, and is amply described
in "Numerical Recipes" as being worthless.

FEATURES OF L'ECUYER'S METHOD:
The method combines two 31-bit MLCGs to generate Å2.1 x 10^9
distinct random 31-bit numbers with a period of 2.30584 x 10^18
(Å2^61) and no detectable correlations.
¥ By combining the outputs of two MLCG generators with different
and suitably chosen M and A values, the period is increased to
2^61 and the correlations are reduced to an undetectable level.
¥ The state of the generator can be uniquely described by two
31-bit seed numbers rather than by a seed number and a large
state table.

C PROGRAM CODE:
Written in the C programming language by Stephen E. Derenzo,
Lawrence Berkeley Laboratory, Berkeley, California, May 1, 1992,
who thanks Orin Dahl of LBL for bringing L'Ecuyer's method to his
attention. L'Ecuyer's method is also used at CERN for Monte Carlo
computations in high energy physics.

USE:
long s1, s2;
double ran;
randome (&s1, &s2, &ran);
The function uses the seeds s1 and s2 to compute two new seeds,
which are combined to form the random number ran.

SEED NUMBERS (INPUT AND OUTPUT):
&s1 and &s2 are pointers to two seed numbers that
uniquely represent the state of the generator.
s1 is between 1 and 2,147,483,562 (2^31 - 86) inclusive.
s2 is between 1 and 2,147,483,398 (2^31 - 250) inclusive.

RANDOM NUMBER (OUTPUT):
ran is a pseudo-random number equal to n/2147483563., where
n is a whole number between 1. and  2147483562. inclusive.

MULTIPLE RUNS:
¥ For multiple debugging runs, where the same sequence of
random numbers is desired, the seeds can be started at the
same values (e.g. s1 = 1 and s2 = 1).
¥ For multiple production runs, where independent sequences of
random numbers are desired, it is suggested that the two seeds
be written to a file at the end of each run and then read at
the start of the next run to continue the sequence.

Note: If the generator is started for example with s1 = 1 and
s2 = 1, then 2.3 x 10^18 numbers are available before the
sequence repeats. By cycling only one of the generators before
starting (e. g. starting with s1 = 1 and s2 = 40692), a different
sequence of 2.3 x 10^18 numbers is available. This is how the
full theoretical period of 2^62 can be realized.  */
void randome(long *, long *, double *);
void randome(long *s1ptr, long *s2ptr, double *ranptr) {
if ( (*s1ptr = 40014*(*s1ptr%53668) - 12211*(*s1ptr/53668) ) < 0)
	*s1ptr += 2147483563;
if ( (*s2ptr = 40692*(*s2ptr%52774) - 3791*(*s2ptr/52774) ) < 0)
	*s2ptr += 2147483399;
if ( (*ranptr = *s1ptr - *s2ptr) < 1.) *ranptr += 2147483562.;
*ranptr = *ranptr/2147483563.;
return; }

