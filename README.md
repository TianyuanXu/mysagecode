# mysagecode

This folder contains data of left and 2-sided Kazhdan-Lusztig cells in all finite Coxeter groups of rank less than or equal to 6, as well as Sage code for computing

1. products of Kazhdan-Lusztig basis elements in Hecke algebras;
2. Lusztig's a-function on fully commutative elements of Coxeter groups;
3. the cacti group action on Coxeter groups;
4. certain products in Lusztig's asymptotic Hecke algebras

The algorithm for 1 uses Sage's implementation of Fokko du Cloux's Coxeter3
package for fast computations of mu coefficients; we expect this makes the
product computation faster than the current algorithm in Sage.

Please email tianyuan@uoregon.edu for comments or questions.

