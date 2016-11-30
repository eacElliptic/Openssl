# Openssl

This repository contains modified versions of Openssl.

These modifications consist of the addition of two algorithms that allow the computation of 'kP' efficiently, where 'k' is a scalar and 'P' a point on an elliptic curve.

The first one is GLV algorithm (see [1]) and the second is an algorithm which is based on euclidian addition chains (EAC), we called EAC-mult.

- Directory 'eac_glv' contains Openssl with EAC-mult and the classical version of GLV (which is not SPA-secure). 

- Directory 'eac_glv_mux' contains Openssl with EAC-mult and a version of GLV which is resistant to data cache timing attacks.

- Directory 'eac_glv_rnaf' contains Openssl with EAC-mult and a version of GLV which is resistant to instruction cache timing attacks.

- Directory 'eac_glv_rnaf_mux' contains Openssl with EAC-mult and a version of GLV which is resistant to instruction/data cache timing attacks, so SPA-secure.


The integration of GLV in Openssl was done by Billy Bob Brumley (see [2, 3]).


Main goal of this repository is to prove that the algorithm based on EAC (EAC-mult) can be an interesting alternative to the SPA-secure version of GLV (see directory 'eac_glv_rnaf_mux'), as EAC-mult is SPA-secure.



[1] R. P. Gallant, R. J. Lambert, and S. A. Vanstone. Faster point multiplication on elliptic curves with efficient endomorphisms. In Advances in Cryptology — CRYPTO, volume 2139 of LNCS, pages 190–200. Springer, 2001.

[2] Billy Bob Brumley. Faster software for fast endomorphisms. 
    (https://eprint.iacr.org/2015/036)

[3] https://github.com/bbbrumley/openssl/tree/ecp_glv
