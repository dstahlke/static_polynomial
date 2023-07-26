# C++ static-length polynomials, with constexpr support

Two classes are provided: `Polynomial<T, N>` and `RationalPolynomial<T, NP, PQ>`.
Template parameter T is the numeric type (e.g., double) and N, NP, and NQ are polynomial degrees.
Degree of polynomial must be specified at compile time.
All data is stored on stack, with no memory allocations.
Basic arithmetic and composition operations are provided.
