# C++ static-length polynomials, with constexpr support

Two classes are provided: `Polynomial<T, N>` and `RationalPolynomial<T, NP, PQ>`.
Template parameter T is the numeric type (e.g., double) and N, NP, and NQ are polynomial degrees.
Degree of polynomial must be specified at compile time.
All data is stored on stack, with no memory allocations.
Basic arithmetic and composition operations are provided.

```
Polynomial p(1, 2, 3);
Polynomial q(4, 5);

// p = (1 + 2*X + 3*X^2)
std::cout << "p = " << p << std::endl;

// p*q = (4 + 13*X + 22*X^2 + 15*X^3)
std::cout << "p*q = " << (p*q) << std::endl;

// p/q = (1 + 2*X + 3*X^2)/(4 + 5*X)
std::cout << "p/q = " << (p/q) << std::endl;

// p.subst(q) = (57 + 130*X + 75*X^2)
std::cout << "p.subst(q) = " << p.subst(q) << std::endl;

// p.derivative = (2 + 6*X)
std::cout << "p.derivative = " << p.derivative() << std::endl;

// p.integral = (0 + 1*X + 1*X^2 + 1*X^3)
std::cout << "p.integral = " << p.integral() << std::endl;
```
