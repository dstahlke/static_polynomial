#include "static_polynomial.h"

using staticpolynomial::Polynomial;

int main() {
    Polynomial p(1, 0, 2);
    Polynomial q(3, 4);
    std::cout << "p = " << p << std::endl;
    std::cout << "p(0.2) = " << p(0.2) << std::endl;
    std::cout << "p(2) = " << p(2) << std::endl;
    std::cout << "p.subst(q) = " << p.subst(q) << std::endl;
    std::cout << "p/q = " << (p/q) << std::endl;
}
