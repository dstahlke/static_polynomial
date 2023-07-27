// g++ -Wall -Wextra tests.cpp -lgtest -o tests
#include "static_polynomial.h"
#include "gtest/gtest.h"

using staticpolynomial::Polynomial;

TEST(staticpolynomial, Basic) {
    Polynomial p(1, 2, 3);
    Polynomial q(4, 5);
    auto const X = staticpolynomial::X;

    std::cout << "p = " << p << std::endl;
    std::cout << "p.subst(q) = " << p.subst(q) << std::endl;
    std::cout << "p*q = " << (p*q) << std::endl;
    std::cout << "p/q = " << (p/q) << std::endl;
    std::cout << "p.derivative = " << p.derivative() << std::endl;
    std::cout << "p.integral = " << p.integral() << std::endl;

    EXPECT_DOUBLE_EQ(p(0.2), 1 + 2*0.2 + 3*0.2*0.2);
    EXPECT_DOUBLE_EQ(p(  7), 1 + 2*7 + 3*7*7);
    EXPECT_DOUBLE_EQ(q(0.2), 4 + 5*0.2);
    EXPECT_DOUBLE_EQ(q(  7), 4 + 5*7);
    EXPECT_DOUBLE_EQ(p.subst(q)(1.23), p(q(1.23)));
    EXPECT_DOUBLE_EQ(p.integral().derivative()(1.23), p(1.23));

    EXPECT_DOUBLE_EQ((p+3)(0.2), p(0.2) + 3);
    EXPECT_DOUBLE_EQ((p-3)(0.2), p(0.2) - 3);
    EXPECT_DOUBLE_EQ((p*3)(0.2), p(0.2) * 3);
    EXPECT_DOUBLE_EQ((p/3)(0.2), p(0.2) / 3);

    EXPECT_DOUBLE_EQ((p+X)(0.2), p(0.2) + 0.2);
    EXPECT_DOUBLE_EQ((p-X)(0.2), p(0.2) - 0.2);
    EXPECT_DOUBLE_EQ((p*X)(0.2), p(0.2) * 0.2);
    EXPECT_DOUBLE_EQ((p/X)(0.2), p(0.2) / 0.2);

    EXPECT_DOUBLE_EQ((p+q)(0.2), p(0.2) + q(0.2));
    EXPECT_DOUBLE_EQ((p-q)(0.2), p(0.2) - q(0.2));
    EXPECT_DOUBLE_EQ((p*q)(0.2), p(0.2) * q(0.2));
    EXPECT_DOUBLE_EQ((p/q)(0.2), p(0.2) / q(0.2));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
