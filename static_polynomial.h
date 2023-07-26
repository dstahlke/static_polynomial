// Copyright 2023 Dan Stahlke
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
// associated documentation files (the “Software”), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge, publish, distribute,
// sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or
// substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
// NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <array>
#include <cassert>
#include <type_traits>
#include <cmath>

namespace staticpolynomial {

// https://stackoverflow.com/a/47563100/1048959
template <class F, std::size_t... Is>
constexpr void for_constexpr(F func, std::index_sequence<Is...>) {
  (func(std::integral_constant<size_t, Is>{}), ...);
}

template <std::size_t N, typename F>
constexpr void for_constexpr(F func)
{
  for_constexpr(func, std::make_index_sequence<N>());
}

template <typename T, size_t NP, size_t NQ>
class RationalPolynomial;

template <typename T, size_t N>
class Polynomial {
public:
    std::array<T, N+1> coeffs;

    static const Polynomial<T, 1> X;

    constexpr Polynomial() : coeffs{} { }

    template <typename... U>
    constexpr Polynomial(U... coeffs_) :
        coeffs{{static_cast<T>(coeffs_)...}}
    { }

    template <size_t M>
    constexpr Polynomial(Polynomial<T, M> const &o) : coeffs{} {
        static_assert(M <= N);
        (*this) += o;
    }

    constexpr Polynomial(T s) : coeffs{} {
        coeffs[0] = s;
    }

    constexpr T const &operator[](size_t i) const { return coeffs[i]; }
    constexpr T       &operator[](size_t i)       { return coeffs[i]; }

    template <size_t M>
    constexpr Polynomial<T, N> &operator+=(Polynomial<T, M> const &o) {
        static_assert(M <= N);
        for(size_t i=0; i<=M; ++i)
            (*this)[i] += o[i];
        return *this;
    }

    template <size_t M>
    constexpr Polynomial<T, std::max(M, N)> operator+(Polynomial<T, M> const &o) const {
        if constexpr(N > M) {
            return Polynomial<T, N>{*this} += o;
        } else {
            return o + *this;
        }
    }

    constexpr Polynomial<T, N> &operator+=(T s) {
        (*this)[0] += s;
        return *this;
    }

    constexpr Polynomial<T, N> operator+(T s) const {
        return Polynomial<T, N>{*this} += s;
    }

    constexpr Polynomial<T, N> &operator-=(Polynomial<T, N> const &o) {
        return *this += -o;
    }

    constexpr Polynomial<T, N> operator-(Polynomial<T, N> const &o) const {
        return Polynomial<T, N>{*this} -= o;
    }

    constexpr Polynomial<T, N> operator-() const {
        auto out = *this;
        for(auto &v : out.coeffs)
            v = -v;
        return out;
    }

    constexpr Polynomial<T, N> &operator*=(T o) {
        for(size_t i=0; i<=N; ++i)
            (*this)[i] *= o;
        return *this;
    }

    constexpr Polynomial<T, N> operator*(T o) const {
        return Polynomial<T, N>{*this} *= o;
    }

    constexpr Polynomial<T, N> &operator/=(T o) {
        for(size_t i=0; i<=N; ++i)
            (*this)[i] /= o;
        return *this;
    }

    constexpr Polynomial<T, N> operator/(T o) const {
        return Polynomial<T, N>{*this} /= o;
    }

    template <size_t M>
    constexpr Polynomial<T, N+M> operator*(Polynomial<T, M> const &o) const {
        Polynomial<T, N+M> out{};
        for(size_t i=0; i<=M; ++i)
            for(size_t j=0; j<=N; ++j)
                out[i+j] += (*this)[j] * o[i];
        return out;
    }

    constexpr size_t leading_term() const {
        for(size_t i=N; i; --i) {
            if((*this)[i])
                return i;
        }
        return 0;
    }

    template <size_t M>
    constexpr Polynomial<T, N> &operator*=(Polynomial<T, M> const &o) {
        size_t n = leading_term();
        size_t m = o.leading_term();
        if(n+m > N)
            throw std::out_of_range{"polynomial order too high"};

        Polynomial<T, N> out{};
        for(size_t i=0; i<=m; ++i)
            for(size_t j=0; j<=n; ++j)
                out[i+j] += (*this)[j] * o[i];
        return *this = out;
    }

    template <size_t E>
    constexpr auto pow() const {
        if constexpr(E == 0) {
            return Polynomial<T, 0>{1};
        } else if constexpr(E == 1) {
            return *this;
        } else if constexpr(E % 2) {
            return *this * pow<E-1>();
        } else {
            return (*this * *this).template pow<E/2>();
        }
    }

    template <size_t M>
    constexpr auto subst(Polynomial<T, M> const &f) const {
        Polynomial<T, N * M> outp{};
        for_constexpr<N+1>([&](auto ic) {
            constexpr size_t e = ic.value;
            outp += f.template pow<e>() * (*this)[e];
        });
        return outp;
    }

    template <size_t MP, size_t MQ>
    constexpr auto subst(RationalPolynomial<T, MP, MQ> const &f) const {
        Polynomial<T, N * std::max(MP, MQ)> outp{};
        for_constexpr<N+1>([&](auto ic) {
            constexpr size_t e = ic.value;
            outp += f.q.template pow<N-e>() * f.p.template pow<e>() * (*this)[e];
        });
        return outp / f.q.template pow<N>();
    }

    // Evaluate polynomial at point.
    template <typename U>
    constexpr auto operator()(U const &x) const {
        // result type
        using R = decltype(std::declval<T&>() * std::declval<U&>());
        R accum{};
        for_constexpr<N+1>([&](auto ic) {
            if(ic)
                accum *= x;
            accum += (*this)[N-ic];
        });
        return accum;
    }
};

template <typename T, size_t N>
constexpr Polynomial<T, 1> Polynomial<T, N>::X{0, 1};

template <typename... U>
Polynomial(U... coeffs) -> Polynomial<double, sizeof...(coeffs)-1>;

template <typename U, typename T, size_t N>
std::enable_if_t<std::is_convertible_v<U, T>, Polynomial<T, N>>
constexpr operator*(U s, Polynomial<T, N> const &p) {
    return p * s;
}

//template <size_t N=1, typename T=double>
//static constexpr Polynomial<T, N> X = Polynomial<T, 1>{T{0}, T{1}}.template pow<N>();
static constexpr Polynomial<double, 1> X = Polynomial<double, 1>{0, 1};

template <typename T, size_t N>
std::ostream &operator<<(std::ostream &os, Polynomial<T, N> const &p) {
    os << "(";
    for(size_t i=0; i<=N; ++i) {
        //if(i) os << ",";
        //os << p[i];
        if(i) os << " + ";
        os << p[i];
        if(i == 1)
            os << "*X";
        else if(i > 1)
            os << "*X^" << i;
    }
    return os << ")";
}

template <typename T, size_t NP, size_t NQ>
class RationalPolynomial {
public:
    Polynomial<T, NP> p;
    Polynomial<T, NQ> q;

    constexpr RationalPolynomial() : p{}, q{T{1}} { }
    constexpr RationalPolynomial(
        Polynomial<T, NP> p_, Polynomial<T, NQ> q_) : p{p_}, q{q_} { }

    template <size_t MP, size_t MQ>
    constexpr auto subst(RationalPolynomial<T, MP, MQ> const &f) const {
        auto outp = p.subst(f).p;
        auto outq = q.subst(f).p;
        if constexpr(NP > NQ) {
            return outp / (outq * f.q.template pow<NP - NQ>());
        } else {
            return (outp * f.q.template pow<NQ - NP>()) / outq;
        }
    }

    template <size_t M>
    constexpr auto subst(Polynomial<T, M> const &f) const {
        auto outp = p.subst(f);
        auto outq = q.subst(f);
        return outp / outq;
    }

    RationalPolynomial<T, NP, NQ> &operator*=(T s) {
        p *= s;
        return *this;
    }

    RationalPolynomial<T, NP, NQ> operator*(T s) const {
        return RationalPolynomial<T, NP, NQ>(*this) *= s;
    }

    RationalPolynomial<T, NP, NQ> &operator/=(T s) {
        p /= s;
        return *this;
    }

    RationalPolynomial<T, NP, NQ> operator/(T s) const {
        return RationalPolynomial<T, NP, NQ>(*this) /= s;
    }

    template <size_t MP, size_t MQ>
    constexpr auto operator*(RationalPolynomial<T, MP, MQ> const &f) const {
        return (p * f.p) / (q * f.q);
    }

    template <size_t MP, size_t MQ>
    constexpr auto operator/(RationalPolynomial<T, MP, MQ> const &f) const {
        return (p * f.q) / (q * f.p);
    }

    template <size_t M>
    constexpr auto operator*(Polynomial<T, M> const &f) const {
        return (p * f) / q;
    }

    template <size_t M>
    constexpr auto operator/(Polynomial<T, M> const &f) const {
        return p / (q * f);
    }

    constexpr auto normalize() const {
        return RationalPolynomial<T, NP, NQ>{p / q[NQ], q / q[NQ]};
    }
};

template <typename U, typename T, size_t NP, size_t NQ>
std::enable_if_t<std::is_convertible_v<U, T>, RationalPolynomial<T, NP, NQ>>
constexpr operator*(U s, RationalPolynomial<T, NP, NQ> const &p) {
    return p * s;
}

template <typename T, size_t NP, size_t NQ>
constexpr RationalPolynomial<T, NP, NQ> operator/(Polynomial<T, NP> const &p, Polynomial<T, NQ> const &q) {
    return {p, q};
}

template <typename U, typename T, size_t N>
std::enable_if_t<std::is_convertible_v<U, T>, RationalPolynomial<T, 0, N>>
constexpr operator/(U s, Polynomial<T, N> const &p) {
    return Polynomial<T, 0>{s} / p;
}

template <typename U, typename T, size_t NP, size_t NQ>
std::enable_if_t<std::is_convertible_v<U, T>, RationalPolynomial<T, NQ, NP>>
constexpr operator/(U s, RationalPolynomial<T, NP, NQ> const &r) {
    return (s * r.q) / r.p;
}

template <typename T, size_t NP, size_t NQ>
std::ostream &operator<<(std::ostream &os, RationalPolynomial<T, NP, NQ> const &r) {
    return os << r.p << "/" << r.q;
}

} // namespace staticpolynomial
