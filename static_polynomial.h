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

#pragma once

#include <iostream>
#include <sstream>
#include <array>
#include <cassert>
#include <type_traits>
#include <cmath>
#include <concepts>
#include <vector>
#include <exception>
#include <numeric>

namespace staticpolynomial {

template <typename T, size_t N>
class Polynomial;

template <typename T, size_t NP, size_t NQ>
class RationalPolynomial;

namespace detail {
    template <class P>
    struct is_scalar : std::true_type { };

    template <typename T, size_t N>
    struct is_scalar<Polynomial<T, N>> : std::false_type { };

    template <typename T, size_t NP, size_t NQ>
    struct is_scalar<RationalPolynomial<T, NP, NQ>> : std::false_type { };

    template <class P>
    inline constexpr bool is_scalar_v = is_scalar<P>::value;

    template <class P>
    concept Scalar = is_scalar_v<P>;

    static_assert(is_scalar_v<int>);
    static_assert(!is_scalar_v<Polynomial<int, 1>>);

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
}

template <typename T, size_t N>
class Polynomial {
public:
    std::array<T, N+1> coeffs;

    static const Polynomial<T, 1> X;

    constexpr Polynomial() : coeffs{} { }

    template <detail::Scalar... U>
    requires(std::convertible_to<U, T> && ...)
    explicit constexpr Polynomial(U... coeffs_) :
        coeffs{{static_cast<T>(coeffs_)...}}
    { }

    template <detail::Scalar U, size_t M>
    requires(std::convertible_to<U, T> && M <= N)
    constexpr Polynomial(Polynomial<U, M> const &o) : coeffs{} {
        for (size_t i=0; i<=M; ++i)
            coeffs[i] = static_cast<T>(o.coeffs[i]);
    }

    constexpr Polynomial(T s) : coeffs{} {
        coeffs[0] = s;
    }

    constexpr size_t order() const { return N; }

    constexpr T const &operator[](size_t i) const { return coeffs[i]; }
    constexpr T       &operator[](size_t i)       { return coeffs[i]; }
    constexpr T coeff_or_zero(size_t i) const { return i <= N ? coeffs[i] : T{0}; }

    template <detail::Scalar U, size_t M>
    constexpr bool operator==(Polynomial<U, M> const &o) const {
        for (size_t i=0; i <= std::max(M,N); ++i) {
            if (coeff_or_zero(i) != o.coeff_or_zero(i))
                return false;
        }
        return true;
    }

    template <detail::Scalar U, size_t M>
    constexpr bool operator!=(Polynomial<U, M> const &o) const {
        return !(*this == o);
    }

    template <detail::Scalar U>
    constexpr bool operator==(U const &s) const {
        return *this == Polynomial<U, 1>(s);
    }

    template <detail::Scalar U>
    constexpr bool operator!=(U const &s) const {
        return *this != Polynomial<U, 1>(s);
    }

    template <detail::Scalar U, size_t M>
    requires(M <= N)
    constexpr Polynomial<T, N> &operator+=(Polynomial<U, M> const &o) {
        for(size_t i=0; i<=M; ++i)
            (*this)[i] += o[i];
        return *this;
    }

    template <detail::Scalar U, size_t M>
    constexpr auto operator+(Polynomial<U, M> const &o) const {
        using R = decltype(coeffs[0] + o.coeffs[0]);
        Polynomial<R, std::max(M,N)> ret{};
        for (size_t i=0; i <= std::max(M,N); ++i) {
            if (i <= std::min(M,N))
                ret.coeffs[i] = coeffs[i] + o.coeffs[i];
            else if (i <= M)
                ret.coeffs[i] = o.coeffs[i];
            else
                ret.coeffs[i] = coeffs[i];
        }
        return ret;
    }

    template <detail::Scalar U>
    constexpr Polynomial<T, N> &operator+=(U s) {
        (*this)[0] += s;
        return *this;
    }

    template <detail::Scalar U>
    constexpr auto operator+(U s) const {
        using R = decltype(coeffs[0] + s);
        return this->cast<R>() += s;
    }

    template <detail::Scalar U, size_t M>
    requires(M <= N)
    constexpr auto &operator-=(Polynomial<U, M> const &o) {
        for(size_t i=0; i<=M; ++i)
            (*this)[i] -= o[i];
        return *this;
    }

    template <detail::Scalar U, size_t M>
    constexpr auto operator-(Polynomial<U, M> const &o) const {
        using R = decltype(coeffs[0] - o.coeffs[0]);
        Polynomial<R, std::max(M,N)> ret{};
        for (size_t i=0; i <= std::max(M,N); ++i) {
            if (i <= std::min(M,N))
                ret.coeffs[i] = coeffs[i] - o.coeffs[i];
            else if (i <= M)
                ret.coeffs[i] = -o.coeffs[i];
            else
                ret.coeffs[i] = coeffs[i];
        }
        return ret;
    }

    template <detail::Scalar U>
    constexpr Polynomial<T, N> &operator-=(U s) {
        (*this)[0] -= s;
        return *this;
    }

    template <detail::Scalar U>
    constexpr auto operator-(U s) const {
        using R = decltype(coeffs[0] - s);
        return this->cast<R>() -= s;
    }

    constexpr Polynomial<T, N> operator-() const {
        auto out = *this;
        for(auto &v : out.coeffs)
            v = -v;
        return out;
    }

    template <detail::Scalar U>
    constexpr Polynomial<T, N> &operator*=(U s) {
        for(size_t i=0; i<=N; ++i)
            (*this)[i] *= s;
        return *this;
    }

    template <detail::Scalar U>
    constexpr auto operator*(U s) const {
        using R = decltype(coeffs[0] * s);
        return this->cast<R>() *= s;
    }

    template <detail::Scalar U>
    constexpr Polynomial<T, N> &operator/=(U s) {
        // Avoid user error of dividing using types that don't support division
        // (i.e., integers).  This can be as simple as X/2, since X is
        // Polynomial<int,1>.
        static_assert((T{1} / U{2}) * U{2} == T{1});
        for(size_t i=0; i<=N; ++i)
            (*this)[i] /= s;
        return *this;
    }

    template <detail::Scalar U>
    constexpr auto operator/(U s) const {
        using R = decltype(coeffs[0] / s);
        return this->cast<R>() /= s;
    }

    template <detail::Scalar U, size_t M>
    constexpr auto operator*(Polynomial<U, M> const &o) const {
        using R = decltype(coeffs[0] * o.coeffs[0]);
        Polynomial<R, N+M> out{};
        for(size_t i=0; i<=M; ++i)
            for(size_t j=0; j<=N; ++j)
                out[i+j] += (*this)[j] * o[i];
        return out;
    }

    template <detail::Scalar U, size_t M>
    constexpr Polynomial<T, N> &operator*=(Polynomial<U, M> const &o) {
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

    // Returns the exponent of the highest order term having non-zero coefficient.
    constexpr size_t leading_term() const {
        for(size_t i=N; i; --i) {
            if((*this)[i] != 0)
                return i;
        }
        return 0;
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

    template <typename U, size_t M>
    constexpr auto subst(Polynomial<U, M> const &f) const {
        using R = decltype(coeffs[0] * f.coeffs[0]);
        Polynomial<R, N * M> outp{};
        detail::for_constexpr<N+1>([&](auto ic) {
            constexpr size_t e = ic.value;
            outp += f.template pow<e>() * (*this)[e];
        });
        return outp;
    }

    template <typename U, size_t MP, size_t MQ>
    constexpr auto subst(RationalPolynomial<U, MP, MQ> const &f) const {
        using R = decltype(coeffs[0] * f.p.coeffs[0]);
        Polynomial<R, N * std::max(MP, MQ)> outp{};
        detail::for_constexpr<N+1>([&](auto ic) {
            constexpr size_t e = ic.value;
            outp += f.q.template pow<N-e>() * f.p.template pow<e>() * (*this)[e];
        });
        return outp / f.q.template pow<N>();
    }

    // Evaluate polynomial at point.
    template <detail::Scalar U>
    constexpr auto operator()(U const &x) const {
        // result type
        using R = decltype(coeffs[0] * x);
        R accum{};
        for(size_t i=0; i<=N; ++i) {
            if(i)
                accum *= x;
            accum += (*this)[N-i];
        }
        return accum;
    }

    constexpr Polynomial<T, std::max(N,size_t(1))-1> derivative() const {
        if constexpr (N == 0)
            return {}; // zero
        Polynomial<T, std::max(N,size_t(1))-1> out;
        for(size_t i=0; i<N; ++i) {
            out[i] = static_cast<T>(i+1) * (*this)[i+1];
        }
        return out;
    }

    constexpr Polynomial<T, N+1> integral() const {
        Polynomial<T, N+1> out;
        for(size_t i=0; i<=N; ++i) {
            out[i+1] = (*this)[i] / static_cast<T>(i+1);
        }
        return out;
    }

    constexpr T integral(T t0, T t1) const {
        auto const q = integral();
        return q(t1) - q(t0);
    }

    // Roots are sorted, and roots with multiplicity are repeated.
    constexpr std::vector<T> real_roots() const;

    constexpr T discriminant() const;

    template <size_t M>
    constexpr Polynomial<T, N+M> times_power_of_x() const;

    constexpr Polynomial<T, N-1> drop_leading_term() const;

    // Divide by the gcd of the coefficents.  Only applicable if T is an integer
    // type.
    constexpr Polynomial<T, N> divide_out_common_factors() const;

    // Returns a polynomial proportional to mod(*this, o).  It differs from the
    // ordinary remainder by a multiplicative constant.  This allows computing
    // remainders of polynomials over the ring of integers, where we don't have
    // division.  T should be a multi-precision integer type, since the
    // coefficients can get very large.
    template <size_t M>
    constexpr Polynomial<T, std::max(M,size_t(1))-1> pseudo_remainder(Polynomial<T, M> const &o) const;

    // Returns a polynomial proportional to gcd(*this, o).  It differs from
    // ordinary gcd by a multiplicative constant.  This allows computing gcd of
    // polynomials over the ring of integers, where we don't have division.  T
    // should be a multi-precision integer type, since the coefficients can get
    // very large.  If do_reduction is set, then it will divide out common
    // integer factors to reduce the magnitude of the numbers involved.
    template <size_t M>
    constexpr Polynomial<T, std::min(N, M)> pseudo_gcd(Polynomial<T, M> const &o, bool do_reduction=true) const;

    template <typename U>
    constexpr Polynomial<U, N> cast() const
    {
        Polynomial<U, N> out;
        for (size_t i=0; i <= N; ++i)
            out.coeffs[i] = static_cast<U>(coeffs[i]);
        return out;
    }

    std::ostream &print(std::ostream &, std::string const &varname) const;

    std::string to_string(std::string const &varname) const;
};

template <typename T, size_t N>
const Polynomial<T, 1> Polynomial<T, N>::X{0, 1};

template <detail::Scalar... U>
Polynomial(U... coeffs) -> Polynomial<std::common_type_t<U...>, sizeof...(coeffs)-1>;

template <detail::Scalar U, typename T, size_t N>
constexpr auto operator*(U s, Polynomial<T, N> const &p) {
    return p * s;
}

template <detail::Scalar U, typename T, size_t N>
constexpr auto operator+(U s, Polynomial<T, N> const &p) {
    return p + s;
}

template <detail::Scalar U, typename T, size_t N>
constexpr auto operator-(U s, Polynomial<T, N> const &p) {
    return -p + s;
}

//template <size_t N=1, typename T=int>
//static constexpr Polynomial<T, N> X = Polynomial<T, 1>{T{0}, T{1}}.template pow<N>();
static constexpr Polynomial<int, 1> X = Polynomial<int, 1>{0, 1};

template <typename T, size_t N>
std::ostream &Polynomial<T, N>::print(std::ostream &os, std::string const &varname) const
{
    os << "(";
    bool first = true;
    for(size_t i=0; i<=N; ++i) {
        if (coeffs[i] == 0)
            continue;
        if(!first) os << " + ";
        os << coeffs[i];
        if(i > 0)
            os << "*" << varname;
        if(i > 1)
            os << "^" << i;
        first = false;
    }
    if (first)
        os << "0";
    return os << ")";
}

template <typename T, size_t N>
std::string Polynomial<T, N>::to_string(std::string const &varname) const
{
    std::ostringstream os;
    print(os, varname);
    return os.str();
}

template <typename T, size_t N>
std::ostream &operator<<(std::ostream &os, Polynomial<T, N> const &p)
{
    return p.print(os, "X");
}

template <typename T>
constexpr T discriminant(Polynomial<T, 2> const &p)
{
    T const a = p[2];
    T const b = p[1];
    T const c = p[0];
    return b*b - 4*a*c;
}

template <typename T>
requires std::floating_point<T>
constexpr std::vector<T> real_roots(Polynomial<T, 1> const &p)
{
    if (p[1] == 0)
        return { };
    else
        return {-p[0] / p[1]};
}

// Numerically stable version of quadratic equation.
// https://math.stackexchange.com/a/2007723/25589
template <typename T>
requires std::floating_point<T>
constexpr std::vector<T> real_roots(Polynomial<T, 2> const &p)
{
    T const a = p[2];
    T const b = p[1];
    T const c = p[0];

    if (a == 0) { // linear polynomial
        return real_roots(Polynomial<T, 1>(c, b));
    } else if (c == 0) { // has zero as a root
        std::vector<T> out = real_roots(Polynomial<T, 1>(b, a));
        out.emplace_back(T(0));
        if (out[0] > out[1])
            std::swap(out[0], out[1]);
        return out;
    } else {
        T const d = discriminant(p);
        if (d < T(0)) // no real roots
            return {};
        T const r1 = -(b + std::copysign(std::sqrt(d), b)) / (2*a);
        T const r2 = c / (a*r1);
        return { std::min(r1, r2), std::max(r1, r2) };
    }
}

template <typename T, size_t N>
constexpr std::vector<T> Polynomial<T, N>::real_roots() const
{
    // Defer to free function since template partial specialization is not
    // allowed for methods.
    return staticpolynomial::real_roots(*this);
}

template <typename T, size_t N>
constexpr T Polynomial<T, N>::discriminant() const
{
    // Defer to free function since template partial specialization is not
    // allowed for methods.
    return staticpolynomial::discriminant(*this);
}

template <typename T, size_t N>
template <size_t M>
constexpr Polynomial<T, N+M> Polynomial<T, N>::times_power_of_x() const
{
    Polynomial<T, N+M> r{};
    std::copy(coeffs.begin(), coeffs.end(), r.coeffs.begin() + M);
    return r;
}

template <typename T, size_t N>
constexpr Polynomial<T, N-1> Polynomial<T, N>::drop_leading_term() const
{
    static_assert(N >= 1);
    Polynomial<T, N-1> r{};
    std::copy(coeffs.begin(), coeffs.end()-1, r.coeffs.begin());
    return r;
}

template <typename T, size_t N>
constexpr Polynomial<T, N> Polynomial<T, N>::divide_out_common_factors() const
{
    using std::gcd;
    T const g = std::reduce(coeffs.begin(), coeffs.end(), T{0},
        [](T const &x, T const &y) { return gcd(x, y); });
    Polynomial<T, N> r = *this;
    for (T &c : r.coeffs)
        c /= g;
    return r;
}

template <typename T, size_t N>
template <size_t M>
constexpr Polynomial<T, std::max(M,size_t(1))-1> Polynomial<T, N>::pseudo_remainder(Polynomial<T, M> const &o) const
{
    if constexpr (M == 0) {
        if (o.coeffs.back() == 0)
            throw std::runtime_error("division by zero");
        return Polynomial<T, 0>{}; // zero
    } else {
        if (o.coeffs.back() == 0) {
            return pseudo_remainder(o.drop_leading_term());
        }
        if constexpr (M > N) {
            return *this;
        } else { // N >= M
            // FIXME should probably factor out common factor of a and b
            T const a = coeffs.back();
            T const b = o.coeffs.back();
            Polynomial<T, N> const r = (*this) * b - o.template times_power_of_x<N-M>() * a;
            Polynomial<T, N-1> s = r.drop_leading_term();
            if constexpr (N == M)
                return s;
            else
                return s.pseudo_remainder(o);
        }
    }
}

template <typename T, size_t NA, size_t NB>
inline Polynomial<T, std::min(NA, NB)>
pseudo_gcd(Polynomial<T, NA> const &a, Polynomial<T, NB> const &b, bool do_reduction=true)
{
    using Ret = Polynomial<T, std::min(NA, NB)>;
    if (a == 0 || b == 0)
        return Ret{}; // zero
    if constexpr (NA < NB) {
        return pseudo_gcd(b, a);
    } else if constexpr (NB == 0) {
        // We already handled the case a==0 or b==0.  So at this point b is a
        // nonzero constant.
        return Ret{1};
    } else {
        auto r = a.pseudo_remainder(b);
        if (r == 0) {
            if (do_reduction)
                return b.divide_out_common_factors();
            else
                return b;
        } else {
            if (do_reduction)
                r = r.divide_out_common_factors();
            return Ret(pseudo_gcd(b, r));
        }
    }
}

template <typename T, size_t N>
template <size_t M>
constexpr Polynomial<T, std::min(N, M)>
Polynomial<T, N>::pseudo_gcd(Polynomial<T, M> const &o, bool do_reduction) const
{
    return staticpolynomial::pseudo_gcd(*this, o, do_reduction);
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

    template <detail::Scalar U>
    RationalPolynomial<T, NP, NQ> &operator*=(U s) {
        p *= s;
        return *this;
    }

    template <detail::Scalar U>
    auto operator*(T s) const {
        return RationalPolynomial{p*s, q};
    }

    template <detail::Scalar U>
    RationalPolynomial<T, NP, NQ> &operator/=(U s) {
        q *= s;
        return *this;
    }

    template <detail::Scalar U>
    auto operator/(U s) const {
        return RationalPolynomial{p, q*s};
    }

    template <typename U, size_t MP, size_t MQ>
    constexpr auto operator*(RationalPolynomial<U, MP, MQ> const &f) const {
        return (p * f.p) / (q * f.q);
    }

    template <typename U, size_t MP, size_t MQ>
    constexpr auto operator/(RationalPolynomial<U, MP, MQ> const &f) const {
        return (p * f.q) / (q * f.p);
    }

    template <typename U, size_t M>
    constexpr auto operator*(Polynomial<U, M> const &f) const {
        return (p * f) / q;
    }

    template <typename U, size_t M>
    constexpr auto operator/(Polynomial<U, M> const &f) const {
        return p / (q * f);
    }

    constexpr auto normalize() const {
        return RationalPolynomial<T, NP, NQ>{p / q[NQ], q / q[NQ]};
    }

    // Evaluate polynomial at point.
    template <detail::Scalar U>
    constexpr auto operator()(U const &x) const {
        return p(x) / q(x);
    }
};

template <detail::Scalar U, typename T, size_t NP, size_t NQ>
constexpr auto operator*(U s, RationalPolynomial<T, NP, NQ> const &p) {
    return p * s;
}

template <typename T, size_t NP, size_t NQ>
constexpr RationalPolynomial<T, NP, NQ> operator/(Polynomial<T, NP> const &p, Polynomial<T, NQ> const &q) {
    return {p, q};
}

template <detail::Scalar U, typename T, size_t N>
constexpr auto operator/(U s, Polynomial<T, N> const &p) {
    return Polynomial<T, 0>{s} / p;
}

template <detail::Scalar U, typename T, size_t NP, size_t NQ>
constexpr auto operator/(U s, RationalPolynomial<T, NP, NQ> const &r) {
    return (s * r.q) / r.p;
}

template <typename T, size_t NP, size_t NQ>
std::ostream &operator<<(std::ostream &os, RationalPolynomial<T, NP, NQ> const &r) {
    return os << r.p << "/" << r.q;
}

} // namespace staticpolynomial
