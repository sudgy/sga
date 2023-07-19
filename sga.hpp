#ifndef SUDGY_GEOMETRIC_ALGEBRA_H
#define SUDGY_GEOMETRIC_ALGEBRA_H

#include <cstdint>
#include <cmath>
#include <array>
#include <vector>
#include <set>
#include <algorithm>
#include <concepts>

constexpr unsigned int GC_max_dimension = 64;

template <typename T>
concept scalar = std::regular<T> and requires(const T& a, const T& b, T c) {
    {a + b} -> std::convertible_to<T>;
    {a - b} -> std::convertible_to<T>;
    {a * b} -> std::convertible_to<T>;
    {a / b} -> std::convertible_to<T>;
    {+b} -> std::convertible_to<T>;
    {-b} -> std::convertible_to<T>;
    {c = 0};
    {c = 1};
    {c = -1};
};

template <
    scalar Scalar,
    std::uint64_t blade1,
    std::uint64_t blade2,
    std::uint64_t N
>
consteval
Scalar basis_blade_product_sign(const std::array<std::int8_t, N>& metric)
{
    static_assert(N <= GC_max_dimension,
            "The metric you are trying to use has too high of a dimension.");
    auto sign = false;
    auto parity = false;
    for (std::uint8_t i = 0; i < N; ++i) {
        if (blade1 & (1UL << i)) parity = !parity;
    }
    for (std::uint8_t i = 0; i < N; ++i) {
        const auto first = static_cast<bool>(blade1 & (1UL << i));
        const auto second = static_cast<bool>(blade2 & (1UL << i));
        if (first) parity = !parity;
        if (second) sign ^= parity;
        if (first and second) {
            if (metric[i] == -1) sign = !sign;
            else if (metric[i] == 0) return 0;
        }
    }
    if (sign) return -1;
    else return 1;
}

template <std::size_t N>
consteval int dual_sign(std::uint64_t b)
{
    auto sign = false;
    auto parity = false;
    for (std::uint8_t i = 0; i < N; ++i) {
        const auto bit = static_cast<bool>(b & (1UL << i));
        if (bit) sign ^= parity;
        else parity = !parity;
    }
    return sign ? -1 : 1;
}

template <
    scalar Scalar,
    auto metric,
    std::uint64_t basis
>
struct BasisBlade {
    Scalar coefficient = 0;
    constexpr BasisBlade(Scalar x) : coefficient(x) {}
    constexpr BasisBlade()=default;
    constexpr bool operator==(const BasisBlade&) const=default;
};

/* Finds the grade of a basis blade */
consteval std::uint64_t basis_grade(std::uint64_t basis)
{
    auto result = 0UL;
    for (std::uint8_t i = 0; i < 64U; ++i) {
        if (basis & (1UL << i)) ++result;
    }
    return result;
}

template <
    scalar Scalar,
    auto metric,
    std::uint64_t... bases
> requires(std::ranges::is_sorted(
            std::array<std::uint64_t, sizeof...(bases)>{bases...}))
class Multivector : public BasisBlade<Scalar, metric, bases>... {
    // This is the one bit of template metaprogramming used, for converting
    // vectors of basis indices to a multivector with those basis multvectors
    template <
        auto arr,
        typename IS = decltype(std::make_index_sequence<arr().size()>())
    >
    struct from_wrapper_type;
    template <auto arr, std::size_t... I>
    struct from_wrapper_type<arr, std::index_sequence<I...>> {
        using type = Multivector<Scalar, metric, arr()[I]...>;
    };
    template <auto arr>
    using from_wrapper = typename from_wrapper_type<arr>::type;

    template <std::uint64_t basis>
    using BB = BasisBlade<Scalar, metric, basis>;
    template <std::uint64_t... bases2>
    using MV = Multivector<Scalar, metric, bases2...>;

    public:
        /* Constructors */
        constexpr Multivector() : BB<bases>()... {}
        constexpr Multivector(BB<bases>... values)
            requires(sizeof...(bases) > 0)
            : BB<bases>(values)... {}

        /* Comparison, note that it only works between multivectors of the same
         * type.
         */
        constexpr bool operator==(const Multivector& other) const=default;
        constexpr bool operator!=(const Multivector& other) const=default;
        template <std::uint64_t... bases2>
        constexpr bool operator==(const MV<bases2...>&) const=delete;
        template <std::uint64_t... bases2>
        constexpr bool operator!=(const MV<bases2...>&) const=delete;
        constexpr bool operator==(Scalar s) const
            requires(sizeof...(bases) == 0)
        {
            return s == 0;
        }
        constexpr bool operator!=(Scalar s) const
            requires(sizeof...(bases) == 0)
        {
            return s != 0;
        }
        constexpr bool operator==(Scalar s) const
            requires(((bases == 0) and ...) and sizeof...(bases) > 0)
        {
            return this->coefficient == s;
        }
        constexpr bool operator!=(Scalar s) const
            requires(((bases == 0) and ...) and sizeof...(bases) > 0)
        {
            return this->coefficient != s;
        }
        constexpr friend bool operator==(Scalar s, const Multivector& m)
        {
            return m == s;
        }
        constexpr friend bool operator!=(Scalar s, const Multivector& m)
        {
            return m != s;
        }

        /* Generic projection */
        template <std::uint64_t... bases2>
        constexpr void generic_project(const MV<bases2...>& b)
        {
            ((BB<bases>::coefficient = b.BB<bases>::coefficient), ...);
        }

        /* Blade projection, note that you must know the binary for your blade
         */
        template <std::uint64_t basis>
        constexpr Scalar blade_project() const
        {
            if constexpr ((... || (bases == basis))) {
                return BB<basis>::coefficient;
            }
            else return 0;
        }

        /* Grade projection */
        template <std::uint64_t grade>
        constexpr auto grade_project() const
        {
            auto result = from_wrapper<[&]{
                auto result = std::vector<std::uint64_t>();
                result.reserve(sizeof...(bases));
                ((basis_grade(bases) == grade ?
                  (result.push_back(bases)) : void()), ...);
                return result;
            }>();
            result.generic_project(*this);
            return result;
        }

        /* Addition */
        constexpr auto operator+() const
        {
            return *this;
        }
        template <std::uint64_t... bases2>
        constexpr auto& operator+=(const MV<bases2...>& other)
        {
            ((BB<bases2>::coefficient +=
                other.BB<bases2>::coefficient), ...);
            return *this;
        }
        // (Hopefully) an optimization for when adding two multivectors of the
        // same type
        constexpr auto operator+(const Multivector& other) const
        {
            auto result = *this;
            result += other;
            return result;
        }
        template <std::uint64_t... bases2>
        constexpr auto operator+(const MV<bases2...>& other) const
        {
            auto result = from_wrapper<[&]{
                auto result = std::vector<std::uint64_t>();
                result.reserve(sizeof...(bases) + sizeof...(bases2));
                auto bs1 = std::array{bases...};
                auto bs2 = std::array<std::uint64_t, sizeof...(bases2)>
                    {bases2...};
                auto i1 = 0UL;
                auto i2 = 0UL;
                while (i1 < bs1.size() and i2 < bs2.size()) {
                    if (bs1[i1] < bs2[i2]) {
                        result.push_back(bs1[i1]);
                        ++i1;
                    }
                    else if (bs2[i2] < bs1[i1]) {
                        result.push_back(bs2[i2]);
                        ++i2;
                    }
                    else {
                        result.push_back(bs1[i1]);
                        ++i1;
                        ++i2;
                    }
                }
                std::copy(
                    bs1.begin() + i1,
                    bs1.end(),
                    std::back_inserter(result)
                );
                std::copy(
                    bs2.begin() + i2,
                    bs2.end(),
                    std::back_inserter(result)
                );
                return result;
            }>();
            result += *this;
            result += other;
            return result;
        }

        constexpr auto operator+=(const Scalar& s)
        {
            *this += MV<0>(s);
        }
        constexpr auto operator+(const Scalar& s) const
        {
            return *this + MV<0>(s);
        }
        constexpr friend auto operator+(const Scalar& s, const Multivector& m)
        {
            return MV<0>(s) + m;
        }
        /* Subtraction */
        constexpr auto operator-() const
        {
            auto result = *this;
            ((result.BB<bases>::coefficient =
                -result.BB<bases>::coefficient), ...);
            return result;
        }
        template <std::uint64_t... bases2>
        constexpr auto& operator-=(const MV<bases2...>& other)
        {
            ((BB<bases2>::coefficient -=
                other.BB<bases2>::coefficient), ...);
            return *this;
        }
        // (Hopefully) an optimization for when adding two multivectors of the
        // same type
        constexpr auto operator-(const Multivector& other) const
        {
            auto result = *this;
            ((result.BB<bases>::coefficient -=
                other.BB<bases>::coefficient), ...);
            return result;
        }
        template <std::uint64_t... bases2>
        constexpr auto operator-(const MV<bases2...>& other) const
        {
            return *this + (-other);
        }

        constexpr auto operator-=(const Scalar& s)
        {
            *this -= MV<0>(s);
        }
        constexpr auto operator-(const Scalar& s) const
        {
            return *this - MV<0>(s);
        }
        constexpr friend auto operator-(const Scalar& s, const Multivector& m)
        {
            return MV<0>(s) - m;
        }

        // Duality
        constexpr auto dual() const
        {
            constexpr auto N = metric.size();
            auto result = from_wrapper<[]{
                auto result = std::vector<std::uint64_t>{bases...};
                for (auto& b : result) {
                    b = b ^ ((1 << N) - 1);
                }
                std::ranges::sort(result);
                return result;
            }>();
            ((result.BB<bases ^ ((1 << N) - 1)>::coefficient =
                BB<bases>::coefficient * dual_sign<N>(bases)), ...);
            return result;
        }
        // I know that there are more efficient versions of this but I don't
        // want to bother
        constexpr auto undual() const
        {
            return dual().dual().dual();
        }
        constexpr auto operator!() const
        {
            return dual();
        }

        /* Involutions */
        constexpr Multivector reverse() const
        {
            auto result = *this;
            ((result.BB<bases>::coefficient *=
              ((basis_grade(bases) / 2) % 2 ? -1 : 1)), ...);
            return result;
        }
        constexpr Multivector operator~() const
        {
            return reverse();
        }
        constexpr Multivector involute() const
        {
            auto result = *this;
            ((result.BB<bases>::coefficient *=
              (basis_grade(bases) % 2 ? -1 : 1)), ...);
            return result;
        }

        /* Norms */
        constexpr Scalar norm2() const
        {
            return (*this * reverse()).template blade_project<0>();
        }
        constexpr Scalar norm() const
        {
            return std::sqrt(std::abs(norm2()));
        }
        constexpr Multivector normalized() const
        {
            return *this / norm();
        }

        /* Multiplication */
        template <std::uint64_t... bases2>
        constexpr auto operator*(const MV<bases2...>& b) const
        {
            return generic_mult<[](auto, auto){return true;}>(b);
        }

        constexpr auto operator*=(const Scalar& s)
        {
            ((BB<bases>::coefficient *= s), ...);
        }
        constexpr auto operator/=(const Scalar& s)
        {
            ((BB<bases>::coefficient /= s), ...);
        }
        constexpr auto operator*(const Scalar& s) const
        {
            auto result = *this;
            result *= s;
            return result;
        }
        constexpr auto operator/(const Scalar& s) const
        {
            auto result = *this;
            result /= s;
            return result;
        }
        constexpr friend auto operator*(const Scalar& s, const Multivector& m)
        {
            return m * s;
        }

        /* Other products */
        template <std::uint64_t... bases2>
        constexpr auto operator^(const MV<bases2...>& b) const
        {
            return generic_mult<[](std::uint64_t m, std::uint64_t n){
                return (m & n) == 0;
            }>(b);
        }
        template <std::uint64_t... bases2>
        constexpr auto operator|(const MV<bases2...>& b) const
        {
            return generic_mult<[](std::uint64_t m, std::uint64_t n){
                return (m | n) == m or (m | n) == n;
            }>(b);
        }
        template <std::uint64_t... bases2>
        constexpr auto operator<<(const MV<bases2...>& b) const
        {
            return generic_mult<[](std::uint64_t m, std::uint64_t n){
                return (m | n) == n;
            }>(b);
        }
        template <std::uint64_t... bases2>
        constexpr auto operator>>(const MV<bases2...>& b) const
        {
            return generic_mult<[](std::uint64_t m, std::uint64_t n){
                return (m | n) == m;
            }>(b);
        }
        template <std::uint64_t... bases2>
        constexpr auto operator&(const MV<bases2...>& b) const
        {
            return (dual() ^ b.dual()).undual();
        }
        constexpr auto operator^(const MV<>& b) const
            {return b;}
        constexpr auto operator|(const MV<>& b) const
            {return b;}
        constexpr auto operator<<(const MV<>& b) const
            {return b;}
        constexpr auto operator>>(const MV<>& b) const
            {return b;}
        constexpr auto operator&(const MV<>& b) const
            {return b;}
        constexpr friend auto operator^(const MV<>& b, const Multivector&)
            {return b;}
        constexpr friend auto operator|(const MV<>& b, const Multivector&)
            {return b;}
        constexpr friend auto operator<<(const MV<>& b, const Multivector&)
            {return b;}
        constexpr friend auto operator>>(const MV<>& b, const Multivector&)
            {return b;}
        constexpr friend auto operator&(const MV<>& b, const Multivector&)
            {return b;}

    private:
        // I know I promised that the above was the only template
        // metaprogramming, but this is basically the same thing
        template <
            auto arr,
            typename IS = decltype(std::make_index_sequence<arr().size()>())
        >
        struct gm_calc;
        template <auto arr, std::size_t... I>
        struct gm_calc<arr, std::index_sequence<I...>> {
            template <
                std::pair<std::uint64_t, std::uint64_t>... bs,
                typename R,
                typename A,
                typename B
            >
            constexpr static void real_calc(R& res, const A& a, const B& b)
            {
                ((res.BB<bs.first ^ bs.second>::coefficient +=
                    basis_blade_product_sign<
                        Scalar, bs.first, bs.second>(metric) *
                    a.BB<bs.first>::coefficient * b.BB<bs.second>::coefficient),
                 ...);
            }
            template <
                typename R,
                typename A,
                typename B
            >
            constexpr static void calc(R& res, const A& a, const B& b)
            {
                real_calc<arr()[I]...>(res, a, b);
            }
        };
        template <auto P, std::uint64_t... bases2>
        constexpr auto generic_mult(const MV<bases2...>& b) const
        {
            constexpr auto pairs = [&]{
                auto result
                    = std::vector<std::pair<std::uint64_t, std::uint64_t>>();
                result.reserve(sizeof...(bases) * sizeof...(bases2));
                for (auto i1 : {bases...}) {
                    for (auto i2 : {bases2...}) {
                        if (P(i1, i2)) {
                            result.emplace_back(i1, i2);
                        }
                    }
                }
                return result;
            };
            auto result = from_wrapper<[&]{
                auto result = std::vector<std::uint64_t>();
                result.reserve(sizeof...(bases) * sizeof...(bases2));
                for (auto i1 : {bases...}) {
                    for (auto i2 : {bases2...}) {
                        if (P(i1, i2)) {
                            result.push_back(i1 ^ i2);
                        }
                    }
                }
                std::ranges::sort(result);
                auto ret = std::ranges::unique(result);
                result.erase(ret.begin(), ret.end());
                return result;
            }>();
            gm_calc<pairs>::calc(result, *this, b);
            return result;
        }
};

// Example usage
constexpr auto metric = std::array<std::int8_t, 4>{0, 1, 1, 1};

template <std::uint64_t... bases>
using PGAMultivector = Multivector<double, metric, bases...>;
using Scalar = PGAMultivector<0>;
using Vector = PGAMultivector<1, 2, 4, 8>;
using Bivector = PGAMultivector<3, 5, 6, 9, 10, 12>;
using Trivector = PGAMultivector<7, 11, 13, 14>;
using Quadvector = PGAMultivector<15>;
constexpr auto e0 = Vector(1, 0, 0, 0);
constexpr auto e1 = Vector(0, 1, 0, 0);
constexpr auto e2 = Vector(0, 0, 1, 0);
constexpr auto e3 = Vector(0, 0, 0, 1);

constexpr std::array<Trivector, 4> plane_points(
    Vector v,
    Trivector center,
    Trivector up
)
{
    const auto axis = v | center;
    const auto center2 = axis ^ v;
    const auto up2 = (v | up) ^ v;
    const auto vert = center2 & up2;
    const auto perp1 = (vert | v).normalized();
    const auto perp2 = (perp1 | axis).normalized();
    return {
        (perp1 + 5*e0) ^ (perp2 + 5*e0) ^ v,
        (perp1 - 5*e0) ^ (perp2 + 5*e0) ^ v,
        (perp1 - 5*e0) ^ (perp2 - 5*e0) ^ v,
        (perp1 + 5*e0) ^ (perp2 - 5*e0) ^ v
    };
}

#endif
