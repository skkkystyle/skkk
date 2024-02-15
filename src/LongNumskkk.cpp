#include "../include/LongNumskkk.hpp"

// includes
#include <complex>
#include <algorithm>
#include <functional>
#include <tuple>

// technical
namespace
{
    // constants
    constexpr unsigned MAXN = (1 << 17);
    const unsigned BASE = 10;
    const double PI = acos(-1);
    const LongNum ZERO(0);
    const double DIGITS_PER_TERM = 14.181647462725477;

    // global parameters
    unsigned gAccuracy = 15;
    unsigned gLimitation = 300;
    bool gInitFlag = false;

    // part of the FFT
    typedef std::complex<double> ftype;
    ftype roots[MAXN];
}

// FFT calculations
void initialization()
{
    for (int i = 0; i < MAXN; i++)
    {
        roots[i] = std::polar(1., 2 * PI / MAXN * i);
    }
    gInitFlag = true;
}

template<typename T>
void fft(T *in, ftype *out, int n, int k = 1)
{
    if (n == 1)
    {
        *out = *in;
        return;
    }
    int t = MAXN / n;
    n >>= 1;
    fft(in, out, n, 2 * k);
    fft(in + k, out + n, n, 2 * k);
    for (int i = 0, j = 0; i < n; i++, j += t)
    {
        ftype p = roots[j] * out[i + n];
        out[i + n] = out[i] - p;
        out[i] += p;
    }
}

std::vector<ftype> evaluation(std::vector<int> p)
{
    while (__builtin_popcount(p.size()) != 1)
    {
        p.push_back(0);
    }
    std::vector<ftype> result(p.size());
    fft(p.data(), result.data(), p.size());
    return result;
}

std::vector<int> interpolation(std::vector<ftype> p)
{
    int n = p.size();
    std::vector<ftype> inv(n);
    fft(p.data(), inv.data(), n);
    std::vector<int> result(n);
    for (int i = 0; i < n; i++)
    {
        result[i] = round(std::real(inv[i]) / n);
    }
    std::reverse(std::begin(result) + 1, std::end(result));
    return result;
}

void aligning(std::vector<int> &a, std::vector<int> &b)
{
    int n = a.size() + b.size() - 1;
    while (a.size() < n)
    {
        a.push_back(0);
    }
    while (b.size() < n)
    {
        b.push_back(0);
    }
}

std::vector<int> normalizating(std::vector<int> c)
{
    int one = 0;
    for (auto &it : c)
    {
        it += one;
        one = it / BASE;
        it %= BASE;
    }
    while (one)
    {
        c.push_back(one % BASE);
        one /= BASE;
    }
    return c;
}

std::vector<int> multiplying(std::vector<int> a, std::vector<int> b)
{
    aligning(a, b);
    auto A = evaluation(a);
    auto B = evaluation(b);
    for (int i = 0; i < A.size(); i++)
    {
        A[i] *= B[i];
    }
    return normalizating(interpolation(A));
}

// constructors
LongNum::LongNum()
{
    number.resize(1, 0);
    fraction.clear();
    depth = 0;
    positivity = true;
}

LongNum::LongNum(long long number)
{
    fraction.clear();
    depth = 0;
    if (number == 0)
    {
        positivity = true;
        this->number.push_back(0);
        return;
    }
    if (number > 0)
    {
        positivity = true;
    }
    else
    {
        positivity = false;
        number = llabs(number);
    }
    while (number)
    {
        this->number.push_back(number % BASE);
        number = number / BASE;
    }
}

LongNum::LongNum(const std::string &s)
{
    positivity = true;
    bool flag = true;
    depth = 0;
    for (auto digit : s)
    {
        if (digit == '-')
        {
            positivity = false;
            continue;
        }
        if (digit == '.' || digit == ',')
        {
            flag = false;
            continue;
        }
        int temp = digit - '0';
        if (flag)
        {
            number.push_back(temp);
        }
        else
        {
            depth++;
            fraction.push_back(temp);
        }
    }
    while (fraction.size() > 1 && !fraction.back())
    {
        fraction.pop_back();
        depth--;
    }
    if (number.size() == 1 && number[0] == 0 && depth == 0)
    {
        positivity = true;
    }
    std::reverse(fraction.begin(), fraction.end());
    std::reverse(number.begin(), number.end());
}

LongNum::LongNum(std::vector<int> number, std::vector<int> fraction, unsigned depth, bool positivity)
{
    this->positivity = positivity;
    this->depth = depth;
    this->number = std::move(number);
    this->fraction = std::move(fraction);
}

// helpers
void display(const LongNum &a)
{
    if (!a.positivity)
    {
        std::cout << '-';
    }
    for (int i = a.number.size() - 1; i >= 0; i--)
    {
        std::cout << a.number[i];
    }
    if (gAccuracy == 0)
    {
        std::cout << std::endl;
        return;
    }
    std::cout << '.';
    int fprinted = 0;
    for (int i = a.fraction.size() - 1; i >= 0 && fprinted < gAccuracy; i--, fprinted++)
    {
        std::cout << a.fraction[i];
    }
    while (fprinted < gAccuracy)
    {
        std::cout << '0';
        fprinted++;
    }
    std::cout << std::endl;
}

std::string to_string(const LongNum &a)
{
    std::string str;
    if (!a.positivity)
    {
        str += '-';
    }
    for (int i = a.number.size() - 1; i >= 0; i--)
    {
        str += a.number[i] + '0';
    }
    if (gAccuracy == 0)
    {
        return str;
    }
    str += '.';
    int fprinted = 0;
    for (int i = a.fraction.size() - 1; i >= 0 && fprinted < gAccuracy; i--, fprinted++)
    {
        str += a.fraction[i] + '0';
    }
    while (fprinted < gAccuracy)
    {
        str += '0';
        fprinted++;
    }
    str += '\n';
    return str;
}

// math functions
LongNum LongNumAbs(const LongNum &a)
{
    LongNum c = a;
    c.positivity = true;
    return c;
}

LongNum LongNumSqrt(const LongNum &a)
{
    if (a < ZERO)
    {
        throw "Square-rooting a negative number error\n";
    }
    LongNum b = a, hunnid(100), ten(10);
    int shift = 0;
    for (int i = 0; i < a.depth / 2 + a.depth % 2 || i < gAccuracy + 2; i++, shift++)
    {
        b *= hunnid;
    }
    if (b.number.size() % 2 == 0)
    {
        b *= hunnid;
        shift++;
    }

    LongNum start, end = b, mid, ans, half("0,5"), one(1);
    while (start <= end)
    {
        mid = (start + end) * half;
        mid.fraction.clear();
        mid.depth = 0;
        if (mid * mid == b)
        {
            ans = mid;
            break;
        }
        if (mid * mid < b)
        {
            ans = start;
            start = mid + one;
        }
        else
        {
            end = mid - one;
        }
    }

    while (shift--)
    {
        ans /= ten;
    }

    return ans;
}

LongNum exponentiation(const LongNum &a, unsigned degree)
{
    LongNum res(1);
    while (degree--)
    {
        res *= a;
    }
    return res;
}

LongNum LongNum::Chudnovskkky(unsigned digits)
{
    unsigned tmp = gAccuracy;
    set_accuracy(0);
    LongNum C(640320), C3_OVER_24 = C * C * C / LongNum(24);
    LongNum one(1), two(2), five(5), six(6);
    LongNum otf("13591409"), fff("545140134");

    std::function<std::tuple<LongNum, LongNum, LongNum>(const int, const int)> binsplit;
    binsplit = [&](const int a, const int b)
    {
        LongNum Pab, Qab, Tab;
        int m;
        if (b - a == 1)
        {
            LongNum aA(a);
            if (a == 0)
            {
                Pab = one, Qab = one;
            }
            else
            {
                Pab = (six * aA - five) * (two * aA - one) * (six * aA - one);
                Qab = aA * aA * aA * C3_OVER_24;
            }
            Tab = Pab * (otf + fff * aA);
            if (a & 1)
            {
                Tab.positivity = !Tab.positivity;
            }
        }
        else
        {
            m = (a + b) / 2;
            LongNum Pam, Qam, Tam, Pmb, Qmb, Tmb;
            std::tie(Pam, Qam, Tam) = binsplit(a, m);
            std::tie(Pmb, Qmb, Tmb) = binsplit(m, b);

            Pab = Pam * Pmb;
            Qab = Qam * Qmb;
            Tab = Qmb * Tam + Pam * Tmb;
        }
        return std::make_tuple(Pab, Qab, Tab);
    };

    int N = (int) (digits / DIGITS_PER_TERM + 1);
    LongNum P, Q, T;
    std::tie(P, Q, T) = binsplit(0, N);
    LongNum one_squared = exponentiation(LongNum(10), 2 * digits);
    LongNum sqrtC = LongNumSqrt(one_squared * LongNum(10005));
    LongNum answer = Q * LongNum(426880) * sqrtC / T;
    set_accuracy(tmp);
    answer.number.pop_back();
    LongNum Pi({3}, answer.number, answer.number.size(), true);

    return Pi;
}

// operators
LongNum operator + (const LongNum &a, const LongNum &b)
{
    if (a.positivity && !b.positivity)
    {
        LongNum c = b;
        c.positivity = true;
        return a - c;
    }
    if (!a.positivity && b.positivity)
    {
        LongNum c = a;
        c.positivity = true;
        return b - c;
    }

    std::vector<int> result, fresult;
    int one = 0;
    bool pushed = false;
    for (int i = 0; i < std::max(a.depth, b.depth) - std::min(a.depth, b.depth); i++)
    {
        pushed = true;
        if (a.depth > b.depth)
        {
            fresult.push_back(a.fraction[i]);
        }
        if (b.depth > a.depth)
        {
            fresult.push_back(b.fraction[i]);
        }
    }
    for (int i = std::max(a.depth, b.depth) - std::min(a.depth, b.depth), j = 0; i < std::max(a.depth, b.depth); i++, j++)
    {
        if (a.depth > b.depth)
        {
            int temp = one + a.fraction[i] + b.fraction[j];
            if (temp % BASE != 0 || pushed)
            {
                fresult.push_back(temp % BASE);
            }
            one = temp / BASE;
        }
        else
        {
            int temp = one + a.fraction[j] + b.fraction[i];
            if (temp % BASE != 0 || pushed)
            {
                fresult.push_back(temp % BASE);
            }
            one = temp / BASE;
        }
    }

    for (int i = 0; i < std::min(a.number.size(), b.number.size()); i++)
    {
        int temp = one + a.number[i] + b.number[i];
        result.push_back(temp % BASE);
        one = temp / BASE;
    }
    for (int i = std::min(a.number.size(), b.number.size()); i < std::max(a.number.size(), b.number.size()); i++)
    {
        if (a.number.size() > b.number.size())
        {
            int temp = one + a.number[i];
            result.push_back(temp % BASE);
            one = temp / BASE;
        }
        else
        {
            int temp = one + b.number[i];
            result.push_back(temp % BASE);
            one = temp / BASE;
        }
    }
    if (one)
    {
        result.push_back(one);
    }
    while (result.size() > 1 && !result.back())
    {
        result.pop_back();
    }

    LongNum c(result, fresult, fresult.size(), a.positivity);
    return c;
}

LongNum operator - (const LongNum &a, const LongNum &b)
{
    if (a.positivity && !b.positivity)
    {
        LongNum c = b;
        c.positivity = true;
        return a + c;
    }
    if (!a.positivity && b.positivity)
    {
        LongNum c = a;
        c.positivity = true;
        LongNum d = b + c;
        d.positivity = false;
        return d;
    }

    std::vector<int> result, fresult;
    int one = 0;
    LongNum c = a, d = b;
    if (LongNumAbs(b) > LongNumAbs(a))
    {
        c = b;
        d = a;
    }

    if (d.depth > c.depth)
    {
        if (!c.fraction.empty())
        {
            c.fraction[0]--;
        }
        else
        {
            for (int & i : c.number)
            {
                if (i > 0)
                {
                    i--;
                    break;
                }
            }
        }
    }
    bool pushed = false;
    for (int i = 0; i < std::max(c.depth, d.depth) - std::min(c.depth, d.depth); i++)
    {
        pushed = true;
        if (d.depth > c.depth)
        {
            if (i == 0)
            {
                fresult.push_back(10 - d.fraction[i]);
            }
            else
            {
                fresult.push_back(9 - d.fraction[i]);
            }
        }
        else
        {
            fresult.push_back(c.fraction[i]);
        }
    }
    for (int i = std::max(c.depth, d.depth) - std::min(c.depth, d.depth), j = 0; i < std::max(c.depth, d.depth); i++, j++)
    {
        int temp = c.fraction[i] - d.fraction[j] - one;
        if (d.depth > c.depth)
        {
            temp = c.fraction[j] - d.fraction[i] - one;
        }
        if (temp < 0)
        {
            temp += BASE;
            one = 1;
        }
        else
        {
            one = 0;
        }
        if (temp != 0 || pushed)
        {
            fresult.push_back(temp);
            pushed = true;
        }
    }

    for (int i = 0; i < d.number.size(); i++)
    {
        int temp = c.number[i] - d.number[i] - one;
        if (temp < 0)
        {
            temp += BASE;
            one = 1;
        }
        else
        {
            one = 0;
        }
        result.push_back(temp);
    }
    for (int i = d.number.size(); i < c.number.size(); i++)
    {
        if (c.number[i] - one < 0)
        {
            result.push_back(c.number[i] + BASE - one);
            one = 1;
        }
        else
        {
            result.push_back(c.number[i] - one);
            one = 0;
        }
    }
    while (result.size() > 1 && !result.back())
    {
        result.pop_back();
    }

    bool positivity;
    if (LongNumAbs(a) > LongNumAbs(b))
    {
        positivity = a.positivity;
    }
    else
    {
        positivity = !a.positivity;
    }
    LongNum res(result, fresult, fresult.size(), positivity);
    return res;
}

LongNum operator * (const LongNum &a, const LongNum &b)
{
    bool positivity;
    if (a.positivity == b.positivity)
    {
        positivity = true;
    }
    else
    {
        positivity = false;
    }
    if (LongNumAbs(a) == LongNum(1))
    {
        LongNum c = b;
        c.positivity = positivity;
        return c;
    }
    if (LongNumAbs(b) == LongNum(1))
    {
        LongNum c = a;
        c.positivity = positivity;
        return c;
    }
    std::vector<int> a_merged(a.number.size() + a.depth), b_merged(b.number.size() + b.depth);
    int x = 0;

    for (int i = a.fraction.size() - a.depth; i < a.fraction.size(); i++)
    {
        a_merged[x] = a.fraction[i];
        x++;
    }
    for (int i : a.number)
    {
        a_merged[x] = i;
        x++;
    }

    x = 0;
    for (int i = b.fraction.size() - b.depth; i < b.fraction.size(); i++)
    {
        b_merged[x] = b.fraction[i];
        x++;
    }
    for (int i : b.number)
    {
        b_merged[x] = i;
        x++;
    }

    if (!gInitFlag)
    {
        initialization();
    }

    std::vector<int> m_result = multiplying(normalizating(a_merged), normalizating(b_merged));
    std::vector<int> result, fresult;
    bool pushed = false;
    for (int i = 0; i < a.depth + b.depth; i++)
    {
        if (m_result[i] != 0 || pushed)
        {
            fresult.push_back(m_result[i]);
            pushed = true;
        }
    }
    for (int i = a.depth + b.depth; i < m_result.size(); i++)
    {
        result.push_back(m_result[i]);
    }
    while (result.size() > 1 && !result.back())
    {
        result.pop_back();
    }

    LongNum c(result, fresult, fresult.size(), positivity);
    return c;
}

LongNum operator / (const LongNum &a, const LongNum &b)
{
    if (b == ZERO)
    {
        throw "Zero division error\n";
    }
    LongNum c = LongNumAbs(a), d = LongNumAbs(b);
    bool positivity;
    if (a.positivity == b.positivity)
    {
        positivity = true;
    }
    else
    {
        positivity = false;
    }
    if (c == d)
    {
        LongNum one({1}, {}, 0, positivity);
        return one;
    }
    if (d == LongNum(1))
    {
        c.positivity = positivity;
        return c;
    }

    unsigned norming = std::max(c.depth, d.depth);
    for (int i = 0; i < norming; i++)
    {
        if (c.depth > 0)
        {
            c.number.insert(c.number.cbegin(), c.fraction.back());
            c.fraction.pop_back();
            c.depth--;
        }
        else
        {
            c.number.insert(c.number.cbegin(), 0);
        }
        if (d.depth > 0)
        {
            d.number.insert(d.number.cbegin(), d.fraction.back());
            d.fraction.pop_back();
            d.depth--;
        }
        else
        {
            d.number.insert(d.number.cbegin(), 0);
        }
    }

    while (c.number.size() > 1 && !c.number.back())
    {
        c.number.pop_back();
    }
    while (d.number.size() > 1 && !c.number.back())
    {
        d.number.pop_back();
    }

    int shift = 0;
    for (int i = 0; i < a.depth + b.depth || i < gAccuracy; i++, shift++)
    {
        c.number.insert(c.number.cbegin(), 0);
    }
    while (c < d)
    {
        c.number.insert(c.number.cbegin(), 0);
        shift++;
    }

    LongNum low(0), high = c, mid, quotient(0);
    LongNum one(1), half("0,5");
    while (low <= high)
    {
        mid = low + (high - low) * half;
        mid.fraction.clear();
        mid.depth = 0;
        if (mid * d > c)
        {
            high = mid - one;
        }
        else
        {
            quotient = mid;
            low = mid + one;
        }
    }

    while (shift > 0 && !quotient.number.empty())
    {
        quotient.fraction.push_back(quotient.number.front());
        quotient.number.erase(quotient.number.cbegin());
        quotient.depth++;
        shift--;
    }
    while (shift)
    {
        quotient.fraction.push_back(0);
        quotient.depth++;
        shift--;
    }
    if (quotient.number.empty())
    {
        quotient.number.push_back(0);
    }

    quotient.positivity = positivity;
    return quotient;
}

LongNum &LongNum::operator = (const LongNum &b)
{
    this->number = b.number;
    this->fraction = b.fraction;
    this->depth = b.depth;
    this->positivity = b.positivity;
    return *this;
}

LongNum &LongNum::operator += (const LongNum &b)
{
    LongNum c = *this + b;
    this->number = c.number;
    this->fraction = c.fraction;
    this->depth = c.depth;
    this->positivity = c.positivity;
    return *this;
}

LongNum &LongNum::operator -= (const LongNum &b)
{
    LongNum c = *this - b;
    this->number = c.number;
    this->fraction = c.fraction;
    this->depth = c.depth;
    this->positivity = c.positivity;
    return *this;
}

LongNum &LongNum::operator *= (const LongNum &b)
{
    LongNum c = *this * b;
    this->number = c.number;
    this->fraction = c.fraction;
    this->depth = c.depth;
    this->positivity = c.positivity;
    return *this;
}

LongNum &LongNum::operator /= (const LongNum &b)
{
    LongNum c = *this / b;
    this->number = c.number;
    this->fraction = c.fraction;
    this->depth = c.depth;
    this->positivity = c.positivity;
    return *this;
}

bool operator == (const LongNum &a, const LongNum &b)
{
    if (&a == &b)
    {
        return true;
    }
    if (a.depth != b.depth || a.number.size() != b.number.size() || a.positivity != b.positivity)
    {
        return false;
    }
    for (int i = 0; i < a.number.size(); i++)
    {
        if (a.number[i] != b.number[i])
        {
            return false;
        }
    }
    for (int i = 0; i < a.depth; i++)
    {
        if (a.fraction[i] != b.fraction[i])
        {
            return false;
        }
    }
    return true;
}

bool operator != (const LongNum &a, const LongNum &b)
{
    return !(a == b);
}

bool operator <  (const LongNum &a, const LongNum &b)
{
    if (a == b)
    {
        return false;
    }
    if (a.positivity != b.positivity)
    {
        return b.positivity;
    }
    if (a.number.size() < b.number.size())
    {
        return a.positivity && b.positivity;
    }
    if (a.number.size() > b.number.size())
    {
        return !(a.positivity && b.positivity);
    }
    for (int i = a.number.size() - 1; i >= 0; i--)
    {
        if (a.number[i] < b.number[i])
        {
            return a.positivity && b.positivity;
        }
        if (a.number[i] > b.number[i])
        {
            return !(a.positivity && b.positivity);
        }
    }
    for (int i = a.depth - 1, j = b.depth - 1; i >= 0 && j >= 0; i--, j--)
    {
        if (a.fraction[i] < b.fraction[j])
        {
            return a.positivity && b.positivity;
        }
        if (a.fraction[i] > b.fraction[j])
        {
            return !(a.positivity && b.positivity);
        }
    }
    if (a.depth < b.depth)
    {
        return a.positivity && b.positivity;
    }
    else
    {
        return !(a.positivity && b.positivity);
    }
}

bool operator >  (const LongNum &a, const LongNum &b)
{
    return !(a < b || a == b);
}

bool operator <= (const LongNum &a, const LongNum &b)
{
    return !(a > b);
}

bool operator >= (const LongNum &a, const LongNum &b)
{
    return !(a < b);
}

// some other helpers
LongNum operator ""_skkk (const char *string)
{
    LongNum c(string);
    return c;
}

void set_accuracy(unsigned value)
{
    if (value > gLimitation)
    {
        std::cout << "The new accuracy value is greater than the current limit, first change the limit value using set_limitation (at your own risk)" << std::endl;
        std::cout << "The current limitation value is: " << gLimitation << std::endl;
        std::cout << "The current accuracy value is: " << gAccuracy << std::endl;
        return;
    }
    gAccuracy = value;
}

void set_limitation(unsigned value)
{
    gLimitation = value;
}

void check_accuracy()
{
    std::cout << gAccuracy << std::endl;
}

void check_limitation()
{
    std::cout << gLimitation << std::endl;
}

unsigned get_accuracy()
{
    return gAccuracy;
}