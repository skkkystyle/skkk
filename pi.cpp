#include "./include/LongNumskkk.hpp"


#include <functional>
#include <chrono>


LongNum Chudnovsky(unsigned digits)
{
    unsigned tmp = get_accuracy();
    set_accuracy(0);
    LongNum C(640320), C3_OVER_24 = C * C * C / LongNum(24);
    LongNum one(1), two(2), five(5), six(6), sone(-1);
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
                Tab *= sone;
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

    int N = digits / 13 + 1;
    LongNum P, Q, T;
    std::tie(P, Q, T) = binsplit(0, N);
    LongNum one_squared = exponentiation(LongNum(10), 2 * digits);
    LongNum sqrtC = LongNumSqrt(one_squared * LongNum(10005));
    LongNum answer = Q * LongNum(426880) * sqrtC / T;
    set_accuracy(digits);
    answer /= exponentiation(LongNum(10), digits);
    set_accuracy(tmp);
    return answer;
}


int main()
{
    std::cout << "please enter the required number of decimal places" << std::endl;
    unsigned len;
    std::cin >> len;
    auto start = std::chrono::steady_clock::now();
    LongNum x = Chudnovsky(len);
    auto end = std::chrono::steady_clock::now();
    std::cout << "the computation time for the first " << len << " digits of pi was: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms" << std::endl;
    set_accuracy(len);
    display(x);
}