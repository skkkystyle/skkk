#include "./include/LongNumskkk.hpp"


#define CHECKER(commands,x, y) do {try {commands; checking(to_string(x), y);} CATCHER
#define CMPRES(x) LongNum(std::__cxx11::to_string(x))
#define OUTPUT printf("::::::TEST_%d::::::", counter); if (counter / 10 == 0) {printf(":");}
#define CATCHER catch (const char* msg){OUTPUT; counter++; std::cout << " " << msg;}} while (0)
const std::string PI200 = "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196\n";


int counter = 0;
void checking(const std::string& skkk, const std::string& ans)
{
    OUTPUT
    counter++;
    if (skkk == ans)
    {
        printf(" TRUE\n");
    }
    else
    {
        printf(" FALSE\n");
    }
}


int main()
{
    // the results for comparison are obtained by the websites https://okcalc.com/ru/ and https://www.wolframalpha.com/
    // constructors testing
    LongNum a, b(0), c(-0), d("0,0"), e("-0.00");
    LongNum f(6606618), g(-9114382), h(992289 + 2052925);
    LongNum i("5350117"), j("-783.0194"), k("0.6186521"), l("-0.205186624682486"), m("229.5711");
    LongNum n = 7154528_skkk, o = 649.9057_skkk, p = 0.6449531_skkk;
    LongNum A, B, C;

    // tests 0 - 2 checking accuracy sync with real output & var initing & to_string func
    CHECKER(set_accuracy(13), l, "-0.2051866246824\n");
    CHECKER(set_accuracy(10), o, "649.9057000000\n");
    CHECKER(set_accuracy(15), l, "-0.205186624682486\n");

    // tests 3 - 5 checking operator +
    CHECKER(, f + a, "6606618.000000000000000\n");
    CHECKER(, g + h, "-6069168.000000000000000\n");
    CHECKER(, k + j, "-782.400747900000000\n");

    // tests 6 - 8 checking operator -
    CHECKER(, i - j, "5350900.019400000000000\n");
    CHECKER(, n - o, "7153878.094300000000000\n");
    CHECKER(, g - m, "-9114611.571100000000000\n");

    // tests 9 - 11 checking operator *
    CHECKER(, i * n, "38277561879776.000000000000000\n");
    CHECKER(, m * o, "149199.566445270000000\n");
    CHECKER(, c * f, "0.000000000000000\n");

    // tests 12 - 14 checking operator /
    CHECKER(, o / p, "1007.679007977479292\n");
    CHECKER(, i / j, "-6832.674899242598586\n");
    CHECKER(, h / n, "0.425634507265888\n");
    CHECKER(, g / a, "NORESULT");

    // tests 16 - 18 checking operator =
    CHECKER(A = f - k, A, to_string(f - k));
    CHECKER(B = h + n, B, to_string(h + n));
    CHECKER(C = o, C, to_string(o));

    // tests 19 - 22 checking +=, -=, *=, /=
    CHECKER(g += m; g += f, g, "-2507534.428900000000000\n");
    CHECKER(i -= j; h -= i, h, "-2305686.019400000000000\n");
    CHECKER(k *= f; i *= k, i, "21870188390429.511152197320000\n");
    CHECKER(g /= o; n /= g, n, "-1854.318917586846777\n");

    // tests 23 - 40 checking ==, !=, <, >, <=, >=
    CHECKER(set_accuracy(0), CMPRES(a == b), "1");
    CHECKER(, CMPRES(j == k), "0");
    CHECKER(, CMPRES(g == g), "1");
    CHECKER(, CMPRES(i != m), "1");
    CHECKER(, CMPRES(d != k), "1");
    CHECKER(, CMPRES(a != c), "0");
    CHECKER(, CMPRES(a < b), "0");
    CHECKER(, CMPRES(i < d), "0");
    CHECKER(, CMPRES(l < k), "1");
    CHECKER(, CMPRES(j > f), "0");
    CHECKER(, CMPRES(e > a), "0");
    CHECKER(, CMPRES(o > f), "0");
    CHECKER(, CMPRES(A <= B), "1");
    CHECKER(, CMPRES(B <= C), "0");
    CHECKER(, CMPRES(C <= A), "1");
    CHECKER(, CMPRES(i >= a), "1");
    CHECKER(, CMPRES(b >= m), "0");
    CHECKER(, CMPRES(a >= a), "1");

    // tests 40 - NUM checking math functions
    CHECKER(set_accuracy(15), LongNumAbs(h), "2305686.019400000000000\n");
    CHECKER(, LongNumAbs(e), "0.000000000000000\n");
    CHECKER(, LongNumAbs(i), "21870188390429.511152197320000\n");
    CHECKER(, LongNumSqrt(i), "4676557.322478738862763\n");
    CHECKER(, LongNumSqrt(o), "25.493248125729290\n");
    CHECKER(, LongNumSqrt(p), "0.803089721015030\n");
    CHECKER(, LongNumSqrt(g), "NORESULT");
    CHECKER(, exponentiation(o, 5), "115944921220582.995856197897841\n");
    CHECKER(, exponentiation(h, 7), "-346418569125010418246863833615170726059384575.689244878704677\n");
    CHECKER(, exponentiation(j, 3), "-480084369.663876941384000\n");

    // test 51 checking first 200 Pi decimal places;
    CHECKER(set_limitation(200); set_accuracy(200), LongNum::Chudnovskkky(200), PI200);

    std::cout << "FINISH" << std::endl;
}