#ifndef LONGNUMSKKK_LONGNUMSKKK_H
#define LONGNUMSKKK_LONGNUMSKKK_H

// includes
#include <vector>
#include <iostream>

// "LongNum" type structure
struct LongNum
{
public:

    // constructors
    explicit LongNum(); // empty (init 0 by default)
    explicit LongNum(long long number);
    explicit LongNum(const std::string &s);

    // helpers
    friend void        display  (const LongNum &a);
    friend std::string to_string(const LongNum &a);

    // math functions
    friend LongNum     LongNumAbs (const LongNum &a); // take the modulus of the number
    friend LongNum     LongNumSqrt(const LongNum &a); // take the square root of the number
    friend LongNum     exponentiation(const LongNum &a, unsigned degree); // power
    static LongNum     Chudnovskkky(unsigned digits); // top secret

    // operators
    friend LongNum operator +  (const LongNum &a, const LongNum &b);
    friend LongNum operator -  (const LongNum &a, const LongNum &b);
    friend LongNum operator *  (const LongNum &a, const LongNum &b);
    friend LongNum operator /  (const LongNum &a, const LongNum &b);
    LongNum&       operator =  (const LongNum &b);
    LongNum&       operator += (const LongNum &b);
    LongNum&       operator -= (const LongNum &b);
    LongNum&       operator *= (const LongNum &b);
    LongNum&       operator /= (const LongNum &b);
    friend bool    operator == (const LongNum &a, const LongNum &b);
    friend bool    operator != (const LongNum &a, const LongNum &b);
    friend bool    operator <  (const LongNum &a, const LongNum &b);
    friend bool    operator >  (const LongNum &a, const LongNum &b);
    friend bool    operator <= (const LongNum &a, const LongNum &b);
    friend bool    operator >= (const LongNum &a, const LongNum &b);


private:
    std::vector<int> number;
    std::vector<int> fraction;
    unsigned depth;
    bool positivity;


    LongNum(std::vector<int> number, std::vector<int> fraction, unsigned depth, bool positivity);
};

// some other helpers
LongNum operator ""_skkk (const char *string);
void set_accuracy(unsigned value); // sets the displayed number of decimal places
void set_limitation(unsigned value); // sets the number of decimal places involved in calculations
void check_accuracy(); // shows the accuracy value
void check_limitation(); // shows the limitation value
unsigned get_accuracy();


#endif //LONGNUMSKKK_LONGNUMSKKK_H