//
// Created by Jonas Wolff on 24/09/2024.
//

#ifndef CODE_NUMBERS_H
#define CODE_NUMBERS_H
#include <stdexcept>
#include "concepts"
#include <cmath>
#include "iostream"

using namespace std;

template<class T1>
class number {

public:

    virtual T1 operator+(const T1 &rhs) const =0 ;

    T1 operator+=(const T1 &rhs){
        *this = *this + rhs ;
    }

    virtual T1 operator-(const T1 &rhs) const =0 ;

    T1 operator-=(const T1 &rhs){
        *this = *this - rhs;
    }

    virtual T1 operator*(const T1 &rhs)  const =0 ;

    T1 operator*=(const T1 &rhs){
        *this = *this * rhs;
    }

    virtual T1 operator/(const T1 &rhs) const =0 ;

    T1 operator/=(const T1 &rhs){
        *this = *this / rhs;
    }

    virtual T1 operator-() const =0;

    virtual T1 operator+() const = 0;

    template<class T2>
    friend T1 operator+(const T2 &lhs, const T1 &rhs){
        T1 myNumber(lhs);
        return myNumber+rhs;
    }

    template<class T2>
    friend T1 operator-(const T2 &lhs, const T1 &rhs){
        T1 myNumber(lhs);
        return myNumber-rhs;
    }

    template<class T2>
    friend T1 operator*(const T2 &lhs, const T1 &rhs){
        T1 myNumber(lhs);
        return myNumber*rhs;
    }

    template<class T2>
    friend T1 operator/(const T2 &lhs, const T1 &rhs){
        T1 myNumber(lhs);
        return myNumber/rhs;
    }

    virtual T1 pow(const T1 &rhs) =0;

    virtual T1 log() =0;

    virtual T1 exp() =0;

    virtual T1 sqrt() =0;

    virtual T1 abs() =0;

    virtual T1 round() =0;

};

template<class T2, class T3>
requires derived_from<T2, number<T2>> and std::is_arithmetic_v<T3>
T2 pow(T2 &myNumber, T3 p){
    T2 pNumber = p;
    return myNumber.pow(pNumber);};

template<class T1>
requires derived_from<T1, number<T1>>
T1 log(const T1 &myNumber){
    T1 myRes = myNumber;
    return myRes.log();};

template<class T1>
requires derived_from<T1, number<T1>>
T1 exp(const T1 &myNumber){
    T1 myRes = myNumber;
    return myRes.exp();};

template<class T1>
requires derived_from<T1, number<T1>>
T1 sqrt(const T1 &myNumber){
    T1 myRes = myNumber;
    return myRes.sqrt();};

template<class T1>
requires derived_from<T1, number<T1>>
T1 abs(const T1 &myNumber){
    T1 myRes = myNumber;
    return myRes.abs();};

template<class T1>
requires derived_from<T1, number<T1>>
T1 round(const T1 &myNumber){
    T1 myRes = myNumber;
    return myRes.round();};

template<class T1, class T2>
requires derived_from<T1, number<T1>> and std::is_arithmetic_v<T2>
T1 max(const T1 &lhs, const T2 &rhs){
    return (lhs > rhs) ? lhs : T1(rhs);};

template<class T1, class T2>
requires derived_from<T1, number<T1>> and std::is_arithmetic_v<T2>
T1 max(const T2 &lhs, const T1 &rhs){
    return (lhs > rhs) ? T1(lhs) : rhs;};

template<class T1, class T2>
requires derived_from<T1, number<T1>> and std::is_arithmetic_v<T2>
T1 min(const T1 &lhs, const T2 &rhs){
    return (lhs < rhs) ? lhs : T1(rhs);};

template<class T1, class T2>
requires derived_from<T1, number<T1>> and std::is_arithmetic_v<T2>
T1 min(const T2 &lhs, const T1 &rhs){
    return (lhs < rhs) ? T1(lhs) : rhs;};

template<class T1>
requires std::is_arithmetic_v<T1>
class dual : public number<dual<T1>>{

    T1 real = 0;
    T1 epsilon = 0;

public:

    // Construct dual numbers manually
    dual<T1>(const T1 &pReal,const T1 &pEpsilon){
        real = pReal;
        epsilon = pEpsilon;
    }

    dual<T1>(T1 &&pReal, T1 &&pEpsilon){
        real = std::move(pReal);
        epsilon = std::move(pEpsilon);
    }

    // Constructor
    dual<T1>(dual<T1> &rhs) = default;

    dual<T1>(dual<T1> &&rhs) {
        real = std::move(rhs.real);
        epsilon = std::move(rhs.epsilon);
    };

    dual<T1>(const dual<T1> &rhs) {
        real = std::move(rhs.real);
        epsilon = std::move(rhs.epsilon);
    };

    // numeric to dual constructors
    dual<T1>(T1 rhs){ real = std::move(rhs);};

    template<class T2>
    requires std::is_arithmetic_v<T2> && (!same_as<T2, T1>)
    dual<T1>(T2 rhs){ real = std::move(rhs);};


    dual<T1>(){};

    ~dual<T1>(){};

    T1 getReal() const {return real;};
    T1 getEpsilon() const {return epsilon;};

    dual<T1> operator=(dual<T1> &rhs){
        epsilon = rhs.epsilon;
        real = rhs.real;
        return *this;
    }

    dual<T1>& operator=(const dual<T1> &rhs){
        epsilon = rhs.epsilon;
        real = rhs.real;
        return *this;
    }

    dual<T1>& operator=(dual<T1> &&rhs) noexcept{
        epsilon = rhs.epsilon;
        real = rhs.real;
        return *this;
    }

    dual<T1> operator-() const{
        dual<T1> myDual(*this);
        myDual.epsilon = -myDual.epsilon;
        myDual.real = -myDual.real;
        return myDual;
    }

    dual<T1> operator+() const{
        dual<T1> myDual(*this);
        myDual.epsilon = +myDual.epsilon;
        myDual.real = +myDual.real;
        return myDual;
    }

    dual<T1> operator+(const dual &rhs) const{
        dual<T1> myDual = *this;
        myDual.epsilon = epsilon + rhs.epsilon;
        myDual.real = real + rhs.real;
        return myDual;
    }

    dual<T1> operator-(const dual<T1> &rhs) const{
        dual<T1> myDual(this->real);
        myDual.epsilon = epsilon - rhs.epsilon;
        myDual.real = real - rhs.real;
        return myDual;
    }

    dual<T1> operator*(const dual<T1> &rhs) const{
        dual<T1> myDual = *this;
        myDual.epsilon = real*rhs.epsilon + rhs.real*epsilon;
        myDual.real = real * rhs.real;
        return myDual;
    }

    dual<T1> operator*(const T1 &rhs) const{
        dual<T1> myDual(rhs);
        return *this*myDual;
    }

    dual<T1> operator/(const dual<T1> &rhs) const{

        if (rhs.real == 0)
            throw logic_error("Division by zero");

        dual<T1> myDual(*this);

        myDual.epsilon = (epsilon * rhs.real - real * rhs.epsilon) / (rhs.real * rhs.real);
        myDual.real = real/rhs.real;
        return myDual;
    }

    // ">" Greater than operators
    bool operator>(const dual<T1> &rhs) const{
        return real > rhs.real;
    }

    bool operator>(dual<T1> &rhs) const{
        return real > rhs.real;
    }

    bool operator>(dual<T1> &rhs){
        return real > rhs.real;
    }

    // Greater than friend operators
    template<class T2>
    requires std::is_arithmetic_v<T2>
    friend bool operator>(const T2 &&lhs, const dual<T1> &&rhs){
        return lhs > rhs.real;
    }

    template<class T2>
    requires std::is_arithmetic_v<T2>
    friend bool operator>(const T2 &lhs, const dual<T1> &rhs){
        return lhs > rhs.real;
    }

    // Less than operators
    bool operator<(const dual<T1> &&rhs) const{
        return real < rhs.real;
    }

    bool operator<(const dual<T1> &rhs) const{
        return real < rhs.real;
    }

    bool operator<(dual<T1> &rhs) const{
        return real < rhs.real;
    }

    bool operator<(dual<T1> &rhs){
        return real < rhs.real;
    }

    // Less than friend operators
    template<class T2>
    requires std::is_arithmetic_v<T2>
    friend bool operator<(const T2 &&lhs, const dual<T1> &&rhs){
        return lhs < rhs.real;
    }

    template<class T2>
    requires std::is_arithmetic_v<T2>
    friend bool operator<(const T2 &lhs, const dual<T1> &rhs){
        return lhs < rhs.real;
    }

    // Less than or equal operators
    bool operator<=(const dual<T1> &rhs) const{
        return real <= rhs.real;
    }

    bool operator<=(dual<T1> &rhs) const{
        return real <= rhs.real;
    }

    // Greater than or equal operators
    bool operator>=(const dual<T1> &rhs) const{
        return real >= rhs.real;
    }

    bool operator>=(dual<T1> &rhs) const{
        return real >= rhs.real;
    }

    // Equality operators
    bool operator==(dual<T1> &rhs) const{
        return real == rhs.real;
    }

    bool operator==(dual<T1> &&rhs) const{
        return real == rhs.real;
    }

    // Equality friend operators
    template<class T2>
    requires std::is_arithmetic_v<T2>
    bool operator==(T2 &rhs) const{
        return real == rhs;
    }


    // Mathematical functions
    dual<T1> pow(const dual<T1> &rhs) {
        *this = ((this->log())*rhs).exp();
        return *this;
    }

    dual<T1> exp() {
        epsilon = std::exp(real)*epsilon;
        real = std::exp(real);
        return *this;
    }

    dual<T1> log() {
        epsilon = 1./real*epsilon;
        real = std::log(real);
        return *this;
    }

    dual<T1> sqrt() {
        epsilon = .5/std::sqrt(real)*epsilon;
        real = std::sqrt(real);
        return *this;
    }

    dual<T1> abs() {
        epsilon = ((real>0) - (real<0))*epsilon;
        real = std::abs(real);
        return *this;
    }

    dual<T1> round() {
        epsilon = std::round(epsilon);
        real = std::round(real);
        return *this;
    }

};

// output stream operator
template<class T1>
requires std::is_arithmetic_v<T1>
std::ostream& operator<<(std::ostream& os, const dual<T1>& myDual) {

    os << myDual.getReal() << ((myDual.getEpsilon() >= 0) ?  "+" : "-") << abs(myDual.getEpsilon()) << "\u03B5";
    return os;
}

// String function
template<class T1>
requires std::is_arithmetic_v<T1>
std::string to_string(const dual<T1>& myDual) {
    string myString = to_string(myDual.getReal()) +
            ((myDual.getEpsilon() >= 0) ?  "+" : "-") + to_string(abs(myDual.getEpsilon())) + "\u03B5";
    return myString;
}

#endif //CODE_NUMBERS_H
