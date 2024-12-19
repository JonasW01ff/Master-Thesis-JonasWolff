//
// Created by Jonas Wolff on 31/08/2024.
//
#ifndef CODE_DATE_H
#define CODE_DATE_H
#include <chrono>
#include <regex>
#include "settings.h"

using namespace std;

class date : public chrono::year_month_day {

public:

    date();

    // Overloaded constructor
    date(const int &year, const int &month, const int &day);

    date(chrono::year cYear, chrono::month cMonth, chrono::day cDay);

    date(chrono::sys_days sDays);

    date(string &myDate);

    date(const num &daysSinceEpoch);

    date& operator=(const date&) = default;

    chrono::days operator- (const date &rhs) const;

    date operator- (const chrono::days &rhs) const;

    date operator+ (const chrono::days &rhs) const;

    date operator+ (const chrono::years &rhs) const{
        return {year()+rhs, month(), day()};
    }

    explicit operator num() const;

    num yearsuntil(const date &rhs) const;

    string tostring() const;


};

static date epochDate(2024y,chrono::January, 1d);

#endif //CODE_DATE_H
