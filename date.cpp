//
// Created by Jonas Wolff on 31/08/2024.
//

#include "date.h"
#include <numeric>

using namespace std;
using namespace chrono;

date::date(){};

date::date(const int &iyear, const int &imonth, const int &iday) : year_month_day(
                chrono::year(iyear),
                chrono::month(imonth),
                chrono::day(iday)
                ){};

date::date(chrono::year cYear, chrono::month cMonth, chrono::day cDay) : year_month_day(
        cYear,
        cMonth,
        cDay){};

date::date(sys_days sysDay) : year_month_day(sysDay) {};

date doubleToDate(const num &daysSinceEpoch){
    int a = floor(daysSinceEpoch);
    days b(a);
    return epochDate + b;
}

date::date(const num &daysSinceEpoch) : date(doubleToDate(daysSinceEpoch)){};


year_month_day dateStringParser(string &myDate){
    if (regex_match(myDate, regex("[0-9][0-9]/[0-9][0-9]/[0-9]{4}"))){
        day myDay(stoi(myDate.substr(0,2)));
        month myMonth(stoi(myDate.substr(3,2)));
        year myYear(stoi(myDate.substr(6,4)));
        return {myYear, myMonth, myDay};
    } else if (myDate.size() < 6) {
        return epochDate;
    }
    throw runtime_error("Date could not be processed");
}

date::date(std::string &myDate) : chrono::year_month_day(dateStringParser(myDate)) {};

days date::operator- (const date &rhs) const{
    return sys_days(*this) - sys_days(rhs);
    //return duration_cast<days>((sys_days(*this) - sys_days(rhs)));
}

date date::operator-(const chrono::days &rhs) const {
    return sys_days(*this) - rhs;
}

date date::operator+(const chrono::days &rhs) const {
    return sys_days(*this) + rhs;
}

int days_in_month(int &&year, int &&month) {

    if (month < 1 or month > 12)
        throw runtime_error("Invalid month");

    switch (month) {
        case 1: case 3: case 5: case 7: case 8: case 10: case 12:
            return 31;
        case 4: case 6: case 9: case 11:
            return 30;
        case 2:
            // Leap year check
            if ((year % 4 == 0 && year % 100 != 0) or (year % 400 == 0))
                return 29;
            else
                return 28;
        default:
            throw runtime_error("Mistake was made in days_in_month");
    }
    throw runtime_error("Mistake was made in days_in_month");
}


int days_in_year(int &&year) {
    return ((year % 4 == 0 and year % 100 != 0) or year % 400 == 0) ? 366 : 365;
}

int days_left_inYear(const date &myDate){

    //return (date(int(myDate.year())+1, 1, 1) - myDate).count();

    int year(myDate.year());
    unsigned int month(myDate.month());
    unsigned int day(myDate.day());

    constexpr int dec = 31;
    constexpr int nov = 31 + 30;
    constexpr int oct = 31 + 30 + 31;
    constexpr int sep = 31 + 30 + 31 + 30;
    constexpr int aug = 31 + 30 + 31 + 30 + 31;
    constexpr int jul = 31 + 30 + 31 + 30 + 31 + 31;
    constexpr int jun = 31 + 30 + 31 + 30 + 31 + 31 + 30;
    constexpr int may = 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31;
    constexpr int apr = 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30;
    constexpr int mar = 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 + 31;
    constexpr int feb = 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 + 31 + 28;
    constexpr int jan = 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 + 31 + 28 + 31;

    switch (month) {
        case 12:
            return dec - day;
        case 11:
            return nov - day;
        case 10:
            return oct - day;
        case 9:
            return sep - day;
        case 8:
            return aug - day;
        case 7:
            return jul - day;
        case 6:
            return jun - day;
        case 5:
            return may - day;
        case 4:
            return apr - day;
        case 3:
            return mar - day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
        case 2:
            return feb - day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
        case 1:
            return jan - day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
    }
}

int days_since_1stJan(const date &myDate) {

    //return (myDate - date(myDate.year(), January, 1d)).count();

    int year(myDate.year());
    unsigned int month(myDate.month());
    unsigned int day(myDate.day());

    constexpr int jan = 31;
    constexpr int feb = 31 + 28;
    constexpr int mar = 31 + 28 + 31;
    constexpr int apr = 31 + 28 + 31 + 30;
    constexpr int may = 31 + 28 + 31 + 30 + 31;
    constexpr int jun = 31 + 28 + 31 + 30 + 31 + 30;
    constexpr int jul = 31 + 28 + 31 + 30 + 31 + 30 + 31;
    constexpr int aug = 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31;
    constexpr int sep = 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30;
    constexpr int oct = 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31;
    constexpr int nov = 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30;

    switch (month) {
        case 1:
            return day;
        case 2:
            return jan + day;
        case 3:
            return feb + day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));;
        case 4:
            return mar + day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
        case 5:
            return apr + day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
        case 6:
            return may + day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
        case 7:
            return jun + day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
        case 8:
            return jul + day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
        case 9:
            return aug + day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
        case 10:
            return sep + day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
        case 11:
            return oct + day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
        case 12:
            return nov + day + (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
    }
}


num date::yearsuntil(const date &rhs) const {

    // If in same year count days
    if (rhs.year() == this->year())
        return num(abs((sys_days(rhs) - sys_days(*this)).count()))/days_in_year(int(year()));// daysinyear.count();
    else {
        const bool thisisfirst = this->year() < rhs.year();
        const date &firstdate = (thisisfirst) ? *this : rhs;
        const date &latestdate = (thisisfirst) ? rhs : *this;

        // Calculate time from first day to 1st January
        //const date next1stJan(firstdate.year() + years(1), January, 1d);
        //const days daysto1stJan = next1stJan - firstdate;
        //date first1stJan(firstdate.year(), January, 1d);
        //num firstYearFrac = days_since_1stJan(firstdate)/days_in_year(int(firstdate.year())); //(next1stJan-first1stJan).count();
        const num firstYearFrac = num(days_left_inYear(firstdate))/days_in_year(int(firstdate.year())); //(next1stJan-first1stJan).count();

        // Calculate time from 1st January to latestdate
        //const date first1stJan = date(latestdate.year() , January, 1d);
        //const days dayssince1stJan = latestdate - first1stJan;
        //next1stJan = date(latestdate.year() + years(1), January, 1d);
        //const num lastYearFrac = num(dayssince1stJan.count())/days_in_year(int(latestdate.year())) ;//(next1stJan - first1stJan).count();
        const num lastYearFrac = num(days_since_1stJan(latestdate))/days_in_year(int(latestdate.year())) ;//(next1stJan - first1stJan).count();

        // Calculate total number of years
        const years yeardiff = abs(firstdate.year() - latestdate.year());
        if (thisisfirst)
            return (firstYearFrac + lastYearFrac) + (yeardiff.count() - 1);
        else
            return -((firstYearFrac + lastYearFrac) + (yeardiff.count() - 1));

    }

}

date::operator num() const {
    return ((*this) - epochDate).count() ;
}

string date::tostring() const {
    string res = "";
    res += to_string(year().operator int());
    res += "-";
    unsigned int myi(month());
    res += ((myi<10) ? "0" : "") + to_string(myi);
    res += "-";
    myi = day().operator unsigned int();
    res += ((myi<10) ? "0" : "") + to_string(myi);
    return res;
}
