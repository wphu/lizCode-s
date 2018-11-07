#include "Timer.h"

#include <iomanip>
#include <string>

#include "Tools.h"

using namespace std;

Timer::Timer()
{
    name_ = "";
    time_acc_=0.0;
}

Timer::~Timer()
{
}

void Timer::init(string name)
{
    last_start_ = clock();
    name_ = name;
}


void Timer::update()
{
    time_acc_ +=  clock()-last_start_;
    last_start_ = clock();
}

void Timer::restart()
{
    last_start_ = clock();
}

void Timer::print(double tot)
{
    if ((time_acc_>0.) && (name_!=""))
        MESSAGE("\t" << setw(12) << name_ << "\t" << time_acc_  << "\t(" << 100.0*time_acc_/tot << "%)");
}

void Timer::print()
{
    MESSAGE("\t" << setw(12) << "The total time: " << "\t" << time_acc_ );
}

string Timer::getDateTime()
{
    long long d, h, m, s;
    long long time_acc_int = time_acc_;
    s = time_acc_int % 60;
    m = (time_acc_int / 60) % 60;
    h = (time_acc_int / 3600) % 24;
    d = time_acc_int / (3600 * 24);
    return to_string(d) + " d  " + to_string(h) + " h  " + to_string(m) + " m  " + to_string(s) + " s";
}