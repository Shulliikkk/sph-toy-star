#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <iostream>

template<typename Duration_Type = std::chrono::microseconds>
class Timer {
private:
    std::chrono::steady_clock::time_point start;
    unsigned* counter;
public:
    Timer();
    Timer(unsigned* counter); 
    ~Timer();
};

#endif //TIMER_H

