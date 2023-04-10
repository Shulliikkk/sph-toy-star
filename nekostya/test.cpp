#include <iostream>
#include "Timer.h"

int main() {
    unsigned counter;
    {
        Timer<std::chrono::nanoseconds> t(&counter);
        for (auto i = 0; i < 1'000'000; i++) {}
    }
    std::cout << counter << std::endl;

    {
        Timer t;
        for (auto i = 0; i < 1'000'000; i++) {}
    }

    return 0;
}
