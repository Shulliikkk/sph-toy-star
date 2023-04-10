template<typename Duration_Type>
Timer<Duration_Type>::Timer() : start(std::chrono::steady_clock::now()),
                 counter(nullptr) {
}

template<typename Duration_Type>
Timer<Duration_Type>::Timer(unsigned* counter) : start(std::chrono::steady_clock::now()),
                                  counter(counter) {
}

template<typename Duration_Type>
Timer<Duration_Type>::~Timer() {
    auto delta_point = std::chrono::steady_clock::now() - start;
    auto delta_time = std::chrono::duration_cast<Duration_Type>(delta_point).count();

    if (counter != nullptr) {
        *counter = delta_time;
    } else {
        std::cout << "Timer counted " << delta_time << std::endl;
    }
}

