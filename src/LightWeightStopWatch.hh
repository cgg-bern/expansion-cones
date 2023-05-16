#pragma once

#include <chrono>


class LightWeightStopWatch
{
public:

    LightWeightStopWatch(bool start = true){
        if(start){
            reset();
        }
    }

    using Clock = std::chrono::high_resolution_clock;

    void reset(){
        start_time_ = Clock::now();
        last_lap_time_ = start_time_;
    }


    template<typename duration_precision = std::chrono::microseconds>
    float lap_duration(){
        auto current_time = Clock::now();
        auto duration = (float)std::chrono::duration_cast<duration_precision>(current_time - last_lap_time_).count() / duration_precision(std::chrono::seconds(1)).count();
        last_lap_time_ = current_time;
        return duration;
    }

    template<typename duration_precision = std::chrono::microseconds>
    float total_duration(){
        auto current_time = Clock::now();
        return (float)std::chrono::duration_cast<duration_precision>(current_time - start_time_).count() / duration_precision(std::chrono::seconds(1)).count();
    }

private:

    std::chrono::time_point<Clock> start_time_;

    std::chrono::time_point<Clock> last_lap_time_;

};
