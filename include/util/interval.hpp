#pragma once

#include <array>

namespace pmondriaan {

struct interval {
    long low;
    long high;
    std::array<long, 2> values = {low, high};
    long length() { return high - low; }
    long operator()(long part) { return values[part]; }
};

} // namespace pmondriaan
