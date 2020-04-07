#pragma once

#include <array>

namespace pmondriaan {

struct interval {
    int low;
    int high;
    std::array<int, 2> values = {low, high};
    int length() { return high - low; }
    int operator()(int part) { return values[part]; }
};

} // namespace pmondriaan
