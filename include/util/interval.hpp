#pragma once

namespace pmondriaan {
	
struct interval {
	int low;
	int high;
	int length() { return high - low; }
};
} // namespace pmondriaan
