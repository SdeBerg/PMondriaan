#pragma once

namespace pmondriaan {

/**
 * A work_item stores the start and end location of the vertices to be split and
 * the highest en lowest label to be assigned.
 */
class work_item {
  public:
    work_item(long start, long end, long label_low, long label_high, long weight)
    : start_(start), end_(end), label_low_(label_low), label_high_(label_high),
      weight_(weight) {}

    long start() { return start_; }
    long end() { return end_; }
    long label_low() { return label_low_; }
    long label_high() { return label_high_; }
    long weight() { return weight_; }

  private:
    long start_;
    long end_;
    long label_low_;
    long label_high_;
    long weight_;
};


} // namespace pmondriaan