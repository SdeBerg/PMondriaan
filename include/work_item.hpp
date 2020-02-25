#pragma once

namespace pmondriaan {

/**
 * A work_item stores the start and end location of the vertices to be split and
 * the highest en lowest label to be assigned.
 */
class work_item {
  public:
    work_item(int start, int end, int label_low, int label_high, long weight)
    : start_(start), end_(end), label_low_(label_low), label_high_(label_high),
      weight_(weight) {}

    int start() { return start_; }
    int end() { return end_; }
    int label_low() { return label_low_; }
    int label_high() { return label_high_; }
    long weight() { return weight_; }

  private:
    int start_;
    int end_;
    int label_low_;
    int label_high_;
    long weight_;
};


} // namespace pmondriaan