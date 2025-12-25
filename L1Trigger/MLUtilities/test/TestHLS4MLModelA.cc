#include <any>
#include <array>

#include "L1Trigger/MLUtilities/interface/HLS4MLModelWrapper.h"

class TestHLS4MLModelA : public l1t::HLS4MLModel {
private:
  using inputs_t = std::array<int, 3>;
  using output_t = int;
  inputs_t inputs_;
  output_t output_;

public:
  TestHLS4MLModelA() : inputs_{0, 0, 0}, output_{0} {}

  virtual ~TestHLS4MLModelA() = default;

  virtual void prepare_input(std::any input) { inputs_ = std::any_cast<inputs_t>(input); }

  virtual void predict() { output_ = inputs_[0] * inputs_[1] + inputs_[2]; }

  virtual void read_result(std::any result) {
    output_t* output = std::any_cast<output_t*>(result);
    *output = output_;
  }
};

extern "C" l1t::HLS4MLModel* create_model() { return new TestHLS4MLModelA; }

extern "C" void destroy_model(l1t::HLS4MLModel* m) { delete m; }
