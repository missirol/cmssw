#include <array>
#include <string>

#include <catch2/catch_test_macros.hpp>

#include "L1Trigger/MLUtilities/interface/HLS4MLModelWrapper.h"

TEST_CASE("HLS4MLModelWrapper basic functionalities", "[HLS4MLModelWrapper]") {
  std::string const libA{"libL1TriggerMLUtilitiesTestHLS4MLModelA"};
  std::string const libB{"libL1TriggerMLUtilitiesTestHLS4MLModelB"};

  l1t::HLS4MLModelWrapper mw1{libA};
  l1t::HLS4MLModelWrapper mw2{libB};
  l1t::HLS4MLModelWrapper mw3{libA};

  std::array<int, 3> mw1_inputs{1, 2, 3};
  std::array<int, 3> mw2_inputs{4, 5, 6};
  std::array<int, 3> mw3_inputs{7, 8, 9};

  int mw1_output{0};
  int mw2_output{0};
  int mw3_output{0};

  SECTION("Run inference for two different models") {
    mw1.prepare_input(mw1_inputs);
    mw1.predict();
    mw1.read_result(&mw1_output);
    REQUIRE(mw1_output == 5);

    mw2.prepare_input(mw2_inputs);
    mw2.predict();
    mw2.read_result(&mw2_output);
    REQUIRE(mw2_output == 34);
  }

  SECTION("Run inference for two instances of the same model") {
    mw1.prepare_input(mw1_inputs);
    mw3.prepare_input(mw3_inputs);

    mw1.predict();
    mw3.predict();

    mw1.read_result(&mw1_output);
    mw3.read_result(&mw3_output);

    REQUIRE(mw1_output == 5);
    REQUIRE(mw3_output == 65);
  }

  SECTION("Reset a model wrapper") { REQUIRE_NOTHROW(mw1.reset()); }

  SECTION("Change model, and rerun inference") {
    mw2.reset(libA);

    mw2.prepare_input(mw2_inputs);
    mw2.predict();
    mw2.read_result(&mw2_output);
    REQUIRE(mw2_output == 26);
  }
}
