#ifndef HLTrigger_HLTfilters_TriggerExpressionL1uGTReader_h
#define HLTrigger_HLTfilters_TriggerExpressionL1uGTReader_h

#include <vector>
#include <string>

#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"

namespace triggerExpression {

  class L1uGTReader : public Evaluator {
  public:
    L1uGTReader(const std::string& pattern) : m_pattern(pattern), m_triggers() {}

    bool operator()(const Data& data) const override;

    void init(const Data& data) override;

    std::vector<std::string> patterns() const override { return std::vector<std::string>{m_pattern}; }

    void dump(std::ostream& out) const override;

    bool can_mask() const override { return true; }

    void mask(Evaluator* eval) override;

    std::vector<std::pair<std::string, unsigned int> > triggers() const { return m_triggers; }

    void maskTriggers(L1uGTReader const&);

  private:
    std::string m_pattern;
    std::vector<std::pair<std::string, unsigned int> > m_triggers;
  };

}  // namespace triggerExpression

#endif  // HLTrigger_HLTfilters_TriggerExpressionL1uGTReader_h
