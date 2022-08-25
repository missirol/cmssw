#ifndef HLTrigger_HLTfilters_TriggerExpressionPathReader_h
#define HLTrigger_HLTfilters_TriggerExpressionPathReader_h

#include <vector>
#include <string>

#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"

namespace triggerExpression {

  class PathReader : public Evaluator {
  public:
    PathReader(const std::string& pattern)
        : m_pattern{pattern},
          m_triggers{},
          m_triggersAfterMasking{},
          m_initialised{false},
          m_useTriggersAfterMasking{false} {}

    bool operator()(const Data& data) const override;

    void init(const Data& data) override;

    std::vector<std::string> patterns() const override { return std::vector<std::string>{m_pattern}; }

    void dump(std::ostream& out) const override;
    void dump(std::ostream& out, bool const) const;

    bool can_mask() const override { return true; }

    void mask(Evaluator* eval) override;

    std::vector<std::pair<std::string, unsigned int>> triggers() const { return m_triggers; }
    std::vector<std::pair<std::string, unsigned int>> triggersAfterMasking() const { return m_triggersAfterMasking; }

    void maskTriggers(PathReader const&);

    bool useTriggersAfterMasking() const { return m_useTriggersAfterMasking; }
    void useTriggersAfterMasking(bool const foo) { m_useTriggersAfterMasking = foo; }

  private:
    std::string m_pattern;
    std::vector<std::pair<std::string, unsigned int>> m_triggers;
    std::vector<std::pair<std::string, unsigned int>> m_triggersAfterMasking;
    bool m_initialised;
    bool m_useTriggersAfterMasking;
  };

}  // namespace triggerExpression

#endif  // HLTrigger_HLTfilters_TriggerExpressionPathReader_h
