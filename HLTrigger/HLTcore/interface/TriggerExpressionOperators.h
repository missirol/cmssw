#ifndef HLTrigger_HLTfilters_TriggerExpressionOperators_h
#define HLTrigger_HLTfilters_TriggerExpressionOperators_h

#include <memory>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"

namespace triggerExpression {

  // abstract unary operator
  class UnaryOperator : public Evaluator {
  public:
    UnaryOperator(Evaluator* arg) : m_arg(arg) {}

    // initialize the depending modules
    void init(const Data& data) override { m_arg->init(data); }

    // mask the depending modules
    void mask(Evaluator* arg) {
      if(arg != nullptr and not arg->can_mask())
        edm::LogWarning("NoOpEvaluatorMask") << "\tEvaluator::mask(arg) call is a no-op:"
          << " arg Evaluator cannot apply masks (arg->dump = \"" << *arg << "\", after masking)";
      m_arg->mask(arg);
    }

    // return the patterns from the depending modules
    std::vector<std::string> patterns() const override { return m_arg->patterns(); }

  protected:
    std::unique_ptr<Evaluator> m_arg;
  };

  // abstract binary operator
  class BinaryOperator : public Evaluator {
  public:
    BinaryOperator(Evaluator* arg1, Evaluator* arg2) : m_arg1(arg1), m_arg2(arg2) {}

    // initialize the depending modules
    void init(const Data& data) override {
      m_arg1->init(data);
      m_arg2->init(data);
    }

    // mask the depending modules
    void mask(Evaluator* arg) {
      if(arg != nullptr and not arg->can_mask())
        edm::LogWarning("NoOpEvaluatorMask") << "\tEvaluator::mask(arg) call is a no-op:"
          << " arg Evaluator cannot apply masks (arg->dump = \"" << *arg << "\", after masking)";
      m_arg1->mask(arg);
      m_arg2->mask(arg);
    }

    // return the patterns from the depending modules
    std::vector<std::string> patterns() const override {
      std::vector<std::string> patterns = m_arg1->patterns();
      auto patterns2 = m_arg2->patterns();
      patterns.insert(patterns.end(), std::make_move_iterator(patterns2.begin()), std::make_move_iterator(patterns2.end()));
      return patterns;
    }

  protected:
    std::unique_ptr<Evaluator> m_arg1;
    std::unique_ptr<Evaluator> m_arg2;
  };

  // concrete operators

  class OperatorNot : public UnaryOperator {
  public:
    OperatorNot(Evaluator* arg) : UnaryOperator(arg) {}

    bool operator()(const Data& data) const override { return not(*m_arg)(data); }

    void dump(std::ostream& out) const override {
      out << "NOT ";
      m_arg->dump(out);
    }
  };

  class OperatorAnd : public BinaryOperator {
  public:
    OperatorAnd(Evaluator* arg1, Evaluator* arg2) : BinaryOperator(arg1, arg2) {}

    bool operator()(const Data& data) const override {
      // force the execution of both arguments, otherwise prescalers won't work properly
      bool r1 = (*m_arg1)(data);
      bool r2 = (*m_arg2)(data);
      return r1 and r2;
    }

    void dump(std::ostream& out) const override {
      m_arg1->dump(out);
      out << " AND ";
      m_arg2->dump(out);
    }
  };

  class OperatorOr : public BinaryOperator {
  public:
    OperatorOr(Evaluator* arg1, Evaluator* arg2) : BinaryOperator(arg1, arg2) {}

    bool operator()(const Data& data) const override {
      // force the execution of both arguments, otherwise prescalers won't work properly
      bool r1 = (*m_arg1)(data);
      bool r2 = (*m_arg2)(data);
      return r1 or r2;
    }

    void dump(std::ostream& out) const override {
      m_arg1->dump(out);
      out << " OR ";
      m_arg2->dump(out);
    }
  };

  class OperatorXor : public BinaryOperator {
  public:
    OperatorXor(Evaluator* arg1, Evaluator* arg2) : BinaryOperator(arg1, arg2) {}

    bool operator()(const Data& data) const override {
      // force the execution of both arguments, otherwise prescalers won't work properly
      bool r1 = (*m_arg1)(data);
      bool r2 = (*m_arg2)(data);
      return r1 xor r2;
    }

    void dump(std::ostream& out) const override {
      m_arg1->dump(out);
      out << " XOR ";
      m_arg2->dump(out);
    }
  };

  class OperatorMasking : public BinaryOperator {
  public:
    OperatorMasking(Evaluator* arg1, Evaluator* arg2) : BinaryOperator(arg1, arg2) {}

    bool operator()(const Data& data) const override {
      // force the execution of both arguments, otherwise prescalers won't work properly
      //!! (*m_arg2)(data); ??
      return (*m_arg1)(data);
    }

    void init(const Data& data) override {
      m_arg1->init(data);
      m_arg2->init(data);
      m_arg1->mask(m_arg2.get());
    }

    void dump(std::ostream& out) const override {
      m_arg1->dump(out);
    }
  };

}  // namespace triggerExpression

#endif  // HLTrigger_HLTfilters_TriggerExpressionOperators_h
