#ifndef L1Trigger_L1TGlobal_AXOL1TLCondition_h
#define L1Trigger_L1TGlobal_AXOL1TLCondition_h

/**
 * \class AXOL1TLCondition
 *
 * Description: evaluation of a CondAXOL1TL condition.
 */

#include <ostream>
#include <string>
#include <utility>

#include "hls4ml/ModelWrapper.h"

#include "L1Trigger/L1TGlobal/interface/ConditionEvaluation.h"

class GlobalCondition;
class AXOL1TLTemplate;

namespace l1t {

  class GlobalBoard;

  // class declaration
  class AXOL1TLCondition : public ConditionEvaluation {
  public:
    // default constructor
    AXOL1TLCondition();

    // constructor from base template condition (from event setup usually)
    AXOL1TLCondition(const GlobalCondition*, const GlobalBoard*);

    // copy constructor
    AXOL1TLCondition(const AXOL1TLCondition&);

    // destructor
    ~AXOL1TLCondition() override = default;

    // assign operator
    AXOL1TLCondition& operator=(const AXOL1TLCondition&);

    // the core function to check if the condition matches
    const bool evaluateCondition(const int bxEval) const override;

    // print condition
    void print(std::ostream& myCout) const override;

    // get/set the pointer to a Condition
    const AXOL1TLTemplate* gtAXOL1TLTemplate() const { return m_gtAXOL1TLTemplate; }

    void setGtAXOL1TLTemplate(const AXOL1TLTemplate* ptr) { m_gtAXOL1TLTemplate = ptr; }

    // get/set the pointer to GTL
    const GlobalBoard* gtGTB() const { return m_gtGTB; }

    void setuGtB(const GlobalBoard* ptr) { m_gtGTB = ptr; }

    // name of the model
    std::string const& model_name() const { return m_model_wrapper.model_name(); }

    // get/set score value
    float getScore() const { return m_saved_score; }

    void setScore(const float scoreval) const { m_saved_score = scoreval; }

  private:
    // copy function for copy constructor and operator=
    void copy(const AXOL1TLCondition& cp);

    // pointer to a AXOL1TLTemplate
    const AXOL1TLTemplate* m_gtAXOL1TLTemplate;

    // pointer to uGt GlobalBoard, to be able to get the trigger objects
    const GlobalBoard* m_gtGTB;

    static constexpr char const* kModelNamePrefix = "GTADModel_";

    hls4mlEmulator::ModelWrapper m_model_wrapper;

    // output score for possible score saving
    mutable float m_saved_score;
  };

}  // namespace l1t

#endif
