#ifndef L1Trigger_L1TGlobal_ADTCondition_h
#define L1Trigger_L1TGlobal_ADTCondition_h

/**
 * \class ADTCondition
 *
 * Description: evaluation of a CondADT condition.
 */

// system include files
#include <iosfwd>
#include <string>

// user include files
//   base classes
#include "L1Trigger/L1TGlobal/interface/ConditionEvaluation.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"

// forward declarations
class GlobalCondition;
class ADTTemplate;

namespace l1t {

  class L1Candidate;
  class GlobalBoard;

  // class declaration
  class ADTCondition : public ConditionEvaluation {
  public:
    /// constructors
    ///     default
    ADTCondition();

    ///     from base template condition (from event setup usually)
    ADTCondition(const GlobalCondition*, const GlobalBoard*);

    // copy constructor
    ADTCondition(const ADTCondition&);
    // destructor
    ~ADTCondition() override;

    // assign operator
    ADTCondition& operator=(const ADTCondition&);

    /// the core function to check if the condition matches
    const bool evaluateCondition(const int bxEval) const override;

    /// print condition
    void print(std::ostream& myCout) const override;

    ///   get / set the pointer to a Condition
    inline const ADTTemplate* gtADTTemplate() const { return m_gtADTTemplate; }

    void setGtADTTemplate(const ADTTemplate*);

    ///   get / set the pointer to GTL
    inline const GlobalBoard* gtGTB() const { return m_gtGTB; }

    void setuGtB(const GlobalBoard*);

    //get / set ADT model version
    inline const std::string gtModelVerion() const { return m_ADTmodelversion; }

    void setModelVersion(const std::string modelversionname);

  private:
    /// copy function for copy constructor and operator=
    void copy(const ADTCondition& cp);

    /// pointer to a ADTTemplate
    const ADTTemplate* m_gtADTTemplate;

    /// pointer to uGt GlobalBoard, to be able to get the trigger objects
    const GlobalBoard* m_gtGTB;

    //to set modelversion from globalboard<-globalproducer<-config
    std::string m_ADTmodelversion = "NONE";
  };

}  // namespace l1t
#endif
