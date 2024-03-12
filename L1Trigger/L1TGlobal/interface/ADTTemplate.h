#ifndef L1Trigger_L1TGlobal_ADTTemplate_h
#define L1Trigger_L1TGlobal_ADTTemplate_h

/**
 * \class ADTTemplate
 *
 *
 * Description: L1 Global Trigger ADT template.
 *
 * \author: Melissa Quinnan (UC San Diego)
 *
 */

// system include files
#include <string>
#include <iosfwd>

// user include files

//   base class
#include "L1Trigger/L1TGlobal/interface/GlobalCondition.h"

// forward declarations

// class declaration
class ADTTemplate : public GlobalCondition {
public:
  // constructor
  ADTTemplate();

  // constructor
  ADTTemplate(const std::string&);

  // constructor
  ADTTemplate(const std::string&, const l1t::GtConditionType&);

  // copy constructor
  ADTTemplate(const ADTTemplate&);

  // destructor
  ~ADTTemplate() override;

  // assign operator
  ADTTemplate& operator=(const ADTTemplate&);

  // typedef for a single object template
  struct ObjectParameter {
    int minADTThreshold;
    int maxADTThreshold;
  };

public:
  inline const std::vector<ObjectParameter>* objectParameter() const { return &m_objectParameter; }

  /// set functions
  void setConditionParameter(const std::vector<ObjectParameter>& objParameter);

  /// print the condition
  void print(std::ostream& myCout) const override;

  /// output stream operator
  friend std::ostream& operator<<(std::ostream&, const ADTTemplate&);

private:
  /// copy function for copy constructor and operator=
  void copy(const ADTTemplate& cp);

  /// variables containing the parameters
  std::vector<ObjectParameter> m_objectParameter;
};

#endif
