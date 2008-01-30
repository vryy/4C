#ifdef CCADISCRET
#ifndef DRT_CONDITIONDEFINITION_H
#define DRT_CONDITIONDEFINITION_H

#include <string>
#include <vector>
#include <iostream>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

#include "drt_condition.H"
#include "drt_inputreader.H"


namespace DRT
{
  class Problem;
  class Discretization;

namespace INPUT
{

  /// pure virtual base class for all condition line parts
  /*!

    The aim of the machinery here is to describe how a valid DRT::Condition
    looks like and to be able to read and write them.

    In general it is assumed that one condition definition in the dat file
    consists of one line. This line might consist of various components, each
    one corresponding to one variable that has to be defined in the
    DRT::Condition object. This class does represent one such component.

    There are different types of variables stored in a DRT::Condition. And
    accordingly there are different subclasses of ConditionComponent that
    implement the reading and writing of these types.

    \see ConditionDefinition

    \author u.kue
    \date 01/08
   */
  class ConditionComponent
  {
  public:

    /// construct with the name of the corresponding variable in the condition
    explicit ConditionComponent(std::string name);

    /// virtual destructor is mandatory
    virtual ~ConditionComponent() {}

    /// write my part of the default (comment) line of the condition
    virtual void DefaultLine(std::ostream& stream) = 0;

    /// write my part of this actual condition
    virtual void Print(std::ostream& stream, const DRT::Condition* cond) = 0;

    /// read my part of an input line and put it into the condition
    virtual void Read(std::istringstream& condline,
                      Teuchos::RCP<DRT::Condition> condition) = 0;

    /// the name of my variable inside a condition
    std::string Name() const { return name_; }

  private:

    /// my condition variable name
    std::string name_;
  };


  /// condition line part that is a string out of a list of valid strings
  /*!
    \see ConditionDefinition

    \author u.kue
    \date 01/08
   */
  class StringConditionComponent : public ConditionComponent
  {
  public:
    StringConditionComponent(std::string name,
                             std::string defaultvalue,
                             const Teuchos::Array<std::string>& datfilevalues,
                             const Teuchos::Array<std::string>& condvalues,
                             bool optional=false);

    virtual void DefaultLine(std::ostream& stream);
    virtual void Print(std::ostream& stream, const DRT::Condition* cond);
    virtual void Read(std::istringstream& condline,
                      Teuchos::RCP<DRT::Condition> condition);

  private:
    std::string defaultvalue_;
    Teuchos::Array<std::string> datfilevalues_;
    Teuchos::Array<std::string> condvalues_;
    bool optional_;
  };


  /// condition line part that is a fixed string without any effect on the DRT::Condition
  /*!

    This is really just a separator at the input line. The DRT::Condition is
    not effected.

    The reason we need this is that the we specify the order of the input line
    part. It might be reasonable to specify names that have to appear in the
    dat file to enhance human readability.

    \see ConditionDefinition

    \author u.kue
    \date 01/08
   */
  class SeparatorConditionComponent : public ConditionComponent
  {
  public:
    SeparatorConditionComponent(std::string separator);

    virtual void DefaultLine(std::ostream& stream);
    virtual void Print(std::ostream& stream, const DRT::Condition* cond);
    virtual void Read(std::istringstream& condline,
                      Teuchos::RCP<DRT::Condition> condition);

  private:
    std::string separator_;
  };


  /// condition line part that is a single integer number
  /*!

    Integer numbers are special as these can be stored in fortran style in the
    dat file (internally we always store C style numbers) and sometimes the
    string "none" is allowed for -1.

    \see ConditionDefinition

    \author u.kue
    \date 01/08
   */
  class IntConditionComponent : public ConditionComponent
  {
  public:
    explicit IntConditionComponent(std::string name, bool fortranstyle=false, bool noneallowed=false);
    virtual void DefaultLine(std::ostream& stream);
    virtual void Print(std::ostream& stream, const DRT::Condition* cond);
    virtual void Read(std::istringstream& condline,
                      Teuchos::RCP<DRT::Condition> condition);

  private:
    bool fortranstyle_;
    bool noneallowed_;
  };


  /// condition line part that is a vector of integer numbers
  /*!

    Integer numbers are special as these can be stored in fortran style in the
    dat file (internally we always store C style numbers) and sometimes the
    string "none" is allowed for -1.

    \see ConditionDefinition

    \author u.kue
    \date 01/08
   */
  class IntVectorConditionComponent : public ConditionComponent
  {
  public:
    IntVectorConditionComponent(std::string name,
                                int length,
                                bool fortranstyle=false,
                                bool noneallowed=false,
                                bool optional=false);
    virtual void DefaultLine(std::ostream& stream);
    virtual void Print(std::ostream& stream, const DRT::Condition* cond);
    virtual void Read(std::istringstream& condline,
                      Teuchos::RCP<DRT::Condition> condition);

  private:
    int length_;
    bool fortranstyle_;
    bool noneallowed_;
    bool optional_;
  };


  /// condition line part that is a single real number
  /*!
    \see ConditionDefinition

    \author u.kue
    \date 01/08
   */
  class RealConditionComponent : public ConditionComponent
  {
  public:
    RealConditionComponent(std::string name);
    virtual void DefaultLine(std::ostream& stream);
    virtual void Print(std::ostream& stream, const DRT::Condition* cond);
    virtual void Read(std::istringstream& condline,
                      Teuchos::RCP<DRT::Condition> condition);
  };


  /// condition line part that is a vector of real numbers
  /*!
    \see ConditionDefinition

    \author u.kue
    \date 01/08
   */
  class RealVectorConditionComponent : public ConditionComponent
  {
  public:
    RealVectorConditionComponent(std::string name, int length);
    virtual void DefaultLine(std::ostream& stream);
    virtual void Print(std::ostream& stream, const DRT::Condition* cond);
    virtual void Read(std::istringstream& condline,
                      Teuchos::RCP<DRT::Condition> condition);

  private:
    int length_;
  };


  /// definition of a valid condition in baci input
  /*!

    A ConditionDefinition is the definition of a condition dat file
    section. This definition includes the knowledge what this section looks
    like, how to read it and how to write it. In particular given a
    ConditionDefinition object it is possible to (a) write an empty dat file
    section that describes this condition, (b) read a dat file and create
    DRT::Condition objects for each line in this section and (c) write the dat
    file section filled with all corresponding conditions from a given
    DRT::Discretization.

    So this is quite sophisticated internal stuff here. If you want to
    introduce a new condition to baci, all you have to do is add an
    appropriate definition in ValidConditions(). This will take care of the
    reading part and you will get your DRT::Discretization filled with proper
    DRT::Condition objects.

    \author u.kue
    \date 01/08
   */
  class ConditionDefinition
  {
  public:

    /// construction of a condition definition
    /*!
      \param sectionname name of dat file section
      \param conditionname name of conditions in DRT::Discretization
      \param description description of condition type
      \param condtype type of conditions to be build
      \param buildgeometry whether we need conditions elements
      \param gtype type of geometry the condition lives on
     */
    ConditionDefinition(std::string sectionname,
                        std::string conditionname,
                        std::string description,
                        Condition::ConditionType condtype,
                        bool buildgeometry,
                        Condition::GeometryType gtype);

    /// add a concrete component to the condition line definition
    /*!
      Add new components to the input line. One at a time. Form left to
      right. The order is important! On reading we try and read component
      after component.
     */
    void AddComponent(Teuchos::RCP<ConditionComponent> c);

    /// read all conditions from my input file section
    /*!
      \param problem (i) global problem instance that manages the input
      \param reader (i) the actual dat file reader that has access to the dat file
      \param cmap (o) the conditions we read here
     */
    void Read(const Problem& problem,
              const DatFileReader& reader,
              std::multimap<int,Teuchos::RCP<DRT::Condition> >& cmap);

    /// print my dat file section and possible conditions from the Discretization
    std::ostream& Print(std::ostream& stream, const Discretization* dis=NULL, bool color=false);

    /// my condition name
    std::string Name() const { return conditionname_; }

    /// my condition description
    std::string Description() const { return description_; }

  private:

    std::string sectionname_;
    std::string conditionname_;
    std::string description_;
    Condition::ConditionType condtype_;
    bool buildgeometry_;
    Condition::GeometryType gtype_;

    ///
    std::vector<Teuchos::RCP<ConditionComponent> > inputline_;
  };


  /// add a separator followed by a single real value
  /*!
    The name on the input line becames the name used to put the value into DRT::Condition.
   */
  inline void AddNamedReal(Teuchos::RCP<ConditionDefinition> conddef, std::string name)
  {
    conddef->AddComponent(Teuchos::rcp(new SeparatorConditionComponent(name)));
    conddef->AddComponent(Teuchos::rcp(new RealConditionComponent     (name)));
  }
}
}

#endif
#endif
