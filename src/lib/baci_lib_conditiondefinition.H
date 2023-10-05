/*---------------------------------------------------------------------*/
/*! \file

\brief Base class for all conditions

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LIB_CONDITIONDEFINITION_H
#define BACI_LIB_CONDITIONDEFINITION_H

#include "baci_lib_condition.H"
#include "baci_lib_inputreader.H"

#include <Teuchos_Array.hpp>
#include <Teuchos_RCP.hpp>

#include <iostream>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>


namespace DRT
{
  class Problem;
  class Discretization;

  namespace INPUT
  {
    class ConditionDefinition;

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
      /* construct with the name of the corresponding variable in the condition
       */
      explicit ConditionComponent(std::string name);

      /// virtual destructor is mandatory
      virtual ~ConditionComponent() = default;

      /// write my part of the default (comment) line of the condition
      virtual void DefaultLine(std::ostream& stream) = 0;

      /* write my part of a default line of the condition
       * as restructuredText for ReadTheDocs
       * For some components it returns the same output as Default line (but as a string).
       * However, for many components the output in ReadTheDocs is more illustrative.
       */
      virtual std::string WriteReadTheDocs() = 0;

      /// write my part of this actual condition
      virtual void Print(std::ostream& stream, const DRT::Condition* cond) = 0;

      /// read my part of an input line and put it into the condition
      virtual Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition) = 0;

      /*!
       * This routine is only used for DRT::INPUT::StringConditionComponent
       *
       * @return list of options as a a Teuchos array of strings
       *
       * \see DRT::INPUT::StringConditionComponent
       */
      //
      virtual Teuchos::Array<std::string> GetOptions() = 0;

      /// the name of my variable inside a condition
      std::string Name() const { return name_; }

     private:
      /// my condition variable name
      std::string name_;

     protected:
      /// push back token into stream so that it can be read again
      Teuchos::RCP<std::stringstream> PushBack(
          const std::string& token, const Teuchos::RCP<std::stringstream>& stream);
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
      /// construct condition component that reads and writes a string variable
      /*!
        A string is read from the input file and (another one is) stored in
        the DRT::Condition.

        It might be a good idea to store the same string in the Condition that
        is read from the dat file. You don't have to. You can make things as
        confusing as you want.

        An optional StringConditionComponent will not complain if does not find
        one of the admissible strings but set the default value instead.

        There is a 1:1 correspondance between datfilevalues and condvalues.

        \param name Name of the corresponding variable in the DRT::Condition
        \param defaultvalue Default value of this variable
        \param datfilevalues Possible values in the dat file
        \param condvalues Corresponding values stored in the DRT::Condition (std::string or enum =^
        int) \param optional Whether this string is optional in the dat file
      */
      StringConditionComponent(std::string name, std::string defaultvalue,
          const Teuchos::Array<std::string>& datfilevalues,
          const Teuchos::Array<std::string>& stringcondvalues, bool optional = false);

      StringConditionComponent(std::string name, std::string defaultvalue,
          const Teuchos::Array<std::string>& datfilevalues,
          const Teuchos::Array<int>& intcondvalues, bool optional = false);

      void DefaultLine(std::ostream& stream) override;

      std::string WriteReadTheDocs() override;

      void Print(std::ostream& stream, const DRT::Condition* cond) override;

      Teuchos::Array<std::string> GetOptions() override;

      Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline,
          Teuchos::RCP<DRT::Condition> condition) override;

     private:
      std::string defaultvalue_;
      Teuchos::Array<std::string> datfilevalues_;
      Teuchos::Array<std::string> stringcondvalues_;
      Teuchos::Array<int> intcondvalues_;
      const bool stringtostring_;
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
      /// construct a separator component that does read a string without setting a variable
      SeparatorConditionComponent(std::string separator, bool optional = false);

      void DefaultLine(std::ostream& stream) override;

      std::string WriteReadTheDocs() override;

      Teuchos::Array<std::string> GetOptions() override;

      void Print(std::ostream& stream, const DRT::Condition* cond) override;
      Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline,
          Teuchos::RCP<DRT::Condition> condition) override;

     private:
      std::string separator_;
      bool optional_;
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
      /// construct a condition component that reads a int variable
      /*!
        A int value is expected in the dat file.

        \param name variable name in DRT::Condition
        \param fortranstyle whether to substract 1 from the value read
        \param noneallowed whether "none" is a legal value (meaning -1)
        \param[in] optional true if the parameter is optional
       */
      explicit IntConditionComponent(std::string name, bool fortranstyle = false,
          bool noneallowed = false, bool optional = false);

      void DefaultLine(std::ostream& stream) override;

      std::string WriteReadTheDocs() override;

      Teuchos::Array<std::string> GetOptions() override;

      void Print(std::ostream& stream, const DRT::Condition* cond) override;
      Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline,
          Teuchos::RCP<DRT::Condition> condition) override;

     private:
      bool fortranstyle_;
      bool noneallowed_;
      bool optional_;
    };


    /**
     * A function type used to determine the length of other components from the
     * @p already_parsed_condition.
     */
    using LengthDefinition = std::function<int(const Condition& already_parsed_condition)>;

    struct LengthFromInt
    {
      /**
       * Determine the length of the vector component at runtime from an IntConditionComponent of
       * given @p name.
       */
      LengthFromInt(std::string name) : name_(std::move(name)) {}
      int operator()(const Condition& already_read_line) { return already_read_line.GetInt(name_); }

     private:
      std::string name_;
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
      /// construct a condition component that reads a vector of int variables
      /*!
        A number of int values are expected in the dat file.

        \param name variable name in DRT::Condition
        \param length number of int values expected
        \param fortranstyle whether to substract 1 from the value read
        \param noneallowed whether "none" is a legal value (meaning -1)
        \param optional whether this condition component can be skipped
       */
      IntVectorConditionComponent(std::string name, int length, bool fortranstyle = false,
          bool noneallowed = false, bool optional = false);

      /**
       * Construct an IntVectorConditionComponent similar to the other constructor but use the
       * given @p length_from_component function to determine the number of vector elements at
       * runtime.
       */
      IntVectorConditionComponent(std::string name, LengthDefinition length_from_component,
          bool fortranstyle = false, bool noneallowed = false, bool optional = false);

      void DefaultLine(std::ostream& stream) override;

      std::string WriteReadTheDocs() override;

      Teuchos::Array<std::string> GetOptions() override;

      void Print(std::ostream& stream, const DRT::Condition* cond) override;
      Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline,
          Teuchos::RCP<DRT::Condition> condition) override;
      virtual void SetLength(int newlength);

     private:
      std::variant<int, LengthDefinition> length_;
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
      /// construct a condition component that reads a double variable
      /*!
        A double value is expected in the dat file.

        \param name variable name in DRT::Condition
       */
      RealConditionComponent(std::string name);

      void DefaultLine(std::ostream& stream) override;

      std::string WriteReadTheDocs() override;

      Teuchos::Array<std::string> GetOptions() override;

      void Print(std::ostream& stream, const DRT::Condition* cond) override;
      Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline,
          Teuchos::RCP<DRT::Condition> condition) override;
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
      /// construct a condition component that reads a vector of double variables
      /*!
        A number of double values are expected in the dat file.

        \param name variable name in DRT::Condition
        \param length number of double values expected
       */
      RealVectorConditionComponent(std::string name, int length, bool optional = false);

      /**
       * Construct an RealVectorConditionComponent similar to the other constructor but use the
       * given @p length_from_component function to determine the number of vector elements at
       * runtime.
       */
      RealVectorConditionComponent(
          std::string name, LengthDefinition length_from_component, bool optional = false);

      void DefaultLine(std::ostream& stream) override;

      std::string WriteReadTheDocs() override;

      Teuchos::Array<std::string> GetOptions() override;

      void Print(std::ostream& stream, const DRT::Condition* cond) override;

      Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline,
          Teuchos::RCP<DRT::Condition> condition) override;
      virtual void SetLength(int newlength);

     private:
      std::variant<int, LengthDefinition> length_;
      bool optional_;
    };

    /**
     * This component contains a series of ConditionComponents that are selected by a key parameter.
     */
    class SwitchConditionComponent : public ConditionComponent
    {
      //! This component only supports integers for keys. Unscoped enums convert automatically to
      //! int and can be used to increase readability.
      using KeyType = int;

     public:
      /**
       * Define a component that selects one of the @p choices at runtime. This component lets the
       * user create composite structures of nested components. Depending on the integer, the
       * corresponding vector of components from the @p choices map is selected and reading
       * continues with these components. By default, the selection is made based on @p default_key.
       */
      SwitchConditionComponent(std::string name, const KeyType& default_key,
          std::map<KeyType, std::pair<std::string, std::vector<Teuchos::RCP<ConditionComponent>>>>
              choices);

      void DefaultLine(std::ostream& stream) override;

      std::string WriteReadTheDocs() override;

      std::vector<std::string> WriteReadTheDocsLines();

      Teuchos::Array<std::string> GetOptions() override;

      void Print(std::ostream& stream, const DRT::Condition* cond) override;

      Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline,
          Teuchos::RCP<DRT::Condition> condition) override;

     private:
      KeyType default_key_;
      std::map<KeyType, std::pair<std::string, std::vector<Teuchos::RCP<ConditionComponent>>>>
          choices_;

      //! Helper component to read the selected key from input.
      std::unique_ptr<StringConditionComponent> component_for_key_;
    };

    //--------------------------------------------------------------

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
      ConditionDefinition(std::string sectionname, std::string conditionname,
          std::string description, Condition::ConditionType condtype, bool buildgeometry,
          Condition::GeometryType gtype);

      /// add a concrete component to the condition line definition
      /*!
        Add new components to the input line. One at a time. Form left to
        right. The order is important! On reading we try and read component
        after component.
       */
      void AddComponent(const Teuchos::RCP<ConditionComponent>& c);

      /// read all conditions from my input file section
      /*!
        \param problem (i) global problem instance that manages the input
        \param reader (i) the actual dat file reader that has access to the dat file
        \param cmap (o) the conditions we read here
       */
      void Read(const Problem& problem, DatFileReader& reader,
          std::multimap<int, Teuchos::RCP<DRT::Condition>>& cmap);

      /// print my dat file section and possible conditions from the Discretization
      std::ostream& Print(std::ostream& stream, const Discretization* dis = nullptr);

      /// name of my section in input file
      std::string SectionName() const { return sectionname_; }

      /// my condition name
      std::string Name() const { return conditionname_; }

      /// my condition description
      std::string Description() const { return description_; }

      /// my condition inputline
      std::vector<Teuchos::RCP<ConditionComponent>> Inputline() const { return inputline_; }

      /// my GeometryType
      Condition::GeometryType GeometryType() const { return gtype_; }

     private:
      std::string sectionname_;
      std::string conditionname_;
      std::string description_;
      Condition::ConditionType condtype_;
      bool buildgeometry_;
      Condition::GeometryType gtype_;

      ///
      std::vector<Teuchos::RCP<ConditionComponent>> inputline_;
    };


    /// add a separator followed by a single integer value
    /*!
      The name on the input line becomes the name used to put the int value into DRT::Condition.
     */
    inline void AddNamedInt(
        const Teuchos::RCP<ConditionDefinition>& conddef, const std::string& name)
    {
      conddef->AddComponent(Teuchos::rcp(new SeparatorConditionComponent(name)));
      conddef->AddComponent(Teuchos::rcp(new IntConditionComponent(name)));
    }

    /// add a separator followed by a single real value
    /*!
      The name on the input line becomes the name used to put the value into DRT::Condition.
     */
    inline void AddNamedReal(
        const Teuchos::RCP<ConditionDefinition>& conddef, const std::string& name)
    {
      conddef->AddComponent(Teuchos::rcp(new SeparatorConditionComponent(name)));
      conddef->AddComponent(Teuchos::rcp(new RealConditionComponent(name)));
    }

    /// add a separator followed by a single real value
    /*!
      The name on the input line becomes the name used to put the value into DRT::Condition.
     */
    inline void AddNamedVector(
        const Teuchos::RCP<ConditionDefinition>& conddef, const std::string& name, const int size)
    {
      conddef->AddComponent(Teuchos::rcp(new SeparatorConditionComponent(name)));
      conddef->AddComponent(Teuchos::rcp(new RealVectorConditionComponent(name, size)));
    }
  }  // namespace INPUT
}  // namespace DRT


#endif
