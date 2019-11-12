/*---------------------------------------------------------------------*/
/*! \file

\brief Base class for all conditions

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

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
      /// construct with the name of the corresponding variable in the condition
      explicit ConditionComponent(std::string name);

      /// virtual destructor is mandatory
      virtual ~ConditionComponent() {}

      /// write my part of the default (comment) line of the condition
      virtual void DefaultLine(std::ostream& stream) = 0;

      /// write my part of this actual condition
      virtual void Print(std::ostream& stream, const DRT::Condition* cond) = 0;

      /// read my part of an input line and put it into the condition
      virtual Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition) = 0;

      /// the name of my variable inside a condition
      std::string Name() const { return name_; }

     protected:
      /// push back token into stream so that it can be read again
      Teuchos::RCP<std::stringstream> PushBack(
          std::string token, Teuchos::RCP<std::stringstream> stream);

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

      virtual void DefaultLine(std::ostream& stream);
      virtual void Print(std::ostream& stream, const DRT::Condition* cond);
      virtual Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition);

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

      virtual void DefaultLine(std::ostream& stream);
      virtual void Print(std::ostream& stream, const DRT::Condition* cond);
      virtual Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition);

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
        \params[in] optional true if the parameter is optional
       */
      explicit IntConditionComponent(std::string name, bool fortranstyle = false,
          bool noneallowed = false, bool optional = false);
      virtual void DefaultLine(std::ostream& stream);
      virtual void Print(std::ostream& stream, const DRT::Condition* cond);
      virtual Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition);

     private:
      bool fortranstyle_;
      bool noneallowed_;
      bool optional_;
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
      virtual void DefaultLine(std::ostream& stream);
      virtual void Print(std::ostream& stream, const DRT::Condition* cond);
      virtual Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition);
      virtual void SetLength(int newlength);

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
      /// construct a condition component that reads a double variable
      /*!
        A double value is expected in the dat file.

        \param name variable name in DRT::Condition
       */
      RealConditionComponent(std::string name);
      virtual void DefaultLine(std::ostream& stream);
      virtual void Print(std::ostream& stream, const DRT::Condition* cond);
      virtual Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition);
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
      virtual void DefaultLine(std::ostream& stream);
      virtual void Print(std::ostream& stream, const DRT::Condition* cond);
      virtual Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition);
      virtual void SetLength(int newlength);

     private:
      int length_;
      bool optional_;
    };

    /// Dirichlet and Neumann boundary condition specific condition component
    /*!
      \see ConditionDefinition

      \author AP and TK
      \date 03/12
     */
    class DirichletNeumannBundle : public ConditionComponent
    {
     public:
      /*! construct a condition component tailored for Dirichlet and Neumann with flexible
       vector length for onoff, val, curve and funct

        The number of dofs and corresponding vectors are expected from the dat-file

        \param name variable name in DRT::Condition
        \param length number of double values expected
       */
      DirichletNeumannBundle(std::string name, Teuchos::RCP<IntConditionComponent> intcomp,
          std::vector<Teuchos::RCP<SeparatorConditionComponent>> intvectsepcomp,
          std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp,
          std::vector<Teuchos::RCP<SeparatorConditionComponent>> realvectsepcomp,
          std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp);

      virtual void DefaultLine(std::ostream& stream);
      virtual void Print(std::ostream& stream, const DRT::Condition* cond);
      virtual Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition);

     private:
      Teuchos::RCP<IntConditionComponent> intcomp_;
      std::vector<Teuchos::RCP<SeparatorConditionComponent>> intvectsepcomp_;
      std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp_;
      std::vector<Teuchos::RCP<SeparatorConditionComponent>> realvectsepcomp_;
      std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp_;
    };

    /// int and real vector boundary condition specific condition component
    /*!
      \see ConditionDefinition

      \author ehrl
      \date 09/12
     */
    class IntRealBundle : public ConditionComponent
    {
     public:
      /*! construct a condition component tailored for a defined number of condition components
       * (int and/or real) separated by an optional separator
       * generalized DirichletNeumannBundle

        The number of condition components and corresponding vectors are expected from the dat-file
        and the drt_validconditions.cpp

        Attention: number of type real is passed to int vector -> numbers after dot are ignored
                   number of type int is passed to real vector -> number is interpreted as a real

        \param name variable name in DRT::Condition
        \param intcomp expected number of condition components per Int/Real-VectorConditionComponent
        \param intsepcomp vector defining an optional separator for the different
       IntVectorConditionComponents (optional: Teuchos::null -> no separator) \param intvectcomp
       vector containing various IntVectorConditionComponents:
               - number of IntVectorConditionComponents in intvectcomp is defined in
       drt_validconditions.cpp
               - number of condition components in a IntVectorConditionComponent is defined by
       intcomp (i.e., the integer called "length" indicating the length of
       IntVectorConditionComponent in the constructor of IntVectorConditionComponent can be
       specified arbitrarily and is eventually overwritten by the integer stored in intcomp)
               - if you do not want to have any IntVectorConditionComponent:
                 pass an empty "std:vector<Teuchos::RCP<IntVectorConditionComponent> > object" to
       the constructor \param realsepcomp vector defining an optional separator for the different
       RealVectorConditionComponents (optional: Teuchos::null -> no separator) \param realvectcomp
       vector containing various RealVectorConditionComponents:
               - number of RealVectorConditionComponents in realvectcomp is defined in
       drt_validconditions.cpp
               - number of condition components in a RealVectorConditionComponent is defined by
       intcomp (i.e., the integer called "length" indicating the length of
       RealVectorConditionComponent in the constructor of RealVectorConditionComponent can be
       specified arbitrarily and is eventually overwritten by the integer stored in intcomp)
               - if you do not want to have any RealVectorConditionComponent:
                 pass an empty "std:vector<Teuchos::RCP<RealVectorConditionComponent> > object" to
       the constructor
       */

      IntRealBundle(std::string name, Teuchos::RCP<IntConditionComponent> intcomp,
          std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp,
          std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp,
          std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepcomp,
          std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp);

      virtual void DefaultLine(std::ostream& stream);
      virtual void Print(std::ostream& stream, const DRT::Condition* cond);
      virtual Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition);

     private:
      //! handles the number of condition components per Int/Real-VectorConditionComponent
      Teuchos::RCP<IntConditionComponent> intcomp_;
      //! vector defining a separator (optional: Teuchos::null -> no separator)
      //! for the different IntVectorConditionComponents
      std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp_;
      //! vector containing IntVectorConditionComponents
      std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp_;
      //! vector defining a separator (optional: Teuchos::null -> no separator)
      //! for the different RealVectorConditionComponents
      std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepcomp_;
      //! vector containing RealVectorConditionComponents
      std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp_;
    };


    class CondCompBundle : public ConditionComponent
    {
     public:
      /*! construct a condition component handling the a condition component bundle associated to a
        special model

        \param name variable name in DRT::Condition
        \param condcomp vector containing condition components associated to a model
        \param model specification of model type (enum->int) as a security check:
               elements of the vector<CondCompBundle> need to have the same order as the enum
        defining your model
       */
      CondCompBundle(
          std::string name, std::vector<Teuchos::RCP<ConditionComponent>> condcomp, int model);

      virtual void DefaultLine(std::ostream& stream);
      virtual void Print(std::ostream& stream, const DRT::Condition* cond);
      virtual Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition);

      int Model() { return model_; };

     private:
      std::vector<Teuchos::RCP<ConditionComponent>> condcomp_;
      const int model_;
    };


    /// Boundary condition for model dependent BC parameters
    /*!
      \see ConditionDefinition

      \author Ehrl
      \date 09/12
     */
    class CondCompBundleSelector : public ConditionComponent
    {
     public:
      /*! construct a selector for a model dependent condition component bundle

        \param name variable name in DRT::Condition
        \param stringcomp type of model
                          model type have to be defined as an enum (strings are not supported)
        \param condcomp vector containing a list of BC parameter assigned to different models
       */
      CondCompBundleSelector(std::string name, Teuchos::RCP<StringConditionComponent> stringcomp,
          std::vector<Teuchos::RCP<CondCompBundle>> condcomp);

      virtual void DefaultLine(std::ostream& stream);
      virtual void Print(std::ostream& stream, const DRT::Condition* cond);
      virtual Teuchos::RCP<std::stringstream> Read(ConditionDefinition* def,
          Teuchos::RCP<std::stringstream> condline, Teuchos::RCP<DRT::Condition> condition);

     private:
      //! handling of the model type
      Teuchos::RCP<StringConditionComponent> stringcomp_;
      //! vector containing a list of BC parameter assigned to different models
      std::vector<Teuchos::RCP<CondCompBundle>> condcomp_;
    };

    //--------------------------------------------------------------------------------------------------


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
      void AddComponent(Teuchos::RCP<ConditionComponent> c);

      /// read all conditions from my input file section
      /*!
        \param problem (i) global problem instance that manages the input
        \param reader (i) the actual dat file reader that has access to the dat file
        \param cmap (o) the conditions we read here
       */
      void Read(const Problem& problem, DatFileReader& reader,
          std::multimap<int, Teuchos::RCP<DRT::Condition>>& cmap);

      /// print my dat file section and possible conditions from the Discretization
      std::ostream& Print(
          std::ostream& stream, const Discretization* dis = NULL, bool color = false);

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
    inline void AddNamedInt(Teuchos::RCP<ConditionDefinition> conddef, std::string name)
    {
      conddef->AddComponent(Teuchos::rcp(new SeparatorConditionComponent(name)));
      conddef->AddComponent(Teuchos::rcp(new IntConditionComponent(name)));
    }

    /// add a separator followed by a single real value
    /*!
      The name on the input line becomes the name used to put the value into DRT::Condition.
     */
    inline void AddNamedReal(Teuchos::RCP<ConditionDefinition> conddef, std::string name)
    {
      conddef->AddComponent(Teuchos::rcp(new SeparatorConditionComponent(name)));
      conddef->AddComponent(Teuchos::rcp(new RealConditionComponent(name)));
    }

    /// add a separator followed by a single real value
    /*!
      The name on the input line becomes the name used to put the value into DRT::Condition.
     */
    inline void AddNamedVector(
        Teuchos::RCP<ConditionDefinition> conddef, const std::string name, const int size)
    {
      conddef->AddComponent(Teuchos::rcp(new SeparatorConditionComponent(name)));
      conddef->AddComponent(Teuchos::rcp(new RealVectorConditionComponent(name, size)));
    }
  }  // namespace INPUT
}  // namespace DRT

#endif
