/*----------------------------------------------------------------------*/
/*!

\brief Definition of one line of an input file.

\level 0

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef DRT_LINEDEFINITION_H
#define DRT_LINEDEFINITION_H

#include <iostream>
#include <string>
#include <vector>

#include "drt_inputreader.H"

namespace DRT
{
  namespace INPUT
  {
    class LineDefinition;
    class Lines;

    /// base class of line components
    /*!
      There are internal subclasses to this class that do the actual work.

      \author u.kue
      \date 08/09
     */
    class LineComponent
    {
     public:
      virtual ~LineComponent() {}

      /// creat an exact copy of this component
      virtual LineComponent* Clone() = 0;

      /// print to a dat file comment
      virtual void Print(std::ostream& stream) = 0;

      /// try to read component from input line
      virtual bool Read(LineDefinition& definition, std::istream& stream) = 0;

      /// try to read component from input line
      virtual bool Read(LineDefinition& definition, std::string& name, std::istream& stream) = 0;

      /// try to read component from input line
      virtual bool Read(std::istream& stream) = 0;

      /// tell if the component has the specified name tag
      virtual bool IsNamed(std::string name) = 0;
    };


    /// definition of one (possible) line in a dat file
    /*!
      LineDefinition defines how one specific line in a dat file looks like. The
      idea is, that each line consists of a list of components.

      Reading a LineDefinition means filling the components with those values
      found at the input line.

      \note This is a low level class. Therefore, we do not use Teuchos::RCP much but
      instead deallocate the internal objects by hand. For this reason there are
      also proper (not exception save) copy operators provided.

      \author u.kue
      \date 08/09
     */
    class LineDefinition
    {
     public:
      LineDefinition();
      ~LineDefinition();

      LineDefinition(const LineDefinition& other);

      LineDefinition& operator=(const LineDefinition& other);

      /// create an exact copy of this LineDefinition
      Teuchos::RCP<DRT::INPUT::LineDefinition> Clone();

      /// @name Add new line components

      /// add a single string definition
      LineDefinition& AddTag(std::string name);

      /// add a single string variable
      LineDefinition& AddString(std::string name);

      /// add a single integer variable
      LineDefinition& AddInt(std::string name);

      /// add a vector of integer variables
      LineDefinition& AddIntVector(std::string name, int length);

      /// add a vector of double variables
      LineDefinition& AddDoubleVector(std::string name, int length);

      /// add a name followed by a variable string
      LineDefinition& AddNamedString(std::string name);

      /// add a name followed by an integer variable
      LineDefinition& AddNamedInt(std::string name);

      /// add a name followed by a vector of integer variables
      LineDefinition& AddNamedIntVector(std::string name, int length);

      /// add a name followed by a double variable
      LineDefinition& AddNamedDouble(std::string name);

      /// add a name followed by a vector of double variables
      LineDefinition& AddNamedDoubleVector(std::string name, int length);

      /// add a name followed by a vector of double variables
      /*!
        The parameter \p lengthdef specifies the name of an integer component
        that gives the length of the vector. The integer component has to
        precede the vector definition on the input line.
       */
      LineDefinition& AddNamedDoubleVector(std::string name, std::string lengthdef);

      /// add a single string definition
      LineDefinition& AddOptionalTag(std::string name);

      /// add a name followed by a variable string
      LineDefinition& AddOptionalNamedString(std::string name);

      /// add a name followed by an integer variable
      LineDefinition& AddOptionalNamedInt(std::string name);

      /// add a name followed by a vector of integer variables
      LineDefinition& AddOptionalNamedIntVector(std::string name, int length);

      /// add a name followed by a double variable
      LineDefinition& AddOptionalNamedDouble(std::string name);

      /// add a name followed by a vector of double variables
      LineDefinition& AddOptionalNamedDoubleVector(std::string name, int length);

      /// add a name followed by a vector of double variables
      /*!
            The parameter \p lengthdef specifies the name of an integer component
            that gives the length of the vector. The integer component has to
            precede the vector definition on the input line.
       */
      LineDefinition& AddOptionalNamedDoubleVector(std::string name, std::string lengthdef);

      /// add a name followed by a vector of string variables
      LineDefinition& AddOptionalNamedStringVector(std::string name, int length);

      /// add a name followed by a vector of string variables
      /*!
            The parameter \p lengthdef specifies the name of an integer component
            that gives the length of the vector. The integer component has to
            precede the vector definition on the input line.
            The space defines the separation between a string and the next one.
       */
      LineDefinition& AddOptionalNamedStringVector(std::string name, std::string lengthdef);

      /// add a name followed by a vector of double variables
      /*!
        The parameter \p lengthdef specifies the name of an integer component
        that gives the length of the vector. The integer component has to
        precede the vector definition on the input line.
       */
      LineDefinition& AddOptionalNamedPairOfStringAndDoubleVector(
          std::string name, std::string lengthdef);

      //@}

      /// print to dat file comment
      void Print(std::ostream& stream);

      /// try to read this line from the stream
      /*!
        If reading succeeds, the LineDefinition contains the values of the
        line. Otherwise, the values are undefined and should not be used.

        @return true on success
       */
      bool Read(std::istream& stream);

      /// try to read this line from the stream, but skip the component named name
      /*!
        If reading succeeds, the LineDefinition contains the values of the
        line. Otherwise, the values are undefined and should not be used.

        @return true on success
       */
      bool Read(std::istream& stream, const std::string* name);

      /// tell if there is a named component with the given name
      bool HaveNamed(std::string name);

      /// @name Extract values from read LineDefinition
      /// There has to be a named component of the given type

      void ExtractString(std::string name, std::string& value);
      bool FindString(std::string name);
      void ExtractInt(std::string name, int& value);
      void ExtractIntVector(std::string name, std::vector<int>& v);
      void ExtractDouble(std::string name, double& value);
      void ExtractDoubleVector(std::string name, std::vector<double>& v);
      void ExtractStringVector(std::string name, std::vector<std::string>& v);
      void ExtractPairOfStringAndDoubleVector(
          std::string name, std::vector<std::pair<std::string, double>>& v);

      //@}

     private:
      /// find named component with given name
      LineComponent* FindNamed(std::string name);

      /// remove all components
      void Clear();

      /// list of components
      std::vector<LineComponent*> components_;

      /// optional tail, only named components here
      std::map<std::string, LineComponent*> optionaltail_;

      /// components from tail read in last line
      std::set<std::string> readtailcomponents_;
    };


    /// a collection of LineDefinitions
    /*!
      The Lines class is used to describe and read input lines from the dat
      file. The point is, the description of the line is used for the actual
      reading. This way the description is never out of sync with the reading
      routines. Furthermore, the reading code and the line definition code are
      separated. All the user needs to do is add a new definition.

      Lines describes a dat file section that consists of any number of specific
      lines. At the time being Lines is used for curves and functions
      sections. Here, each line can be of one of a list of given formats.

      Reading means finding the LineDefinition that matches the input line
      exactly and filling its internal values. After reading the user has to
      extract the values from the LineDefinition.

      \author u.kue
      \date 08/09
     */
    class Lines
    {
     public:
      /// construct to read a given section
      explicit Lines(std::string sectionname) : sectionname_(sectionname) {}

      /// add a new LineDefinition
      void Add(const LineDefinition& ld) { definitions_.push_back(ld); }

      /// print all LineDefinitions as dat file comment
      void Print(std::ostream& stream);

      /// read section and return one LineDefinition per line read
      std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> Read(
          DatFileReader& reader, int suffix = -1);

     private:
      /// read given line and return LineDefinition
      Teuchos::RCP<DRT::INPUT::LineDefinition> Read(const char* line);

      /// name of the section we read
      std::string sectionname_;

      /// all possible LineDefinitions
      std::vector<LineDefinition> definitions_;
    };

  }  // namespace INPUT
}  // namespace DRT

#endif
