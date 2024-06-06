/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of one line of an input file.

\level 0


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_IO_LINEDEFINITION_HPP
#define FOUR_C_IO_LINEDEFINITION_HPP

#include "4C_config.hpp"

#include "4C_io_input_parameter_container.hpp"
#include "4C_io_inputreader.hpp"

#include <functional>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Input
{
  namespace INTERNAL
  {
    class LineDefinitionImplementation;
  }

  /**
   * LineDefinition defines how one specific line in a dat file looks like. The
   * idea is, that each line consists of a list of components.
   *
   * Reading a LineDefinition means filling the components with those values
   * found at the input line.
   *
   * @note This class has value semantics.
   */
  class LineDefinition
  {
   public:
    /**
     * An empty LineDefinition without any components.
     */
    LineDefinition();

    /**
     * Destructor.
     */
    ~LineDefinition();

    /**
     * Copy constructor.
     */
    LineDefinition(const LineDefinition& other);

    /**
     * Copy assignment.
     */
    LineDefinition& operator=(const LineDefinition& other);

    /**
     * Move constructor.
     */
    LineDefinition(LineDefinition&& other) noexcept;

    /**
     * Move assignment.
     */
    LineDefinition& operator=(LineDefinition&& other) noexcept;

    /**
     * Builder class to incrementally add components and then build a LineDefinition from them.
     * Usage example:
     *
     * @code
     *   const LineDefinition line =
     *     LineDefinition::Builder().AddTag("a").AddInt("i).Build();
     * @endcode
     */
    class Builder
    {
     public:
      /**
       * A function type that may be supplied to some of the Add... functions in this class to
       * define how the number of vector entries should be determined from the @p already_read_line.
       */
      using LengthDefinition =
          std::function<std::size_t(const Core::IO::InputParameterContainer& already_read_line)>;

      /**
       * Create a new Builder.
       */
      Builder();

      /**
       * Destructor.
       */
      ~Builder();
      /**
       * Copy constructor.
       */
      Builder(const Builder& other);

      /**
       * Copy assignment.
       */
      Builder& operator=(const Builder& other);

      /**
       * Move constructor.
       */
      Builder(Builder&& other) noexcept;

      /**
       * Move assignment.
       */
      Builder& operator=(Builder&& other) noexcept;

      /**
       * Initialize a Builder with all components from an existing LineDefinition.
       */
      Builder(const LineDefinition& line_definition);

      /**
       * Convert the gathered components into a LineDefinition.
       *
       * @note This overload is chosen for rvalues and can steal the data gathered in the Builder.
       * After calling this function, the builder is in a moved-from state.
       */
      [[nodiscard]] LineDefinition Build() &&;

      /**
       * Convert the gathered components into a LineDefinition.
       *
       * @note This overload copies the gathered data to the newly created LineDefinition.
       */
      [[nodiscard]] LineDefinition Build() const&;

      /// Add a single string definition without a value.
      Builder& AddOptionalTag(const std::string& name);

      /// Add a single string definition without a value.
      Builder& AddTag(std::string name);

      /// Add a single string variable
      Builder& AddString(std::string name);

      /// Add a single integer variable
      Builder& AddInt(std::string name);

      /// Add a vector of integer variables
      Builder& AddIntVector(std::string name, int length);

      /// Add a vector of double variables
      Builder& AddDoubleVector(std::string name, int length);

      /// Add a name followed by a variable string
      Builder& AddNamedString(std::string name);

      /// Add a name followed by an integer variable
      Builder& AddNamedInt(std::string name);

      /// Add a name followed by a vector of integer variables
      Builder& AddNamedIntVector(std::string name, int length);

      /// Add a name followed by a double variable
      Builder& AddNamedDouble(std::string name);

      /// Add a name followed by a vector of double variables
      Builder& add_named_double_vector(std::string name, int length);

      /*!
       * Add a name followed by a vector of double variables.
       *
       * The function @p length_definition specifies how to obtain the number of to be read vector
       * entries from the previously added components. See LengthDefinition for details.
       */
      Builder& add_named_double_vector(std::string name, LengthDefinition length_definition);

      /// Add a name followed by a variable string
      Builder& add_optional_named_string(const std::string& name);

      /// Add a name followed by an integer variable
      Builder& AddOptionalNamedInt(const std::string& name);

      /// Add a name followed by a vector of integer variables
      Builder& add_optional_named_int_vector(const std::string& name, int length);

      /// Add a name followed by a double variable
      Builder& add_optional_named_double(const std::string& name);

      /// Add a name followed by a vector of double variables
      Builder& add_optional_named_double_vector(const std::string& name, int length);

      /*!
       * Add a name followed by a vector of double variables.
       *
       * The parameter \p lengthdef specifies the name of an integer component
       * that gives the length of the vector. The integer component has to
       * precede the vector definition on the input line.
       */
      Builder& add_optional_named_double_vector(
          const std::string& name, LengthDefinition lengthdef);

      /// Add a name followed by a vector of string variables.
      Builder& add_optional_named_string_vector(const std::string& name, int length);

      /**
       * Add a name followed by a vector of string variables.
       *
       * The parameter \p lengthdef specifies the name of an integer component
       * that gives the length of the vector. The integer component has to
       * precede the vector definition on the input line.
       * The space defines the separation between a string and the next one.
       */
      Builder& add_optional_named_string_vector(
          const std::string& name, LengthDefinition lengthdef);

      /*!
       * Add a name followed by a vector of double variables.
       *
       * The parameter \p lengthdef specifies the name of an integer component
       * that gives the length of the vector. The integer component has to
       * precede the vector definition on the input line.
       */
      Builder& add_optional_named_pair_of_string_and_double_vector(
          const std::string& name, LengthDefinition lengthdef);

     private:
      /// Implementation details are hidden behind the PIMPL idiom.
      std::unique_ptr<INTERNAL::LineDefinitionImplementation> pimpl_;
    };

    /// print to dat file comment
    void Print(std::ostream& stream) const;

    /**
     * If reading succeeds, returns the data. Otherwise, returns an empty std::optional.
     */
    std::optional<Core::IO::InputParameterContainer> Read(std::istream& stream);

    /**
     * If reading succeeds, returns the data. Otherwise, returns an empty std::optional.
     */
    std::optional<Core::IO::InputParameterContainer> Read(
        std::istream& stream, const std::string* name);

    /// tell if there is a named component with the given name
    [[nodiscard]] bool HaveNamed(const std::string& name) const;

    /// @name Extract values from read LineDefinition
    /// There has to be a named component of the given type

    void ExtractString(const std::string& name, std::string& value) const;
    [[nodiscard]] bool HasString(const std::string& name) const;
    void ExtractInt(const std::string& name, int& value) const;
    void ExtractIntVector(const std::string& name, std::vector<int>& v) const;
    void ExtractDouble(const std::string& name, double& value) const;
    void ExtractDoubleVector(const std::string& name, std::vector<double>& v) const;
    void ExtractStringVector(const std::string& name, std::vector<std::string>& v) const;
    void extract_pair_of_string_and_double_vector(
        const std::string& name, std::vector<std::pair<std::string, double>>& v) const;

    //@}

   private:
    /// Constructor called by the Builder to directly pass on implementation.
    explicit LineDefinition(std::unique_ptr<INTERNAL::LineDefinitionImplementation>&& pimpl);

    /// Implementation details are hidden behind the PIMPL idiom.
    std::unique_ptr<INTERNAL::LineDefinitionImplementation> pimpl_;
  };


  /**
   * Helper functor to parse the length of a vector input from an integer component. This
   * functor is compatible with LengthDefinition.
   * Example:
   *
   * @code
   *    [...].AddDoubleVector("values", FromIntNamed("NUMVALUES")) [...]
   * @endcode
   */
  struct LengthFromIntNamed
  {
    LengthFromIntNamed(std::string definition_name);

    std::size_t operator()(const Core::IO::InputParameterContainer& already_read_line);

   private:
    std::string definition_name_;
  };

}  // namespace Input

FOUR_C_NAMESPACE_CLOSE

#endif
