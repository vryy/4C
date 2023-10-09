/*----------------------------------------------------------------------*/
/*! \file
\brief Base class for components of input lines
\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_LIB_LINECOMPONENT_H
#define BACI_LIB_LINECOMPONENT_H


#include "baci_lib_container.H"

#include <Teuchos_RCP.hpp>

namespace INPUT
{

  /**
   * Interface for components in an input line.
   *
   * TODO: this class unifies duplicated interfaces and is currently a collection of various
   * methods. This is not the final interface.
   *
   */
  class LineComponent
  {
   public:
    /// construct with the name of the corresponding variable in the material
    explicit LineComponent(std::string name, bool optional = false);

    /// virtual destructor is mandatory
    virtual ~LineComponent() = default;

    /// write my part of the default (comment) line of the condition
    virtual void DefaultLine(std::ostream& stream) = 0;

    /// Write whatever this LineComponent owns in the given @p container.
    virtual void Print(std::ostream& stream, const DRT::Container& container) = 0;

    /// A human-readable description of this component used in help messages.
    virtual void Describe(std::ostream& stream) {}

    /// TODO: ConditionDefinition* def is only used to get the SectionName() could be a string or
    /// droppped. Same is true for MaterialDefinition* def!
    /// read my part of an input line and put it into the condition
    virtual Teuchos::RCP<std::stringstream> Read(const std::string& section_name,
        Teuchos::RCP<std::stringstream> condline, DRT::Container& container) = 0;

    /* write my part of a default line of the condition
     * as restructuredText for ReadTheDocs
     * For some components it returns the same output as Default line (but as a string).
     * However, for many components the output in ReadTheDocs is more illustrative.
     */
    virtual std::string WriteReadTheDocs() { return ""; }


    virtual Teuchos::Array<std::string> GetOptions() { return {}; }

    /// the name of my variable inside a material
    [[nodiscard]] std::string Name() const { return name_; }

   protected:
    /// push back token into stream so that it can be read again
    Teuchos::RCP<std::stringstream> PushBack(
        const std::string& token, const Teuchos::RCP<std::stringstream>& stream);

    /// for optional components
    bool optional_;

   private:
    /// my material variable name
    std::string name_;
  };

}  // namespace INPUT

#endif
