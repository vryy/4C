/*----------------------------------------------------------------------*/
/*! \file
 *
\brief Definition of Constitutive contact Laws

\level 3

\maintainer Nora Hagmeyer
*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_COCONSTLAW_DEFINITION_H
#define CONTACT_COCONSTLAW_DEFINITION_H

#include <Teuchos_RCP.hpp>
#include <iostream>
#include <list>
#include <string>
#include <vector>

namespace DRT
{
  class Problem;
  class Discretization;
  namespace INPUT
  {
    class DatFileReader;
  }  // namespace INPUT

}  // namespace DRT

namespace INPAR
{
  namespace CONTACT
  {
    enum class ConstitutiveLawType;
  }
}  // namespace INPAR

namespace CONTACT
{
  class CoInterface;
  namespace CONSTITUTIVELAW
  {
    class Container;
    class ContactConstitutiveLawComponent;
    class Bundle;
    /**
     * \brief Definition of a valid contact constitutive law in BACI input
     *
     * This is basically a clone of \see MaterialDefinition
     * A ConstitutivelawDefinition is the definition of a --Contact Constitutive Law in the DAT
     * file section. This definition includes the knowledge what this section looks like, how to
     * read it and how to write it. In particular given a ConstitutivelawDefinition object it is
     * possible to read a DAT file and create Constitutivelaw::Container objects for each line in
     * this section.
     */
    class LawDefinition
    {
     public:
      /**
       * \brief default constructor
       *
       * \param[in] name name of contact constitutive law
       * \param[in] description of contact constitutive law type
       * \param[in] type of contact constitutive law to be built
       */
      LawDefinition(std::string name, std::string description,
          INPAR::CONTACT::ConstitutiveLawType contactconstitutivelawlawtype);

      // destructor
      virtual ~LawDefinition() = default;

      /** add a concrete component to the contact constitutive law line definition
       *
       * \param[in] c new component needed as parameter for the contact constitutive law
       */
      void AddComponent(Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ContactConstitutiveLawComponent> c);

      /** \brief read all materials from my input file section
       *
       * \param[in] problem global problem instance that manages the input
       * \param[in] reader the actual dat file reader that has access to the dat file
       * \params[inout] mmap map mapping IDs to parameters of the contact constitutive model
       */
      void Read(const DRT::Problem& problem, DRT::INPUT::DatFileReader& reader,
          Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Bundle> mmap);

      /// print my DAT file section and possible contact constitutive laws
      std::ostream& Print(std::ostream& stream,  ///< the output stream
          const DRT::Discretization* dis = NULL, const bool color = false);

      /// get contact constitutive law name
      std::string Name() const { return coconstlawname_; }

      // get contact constitutive law type
      INPAR::CONTACT::ConstitutiveLawType Type() const { return coconstlawtype_; }

      /// get contact constitutive law description
      std::string Description() const { return description_; }

      /// get contact constitutive law inputline
      std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ContactConstitutiveLawComponent>>
      Inputline() const
      {
        return inputline_;
      }

     private:
      /// name of ontact constitutive law
      std::string coconstlawname_;
      /// description of ontact constitutive law type
      std::string description_;
      /// type of ontact constitutive law to be build
      INPAR::CONTACT::ConstitutiveLawType coconstlawtype_;
      /// the list of valid components
      std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ContactConstitutiveLawComponent>>
          inputline_;
    };

    /**
     * \brief Base class for a component in the input definition of a contact constitutive law
     * component
     */
    class ContactConstitutiveLawComponent
    {
     public:
      /**
       * \brief Constructor
       * \param[in] name  name of the parameter
       * \param[in] optional true if the parameter is optional
       */
      explicit ContactConstitutiveLawComponent(std::string name, bool optional = false);

      /// virtual destructor
      virtual ~ContactConstitutiveLawComponent() {}

      /**
       * \brief write my part of the default line of the contact constitutive law
       * definition
       */
      virtual void DefaultLine(std::ostream& stream) = 0;

      /**
       * \brief write my part of the actual line of the contact constitutive law
       * definition
       */
      virtual void Print(std::ostream& stream, const CONTACT::CONSTITUTIVELAW::Container* cond) = 0;

      /**
       * \brief describe my part of the default (comment) line of the contact constitutive law
       * definition
       */
      virtual void Describe(std::ostream& stream) = 0;

      /**
       * \brief read my part of the the contact constitutive law
       * definition
       *
       *\param[in] def Contact Constitutive Law definition
       *\param[in] condline read input line
       *\params[inout] container container to add the read parameter to
       */
      virtual Teuchos::RCP<std::stringstream> Read(LawDefinition* def,
          Teuchos::RCP<std::stringstream> condline,
          Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container) = 0;

      /// get name of parameter inside the contact constitutive law
      std::string Name() const { return name_; }

     protected:
      /// for optional components
      bool optional_;

     private:
      /// my material variable name
      std::string name_;
    };


    /**
     * \brief This is a seperator component within a contact constitutive law definition
     *
     */
    class SeparatorContactConstitutiveLawComponent : public ContactConstitutiveLawComponent
    {
     public:
      /**
       * \brief default constructor
       *
       * \param[in] name name of contact constitutive law
       * \param[in] description of contact constitutive law type
       * \param[in] type of contact constitutive law to be built
       */
      SeparatorContactConstitutiveLawComponent(
          std::string separator, std::string description, bool optional = false);

      /**
       * \brief write my part of the default line of the contact constitutive law
       * definition
       */
      void DefaultLine(std::ostream& stream) override;

      /**
       * \brief write my part of the actual line of the contact constitutive law
       * definition
       */
      void Print(std::ostream& stream, const CONTACT::CONSTITUTIVELAW::Container* cond) override;

      /**
       * \brief describe my part of the default (comment) line of the contact constitutive law
       * definition
       */
      void Describe(std::ostream& stream) override;

      /**
       * \brief read my part of the the contact constitutive law
       * definition
       *
       *\param[in] def Contact Constitutive Law definition
       *\param[in] condline read input line
       *\params[inout] container container to add the read parameter to
       */
      Teuchos::RCP<std::stringstream> Read(LawDefinition* def,
          Teuchos::RCP<std::stringstream> condline,
          Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container) override;

     private:
      /// separator string, i.e. the NAME of variable in DAT input file
      std::string separator_;

      /// description attached to the field separator
      std::string description_;
    };

    /**
     * \brief This is a floating point number component within a contact constitutive law definition
     *
     */
    class RealContactConstitutiveLawComponent : public ContactConstitutiveLawComponent
    {
     public:
      /**
       * \brief default constructor
       *
       * \param[in] name name of contact constitutive law
       * \param[in] description of contact constitutive law type
       * \param[in] type of contact constitutive law to be built
       */
      RealContactConstitutiveLawComponent(std::string name,  ///< variable name in DRT::Material
          const double defaultvalue = 0,                     ///< a default value
          bool optional = false);

      /**
       * \brief write my part of the default line of the contact constitutive law
       * definition
       */
      void DefaultLine(std::ostream& stream) override;

      /**
       * \brief write my part of the actual line of the contact constitutive law
       * definition
       */
      void Print(std::ostream& stream, const CONTACT::CONSTITUTIVELAW::Container* cond) override;

      /**
       * \brief describe my part of the default (comment) line of the contact constitutive law
       * definition
       */
      void Describe(std::ostream& stream) override;

      /**
       * \brief read my part of the the contact constitutive law
       * definition
       *
       *\param[in] def Contact Constitutive Law definition
       *\param[in] condline read input line
       *\params[inout] container container to add the read parameter to
       */
      Teuchos::RCP<std::stringstream> Read(LawDefinition* def,
          Teuchos::RCP<std::stringstream> condline,
          Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Container> container) override;

     private:
      /// a default value
      int defaultvalue_;
    };

    /** Add a separator followed by a single real value
     *
     * The name on the input line becomes the name used to put the value into the contact
     * constitutive law container container
     *
     * \params[inout] def contact constitutive law definition to which the double should be added
     * \param[in] name name of the double parameter
     * \param[in] description of the parameter
     * \param[in] defaultvalue default value to be used when parameter is optional and not given
     * \param[in] optional true if parameter is optional
     *
     */
    inline void AddNamedReal(Teuchos::RCP<LawDefinition> def, const std::string name,
        const std::string description, const double defaultvalue = 0.0, const bool optional = false)
    {
      def->AddComponent(
          Teuchos::rcp(new SeparatorContactConstitutiveLawComponent(name, description, optional)));
      def->AddComponent(
          Teuchos::rcp(new RealContactConstitutiveLawComponent(name, defaultvalue, optional)));
    }

    /** Add a separator
     *
     * The name on the input line becomes the name used to put the value into the contact
     * constitutive law container container
     *
     * \params[inout] def contact constitutive law definition to which the double should be added
     * \param[in] name name of the double parameter
     * \param[in] description of the parameter
     * \param[in] optional true if parameter is optional
     *
     */
    inline void AddNamedSeparator(Teuchos::RCP<LawDefinition> def, const std::string name,
        const std::string description, const bool optional = false)
    {
      def->AddComponent(
          Teuchos::rcp(new SeparatorContactConstitutiveLawComponent(name, description, optional)));
    }

    /** Add contact constitutive law definition to list of defined definitions
     *
     * \params[inout] list list of defined contact constitutive law definitions
     * \param[in] def contact constitutive law definition to add to the list
     *
     */
    void AppendCoConstLawComponentDefinition(
        std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>& list,
        Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition> def);

  }  // namespace CONSTITUTIVELAW
}  // namespace CONTACT

#endif  // #ifndef CONTACT_COCONSTLAW_INTERFACE_H
