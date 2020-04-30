/*----------------------------------------------------------------------*/
/*! \file
 *
  \brief Interface class for contact constitutive law parameters, i.e. paramters for laws that
relate the contact gap to the contact pressure based on micro interactions

\level 3

\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/
/* macros */
#ifndef COCONSTLAW_PARAMETER_H
#define COCONSTLAW_PARAMETER_H


/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include "Epetra_Vector.h"
#include "../drt_lib/drt_container.H"

/*----------------------------------------------------------------------*/
/* forward declarations */
namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    class ConstitutiveLaw;
    class Container;
  }  // namespace CONSTITUTIVELAW
}  // namespace CONTACT

namespace INPAR
{
  namespace CONTACT
  {
    /// Type of contact constitutive law
    enum class ConstitutiveLawType
    {
      colaw_none,            ///< undefined
      colaw_brokenrational,  ///< brokenrational constitutive law
      colaw_linear,          ///< linear constitutive law
      colaw_cubic,           ///< cubic constitutive law
      colaw_power            ///< simple power law as constitutive law
    };
  }  // namespace CONTACT
}  // namespace INPAR

/*----------------------------------------------------------------------*/
/* declarations */

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    /**
     * \brief Base object to hold 'quick' access contact constitutive law parameters
     */
    class Parameter
    {
     public:
      Parameter() = delete;
      /** construct the contact constitutive law object given the parameters
       * params[in] contactconstitutivelawdata "dumb" container containing the
       * contactconstitutivelaw data from the input file
       */
      Parameter(
          const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> contactconstitutivelawdata);

      /// destructor
      virtual ~Parameter() = default;

      /// create CoConstLaw instance of matching type with my parameters
      virtual Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> CreateConstitutiveLaw() = 0;

      // Access offset of the function
      double GetOffset() { return offset_; }

      /**
       * \brief Offset from the edge (gap==0) from where the constitutive law will be used
       *
       * When regarding different smoothness patches, the maximum peaks of the patches are in
       * general not aligned. To model this phenomenon, an offset is introduced into the
       * constitutive laws
       */
      const double offset_;
    };  // class Parameter

    /// For every DRT::ParObject a corresponding DRT::ParObjectType has to be implemented
    class ConstitutiveLawType : public DRT::ParObjectType
    {
     public:
      virtual std::string Name() const { return "ParConstitutiveLawType"; }

      static ConstitutiveLawType& Instance() { return instance_; };

     private:
      static ConstitutiveLawType instance_;
    };

    /**
     * \brief Container to pass Contact Constitutive Law parameters around
     */

    class Container : public DRT::Container
    {
     public:
      /// @name life span
      //@{

      /// standard constructor
      Container(const int id,                              ///< unique contact constitutivelaw ID
          const INPAR::CONTACT::ConstitutiveLawType type,  ///< type of contact constitutivelaw
          const std::string name                           ///< name of contact constitutivelaw
      );

      /// Destructor
      virtual ~Container() = default;

      /**
       * \brief Return unique ParObject id
       *
       * every class implementing ParObject needs a unique id defined at the
       * top of this file.
       *
       */
      int UniqueParObjectId() const override
      {
        return ConstitutiveLawType::Instance().UniqueParObjectId();
      }

      /**
       * \brief Pack this class so it can be communicated
       *
       * \ref Pack and \ref Unpack are used to communicate this class
       *
       */
      void Pack(DRT::PackBuffer& data) const override;

      /** \brief Unpack data from a char vector into this class
       *
       * \ref Pack and \ref Unpack are used to communicate this class
       */
      void Unpack(const std::vector<char>& data) override;

      //@}

      /// @name Query methods
      //@{

      /// Return material id
      inline virtual int Id() const
      {
        return id_;
      }  // todo does not override anything.. is it supposed to be this way?

      /// Return material name
      inline virtual std::string Name() const
      {
        return name_;
      }  // todo does not override anything.. is it supposed to be this way?

      /// Print this ConstitutiveLaw
      void Print(std::ostream& os) const override;

      /// Return type of constitutivelaw
      inline virtual INPAR::CONTACT::ConstitutiveLawType Type() const { return type_; }

      /**
       * \brief Return quickly accessible material parameter data
       *
       * These quick access parameters are stored in separate member #params_;
       * whereas the originally read ones are stored in DRT::Container base
       */
      inline CONTACT::CONSTITUTIVELAW::Parameter* Parameter() const { return params_.get(); }

      //@}

     protected:
      /// Unique ID of this ConstitutiveLaw, no second ConstitutiveLaw of same ID may exist
      int id_;

      /// Type of this condition
      INPAR::CONTACT::ConstitutiveLawType type_;

      /// Name
      std::string name_;

      /// Unwrapped constitutivelaw data for 'quick' access
      Teuchos::RCP<CONTACT::CONSTITUTIVELAW::Parameter> params_;
    };
  }  // namespace CONSTITUTIVELAW
}  // namespace CONTACT
#endif  // #ifndef COCONSTLAW_PARAMETER_H
