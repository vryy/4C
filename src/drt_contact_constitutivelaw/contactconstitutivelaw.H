/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for contact constitutive laws, i.e. laws that relate the contact gap to the
contact pressure based on micro interactions

\level 1

\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/

#ifndef CONSTITUTIVELAW_H
#define CONSTITUTIVELAW_H



#include <Teuchos_RCP.hpp>
#include "../drt_lib/drt_parobject.H"
#include "contactconstitutivelaw_parameter.H"

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    class Parameter;

    /**
     * \brief The ConstitutiveLaw class provides a framework in order to relate the contact gap to
     * the contact pressure using information like micro roughness for contact problems.
     *
     * In order to support storage of material data (for restart, post
     * processing, ...) this class implements DRT::ParObject.
     */
    class ConstitutiveLaw : public DRT::ParObject
    {
     public:
      /// return type of this constitutive law
      virtual INPAR::CONTACT::ConstitutiveLawType GetConstitutiveLawType() const
      {
        return INPAR::CONTACT::ConstitutiveLawType::colaw_none;
      };

      /// return copy of this Contact ConstitutiveLaw object
      virtual Teuchos::RCP<ConstitutiveLaw> Clone() const = 0;

      /// Return quick accessible Contact Constitutive Law parameter data
      virtual CONTACT::CONSTITUTIVELAW::Parameter* Parameter() const = 0;

      virtual double Evaluate(double gap) = 0;
      virtual double EvaluateDeriv(double gap) = 0;

      /* \brief create Contact ConstitutiveLaw object given the id of the constitutive law in the
       * input file
       */
      static Teuchos::RCP<ConstitutiveLaw> Factory(const int id);

      /* \brief create Contact Constitutivelaw object given input information for the constitutive
       * law \param[in] Container holding the Coefficients for the Contact ConstitutiveLaw
       */
      static Teuchos::RCP<ConstitutiveLaw> Factory(
          const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> contactconstitutivelawdata);

      virtual ~ConstitutiveLaw() = default;
    };
  }  // namespace CONSTITUTIVELAW
}  // namespace CONTACT

#endif
