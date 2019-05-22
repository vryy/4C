/*---------------------------------------------------------------------*/
/*!
\brief Factory to create the desired integrator object.

\level 2

\maintainer Matthias Mayr

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_CONTACT_INTEGRATOR_FACTORY_H_
#define SRC_DRT_CONTACT_CONTACT_INTEGRATOR_FACTORY_H_

#include "../drt_contact_xcontact/xcontact_integrator.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_element.H"

// forward declaration
class Epetra_Comm;
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace CONTACT
{
  class CoIntegrator;
  namespace INTEGRATOR
  {
    class Factory
    {
     public:
      //! constructor
      Factory();

      //! destructor
      virtual ~Factory(){};

      /*! \brief Build the desired contact integrator
       *
       * \date 04/16
       * \author hiermeier */
      Teuchos::RCP<CONTACT::CoIntegrator> BuildIntegrator(
          const INPAR::CONTACT::SolvingStrategy& stype, Teuchos::ParameterList& p_mortar,
          const DRT::Element::DiscretizationType& slave_type, const Epetra_Comm& comm) const;

     private:
      /*! \brief Build the eXtended integrator considering the slave element type
       *
       *  \author hiermeier
       *  \date 08/16 */
      CONTACT::CoIntegrator* BuildXIntegrator(Teuchos::ParameterList& p_mortar,
          const DRT::Element::DiscretizationType& slave_type, const Epetra_Comm& comm) const;

      /*! \brief Build the eXtended integrator considering the problem dimension
       *
       *  \author hiermeier
       *  \date 08/16 */
      template <DRT::Element::DiscretizationType slave_type>
      CONTACT::CoIntegrator* BuildConcreteXIntegrator(
          Teuchos::ParameterList& p_mortar, const Epetra_Comm& comm) const
      {
        const int probdim = p_mortar.get<int>("DIMENSION");
        const int dim = DRT::UTILS::DisTypeToDim<slave_type>::dim;
        CONTACT::CoIntegrator* xintegrator = NULL;

        if (probdim <= dim)
        {
          dserror(
              "DIMENSION MISMATCH -- Problem dimension must be"
              "larger than the element dimension! (probdim <= dim | "
              "%d <= %d)",
              probdim, dim);
        }

        switch (probdim)
        {
          case 2:
          {
            xintegrator = new XCONTACT::Integrator<2, slave_type>(p_mortar, comm);
            break;
          }
          default:
          {
            dserror("Currently unsupported problem dimension! (probdim = %d)", probdim);
            break;
          }
        }

        return xintegrator;
      }
    };  // class Factory

    // non-member function, please call this one from outside!
    Teuchos::RCP<CONTACT::CoIntegrator> BuildIntegrator(
        const INPAR::CONTACT::SolvingStrategy& stype, Teuchos::ParameterList& p_mortar,
        const DRT::Element::DiscretizationType& eletype, const Epetra_Comm& comm);
  }  // namespace INTEGRATOR
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_CONTACT_INTEGRATOR_FACTORY_H_ */
