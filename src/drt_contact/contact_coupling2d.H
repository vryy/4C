/*---------------------------------------------------------------------*/
/*!

\file contact_coupling2d.H

\brief Classes for mortar contact coupling in 2D.

\level 2

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_COUPLING2D_H
#define CONTACT_COUPLING2D_H

#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_wear.H"
#include "../drt_mortar/mortar_coupling2d.H"

/*!
\brief CONTACT: namespace of the BACI contact and meshtying algorithms

*/
namespace CONTACT
{
  /*!
  \brief A class representing the framework for mortar coupling of ONE
         slave element and ONE master element of a contact interface in
         2D. This is a derived class from MORTAR::Coupling2d which does
         the contact-specific stuff for 2d mortar coupling.

  \author popp (popp@lnm.mw.tum.de)
  */

  class CoCoupling2d : public MORTAR::Coupling2d
  {
   public:
    /*!
    \brief Constructor with shape function specification

    Constructs an instance of this class and enables custom shape function types.<br>
    Note that this is \b not a collective call as coupling is
    performed in parallel by individual processes.

    */
    CoCoupling2d(DRT::Discretization& idiscret, int dim, bool quad, Teuchos::ParameterList& params,
        MORTAR::MortarElement& sele, MORTAR::MortarElement& mele);

    /*!
    \brief Destructor

    */
    virtual ~CoCoupling2d() {}

    //! @name Evlauation methods

    /*!
    \brief Integrate overlap of slave / master pair (2D)

    Derived version! Most importantly, in this derived version
    a CONTACT::CoIntegrator instance is created, which also
    does integration of the mortar quantity linearizations

    This method integrates the overlap of the current MortarElement
    pair sele_ / mele_ based on the integration limits (xiproj). The
    integration includes the Mortar matrices D/M and the gap g.

    */
    bool IntegrateOverlap(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    /*!
    \brief Return type of wear surface definition

    */
    INPAR::WEAR::WearSide WearSide()
    {
      return DRT::INPUT::IntegralValue<INPAR::WEAR::WearSide>(imortar_, "BOTH_SIDED_WEAR");
    }

    /*!
    \brief Return type of wear surface definition

    */
    INPAR::WEAR::WearType WearType()
    {
      return DRT::INPUT::IntegralValue<INPAR::WEAR::WearType>(imortar_, "WEARTYPE");
    }

    //@}


   protected:
    // don't want = operator and cctor
    CoCoupling2d operator=(const CoCoupling2d& old);
    CoCoupling2d(const CoCoupling2d& old);

    // new variables as compared to base class
    INPAR::CONTACT::SolvingStrategy stype_;

  };  // class CoCoupling2d

  /*!
  \brief A class representing the framework for mortar coupling of ONE
         slave element and SEVERAL master elements of a mortar interface in
         2D. Concretely, this class simply stores several Coupling2d objects.

  \author popp (popp@lnm.mw.tum.de)
  */

  class CoCoupling2dManager : public MORTAR::Coupling2dManager
  {
   public:
    /*!
    \brief Constructor with shape function specification

    Constructs an instance of this class and enables custom shape function types.<br>
    Note that this is \b not a collective call as coupling is
    performed in parallel by individual processes.

    */
    CoCoupling2dManager(DRT::Discretization& idiscret, int dim, bool quad,
        Teuchos::ParameterList& params, MORTAR::MortarElement* sele,
        std::vector<MORTAR::MortarElement*> mele);

    /*!
    \brief Destructor

    */
    virtual ~CoCoupling2dManager() {}

    /*!
    \brief Get communicator

    */
    virtual const Epetra_Comm& Comm() const;

    /*!
    \brief Get problem dimension

    */
    virtual const int& Dim() const { return dim_; }

    /*!
    \brief Return the LM shape fcn type

    */
    INPAR::MORTAR::ShapeFcn ShapeFcn()
    {
      return DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(imortar_, "LM_SHAPEFCN");
    }

    /*!
    \brief Evaluate mortar coupling

    */
    virtual void IntegrateCoupling(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    virtual bool EvaluateCoupling(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);
    //@}
   private:
    /*!
    \brief Calculate consistent dual shape functions in boundary elements

    */
    virtual void ConsistDualShape();

   protected:
    // don't want = operator and cctor
    CoCoupling2dManager operator=(const CoCoupling2dManager& old);
    CoCoupling2dManager(const CoCoupling2dManager& old);

    INPAR::CONTACT::SolvingStrategy stype_;  // solving strategy

  };  // class CoCoupling2dManager

}  // namespace CONTACT


#endif  // #ifndef CONTACT_COUPLING2D_H
