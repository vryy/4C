/*---------------------------------------------------------------------*/
/*! \file
\brief Classes for mortar contact coupling in 2D.

\level 2


*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_COUPLING2D_HPP
#define FOUR_C_CONTACT_COUPLING2D_HPP

#include "4C_config.hpp"

#include "4C_inpar_contact.hpp"
#include "4C_inpar_wear.hpp"
#include "4C_mortar_coupling2d.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  /*!
  \brief A class representing the framework for mortar coupling of ONE
         slave element and ONE master element of a contact interface in
         2D. This is a derived class from Mortar::Coupling2d which does
         the contact-specific stuff for 2d mortar coupling.

  */

  class Coupling2d : public Mortar::Coupling2d
  {
   public:
    /*!
    \brief Constructor with shape function specification

    Constructs an instance of this class and enables custom shape function types.<br>
    Note that this is \b not a collective call as coupling is
    performed in parallel by individual processes.

    */
    Coupling2d(Core::FE::Discretization& idiscret, int dim, bool quad,
        Teuchos::ParameterList& params, Mortar::Element& sele, Mortar::Element& mele);

    //! @name Evlauation methods

    /*!
    \brief Integrate overlap of slave / master pair (2D)

    Derived version! Most importantly, in this derived version
    a CONTACT::Integrator instance is created, which also
    does integration of the mortar quantity linearizations

    This method integrates the overlap of the current Mortar::Element
    pair sele_ / mele_ based on the integration limits (xiproj). The
    integration includes the Mortar matrices D/M and the gap g.

    */
    bool IntegrateOverlap(const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr) override;

    /*!
    \brief Return type of wear surface definition

    */
    Inpar::Wear::WearSide WearSide()
    {
      return Core::UTILS::IntegralValue<Inpar::Wear::WearSide>(imortar_, "BOTH_SIDED_WEAR");
    }

    /*!
    \brief Return type of wear surface definition

    */
    Inpar::Wear::WearType WearType()
    {
      return Core::UTILS::IntegralValue<Inpar::Wear::WearType>(imortar_, "WEARTYPE");
    }

    //@}


   protected:
    // don't want = operator and cctor
    Coupling2d operator=(const Coupling2d& old) = delete;
    Coupling2d(const Coupling2d& old) = delete;

    // new variables as compared to base class
    Inpar::CONTACT::SolvingStrategy stype_;

  };  // class Coupling2d

  /*!
  \brief A class representing the framework for mortar coupling of ONE
         slave element and SEVERAL master elements of a mortar interface in
         2D. Concretely, this class simply stores several Coupling2d objects.

  */

  class Coupling2dManager : public Mortar::Coupling2dManager
  {
   public:
    /*!
    \brief Constructor with shape function specification

    Constructs an instance of this class and enables custom shape function types.<br>
    Note that this is \b not a collective call as coupling is
    performed in parallel by individual processes.

    */
    Coupling2dManager(Core::FE::Discretization& idiscret, int dim, bool quad,
        Teuchos::ParameterList& params, Mortar::Element* sele, std::vector<Mortar::Element*> mele);


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
    Inpar::Mortar::ShapeFcn ShapeFcn()
    {
      return Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(imortar_, "LM_SHAPEFCN");
    }

    /*!
    \brief Evaluate mortar coupling

    */
    void integrate_coupling(const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr) override;

    bool evaluate_coupling(const Teuchos::RCP<Mortar::ParamsInterface>& mparams_ptr) override;
    //@}
   private:
    /*!
    \brief Calculate consistent dual shape functions in boundary elements

    */
    void consist_dual_shape() override;

   protected:
    // don't want = operator and cctor
    Coupling2dManager operator=(const Coupling2dManager& old) = delete;
    Coupling2dManager(const Coupling2dManager& old) = delete;

    Inpar::CONTACT::SolvingStrategy stype_;  // solving strategy

  };  // class Coupling2dManager

}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
