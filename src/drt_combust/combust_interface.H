/*!-----------------------------------------------------------------------------------------------*
 \file combust_interface.H

 \brief interface handle that transports the intersection related things around for combustion problems

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef COMBUST_INTERFACE_H
#define COMBUST_INTERFACE_H

#include "../drt_xfem/interface.H"
#include "combust_flamefront.H"

namespace COMBUST
{

  class DRT::Discretization;

  /// COMBUST Interface class
  /*!

    detailed description

    \author henke
    \date 10/08
   */
  class InterfaceHandleCombust : public XFEM::InterfaceHandle
  {
  public:
    /// constructor
    explicit InterfaceHandleCombust(
             const Teuchos::RCP<DRT::Discretization> fluiddis,  ///< discretization (combust3 elements) which will be intersected
             const Teuchos::RCP<const DRT::Discretization> gfuncdis,  ///< discretization (condif3 elements) with level set function (G-function)
             const Teuchos::RCP<const COMBUST::FlameFront> flamefront  ///< geometrical interface (flame front) information
             );

    /// destructor
    ~InterfaceHandleCombust();

    //! return pointer to the intersecting discretization (does anyone need this?)
    const Teuchos::RCP<const DRT::Discretization>& GfuncDis() const { return gfuncdis_; }

    //! return pointer to flame front (does anyone need this?)
    const Teuchos::RCP<const COMBUST::FlameFront>& FlameFront() const { return flamefront_; }

    //! update interfacehandle according to current flame front
    void UpdateInterfaceHandle();

    //! implement this function if needed for combustion!
    void toGmsh(const int step) const;

    //! implement this function if needed for combustion!
    int PositionWithinConditionNP(const LINALG::Matrix<3,1>& x_in) const;

    //! implement this function if needed for combustion!
    int PositionWithinConditionN(const LINALG::Matrix<3,1>& x_in) const;

    //! implement this function if needed for combustion!
    int PositionWithinConditionNP(const LINALG::Matrix<3,1>&     x_in,
                                  GEO::NearestObject&  nearestobject) const;

    //! implement this function if needed for combustion!
    int PositionWithinConditionN(const LINALG::Matrix<3,1>&     x_in,
                                 GEO::NearestObject&  nearestobject) const;

  private:

    /// don't want default constructor, therefore declared private!
    explicit InterfaceHandleCombust();

    /// don't want copy constructor, therefore declared private!
    explicit InterfaceHandleCombust(const InterfaceHandleCombust&);

    /// pointer to level set (G-function) discretization that intersects the fluiddis_
    const Teuchos::RCP<const DRT::Discretization> gfuncdis_;

    /// pointer to the flame front
    const Teuchos::RCP<const COMBUST::FlameFront> flamefront_;

  }; // class InterfaceHandleCombust
} // namespace COMBUST

#endif // COMBUST_INTERFACE_H
#endif // CCADISCRET
