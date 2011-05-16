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
#include "combust_defines.H"

namespace COMBUST
{

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

    //! compute volume of domain minus
    const double ComputeVolumeMinus();

    //! output of domain and boundary integration cells to Gmsh!
    void toGmsh(const int step) const;

    //! tell me, whether the element has a whole touched face and lies in the plus domain or not
    bool ElementTouched(const int xfemeleid) const;

    //! tell me, whether the element is bisected or not
    bool ElementBisected(const int xfemeleid) const;

    //! tell me, whether the element is trisected or not
    bool ElementTrisected(const int xfemeleid) const;

    //! tell me, whether the element is trisected or not
    COMBUST::FlameFront::CutStatus ElementCutStatus(const int xfemeleid) const;


  private:

    /// don't want default constructor, therefore declared private!
    explicit InterfaceHandleCombust();

    /// don't want copy constructor, therefore declared private!
    explicit InterfaceHandleCombust(const InterfaceHandleCombust&);

    /// pointer to level set (G-function) discretization that intersects the fluiddis_
    const Teuchos::RCP<const DRT::Discretization> gfuncdis_;

    /// pointer to the flame front
    const Teuchos::RCP<const COMBUST::FlameFront> flamefront_;

    /// map holding status of elements with respect to the interface
    std::map<int,COMBUST::FlameFront::CutStatus > elementcutstatus_;

  }; // class InterfaceHandleCombust
} // namespace COMBUST

#endif // COMBUST_INTERFACE_H
#endif // CCADISCRET
