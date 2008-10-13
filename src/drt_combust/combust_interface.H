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
             const Teuchos::RCP<DRT::Discretization>  xfemdis,   ///< discretization (combust3 elements) which will be intersected
             const Teuchos::RCP<DRT::Discretization>  gfuncdis   ///< discretization (condif3 elements) with level set function (G-function)
             );
    /// destructor
    ~InterfaceHandleCombust();

    //! return pointer to the intersecting discretization
    Teuchos::RCP<DRT::Discretization> gfuncdis() const {return gfuncdis_;};

    //! implement this function if needed for combustion!
    void toGmsh(const int step) const;

    //! implement this function if needed for combustion!
    int PositionWithinConditionNP(const BlitzVec3& x_in) const;

    //! implement this function if needed for combustion!
    int PositionWithinConditionN(const BlitzVec3& x_in) const;

    //! implement this function if needed for combustion!
    int PositionWithinConditionNP(const BlitzVec3&     x_in,
                                  GEO::NearestObject&  nearestobject) const;

    //! implement this function if needed for combustion!
    int PositionWithinConditionN(const BlitzVec3&     x_in,
                                 GEO::NearestObject&  nearestobject) const;

  private:

    /// don't want default constructor
    explicit InterfaceHandleCombust();
    
    /// don't want copy constructor
    explicit InterfaceHandleCombust(const InterfaceHandleCombust&);

    /// discretization that intersects the fluiddis_
    const Teuchos::RCP<DRT::Discretization> gfuncdis_;

  }; // class InterfaceHandleCombust
} // namespace COMBUST

#endif // COMBUST_INTERFACE_H
#endif // CCADISCRET
