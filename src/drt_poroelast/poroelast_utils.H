/*----------------------------------------------------------------------*/
/*!
 \file poroelast_utils.cpp

 \brief utility functions for poroelasticity problems

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef POROELAST_UTILS_H_
#define POROELAST_UTILS_H_

#include <set>
#include <string>
#include <functional>
#include <iostream>

/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_condition_utils.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace POROELAST
{
  /// POROELAST::UTILS: Random stuff that might be helpful when dealing with poroelasticity problems
  namespace UTILS
  {
    //! \brief implementation of special clone strategy for automatic generation
    //!        of fluid discretization from a given structure discretization
    class PoroelastCloneStrategy {
    public:

      /// constructor
      explicit PoroelastCloneStrategy() {}
      /// destructor
      virtual ~PoroelastCloneStrategy() {}

    protected:

      /// determine element type string and whether element is copied or not
      bool DetermineEleType(DRT::Element* actele, const bool ismyele,
          vector<string>& eletype);
      /// set element-specific data (material etc.)
      void SetElementData(Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele,
          const int matid, const bool isnurbs);

      /// returns conditions names to be copied (source and target name)
      map<string, string> ConditionsToCopy();

      /// check for correct material
      void CheckMaterialType(const int matid);

    private:

    }; // class PoroelastCloneStrategy

    /// setup Poroelasticity, clone the structural discretization
    void SetupPoro(const Epetra_Comm& comm);

  } // namespace UTILS

} // namespace POROELAST

/*----------------------------------------------------------------------*/
#endif // POROELAST_UTILS_H
#endif // CCADISCRET
