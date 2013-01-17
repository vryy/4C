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

#ifndef POROELAST_UTILS_H_
#define POROELAST_UTILS_H_

/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "../drt_inpar/inpar_poroelast.H"

namespace DRT
{
  class Element;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace POROELAST
{
  class PoroBase;

  enum coupltype
  {
    fluidfluid,
    fluidstructure
  };

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
          std::vector<std::string>& eletype);
      /// set element-specific data (material etc.)
      void SetElementData(Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele,
          const int matid, const bool isnurbs);

      /// returns conditions names to be copied (source and target name)
      std::map<std::string, std::string> ConditionsToCopy();

      /// check for correct material
      void CheckMaterialType(const int matid);

    private:
      ///check if element is a poro-element
      bool CheckPoro(DRT::Element* actele);

    }; // class PoroelastCloneStrategy

    /// setup poro discretization,i.e. clone the structural discretization
    void SetupPoro();

    /// setup poro discretization,i.e. clone the structural discretization
    Teuchos::RCP<POROELAST::PoroBase> CreatePoroAlgorithm(
                             const Teuchos::ParameterList& timeparams, //!< coupling scheme (i)
                             const Epetra_Comm& comm //!< communicator(i)
                             );

  } // namespace UTILS

  void PrintLogo();

} // namespace POROELAST

/*----------------------------------------------------------------------*/
#endif // POROELAST_UTILS_H
