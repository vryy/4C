/*----------------------------------------------------------------------*/
/*!
 \file poroelast_utils.cpp

 \brief utility functions for poroelasticity problems

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

#ifndef POROELAST_UTILS_H_
#define POROELAST_UTILS_H_

/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "../drt_scatra/scatra_utils_clonestrategy.H"

/*----------------------------------------------------------------------*
 |  forward declaration                                                 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Element;
  class Discretization;
}

namespace LINALG
{
  class MapExtractor;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace POROELAST
{
  class PoroBase;
  class PORO_SCATRA_Base;

  //! coupling types of porous media problem
  enum coupltype
  {
    fluidfluid,
    fluidstructure,
    undefined
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
      bool DetermineEleType(
          DRT::Element*             actele,
          const bool                ismyele,
          std::vector<std::string>& eletype);
      /// set element-specific data (material etc.)
      void SetElementData(
          Teuchos::RCP<DRT::Element> newele,
          DRT::Element*              oldele,
          const int                  matid,
          const bool                 isnurbs);

      /// returns conditions names to be copied (source and target name)
      std::map<std::string, std::string> ConditionsToCopy();

      /// check for correct material
      void CheckMaterialType(const int matid);

    }; // class PoroelastCloneStrategy

    //! \brief implementation of special clone strategy for automatic generation
    //!        of scatra discretization from a given structure discretization for porous media
    class PoroScatraCloneStrategy: public SCATRA::ScatraFluidCloneStrategy {
    public:

      /// constructor
      explicit PoroScatraCloneStrategy() {}
      /// destructor
      virtual ~PoroScatraCloneStrategy() {}

    protected:

      /// determine element type string and whether element is copied or not
      virtual bool DetermineEleType(
          DRT::Element* actele,
          const bool ismyele,
          std::vector<std::string>& eletype);
    }; // class PoroScatraCloneStrategy

    /// setup poro discretization,i.e. clone the structural discretization
    void SetupPoro();

    ///check if element is a poro-element
    bool CheckPoro(DRT::Element* actele);

    ///check if element is a poro-p1-element
    bool CheckPoroP1(DRT::Element* actele);

    Teuchos::RCP<LINALG::MapExtractor> BuildPoroSplitter(Teuchos::RCP<DRT::Discretization> dis);

    /// create solution algorithm depending on input file
    Teuchos::RCP<POROELAST::PoroBase> CreatePoroAlgorithm(
                             const Teuchos::ParameterList& timeparams, //!< problem parameters (i)
                             const Epetra_Comm& comm //!< communicator(i)
                             );

    /// create solution algorithm depending on input file
    Teuchos::RCP<POROELAST::PORO_SCATRA_Base> CreatePoroScatraAlgorithm(
                             const Teuchos::ParameterList& timeparams, //!< problem parameters (i)
                             const Epetra_Comm& comm //!< communicator(i)
                             );

    void SetMaterialPointersMatchingGrid(
        Teuchos::RCP<const DRT::Discretization> sourcedis,
        Teuchos::RCP<const DRT::Discretization> targetdis);

  } // namespace UTILS

  void PrintLogo();

} // namespace POROELAST

/*----------------------------------------------------------------------*/
#endif // POROELAST_UTILS_H
