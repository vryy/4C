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
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "../drt_inpar/inpar_poroelast.H"

#include "../drt_volmortar/volmortar_utils.H"

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

    //! Determine norm of vector
    double CalculateVectorNorm(
      const enum INPAR::POROELAST::VectorNorm norm,  //!< norm to use
      const Teuchos::RCP<const Epetra_Vector> vect  //!< the vector of interest
      );

    //! strategy for material assignment for non matching meshes with TSI
    class PoroMaterialStrategy: public VOLMORTAR::UTILS::DefaultMaterialStrategy
    {
    public:
      //! constructor
      PoroMaterialStrategy(){};

      //! destructor
      virtual ~PoroMaterialStrategy(){};

      //! assignment of thermo material to structure material
      virtual void AssignMaterialBToA(const VOLMORTAR::VolMortarCoupl* volmortar,
                                  DRT::Element* Aele,
                                  const std::vector<int>& Bids,
                                  Teuchos::RCP<DRT::Discretization> disA,
                                  Teuchos::RCP<DRT::Discretization> disB);

      //! assignment of structure material to thermo material
      virtual void AssignMaterialAToB(const VOLMORTAR::VolMortarCoupl* volmortar,
                                  DRT::Element* Bele,
                                  const std::vector<int>& Aids,
                                  Teuchos::RCP<DRT::Discretization> disA,
                                  Teuchos::RCP<DRT::Discretization> disB);
    };
  } // namespace UTILS

  void PrintLogo();

} // namespace POROELAST

/*----------------------------------------------------------------------*/
#endif // POROELAST_UTILS_H
