/*----------------------------------------------------------------------*/
/*!
 \file poroelast_utils.H

 \brief utility functions for poroelasticity problems

\level 2

\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
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

namespace MAT
{
class Material;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace POROELAST
{
class PoroBase;
class PoroScatraBase;

//! coupling types of porous media problem
enum coupltype
{
  fluidfluid, fluidstructure, undefined
};

/// POROELAST::UTILS: Random stuff that might be helpful when dealing with poroelasticity problems
namespace UTILS
{

///check if element is a poro-element
bool CheckPoro(const DRT::Element* actele);

///check if element is a poro-p1-element
bool CheckPoroP1(const DRT::Element* actele);

///check if element is a poro-material
bool CheckPoroMaterial(Teuchos::RCP<const MAT::Material> material);

Teuchos::RCP<LINALG::MapExtractor> BuildPoroSplitter(
    Teuchos::RCP<DRT::Discretization> dis);

/// create solution algorithm depending on input file
Teuchos::RCP<POROELAST::PoroBase> CreatePoroAlgorithm(
    const Teuchos::ParameterList& timeparams, //!< problem parameters (i)
    const Epetra_Comm& comm, //!< communicator(i)
    bool setup_solve = true //!< setup linear solver for Poroelastic problem (only required if Solve() is called) (i)
    );

/// create solution algorithm depending on input file
Teuchos::RCP<POROELAST::PoroScatraBase> CreatePoroScatraAlgorithm(
    const Teuchos::ParameterList& timeparams, //!< problem parameters (i)
    const Epetra_Comm& comm //!< communicator(i)
    );

void SetMaterialPointersMatchingGrid(
    Teuchos::RCP<const DRT::Discretization> sourcedis,
    Teuchos::RCP<const DRT::Discretization> targetdis);

/*!
 Create volume ghosting:
 To ease the contact search algorithms we'll afford the luxury to ghost all nodes
 on all processors in the general mortar coupling framework.
 Here we additionally ghost volume elements and their nodes
 (required if an evaluation of gradients is required)!
 */
void CreateVolumeGhosting(DRT::Discretization& idiscret); //redistributed interface discretization of contact!

/*! Reconnect Face Element - Parent Element Pointers!
 Parent Element need to be ghosted on the processors where Face Elements
 exist already.
 */
 void ReconnectParentPointers(DRT::Discretization& idiscret,
                             DRT::Discretization& voldiscret);

//! Determine norm of vector
double CalculateVectorNorm(const enum INPAR::POROELAST::VectorNorm norm, //!< norm to use
    const Teuchos::RCP<const Epetra_Vector> vect //!< the vector of interest
    );

//! strategy for material assignment for non matching meshes with poro

/// Helper class for assigning materials for volumetric coupling of non conforming meshes (poro)
/*!
 When coupling two overlapping discretizations, most often one discretization needs access
 to the corresponding element/material on the other side. For conforming meshes this is straight forward
 as there is one unique element on the other side and therefore one unique material, which can
 be accessed. However, for non conforming meshes there are potentially several elements overlapping.
 Therefore, some rule for assigning materials is needed. This class is meant to do that.
 It gets the element to which it shall assign a material and a vector of IDs of the overlapping elements
 of the other discretization.

 In case of poro, we also need the initial porosity type of the structural element to be known in the fluid
 element, which is why there is a special strategy for poro. Note that this is not yet working for
 inhomogeneous material properties.

 \author vuong 10/14
 */
class PoroMaterialStrategy: public VOLMORTAR::UTILS::DefaultMaterialStrategy
{
public:
  //! constructor
  PoroMaterialStrategy(){};

  //! destructor
  virtual ~PoroMaterialStrategy(){};

  //! assignment of fluid material to structure material
  virtual void AssignMaterial2To1(
      const VOLMORTAR::VolMortarCoupl* volmortar,
      DRT::Element* ele1, const std::vector<int>& ids_2,
      Teuchos::RCP<DRT::Discretization> dis1,
      Teuchos::RCP<DRT::Discretization> dis2);

  //! assignment of structure material to fluid material
  virtual void AssignMaterial1To2(
      const VOLMORTAR::VolMortarCoupl* volmortar,
      DRT::Element* ele2, const std::vector<int>& ids_1,
      Teuchos::RCP<DRT::Discretization> dis1,
      Teuchos::RCP<DRT::Discretization> dis2);
};
} // namespace UTILS

void PrintLogo();

} // namespace POROELAST

/*----------------------------------------------------------------------*/
#endif // POROELAST_UTILS_H
