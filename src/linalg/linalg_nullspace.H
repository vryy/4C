/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of calculation routines related to Nullspaces

\level 0

*----------------------------------------------------------------------*/
#ifndef LINALG_NULLSPACE_H
#define LINALG_NULLSPACE_H

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_elementtype.H"

namespace LINALG
{
  namespace NULLSPACE
  {
    /*!
       \brief Calculate the nullspace based on a given discretization

       The nullspace is build by looping over all nodes of a discretization and stored
       in the respective variable.

       \param dis (in): discretization
       \param nullspace (in): nullspace
       \param numdf (in): number of degrees of freedom
       \param dimns (in): nullspace dimension
    */
    Teuchos::RCP<Epetra_MultiVector> ComputeNullSpace(
        const DRT::Discretization& dis, const int numdf, const int dimns);

    /*!
       \brief Fix the nullspace to match a new given map

       The nullspace is looked for in the parameter list. If found, it is assumed that
       it matches the oldmap. Then it is fixed to match the new map.

       \param field (in): field name (just used for output)
       \param oldmap (in): row map of nullspace
       \param newmap (in): row map of nullspace upon exit
       \param solveparams (in): parameterlist including nullspace vector
    */
    void FixNullSpace(std::string field, const Epetra_Map& oldmap, const Epetra_Map& newmap,
        Teuchos::ParameterList& solveparams);

  }  // namespace NULLSPACE
}  // namespace LINALG

#endif
