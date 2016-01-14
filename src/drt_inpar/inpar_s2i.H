/*----------------------------------------------------------------------*/
/*!
\file inpar_s2i.H

\brief input quantities and globally accessible enumerations for scatra-scatra interface coupling

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_S2I_H
#define INPAR_S2I_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}

namespace INPAR
{
  namespace S2I
  {
    //! type of interface side
    enum InterfaceSides
    {
      side_undefined,
      side_slave,
      side_master
    };

    //! type of mortar meshtying
    enum MortarType
    {
      mortar_undefined,
      mortar_none,
      mortar_standard,
      mortar_saddlepoint,
      mortar_condensed
    };

    //! models for interface kinetics
    enum KineticModels
    {
      kinetics_constperm,
      kinetics_butlervolmer,
      kinetics_butlervolmerpeltier
    };

    //! type of global system matrix in global system of equations
    enum MatrixType
    {
      matrix_sparse,
      matrix_block_geometry,
      matrix_block_condition,
      matrix_block_condition_dof
    };

    //! equilibration methods for global system of equations
    enum EquilibrationMethods
    {
      equilibration_none,
      equilibration_rows,
      equilibration_columns,
      equilibration_full
    };

    //! set valid parameters for scatra-scatra interface coupling
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    //! set valid conditions for scatra-scatra interface coupling
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);
  } // namespace S2I
} // namespace INPAR
#endif // #ifndef INPAR_S2I_H
