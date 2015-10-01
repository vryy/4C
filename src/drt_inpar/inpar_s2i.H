/*----------------------------------------------------------------------*/
/*!
\file inpar_s2i.H

\brief globally accessible enumerations for scatra-scatra interface coupling

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

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace S2I
  {
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
      kinetics_butlervolmer
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


    /// set the s2i parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  } // namespace S2I
} // namespace INPAR
#endif // #ifndef INPAR_S2I_H
