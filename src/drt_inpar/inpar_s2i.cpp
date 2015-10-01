/*----------------------------------------------------------------------*/
/*!
\file inpar_s2i.cpp

\brief Input parameters for s2i

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_s2i.H"



void INPAR::S2I::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& s2idyn = list->sublist(
      "S2I DYNAMIC",
      false,
      "control parameters for scatra-scatra interface coupling"
      );

  // type of global system matrix in global system of equations
  setStringToIntegralParameter<int>(
      "MATRIXTYPE",
      "sparse",
      "type of global system matrix in global system of equations",
      tuple<std::string>(
          "sparse",
          "block_geometry",
          "block_condition",
          "block_condition_dof"
          ),
      tuple<int>(
          matrix_sparse,
          matrix_block_geometry,
          matrix_block_condition,
          matrix_block_condition_dof
          ),
      &s2idyn
      );

  // type of mortar meshtying
  setStringToIntegralParameter<int>(
      "MORTARTYPE",
      "Undefined",
      "type of mortar meshtying",
      tuple<std::string>(
          "Undefined",
          "NoMortar",
          "StandardMortar",
          "SaddlePointMortar",
          "CondensedMortar"
          ),
      tuple<int>(
          mortar_undefined,
          mortar_none,
          mortar_standard,
          mortar_saddlepoint,
          mortar_condensed
          ),
      &s2idyn
      );

  // flag for equilibration of global system of equations
  setStringToIntegralParameter<int>(
      "EQUILIBRATION",
      "none",
      "flag for equilibration of global system of equations",
      tuple<std::string>(
          "none",
          "rows",
          "columns",
          "full"
          ),
      tuple<int>(
          equilibration_none,
          equilibration_rows,
          equilibration_columns,
          equilibration_full
          ),
      &s2idyn
      );
}
