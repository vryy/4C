/*-----------------------------------------------------------*/
/*!
\file inpar_biopolynet_dbc.H

\brief input parameter for statistical mechanic problem

\maintainer Jonas Eichinger, Maximilian Grill

\level 2

*/
/*-----------------------------------------------------------*/
#ifndef INPAR_BEAMINTERACTION_DBC_H_
#define INPAR_BEAMINTERACTION_DBC_H_

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters used in statistical mechanics
  namespace BEAMINTERACTIONDBC
  {

    /// (this enum represents the input file parameter DBCTYPE)
    enum DBCType
    {
      dbctype_none,        ///< no Dirichlet Boundary conditions
      dbctype_std,          ///< standard Dirichlet Boundary conditions
      dbctype_shearfixed,  ///< special DBC application in case of rheological measurements with fixed set of nodes
      dbctype_shearfixeddel,  ///< special DBC application in case of rheological measurements with fixed set of nodes, delete cut elements
      dbctype_sheartrans,  ///< special DBC application in case of rheological measurements with transient set of nodes
      dbctype_pinnodes,    ///< pin down individual nodes when a linker is attached
      dbctype_affineshear,  ///< apply affine shear displacement to all nodes
      dbctype_affinesheardel,  ///< apply affine shear displacement to all nodes, delete cut elements
      dbctype_movablesupport1d,  ///< constrain all displacement except for one direction
    };

    /// set the statmech parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }

}

/*----------------------------------------------------------------------*/
#endif





#endif /* SRC_DRT_INPAR_INPAR_BEAMINTERACTION_DBC_H_ */
