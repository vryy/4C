#ifdef TRILINOS_PACKAGE
#ifndef SOLVER_TRILINOS_ML_H
#define SOLVER_TRILINOS_ML_H

#ifdef __cplusplus
extern "C" 
{
#endif

#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#ifdef D_SHELL8
#include "../shell8/shell8.h"
#endif

#ifdef __cplusplus
}
#endif

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
using namespace Teuchos;

void create_ml_parameterlist(struct _SOLVAR         *actsolv,
                             ParameterList&          mllist,
                             FIELD*                  actfield,
                             int                     disnum,
                             const Epetra_CrsMatrix& matrix,
                             double**                nullspace);


#endif /* SOLVER_TRILINOS_ML_H */
#endif /* TRILINOS_PACKAGE */
