/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu repartition algorithm class for contact
\level 2

*----------------------------------------------------------------------*/

#ifdef TRILINOS_Q1_2015

#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_ContactSPRepartitionInterface_def.hpp"
#include "MueLu_NodeDefinition.hpp"

#ifdef HAVE_MUELU_INST_DOUBLE_INT_INT
template class MueLu::ContactSPRepartitionInterface<int, int, KokkosSerialNode,
    KokkosClassic::DefaultKernels<void, int, KokkosSerialNode>::SparseOps>;
#endif

#ifdef HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
template class MueLu::ContactSPRepartitionInterface<int, long long int, KokkosSerialNode,
    KokkosClassic::DefaultKernels<void, int, KokkosSerialNode>::SparseOps>;
#else
#warning To compile MueLu with 'long long int' support, please turn on Teuchos_ENABLE_LONG_LONG_INT
#endif
#endif

#endif  // TRILINOS_Q1_2015
