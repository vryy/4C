/*
 * MueLu_ContactSPRepartitionInterface.cpp
 *
 *  Created on: 11 Sep 2013
 *      Author: wiesner
 */


#include "MueLu_ExplicitInstantiation.hpp"

#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q3_2013

#include "MueLu_ContactSPRepartitionInterface_def.hpp"

#ifdef HAVE_MUELU_INST_DOUBLE_INT_INT
template class MueLu::ContactSPRepartitionInterface<int, int, KokkosClassic::DefaultNode::DefaultNodeType, KokkosClassic::DefaultKernels<void, int, KokkosClassic::DefaultNode::DefaultNodeType>::SparseOps>;
#endif

#ifdef HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT
# ifdef HAVE_TEUCHOS_LONG_LONG_INT
template class MueLu::ContactSPRepartitionInterface<int, long long int, KokkosClassic::DefaultNode::DefaultNodeType, KokkosClassic::DefaultKernels<void, int, KokkosClassic::DefaultNode::DefaultNodeType>::SparseOps>;
# else
# warning To compile MueLu with 'long long int' support, please turn on Teuchos_ENABLE_LONG_LONG_INT
# endif
#endif

#endif
#endif
