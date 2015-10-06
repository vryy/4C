/*
 * MueLu_IterationAFactory.cpp
 *
 *  Created on: Jan 10, 2013
 *      Author: tobias
 */

#ifdef HAVE_MueLu

#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_IterationAFactory_def.hpp"

template class MueLu::IterationAFactory<double, int, int, KokkosClassic::DefaultNode::DefaultNodeType>;

#ifdef HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT
# ifdef HAVE_TEUCHOS_LONG_LONG_INT
template class MueLu::IterationAFactory<double, int, long long int, KokkosClassic::DefaultNode::DefaultNodeType>;
# else
# warning To compile MueLu with 'long long int' support, please turn on Teuchos_ENABLE_LONG_LONG_INT
# endif
#endif

#ifdef HAVE_MUELU_INST_COMPLEX_INT_INT
# ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
template class MueLu::IterationAFactory<std::complex<double>, int, int, KokkosClassic::DefaultNode::DefaultNodeType>;
# else
# warning To compile MueLu with 'complex' support, please turn on Teuchos_ENABLE_COMPLEX
# endif
#endif
#endif // HAVE_MueLu


