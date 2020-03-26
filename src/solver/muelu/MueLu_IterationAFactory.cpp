/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu iteration factory class
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifdef TRILINOS_Q1_2015

#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_NodeDefinition.hpp"

#include "MueLu_IterationAFactory_def.hpp"
#include "MueLu_NodeDefinition.hpp"

template class MueLu::IterationAFactory<double, int, int, KokkosSerialNode>;

#ifdef HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
template class MueLu::IterationAFactory<double, int, long long int, KokkosSerialNode>;
#else
#warning To compile MueLu with 'long long int' support, please turn on Teuchos_ENABLE_LONG_LONG_INT
#endif
#endif

#ifdef HAVE_MUELU_INST_COMPLEX_INT_INT
#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
template class MueLu::IterationAFactory<std::complex<double>, int, int, KokkosSerialNode>;
#else
#warning To compile MueLu with 'complex' support, please turn on Teuchos_ENABLE_COMPLEX
#endif
#endif
#endif  // TRILINOS_Q1_2015
