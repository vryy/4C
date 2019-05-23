/*!----------------------------------------------------------------------

\brief MueLu contact aggregation factory class for saddle point formulations
\level 2
\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/
#ifdef HAVE_MueLu

#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_NodeDefinition.hpp"
#include "MueLu_ContactSPAggregationFactory_def.hpp"

template class MueLu::ContactSPAggregationFactory<double, int, int, KokkosSerialNode>;

#ifdef HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
template class MueLu::ContactSPAggregationFactory<double, int, long long int, KokkosSerialNode>;
#else
#warning To compile MueLu with 'long long int' support, please turn on Teuchos_ENABLE_LONG_LONG_INT
#endif
#endif

//#ifdef HAVE_MUELU_INST_COMPLEX_INT_INT
//# ifdef HAVE_TEUCHOS_COMPLEX
//#include <complex>
// template class MueLu::ContactSPAggregationFactory<std::complex<double>, int, int,
// KokkosSerialNode>; # else # warning To compile MueLu with 'complex' support, please turn on
// Teuchos_ENABLE_COMPLEX # endif #endif
#endif  // HAVE_MueLu
