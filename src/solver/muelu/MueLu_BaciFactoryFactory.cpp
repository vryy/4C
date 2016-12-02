/*!----------------------------------------------------------------------
\file MueLu_BaciFactoryFactory.cpp

\brief MueLu factory class for BACI
\level 2
\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#ifdef HAVE_MueLu

#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_BaciFactoryFactory_def.hpp"

#ifdef HAVE_MUELU_INST_DOUBLE_INT_INT
template class MueLu::BaciFactoryFactory<double, int, int>;
#endif

#ifdef HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT
# ifdef HAVE_TEUCHOS_LONG_LONG_INT
template class MueLu::BaciFactoryFactory<double, int, long long int>;
# else
# warning To compile MueLu with 'long long int' support, please turn on Teuchos_ENABLE_LONG_LONG_INT
# endif
#endif

#ifdef HAVE_MUELU_INST_COMPLEX_INT_INT
# ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
template class MueLu::BaciFactoryFactory<std::complex<double>, int, int>;
# else
# warning To compile MueLu with 'complex' support, please turn on Teuchos_ENABLE_COMPLEX
# endif
#endif

#endif /* HAVE_MueLu */
