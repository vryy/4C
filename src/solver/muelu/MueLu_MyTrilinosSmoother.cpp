/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu smoother interface
\level 2

*----------------------------------------------------------------------*/

#ifdef TRILINOS_Q1_2015

#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_MyTrilinosSmoother_def.hpp"
#include "MueLu_NodeDefinition.hpp"

#ifdef HAVE_MUELU_INST_DOUBLE_INT_INT
template class MueLu::MyTrilinosSmoother<double, int, int, KokkosSerialNode>;
#endif

#ifdef HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
template class MueLu::MyTrilinosSmoother<double, int, long long int, KokkosSerialNode>;
#else
#warning To compile MueLu with 'long long int' support, please turn on Teuchos_ENABLE_LONG_LONG_INT
#endif
#endif

#ifdef HAVE_MUELU_INST_COMPLEX_INT_INT
#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
template class MueLu::MyTrilinosSmoother<std::complex<double>, int, int, KokkosSerialNode>;
#else
#warning To compile MueLu with 'complex' support, please turn on Teuchos_ENABLE_COMPLEX
#endif
#endif

#ifdef HAVE_MUELU_INST_PCE_INT_INT
#include "Stokhos_Sacado.hpp"
typedef Stokhos::StandardStorage<int, double> Storage;
typedef Sacado::PCE::OrthogPoly<double, Storage> pce_type;
template class MueLu::MyTrilinosSmoother<pce_type, int, int, KokkosSerialNode>;
#endif

#endif  // TRILINOS_Q1_2015
