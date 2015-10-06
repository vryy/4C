/*
 * MueLu_SelectiveSaPFactory.cpp
 *
 *  Created on: Jan 17, 2013
 *      Author: wiesner
 */


#ifdef HAVE_MueLu

#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_SelectiveSaPFactory_def.hpp"

template class MueLu::SelectiveSaPFactory<double, int, int, KokkosClassic::DefaultNode::DefaultNodeType>;

#ifdef HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT
# ifdef HAVE_TEUCHOS_LONG_LONG_INT
template class MueLu::SelectiveSaPFactory<double, int, long long int, KokkosClassic::DefaultNode::DefaultNodeType>;
# else
# warning To compile MueLu with 'long long int' support, please turn on Teuchos_ENABLE_LONG_LONG_INT
# endif
#endif

//#ifdef HAVE_MUELU_INST_COMPLEX_INT_INT
//# ifdef HAVE_TEUCHOS_COMPLEX
//#include <complex>
//template class MueLu::SelectiveSaPFactory<std::complex<double>, int, int, KokkosClassic::DefaultNode::DefaultNodeType>;
//# else
//# warning To compile MueLu with 'complex' support, please turn on Teuchos_ENABLE_COMPLEX
//# endif
//#endif

#endif // HAVE_MueLu

