#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q1_2013

#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_MeshtyingSPAmalgamationFactory_def.hpp"

template class MueLu::MeshtyingSPAmalgamationFactory<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void, int, Kokkos::DefaultNode::DefaultNodeType>::SparseOps>;

#ifdef HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT
# ifdef HAVE_TEUCHOS_LONG_LONG_INT
template class MueLu::MeshtyingSPAmalgamationFactory<double, int, long long int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void, int, Kokkos::DefaultNode::DefaultNodeType>::SparseOps>;
# else
# warning To compile MueLu with 'long long int' support, please turn on Teuchos_ENABLE_LONG_LONG_INT
# endif
#endif

#ifdef HAVE_MUELU_INST_COMPLEX_INT_INT
# ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
template class MueLu::MeshtyingSPAmalgamationFactory<std::complex<double>, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<void, int, Kokkos::DefaultNode::DefaultNodeType>::SparseOps>;
# else
# warning To compile MueLu with 'complex' support, please turn on Teuchos_ENABLE_COMPLEX
# endif
#endif
#endif // Q1/2013
#endif // HAVE_MueLu
