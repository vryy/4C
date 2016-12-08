#ifdef HAVE_Kokkos
/*----------------------------------------------------------------------*/
/*!
\file objective_funct_surfcurr.H
\brief Surface current based objective function
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
*/
/*----------------------------------------------------------------------*/

#ifndef OBJECTIVE_FUNCT_SURFCURR_H
#define OBJECTIVE_FUNCT_SURFCURR_H

#include "objective_funct.H"
//#include "../drt_lib/drt_globalproblem.H"

#include "Epetra_Map.h"

#include <Kokkos_Core.hpp>

// forward declarations
namespace DRT
{
  class Discretization;
  class Condition;
  class Problem;

  namespace INPUT
  {
    class DatFileReader;
  }
}

namespace LINALG
{
  class SerialDenseMatrix;
}

namespace INVANA
{

class SurfCurrentPair;
class Triangulation;

// the kernel
inline double kernel(double x1, double x2, double x3,
                     double y1, double y2, double y3,
                     double lscale) {
  return exp(-((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2)+(x3-y3)*(x3-y3))/lscale);}


/*!
  \brief An objective function comparing sets of surfaces

  This is a handler to have easy access to currents based on different conditions,
  i.e every source surface current condition should be evaluated against the corresponding
  target surface current condition.

  This object just handles correspondences between surfaces to be compared against each other.
  The necessary information comes from the input files of meshes to be compared with each other.

  The actual computation of a norm on surfaces and the gradient thereof is forwarded to
  INVANA::SurfCurrentPair.

  The whole objective function is made up as sums of pair-wise norm.

  Eventually it is expected that this class might also manage the tracking of surfaces warped by
  the solution of the timestepping of a problem. But sofar only one timestep, which is considered the
  "final" timestep is taken into account.
  */
class SurfCurrentGroup : public ObjectiveFunct
{

private:

  //! Discretization of the source model
  Teuchos::RCP<DRT::Discretization> sourcedis_;

  //! Discretization of the target model
  Teuchos::RCP<DRT::Discretization> targetdis_;

  //! Measured timesteps
  std::vector<double> timesteps_;

  /*! \brief Vector of single surface objective functions
   * 'std::vector wrap' to account for measurement
   * for more than the last time step in the future
   */
  std::vector< std::vector<Teuchos::RCP<INVANA::SurfCurrentPair> > > currents_;


public:

  //! Constructor
  SurfCurrentGroup(Teuchos::RCP<DRT::Discretization> discret);

  //! Standard destructor
  virtual ~SurfCurrentGroup(){Kokkos::finalize();}

  //! read target discretization
  DRT::Problem* ReadReferenceDiscretization();

  virtual double GetScaleFac();

  //! Evaluate value of the objective function
  virtual void Evaluate(Teuchos::RCP<Epetra_Vector> state, double time, double& val);

  //! Evaluate gradient of objective function w.r.t to displacements
  virtual void EvaluateGradient(Teuchos::RCP<Epetra_Vector> state, double time, Teuchos::RCP<Epetra_Vector> gradient);

  //! Get measured time steps
  virtual std::vector<double> MeasuredTime() {return timesteps_;}

  //! Query for step of measurements of this time step
  int FindStep(double time);

}; //SurfCurrentGroup


//! @name Typedefs used by INVANA::SurfCurrentPair
//@{

//! Kokkos::View to hold data of type centers and normals
typedef Kokkos::View<double*[3]> currents_type;
typedef currents_type::HostMirror currents_host_type;

//! Kokkos::View to hold data of type points (and associated data) of a triangulation
typedef Kokkos::View<double*[9]> trimesh_type;
typedef trimesh_type::HostMirror trimesh_host_type;
typedef Kokkos::View<double*[9],Kokkos::MemoryTraits<Kokkos::Atomic>> trimesh_atomic_type;

//! Kokkos::View to hold data of type derivative of currents_type wrt nodal dofs
typedef Kokkos::View<double*[27]> trimesh_vecdata_type;

//! Strided Layout to let Kokkos process only blocks of data
typedef Kokkos::View<double**, Kokkos::LayoutStride> ViewStride;

// Map keeping trimesh data used by the Triangulation
typedef std::map<int, std::vector<double> > extract_type;

//@}



/*! @name Functors used by INVANA::SurfCurrentPair
 *
 * These functors are used by the Kokkos-based evaluation of the surface
 * current norm and its gradient wrt. to some global dofs
 */

//! Normal functor
struct Normal {
  currents_type n;
  trimesh_type p;

  // Constructor takes View by "value"; this does a shallow copy.
  Normal (currents_type n_, trimesh_type p_) : n (n_),p (p_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    n(i,0)=0.5*((p(i,7)-p(i,1))*(p(i,2)-p(i,5))-(p(i,8)-p(i,2))*(p(i,1)-p(i,4)));
    n(i,1)=0.5*((p(i,8)-p(i,2))*(p(i,0)-p(i,3))-(p(i,6)-p(i,0))*(p(i,2)-p(i,5)));
    n(i,2)=0.5*((p(i,6)-p(i,0))*(p(i,1)-p(i,4))-(p(i,7)-p(i,1))*(p(i,0)-p(i,3)));
  }
};

//! Derivative of normals wrt to triangular coordinates functor
struct DNormal {
  trimesh_type p;
  trimesh_vecdata_type dn;

  // Constructor takes View by "value"; this does a shallow copy.
  DNormal (trimesh_type p_, trimesh_vecdata_type dn_) : p (p_),dn (dn_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    // derivative of normal (see Normal functor) wrt the single dofs
    // -> a 3by9 matrix [dN/du1 dN/du2 ... dN/du9]

    // dofs of the first node
    // dN/du1 -> dN/p(i,0)
    dn(i,0)=0.0;
    dn(i,1)=0.5*((p(i,8)-p(i,2))+(p(i,2)-p(i,5)));
    dn(i,2)=0.5*(-(p(i,1)-p(i,4))-(p(i,7)-p(i,1)));

    // dN/du2 -> dN/x(0,1)
    dn(i,3)=0.5*(-(p(i,2)-p(i,5))-(p(i,8)-p(i,2)));
    dn(i,4)=0.0;
    dn(i,5)=0.5*((p(i,6)-p(i,0))+(p(i,0)-p(i,3)));

    // dN/du3 -> dN/dx(0,2)
    dn(i,6)=0.5*((p(i,7)-p(i,1))+(p(i,1)-p(i,4)));
    dn(i,7)=0.5*(-(p(i,0)-p(i,3))-(p(i,6)-p(i,0)));
    dn(i,8)=0.0;

    // dofs of the second node
    // dN/du4 -> dN/dx(1,0)
    dn(i,9)=0.0;
    dn(i,10)=0.5*(-(p(i,8)-p(i,2)));
    dn(i,11)=0.5*((p(i,7)-p(i,1)));

    // dN/du5 -> dN/dp(i,4)
    dn(i,12)=0.5*((p(i,8)-p(i,2)));
    dn(i,13)=0.0;
    dn(i,14)=0.5*(-(p(i,6)-p(i,0)));

    // dN/du6 -> dN/dx(1,2)
    dn(i,15)=0.5*(-(p(i,7)-p(i,1)));
    dn(i,16)=0.5*((p(i,6)-p(i,0)));
    dn(i,17)=0.0;

    // dofs of the third node
    // dN/du7 -> dN/dx(2,0)
    dn(i,18)=0.0;
    dn(i,19)=0.5*(-(p(i,2)-p(i,5)));
    dn(i,20)=0.5*((p(i,1)-p(i,4)));

    // dN/du8 -> dN/dx(2,1)
    dn(i,21)=0.5*((p(i,2)-p(i,5)));
    dn(i,22)=0.0;
    dn(i,23)=0.5*(-(p(i,0)-p(i,3)));

    // dN/du9 -> dN/dx(2,2)
    dn(i,24)=0.5*(-(p(i,1)-p(i,4)));
    dn(i,25)=0.5*((p(i,0)-p(i,3)));
    dn(i,26)=0.0;
  }
};

//! Center functor
struct Centers {
  currents_type c;
  trimesh_type p;

  // Constructor takes View by "value"; this does a shallow copy.
  Centers (currents_type c_, trimesh_type p_) : c (c_),p (p_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    c(i,0)=(p(i,0)+p(i,3)+p(i,6))/3;
    c(i,1)=(p(i,1)+p(i,4)+p(i,7))/3;
    c(i,2)=(p(i,2)+p(i,5)+p(i,8))/3;
  }
};

//! Convolution of dot-product of normals with a kernel functor
template<class ViewType>
struct Convolute {
  ViewType c1;
  ViewType n1;
  currents_type c2;
  currents_type n2;
  double lscale;

  // Constructor takes View by "value"; this does a shallow copy.
  Convolute (ViewType c1_, ViewType n1_,
             currents_type c2_, currents_type n2_,
             double lscale_)
  : c1(c1_), n1(n1_) ,c2(c2_), n2(n2_), lscale(lscale_)  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const typename ViewType::size_type i,double& sum) const {
    for (currents_type::size_type j=0; j<c2.dimension_0(); j++)
    {
      sum += (n1(i,0)*n2(j,0)+n1(i,1)*n2(j,1)+n1(i,2)*n2(j,2))*
          kernel(c1(i,0),c1(i,1),c1(i,2),c2(j,0),c2(j,1),c2(j,2),lscale);
    }
  }
};

//! Update points functor
struct Update {
  trimesh_type in;
  trimesh_type add;

  // Constructor takes View by "value"; this does a shallow copy.
  Update (trimesh_type in_, trimesh_type add_) : in (in_),add (add_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    in(i,0)+=add(i,0);
    in(i,1)+=add(i,1);
    in(i,2)+=add(i,2);
    in(i,3)+=add(i,3);
    in(i,4)+=add(i,4);
    in(i,5)+=add(i,5);
    in(i,6)+=add(i,6);
    in(i,7)+=add(i,7);
    in(i,8)+=add(i,8);
  }
};

//! Init View of type trimesh_type with zeros
struct Init9{
  trimesh_type in;

  Init9(trimesh_type in_): in (in_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i) const {
    for (trimesh_type::size_type j=0; j<in.dimension_1(); j++)
      in(i,j)=0.0;
  }
};

//! Convolution of "contraction" of a derivative of normals and normals with a kernel
template <class ViewType1, class ViewType2>
struct ConvoluteDN {
  ViewType1 c1;
  ViewType2 dn1;
  currents_type c2;
  currents_type n2;
  trimesh_atomic_type g;
  int offset;
  double lscale;
  double fac;

  // Constructor takes View by "value"; this does a shallow copy.
  ConvoluteDN (ViewType1 c1_, ViewType2 dn1_,
               currents_type c2_, currents_type n2_, trimesh_atomic_type g_,
               int offset_, double lscale_,double fac_)
  : c1(c1_), dn1(dn1_), c2(c2_), n2(n2_), g(g_), offset(offset_),
    lscale(lscale_), fac (fac_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const typename ViewType1::size_type i) const {
    for (currents_type::size_type j=0; j<c2.dimension_0(); j++)
    {
      double ker = fac*kernel(c1(i,0),c1(i,1),c1(i,2),c2(j,0),c2(j,1),c2(j,2),lscale);

      int ii=i+offset; // the global row to assemble
      g(ii,0) += (dn1(i,0)*n2(j,0)+dn1(i,1)*n2(j,1)+dn1(i,2)*n2(j,2))*ker;
      g(ii,1) += (dn1(i,3)*n2(j,0)+dn1(i,4)*n2(j,1)+dn1(i,5)*n2(j,2))*ker;
      g(ii,2) += (dn1(i,6)*n2(j,0)+dn1(i,7)*n2(j,1)+dn1(i,8)*n2(j,2))*ker;
      g(ii,3) += (dn1(i,9)*n2(j,0)+dn1(i,10)*n2(j,1)+dn1(i,11)*n2(j,2))*ker;
      g(ii,4) += (dn1(i,12)*n2(j,0)+dn1(i,13)*n2(j,1)+dn1(i,14)*n2(j,2))*ker;
      g(ii,5) += (dn1(i,15)*n2(j,0)+dn1(i,16)*n2(j,1)+dn1(i,17)*n2(j,2))*ker;
      g(ii,6) += (dn1(i,18)*n2(j,0)+dn1(i,19)*n2(j,1)+dn1(i,20)*n2(j,2))*ker;
      g(ii,7) += (dn1(i,21)*n2(j,0)+dn1(i,22)*n2(j,1)+dn1(i,23)*n2(j,2))*ker;
      g(ii,8) += (dn1(i,24)*n2(j,0)+dn1(i,25)*n2(j,1)+dn1(i,26)*n2(j,2))*ker;
    }
  }
};

//! Convolution of dot-product normals with a derivative of a kernel
template <class ViewType>
struct ConvoluteDk {
  ViewType c1;
  ViewType n1;
  currents_type c2;
  currents_type n2;
  trimesh_atomic_type g;
  int of;
  double lscale;
  double fac;
  bool sym;

  // Constructor takes View by "value"; this does a shallow copy.
  ConvoluteDk (ViewType c1_, ViewType n1_,
      currents_type c2_, currents_type n2_, trimesh_atomic_type g_,
      int of_, double lscale_, double fac_, bool sym_)
  : c1(c1_), n1(n1_), c2(c2_), n2(n2_), g(g_),
    of(of_), lscale(lscale_), fac(fac_), sym(sym_) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const typename ViewType::size_type i) const {
    for(currents_type::size_type j=0;j<c2.dimension_0();j++)
    {
      double ker[3];
      double arg = fac*kernel(c1(i,0),c1(i,1),c1(i,2),c2(j,0),c2(j,1),c2(j,2),
          lscale)*(-2.0/3.0)/lscale;
      ker[0] =  arg*(c1(i,0)-c2(j,0));
      ker[1] =  arg*(c1(i,1)-c2(j,1));
      ker[2] =  arg*(c1(i,2)-c2(j,2));

      double dot= (n1(i,0)*n2(j,0)+n1(i,1)*n2(j,1)+n1(i,2)*n2(j,2));

      int ii = i+of;
      g(ii,0) += dot*ker[0];
      g(ii,1) += dot*ker[1];
      g(ii,2) += dot*ker[2];
      g(ii,3) += dot*ker[0];
      g(ii,4) += dot*ker[1];
      g(ii,5) += dot*ker[2];
      g(ii,6) += dot*ker[0];
      g(ii,7) += dot*ker[1];
      g(ii,8) += dot*ker[2];

      if (sym){
        g(j,0) -= dot*ker[0];
        g(j,1) -= dot*ker[1];
        g(j,2) -= dot*ker[2];
        g(j,3) -= dot*ker[0];
        g(j,4) -= dot*ker[1];
        g(j,5) -= dot*ker[2];
        g(j,6) -= dot*ker[0];
        g(j,7) -= dot*ker[1];
        g(j,8) -= dot*ker[2];
      }
    }
  }
};


  /*!
  \brief Comparison of 2 surfaces in the space of surface current W.

  This class computes the difference between two surfaces in the space of surface currents.
  In here surfaces are already discretized objects in the sense of Triangulations.
  These must be provided by INVANA:: Triangulation.

  Further it provides gradient information of the above mentioned norm wrt. degrees of freedom
  of the triangulation. Track of how these degrees of freedom are matched to the degrees of
  freedom of some global problem this norm is part of is kept by the INVANA::Triangulation class.

  This class makes heavy use of OpenMP(if available) through Kokkos. Otherwise one will at least
  benefit from the Kokkos:Serial mode which is anyhow way faster than plain mpi parallelization.

  */
class SurfCurrentPair
{
private:

  //! triangulation of the source surface
  Teuchos::RCP<INVANA::Triangulation> tri_source_;

  //! triangulation of the target surface
  Teuchos::RCP<INVANA::Triangulation> tri_target_;

  //! Kokkos view on the target centers
  currents_type c_target_;

  //! Kokkos view on the target normals
  currents_type n_target_;

  //! spatial scale of the kernel
  double sigmaW_;

  //! estimation of the measurement noise
  double var_estim_;

  //! want scaling of the optimization problem?
  bool scaling_;

  //! scale factor
  double scalefac_;

  //! pre computed convolution of target with itself; constant since construction
  double preconvtarget_;

  //! extract_type to h_view_type conversion
  template <typename host_data_type>
  void ExtractToHView(const extract_type& in, host_data_type& out, Teuchos::RCP<Epetra_Map> map);
  template <typename host_data_type>
  void ExtractToHView(const extract_type& in, host_data_type& out);

  //! trimesh_host_type to extract_type conversion
  template <typename host_data_type>
  void HViewToExtract(const host_data_type& in, const Epetra_Map& inmap, extract_type& out);
  template <typename host_data_type>
  void HViewToExtract(const host_data_type& in, extract_type& out);

  //! subdivide a range of indices (magnitude "size") into nrnk number of chunks
  std::pair<int, int> MyIndices(int nrnk, int myrnk, int size);

  /*! \brief Apply noise for testing
   *
   *  Apply normally distributed noise in the space of currents.
   *
   *  Dont use this function on large meshes! The purpose is set up
   *  proper measurement noise fitting in the framework to avoid the
   *  inverse crime and have a well-posed problem.
   *
   *  \author kehl \date 06/2016
   */
  void ApplyNoise(int seed);

public:

  //! constructor
  SurfCurrentPair(Teuchos::RCP<DRT::Discretization> sourcedis,
                  Teuchos::RCP<DRT::Discretization> targetdis,
                  DRT::Condition* sourcecond,
                  DRT::Condition* targetcond);

  //! standard destructor
  ~SurfCurrentPair(){}

  double GetScaleFac() {return scalefac_;}

  //! Compute Norm on current difference in W-space
  double WSpaceNorm();

  //! Compute gradient of Norm on current difference in W-space
  void GradientWSpaceNorm(Teuchos::RCP<Epetra_MultiVector> gradient);
}; // class SurfCurrentPair


/*!
\brief A triangulation object

  * A Triangulation of a surface condition allowing to access points, and data associated
  * to these points, to be used independently.
  *
  * If there will be more ways to obain triangulations, e.g. via input of STL, this
  * should inherit from a proper base class. Since the main focus - at the moment -
  * is that this class is used by the INVANA::SurfCurrrentPair class to have it easily
  * being reused for any sort of opmtimization problem.

*/
class Triangulation {

public:

  //! Constructor
  Triangulation(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<DRT::Condition> surface)
  : discret_ (discret), surface_ (surface), trimap_(Teuchos::null) {
    EvaluatePoints();
    EvaluateDofs();
  }

  //! Destructor
  ~Triangulation(){};

  //! Definition for maptype describing tri to global dof layout
  typedef std::map<int, std::vector<int> > facetdofstype;

  //! @name Global Data Extraction
  //@{

  //! Get points reduced to mpi-rank rank
  extract_type Points(int rank);

  //! Get Points reduced to all ranks
  extract_type Points();

  //! Get Data reduced to mpi-rank rank
  extract_type Data(int rank);

  //! Get Data reduced to all mpi-ranks
  extract_type Data();
  //@}

  //! @name Local Data Extraction
  //@{
  //! get Points from myrank only
  extract_type MyPoints();

  //! Get Data from myrank only
  extract_type MyData();
  //@}

  //! @name Query
  //@{

  //! global number of tris in this triangulation
  int NumTris();
  //@}

  //! Evaluate the warp of the triangulation according to a state being set in the discretization
  void EvaluateWarp();

  //! Set data back to global system vector
  void SetDataGlobally(const extract_type& data,const Epetra_Map&, Teuchos::RCP<Epetra_MultiVector> vector);

  //! The Triangulations Communicator to query
  const Teuchos::RCP<Epetra_Comm> Comm();

  //! Apply noise to the triangulation
  void ApplyNoise(const extract_type& normals, const extract_type& centers,
      extract_type& blurred, double variance, double lengthscale, int seed) const;

private:

  //! @name Basis for this triangulation
  //{@

  //! The underlying discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! The surface condition which is used for this triangulation
  Teuchos::RCP<DRT::Condition> surface_;
  //@}

  //! @name Triangulation data
  //{@

  //! Points making up the tris
  INVANA::extract_type points_;

  //! Data associated to the tris
  INVANA::extract_type data_;

  //! Map of GIDs for each tri (inherited from the associted surface elements)
  Teuchos::RCP<Epetra_Map> trimap_;

  //! Dofs associated to the tris
  facetdofstype facetmap_;
  //@}

  //! @name Initialization bools
  //{@
  bool haspoints_=false;
  bool hasdofs_=false;
  bool hasdata_=false;
  //@}

  //! @name Evaluation Methods
  //{@

  //! Evaluate the Points according to the surface condition
  void EvaluatePoints();
  //! Evalaute the associated global Dofs
  void EvaluateDofs();
  //!@}

  //! @name Communicate distributed data around
  //{@

  //! Bring data to mpi-rank rank
  void CommunicateData(extract_type& data, int rank);

  //! Bring data to all mpi-ranks
  void CommunicateData(extract_type& data);
  //@}
};

}  // namespace INVANA

#endif // OBJECTIVE_FUNCT_SURFCURR_H
#endif
