/*!-----------------------------------------------------------------------------------------------*
\file xfluid_timeinnt_std_SemiLagrange.H

\brief provides the SemiLagrangean class

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241

</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef XFLUID_TIMEINT_STD_SEMILAGRANGE_H_
#define XFLUID_TIMEINT_STD_SEMILAGRANGE_H_


#include "xfluid_timeInt_base.H"


namespace DRT
{
  class Element;
}


namespace XFEM
{

  class XFLUID_TIMEINT_BASE;
  class TimeIntData;


/*!
  \brief this class is used in XFEM to compute new values for standard degrees of freedom
  for nodes which change their interface side between two consecutive time steps.
  It bases on the Semi-Lagrangian approach described in
  "/intern/arbeiten/diplomarbeiten/WinklmaierMartin2010.pdf"
 */
class XFLUID_SemiLagrange : public XFLUID_STD
{
public:

  //! constructor
  explicit XFLUID_SemiLagrange(
      XFEM::XFLUID_TIMEINT_BASE&                                     timeInt,          /// time integration base class object
      const std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method,  /// reconstruction map for nodes and its dofsets
      INPAR::XFEM::XFluidTimeInt&                                    timeIntType,      /// type of time integration
      const RCP<Epetra_Vector>                                       veln,             /// velocity at time t^n
      const double&                                                  dt,               /// time step size
      const double&                                                  theta,            /// OST theta
      bool                                                           initialize        /// is initialization?

  );

  //! destructor
  ~XFLUID_SemiLagrange()
  {return;}

  //! get startvalues in time step for nodes which changed interface-side
  void compute(
      vector<RCP<Epetra_Vector> >& newRowVectorsn
  );

private:

  //! disabled copy constructor
  explicit XFLUID_SemiLagrange(const XFEM::XFLUID_SemiLagrange& slvalues);

  //! disabled default constructor
  explicit XFLUID_SemiLagrange();

  //! determine point's dofset in element ele w.r.t old or new interface position
  void getNodalDofSet(
      DRT::Element*           ele,    /// pointer to element
      LINALG::Matrix<3,1>&    x ,     /// global coordinates of point
      std::vector<int>&       nds,    /// determine the points dofset w.r.t old/new interface position
      bool                    step_np /// computation w.r.t old or new interface position?
      );

  //! run a Newton loop in order to compute the exact Lagrangian origin for a node which changed interface side
  void NewtonLoop(
      DRT::Element*&          ele,              /// pointer to element
      TimeIntData*            data,             /// current data
      LINALG::Matrix<3,1>&    xi,               /// local coordinates of point
      LINALG::Matrix<3,1>&    vel,              /// velocity at current point
      bool&                   elefound          /// is element found ?
  );

  //! perform one Newton iteration in order to compute the exact Lagrangian origin for a node which changed its interface side
  template<const int numnode, DRT::Element::DiscretizationType DISTYPE>
  void NewtonIter(
      DRT::Element*&          ele,              /// pointer to element to be updated
      TimeIntData*            data,             /// current data to be updated
      LINALG::Matrix<3,1>&    xi,               /// local coordinates w.r.t ele to be updated
      LINALG::Matrix<3,1>&    residuum,         /// residual for semilagrangean backtracking to be updated
      LINALG::Matrix<3,1>&    incr,             /// computed increment for lagrangean origin to be updated
      bool&                   elefound          /// element found ?
  );

  //! determine velocity and pressure for nodes where the "normal" semi-lagrange startfinder failed
  void getDataForNotConvergedNodes();

  //! call the back tracking which computes the final values
  void callBackTracking(
      DRT::Element*&        ele,                 /// pointer to element
      TimeIntData*          data,                /// data
      LINALG::Matrix<3,1>&  xi,                  /// local coordinates
      const char*           backTrackingType     /// type of backTracking
  );

  template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
  void extractNodalValuesFromVector(
      LINALG::Matrix<3,numnode>& evel,
      LINALG::Matrix<numnode,1>& epre,
      RCP<Epetra_Vector> vel_vec,
      std::vector<int>& lm);

  //! track back the Lagrangian origin to get final values
  template<const int numnode, DRT::Element::DiscretizationType DISTYPE>
  void backTracking(
      DRT::Element*&        fittingele,          /// pointer to element
      TimeIntData*          data,                /// data
      LINALG::Matrix<3,1>&  xi,                  /// local coordinates
      const char*           backTrackingType     /// type of backTrackingwVectors
  );

  void newIteration_prepare(
      vector<RCP<Epetra_Vector> > newRowVectors
      );

  //! gradients at a node
  void newIteration_nodalData(
      vector<RCP<Epetra_Vector> > newColVectors
  );

  //! reinitialize some data for new computations, f.e. at a new FGI
  void reinitializeData();

  //! compute the nodal gradient
  template<const int numnode, DRT::Element::DiscretizationType DISTYPE>
  void computeNodalGradient(
      vector<RCP<Epetra_Vector> >& newColVectors,
      const Epetra_Map& newdofcolmap,
      map<XFEM::DofKey,XFEM::DofGID>& newNodalDofColDistrib,
      const DRT::Element* ele,
      DRT::Node* node,
      vector<LINALG::Matrix<3,3> >& velnpDeriv1,
      vector<LINALG::Matrix<1,3> >& presnpDeriv1
  ) const;

  //! compute the theta which has to be used for computation
  double Theta(TimeIntData* data) const;

  //! check if newton iteration has finished
  bool globalNewtonFinished(
      int counter = 0
  ) const;

# ifdef PARALLEL
  //! export data to startpoint processor when Semi-Lagrange algorithm failed
  void exportAlternativAlgoData(
  );

  //! export data to neighbour proc in Newton loop
  void exportIterData(
      bool& procfinished
  );
# endif

  // constants
  const double theta_default_; // factor of one-step theta scheme
//  double flamespeed_; // absolute flame velocity without fluid influence
}; // class Startvalues
} // namespace XFEM
#endif /*XFLUID_TIMEINT_SEMILAGRANGE_H_*/
