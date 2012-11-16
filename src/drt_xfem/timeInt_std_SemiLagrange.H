/*!-----------------------------------------------------------------------------------------------*
\file timeInt_std_SemiLagrange.H

\brief provides the SemiLagrange class

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef TIMEINT_STD_SEMILAGRANGE_H_
#define TIMEINT_STD_SEMILAGRANGE_H_


#include "timeInt.H"


namespace XFEM
{
/*!
  \brief this class is used in XFEM to compute new values for standard degrees of freedom
  for nodes which change their interface side between two consecutive time steps.
  It bases on the Semi-Lagrangian approach described in
  "/intern/arbeiten/diplomarbeiten/WinklmaierMartin2010.pdf"
 */
class SemiLagrange : public STD
{
public:

  //! constructor
  explicit SemiLagrange(
      XFEM::TIMEINT& timeInt,
      INPAR::COMBUST::XFEMTimeIntegration& timeIntType,
      const RCP<Epetra_Vector> veln,
      const double& dt,
      const double& theta,
      const RCP<COMBUST::FlameFront> flamefront,
      const double& veljump,
      bool initialize
  );

  //! destructor
  ~SemiLagrange()
  {return;}

  //! get startvalues in time step for nodes which changed interface-side
  void compute(
      vector<RCP<Epetra_Vector> > newRowVectorsn,
      vector<RCP<Epetra_Vector> > newRowVectorsnp
  );

private:

  //! disabled copy constructor
  explicit SemiLagrange(const XFEM::SemiLagrange& slvalues);

  //! disabled default constructor
  explicit SemiLagrange();

  //! run a Newton loop in order to compute the exact Lagrangian origin for a node which changed interface side
  void NewtonLoop(
      DRT::Element*& ele,
      TimeIntData* data,
      LINALG::Matrix<3,1>& xi,
      LINALG::Matrix<3,1>& vel,
      double& phi,
      bool& eleFound,
      bool& stdBackTracking
  );

  //! perform one Newton iteration in order to compute the exact Lagrangian origin for a node which changed its interface side
  template<const int numnode, DRT::Element::DiscretizationType DISTYPE>
  void NewtonIter(
      DRT::Element*& ele,
      TimeIntData* data,
      LINALG::Matrix<3,1>& xi,
      LINALG::Matrix<3,1>& vel,
      LINALG::Matrix<3,1>& residuum,
      LINALG::Matrix<3,1>& incr,
      double& phi,
      bool& elefound
  );

  //! determine velocity and pressure for nodes where the "normal" semi-lagrange startfinder failed
  void getDataForNotConvergedNodes(
  );

  //! call the back tracking which computes the final values
  void callBackTracking(
      DRT::Element*& ele,
      TimeIntData* data,
      LINALG::Matrix<3,1>& xi,
      const char* backTrackingType
  );

  //! track back the Lagrangian origin to get final values
  template<const int numnode, DRT::Element::DiscretizationType DISTYPE>
  void backTracking(
      DRT::Element*& fittingele,
      TimeIntData* data,
      LINALG::Matrix<3,1>& xi,
      const char* backTrackingType
  );

  //! preparation of new iteration
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
      std::map<XFEM::DofKey,XFEM::DofGID>& newNodalDofColDistrib,
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
  void exportDataToStartpointProc(
  );

  //! export data to neighbour proc in Newton loop
  void exportIterData(
      bool& procfinished
  );
# endif

  // constants
  const double theta_default_; // factor of one-step theta scheme
  const double veljump_; // absolute flame velocity without fluid influence
}; // class Startvalues
} // namespace XFEM
#endif /*TIMEINT_SEMILAGRANGE_H_*/
