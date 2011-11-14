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


#ifdef CCADISCRET
#ifndef TIMEINT_STD_SEMILAGRANGE_H_
#define TIMEINT_STD_SEMILAGRANGE_H_


#include "timeInt_std_extrapolation.H"


namespace XFEM
{

  /*!
  \brief this class is used to compute startvalues for critical interface-near nodes

  */
  class SemiLagrange : public STD
  {
  public:

    //! constructor
    explicit SemiLagrange(
    	XFEM::TIMEINT& timeInt,
    	INPAR::COMBUST::XFEMTimeIntegration& timeIntType,
        const RCP<Epetra_Vector> veln,
        const RCP<COMBUST::FlameFront> flamefront,
        const double& dt,
        const double& theta,
        const double& flamespeed,
        bool initialize
    );

    //! destructor
    ~SemiLagrange()
    {
      return;
    }

    //! get startvalues in time step for nodes which changed interface-side
    void compute(
        vector<RCP<Epetra_Vector> > newColVectorsn,
    	vector<RCP<Epetra_Vector> > newRowVectorsnp
    );

  private:

    //! disabled copy constructor
    explicit SemiLagrange(const XFEM::SemiLagrange& slvalues);

    //! disabled default constructor
    explicit SemiLagrange();

    //! run a Newton loop in order to compute the exact lagrangian
    // origin for a node which changed interface side
    void NewtonLoop(
        DRT::Element*& ele,
        TimeIntData* data,
        LINALG::Matrix<3,1>& xi,
        LINALG::Matrix<3,1>& vel,
        double& phi,
        bool& eleFound,
        bool& stdBackTracking
    );

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

    //! determine velocity and pressure for nodes where the
    // "normal" semi-lagrange startfinder failed
    void getDataForNotConvergedNodes(
    );

    void callBackTracking(
        DRT::Element*& ele,
        TimeIntData* data,
        LINALG::Matrix<3,1>& xi,
        const char* backTrackingType
    );

    //! track back the lagrangian origin to get final startvalues
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

    void reinitializeData();

    template<const int numnode, DRT::Element::DiscretizationType DISTYPE>
    void computeNodalGradient(
        vector<RCP<Epetra_Vector> >& newColVectors,
        const Epetra_Map& newdofcolmap,
        map<XFEM::DofKey<XFEM::onNode>,XFEM::DofGID>& newNodalDofColDistrib,
        DRT::Element*& ele,
        LINALG::Matrix<3,1> coords,
        vector<LINALG::Matrix<3,3> >& velnpDeriv1,
        vector<LINALG::Matrix<1,3> >& presnpDeriv1
    ) const;

    double Theta(TimeIntData* data) const;

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

    // data
    const RCP<COMBUST::FlameFront> flamefront_;

    // constants
    const double theta_default_; // factor of one-step theta scheme
    double flamespeed_;

  }; // class Startvalues
} // namespace XFEM
#endif /*TIMEINT_SEMILAGRANGE_H_*/
#endif  // CCADISCRET
