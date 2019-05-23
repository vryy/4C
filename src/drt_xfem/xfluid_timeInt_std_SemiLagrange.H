/*----------------------------------------------------------------------*/
/*!
\file xfluid_timeInt_std_SemiLagrange.H

\brief provides the SemiLagrangean class

\level 3

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

*/
/*----------------------------------------------------------------------*/


#ifndef XFLUID_TIMEINT_STD_SEMILAGRANGE_H_
#define XFLUID_TIMEINT_STD_SEMILAGRANGE_H_

#include "xfluid_timeInt_base.H"


namespace DRT
{
  class Element;
}


namespace XFEM
{
  /*========================================================================*/
  // forward declarations
  /*========================================================================*/

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
    /*========================================================================*/
    //! constructor/destructor
    /*========================================================================*/

    //! constructor
    explicit XFLUID_SemiLagrange(
        XFEM::XFLUID_TIMEINT_BASE& timeInt,  ///< time integration base class object
        const std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt>>&
            reconstr_method,                      ///< reconstruction map for nodes and its dofsets
        INPAR::XFEM::XFluidTimeInt& timeIntType,  ///< type of time integration
        const Teuchos::RCP<Epetra_Vector> veln,   ///< velocity at time t^n
        const double& dt,                         ///< time step size
        const double& theta,                      ///< OST theta
        bool initialize                           ///< is initialization?

    );

    //! destructor
    ~XFLUID_SemiLagrange() { return; }

    /*========================================================================*/
    //! compute routines
    /*========================================================================*/

    //! get startvalues in time step for nodes which changed interface-side
    void compute(std::vector<Teuchos::RCP<Epetra_Vector>>& newRowVectorsn);

   private:
    /*========================================================================*/
    //! control routines for Semi-Lagrangean Newton-loop
    /*========================================================================*/

    //! run a Newton loop in order to compute the exact Lagrangian origin for a node which changed
    //! interface side
    void NewtonLoop(DRT::Element*& ele,  ///< pointer to element
        TimeIntData* data,               ///< current data
        LINALG::Matrix<3, 1>& xi,        ///< local coordinates of point
        LINALG::Matrix<3, 1>& vel,       ///< velocity at current point
        bool& elefound                   ///< is element found ?
    );

    //! perform one Newton iteration in order to compute the exact Lagrangian origin for a node
    //! which changed its interface side
    void NewtonIter(DRT::Element*& ele,  ///< pointer to element to be updated
        TimeIntData* data,               ///< current data to be updated
        LINALG::Matrix<3, 1>& xi,        ///< local coordinates w.r.t ele to be updated
        LINALG::Matrix<3, 1>& residuum,  ///< residual for semilagrangean backtracking to be updated
        LINALG::Matrix<3, 1>& incr,      ///< computed increment for lagrangean origin to be updated
        bool& elefound                   ///< element found ?
    );

    //! check if newton iteration has finished
    bool globalNewtonFinished(int counter = 0) const;

    //! Decide how or if to continue when the startpoint approximation changed the side
    bool continueForChangingSide(TimeIntData* data,  ///< current data to be updated
        DRT::Element* ele,          ///< pointer to element the current point lies in
        std::vector<int>& nds_curr  ///< nds-vector of current volumecell the current startpoint
                                    ///< approximation lies in
    );

    //! determine velocity and pressure for nodes where the "normal" semi-lagrange startfinder
    //! failed
    void getDataForNotConvergedNodes();

    //! prepare new iteration
    void newIteration_prepare(std::vector<Teuchos::RCP<Epetra_Vector>> newRowVectors);

    //! gradients at a node
    void newIteration_nodalData(std::vector<Teuchos::RCP<Epetra_Vector>> newColVectors);

    //! reinitialize some data for new computations, e.g. at a new FGI
    void reinitializeData();


    /*========================================================================*/
    //! Semi-Lagrangean backtracking routines
    /*========================================================================*/

    //! call the back tracking which computes the final values
    void callBackTracking(DRT::Element*& ele,  ///< pointer to element
        TimeIntData* data,                     ///< data
        LINALG::Matrix<3, 1>& xi,              ///< local coordinates
        const char* backTrackingType           ///< type of backTracking
    );

    //! track back the Lagrangian origin to get final values
    template <const int numnode, DRT::Element::DiscretizationType DISTYPE>
    void backTracking(DRT::Element*& fittingele,  ///< pointer to element
        TimeIntData* data,                        ///< data
        LINALG::Matrix<3, 1>& xi,                 ///< local coordinates
        const char* backTrackingType              ///< type of backTrackingwVectors
    );

    /*========================================================================*/
    //! element/dofset based routines
    /*========================================================================*/

    //! determine point's dofset in element ele w.r.t old or new interface position
    void getNodalDofSet(DRT::Element* ele,  ///< pointer to element
        LINALG::Matrix<3, 1>& x,            ///< global coordinates of point
        std::vector<int>& nds,  ///< determine the points dofset w.r.t old/new interface position
        GEO::CUT::VolumeCell*& vc,  ///< valid fluid volumecell the point x lies in
        bool step_np                ///< computation w.r.t old or new interface position?
    );

    //! compute the nodal gradient
    void computeNodalGradient(
        const std::vector<Teuchos::RCP<Epetra_Vector>>&
            colVectors,   ///< all vectors for that we reconstruct the their gradients
        DRT::Node* node,  ///< node at which we reconstruct the gradients
        std::vector<DRT::Element*>& eles,  ///< elements around node used for the reconstruction
        std::vector<std::vector<int>>& ele_nds,  ///< corresonding elements nodal dofset information
        XFEM::XFEMDofSet& dofset,                ///< XFEM dofset
        std::vector<LINALG::Matrix<3, 3>>&
            velDeriv_avg,  ///< velocity/acc component derivatives for several vectors
        std::vector<LINALG::Matrix<1, 3>>&
            preDeriv_avg  ///< pressure-component derivatives for several vectors
        ) const;


    /*========================================================================*/
    //! others
    /*========================================================================*/

    //! compute the theta which has to be used for computation
    double Theta(TimeIntData* data) const;



    /*========================================================================*/
    //! parallel routines
    /*========================================================================*/

    //! export data to startpoint processor when Semi-Lagrange algorithm failed
    void exportAlternativAlgoData();

    //! export data to neighbour proc in Newton loop
    void exportIterData(bool& procfinished);


    /*========================================================================*/
    //! constants
    /*========================================================================*/

    const double theta_default_;  //! factor of one-step theta scheme

    const double relTolIncr_;  //! tolerance for the increment
    const double relTolRes_;   //! tolerance for the residual

  };  // class XFLUID_SemiLagrange
}  // namespace XFEM


#endif /*XFLUID_TIMEINT_SEMILAGRANGE_H_*/
