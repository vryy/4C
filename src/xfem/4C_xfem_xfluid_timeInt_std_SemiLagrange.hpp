/*----------------------------------------------------------------------*/
/*! \file

\brief provides the SemiLagrangean class

\level 3


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_XFEM_XFLUID_TIMEINT_STD_SEMILAGRANGE_HPP
#define FOUR_C_XFEM_XFLUID_TIMEINT_STD_SEMILAGRANGE_HPP

#include "4C_config.hpp"

#include "4C_xfem_xfluid_timeInt_base.hpp"

FOUR_C_NAMESPACE_OPEN


namespace DRT
{
  class Element;
}


namespace XFEM
{
  /*========================================================================*/
  // forward declarations
  /*========================================================================*/

  class XfluidTimeintBase;
  class TimeIntData;


  /*!
    \brief this class is used in XFEM to compute new values for standard degrees of freedom
    for nodes which change their interface side between two consecutive time steps.
    It bases on the Semi-Lagrangian approach described in
    "/intern/arbeiten/diplomarbeiten/WinklmaierMartin2010.pdf"
   */
  class XfluidSemiLagrange : public XfluidStd
  {
   public:
    /*========================================================================*/
    //! constructor/destructor
    /*========================================================================*/

    //! constructor
    explicit XfluidSemiLagrange(
        XFEM::XfluidTimeintBase& timeInt,  ///< time integration base class object
        const std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt>>&
            reconstr_method,                      ///< reconstruction map for nodes and its dofsets
        INPAR::XFEM::XFluidTimeInt& timeIntType,  ///< type of time integration
        const Teuchos::RCP<Epetra_Vector> veln,   ///< velocity at time t^n
        const double& dt,                         ///< time step size
        const double& theta,                      ///< OST theta
        bool initialize                           ///< is initialization?

    );

    /*========================================================================*/
    //! compute routines
    /*========================================================================*/

    //! get startvalues in time step for nodes which changed interface-side
    void compute(std::vector<Teuchos::RCP<Epetra_Vector>>& newRowVectorsn) override;

   private:
    /*========================================================================*/
    //! control routines for Semi-Lagrangean Newton-loop
    /*========================================================================*/

    //! run a Newton loop in order to compute the exact Lagrangian origin for a node which changed
    //! interface side
    void NewtonLoop(DRT::Element*& ele,   ///< pointer to element
        TimeIntData* data,                ///< current data
        CORE::LINALG::Matrix<3, 1>& xi,   ///< local coordinates of point
        CORE::LINALG::Matrix<3, 1>& vel,  ///< velocity at current point
        bool& elefound                    ///< is element found ?
    );

    //! perform one Newton iteration in order to compute the exact Lagrangian origin for a node
    //! which changed its interface side
    void NewtonIter(DRT::Element*& ele,  ///< pointer to element to be updated
        TimeIntData* data,               ///< current data to be updated
        CORE::LINALG::Matrix<3, 1>& xi,  ///< local coordinates w.r.t ele to be updated
        CORE::LINALG::Matrix<3, 1>&
            residuum,  ///< residual for semilagrangean backtracking to be updated
        CORE::LINALG::Matrix<3, 1>&
            incr,       ///< computed increment for lagrangean origin to be updated
        bool& elefound  ///< element found ?
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
    void newIteration_nodalData(std::vector<Teuchos::RCP<Epetra_Vector>> newRowVectors);

    //! reinitialize some data for new computations, e.g. at a new FGI
    void reinitializeData();


    /*========================================================================*/
    //! Semi-Lagrangean backtracking routines
    /*========================================================================*/

    //! call the back tracking which computes the final values
    void callBackTracking(DRT::Element*& ele,  ///< pointer to element
        TimeIntData* data,                     ///< data
        CORE::LINALG::Matrix<3, 1>& xi,        ///< local coordinates
        const char* backTrackingType           ///< type of backTracking
    );

    //! track back the Lagrangian origin to get final values
    template <const int numnode, CORE::FE::CellType DISTYPE>
    void backTracking(DRT::Element*& fittingele,  ///< pointer to element
        TimeIntData* data,                        ///< data
        CORE::LINALG::Matrix<3, 1>& xi,           ///< local coordinates
        const char* backTrackingType              ///< type of backTrackingwVectors
    );

    /*========================================================================*/
    //! element/dofset based routines
    /*========================================================================*/

    //! determine point's dofset in element ele w.r.t old or new interface position
    void getNodalDofSet(DRT::Element* ele,  ///< pointer to element
        CORE::LINALG::Matrix<3, 1>& x,      ///< global coordinates of point
        std::vector<int>& nds,  ///< determine the points dofset w.r.t old/new interface position
        CORE::GEO::CUT::VolumeCell*& vc,  ///< valid fluid volumecell the point x lies in
        bool step_np                      ///< computation w.r.t old or new interface position?
    );

    //! compute the nodal gradient
    void computeNodalGradient(
        const std::vector<Teuchos::RCP<Epetra_Vector>>&
            colVectors,   ///< all vectors for that we reconstruct the their gradients
        DRT::Node* node,  ///< node at which we reconstruct the gradients
        std::vector<DRT::Element*>& eles,  ///< elements around node used for the reconstruction
        std::vector<std::vector<int>>& ele_nds,  ///< corresonding elements nodal dofset information
        XFEM::XFEMDofSet& dofset,                ///< XFEM dofset
        std::vector<CORE::LINALG::Matrix<3, 3>>&
            velDeriv_avg,  ///< velocity/acc component derivatives for several vectors
        std::vector<CORE::LINALG::Matrix<1, 3>>&
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
    void exportIterData(bool& procDone);


    /*========================================================================*/
    //! constants
    /*========================================================================*/

    const double theta_default_;  //! factor of one-step theta scheme

    const double rel_tol_incr_;  //! tolerance for the increment
    const double rel_tol_res_;   //! tolerance for the residual

  };  // class XFLUID_SemiLagrange
}  // namespace XFEM


FOUR_C_NAMESPACE_CLOSE

#endif
