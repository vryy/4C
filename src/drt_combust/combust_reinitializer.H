/*!-----------------------------------------------------------------------------------------------*
\file combust_reinitializer.H

  \brief reinitialization of G-function (level set function)

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef COMBUST_REINITIALIZER_H
#define COMBUST_REINITIALIZER_H


#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/geo_utils.H"
#include "../drt_geometry/integrationcell.H"
#include "../drt_inpar/inpar_combust.H"


namespace SCATRA
{
  class ScaTraTimIntImpl;
}

namespace COMBUST
{
  /// COMBUST Reinitializer class
  /*!
    \brief not documented

    \detailed not documented

    \author henke
    \date 12/09
   */
  class Reinitializer
  {
  public:

    /// constructor (overloaded for restart with level-set field from input file)
    explicit Reinitializer(
        const Teuchos::ParameterList&              combustdyn,
        SCATRA::ScaTraTimIntImpl&                  scatra,
        const std::map<int,GEO::BoundaryIntCells>& boundaryintcells,
        Teuchos::RCP<Epetra_Vector>                phivector,
        const bool                                 compdist=true  // flag for initial computation of signed distance function
    );

    /// destructor
    ~Reinitializer();

  private:

    /// don't allow default constructor
    explicit Reinitializer();

    /// don't allow copy constructor
    explicit Reinitializer(const Reinitializer&);

    /*!
      This member function constructs a G-function field that meets the signed distance property.
      The algorithm assigns the value of the distance between each node and the surface defined by
      G=0 as a scalar value to every node in the G-function discretization.
     */
    /// build signed distance function
    void SignedDistanceFunction(Teuchos::RCP<Epetra_Vector> phivector);

    /*!
      This member function does the same as SignedDistanceFunction, only faster.
     */
    /// build signed distance function faster
    void FastSignedDistanceFunction(Teuchos::RCP<Epetra_Vector> phivector);

    /// compute normal vector of flame front patch
    void ComputeNormalVectorToFlameFront(
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        LINALG::Matrix<3,1>&             normal
     );

    /// compute distance to edge of patch
    void ComputeDistanceToEdge(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        double&                          edgedist
    );


    /// compute distance to vertex of patch
    void ComputeDistanceToPatch(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        double&                          vertexdist
    );

    /// find a facing flame front patch by projecton of node on boundary cell edges
    void FindFacingPatchProjCellEdges(
        const LINALG::Matrix<3,1>&       node,
        const LINALG::SerialDenseMatrix& patch,
        bool&                            facenode,
        double&                          vertexdist
    );

    /// find a facing flame front patch by projecton of node into boundary cell space
    void FindFacingPatchProjCellSpace(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        const LINALG::Matrix<3,1>&       normal,
        bool&                            facenode,
        double&                          patchdist
    );

    /// project node into the boundary cell space (2D)
    template<DRT::Element::DiscretizationType DISTYPE>
    bool ProjectNodeOnPatch(
        const LINALG::Matrix<3,1>&       node,
        const GEO::BoundaryIntCell&      patch,
        const LINALG::SerialDenseMatrix& patchcoord,
        const LINALG::Matrix<3,1>&       normal,
        LINALG::Matrix<2,1>&             eta,
        double&                          alpha)
    {
      // indicator for convergence of Newton-Raphson scheme
      bool converged = false;
      // number space dimensions for 3d combustion problems
      const size_t nsd = 3;
      // here, a triangular boundary integration cell is assumed (numvertices = 3)
      const size_t numvertices = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

      // get coordinates of vertices of flame front patch
      // remark: here we only get a view (bool true) on the SerialDenseMatrix returned by CellNodalPosXYZ()
      LINALG::Matrix<nsd,numvertices> patchcoordfix(patchcoord.A(),true);

      static LINALG::Matrix<numvertices,1> funct(true);
      static LINALG::Matrix<2,numvertices> deriv(true);
      static LINALG::Matrix<nsd,1> projX(true);
      static LINALG::Matrix<nsd,2> gradprojX(true);

      //----------------------------------
      // start values for iterative scheme
      //----------------------------------
      // start position (barycenter of triangular boundary cell)
      eta(0) = 1.0/3.0;
      eta(1) = 1.0/3.0;
      // auxiliary variable
      // remark: third unknown to close system of equations; arbitrary value
      alpha = 0.0;

      // function F (system of equations)
      static LINALG::Matrix<nsd,1> f(true);
      // gradient of function F (dF/deta(0), dF/deta(1), dF/dalpha)
      static LINALG::Matrix<nsd,nsd> gradf(true);
      // increment in Newton iteration (unknown to be solved for)
      static LINALG::Matrix<nsd,1> incr(true);

      // maximum number Newton iterations
      size_t maxiter = 3;
      // convergence tolerance
      double conv = 0.0;

      //------------------------------------------------------
      // Newton-Raphson loop for non-linear projection problem
      //------------------------------------------------------
      for (size_t iter=0;iter<maxiter;++iter)
      {
        // evaluate shape functions in boundary cell space at current position \eta_1,\eta_2 on the patch
        funct.Clear();
        DRT::UTILS::shape_function_2D(funct,eta(0),eta(1),patch.Shape());
        // evaluate derivatives of shape functions in boundary cell space at current position \eta_1,\eta_2 on the patch
        deriv.Clear();
        DRT::UTILS::shape_function_2D_deriv1(deriv,eta(0),eta(1),patch.Shape());

        // evaluate projection X of node P at current position \eta_1,\eta_2 on the patch
        // projX(i,j) = patchcoord(i,k)*funct(k,1)
        projX.Clear();
        projX.MultiplyNN(patchcoordfix,funct);

        // evaluate gradient of projection X of node P at current position \eta_1,\eta_2 on the patch
        // gradprojX(i,j) = patchcoord(i,k)*deriv(j,k)
        gradprojX.Clear();
        gradprojX.MultiplyNT(patchcoordfix,deriv);

        //---------------------------------------------------
        // build system of equations F and its gradient gradF
        //---------------------------------------------------
        // TODO documentaton missing
        f.Clear();
        gradf.Clear();
        incr.Clear();
        for (size_t icoord=0;icoord<nsd;++icoord)
        {
          // evaluate function f
          f(icoord) = projX(icoord) + alpha*normal(icoord) - node(icoord);
          // evaluate gradient of function at current position on patch in boundary cell space
          gradf(icoord,0) = gradprojX(icoord,0);
          gradf(icoord,1) = gradprojX(icoord,1);
          gradf(icoord,2) = normal(icoord);
        }

        // check convergence
        conv = sqrt(f(0)*f(0)+f(1)*f(1)+f(2)*f(2));
        //cout << "iteration " << iter << ": -> |f|=" << conv << endl;
        if (conv <= 1.0E-12) break;

        //----------------------------------------------------
        // solve linear system of equations: gradF * incr = -F
        //----------------------------------------------------
        // F = F*-1.0
        f.Scale(-1.0);
        // solve A.X=B
        LINALG::FixedSizeSerialDenseSolver<nsd,nsd,1> solver;
        solver.SetMatrix(gradf);              // set A=gradF
        solver.SetVectors(incr, f);           // set X=incr, B=F
        solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
        int err2 = solver.Factor();           // ?
        int err = solver.Solve();             // incr = gradF^-1.F
        if ((err != 0) || (err2!=0))
          dserror("solving linear system in Newton-Raphson method for projection failed");

        // update eta and alpha
        eta(0) += incr(0);
        eta(1) += incr(1);
        alpha  += incr(2);
        //cout << "solution vector: component 1: " << eta(0) << " component 2: " << eta(1) << " alpha: " << alpha << endl;
      }
      // change sign to preserve sign of G-function
      alpha = -alpha;

      // Newton iteration unconverged
      if (conv > 1.0E-12)
      {
        alpha = 7777.7;
//        cout << "projection did not converge" << endl;
        //dserror("projection did not converge!");
      }
      else
      {
        converged = true;
        //cout << "convergence criterion " << conv << endl;
        //cout << "solution vector: component 1: " << eta(0) << " component 2: " << eta(1) << " alpha: " << alpha << endl;
      }

      return converged;
    }

    /// compares the second entry of a pair<int,double>. To be passed to the sorting algo
    static bool MyComparePairs (const std::pair<int,double>& first, const std::pair<int,double>& second)
    {
      if (fabs(first.second) < fabs(second.second))
        return true;
      else
        return false;
    }

    /// parameter list with combustion input parameters
    const Teuchos::ParameterList& combustdyn_;
    /// option for reinitialization of G-function;
    const INPAR::COMBUST::ReInitialActionGfunc reinitaction_;
    /// switch for reinitialization only within a band around the interface
    const bool reinitband_;
    /// band width for reinitialization (maximum G-function value)
    const double reinitbandwidth_;

    /// access to the scalar transport field time integration scheme
    SCATRA::ScaTraTimIntImpl& scatra_;

    /// boundary integration cells (flame front patches)
    const std::map<int,GEO::BoundaryIntCells>& flamefront_;

  }; // class Reinitializer
} // namespace COMBUST

#endif // COMBUST_REINITIALIZER_H
