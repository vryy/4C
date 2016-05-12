/*!--------------------------------------------------------------------------
\file drt_utils_maxent_basisfunctions.H

\brief Specification of max entropy basis functions

<pre>
\level 2
\maintainer Keijo Nissen
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/

#ifndef DRT_UTILS_MAXENT_BASISFUNCTION_H
#define DRT_UTILS_MAXENT_BASISFUNCTION_H

#include <vector>
#include "../linalg/linalg_fixedsizematrix.H"
#include <Teuchos_ParameterList.hpp>

/*==========================================================================*/
// forward declarations
/*==========================================================================*/
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}

namespace Teuchos
{
  class ParameterList;
}

/*==========================================================================*/

namespace DRT
{
namespace MESHFREE
{

/*----------------------------------------------------------------------------*/
/*!
 * \brief abstract meshfree problem class
 *///                                                                nis Mar12
/*----------------------------------------------------------------------------*/
class MeshfreeApprox
{

public:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Returns range of prior function
   *///                                                              nis Mar12
  /*--------------------------------------------------------------------------*/
  double GetRange(){return range_;};

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Evaluating meshfree basis functions and derivatives
   *        Bubnov-Galerkin-wrapper
   *///                                                           nis Mar12
  /*--------------------------------------------------------------------------*/
  virtual int GetMeshfreeBasisFunction(
    const int                         dim,                         //!< (i) dimension
    Teuchos::RCP<LINALG::SerialDenseMatrix>  diffx,                //!< (i) distance vector between node and integration point
    Teuchos::RCP<LINALG::SerialDenseVector>& funct                 //!< (o) basis functions values
  ) const
  {
    Teuchos::RCP<LINALG::SerialDenseMatrix> deriv = Teuchos::null;
    return GetMeshfreeBasisFunction(dim,diffx,funct,deriv);
  };

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Evaluating meshfree basis functions and derivatives
   *///                                                           nis Mar12
  /*--------------------------------------------------------------------------*/
  virtual int GetMeshfreeBasisFunction(
    const int                         dim,              //!< (i) dimension
    Teuchos::RCP<LINALG::SerialDenseMatrix>  diffx,     //!< (i) distance vector between node and integration point
    Teuchos::RCP<LINALG::SerialDenseVector>& funct,     //!< (o) basis functions values
    Teuchos::RCP<LINALG::SerialDenseMatrix>& deriv,     //!< (o) spatial derivatives of basis functions
    Teuchos::RCP<const LINALG::SerialDenseVector> cconv = Teuchos::null, //!< (i) convection at cell points for inf-flux weighting basis functions
    Teuchos::RCP<LINALG::SerialDenseVector> sfunct = Teuchos::null,      //!< (i) basis functions values
    Teuchos::RCP<LINALG::SerialDenseMatrix> sderiv = Teuchos::null       //!< (i) spatial derivatives of basis functions
  ) const = 0;

protected:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Determines/sets range of prior function for given cut-off tolerance
   *///                                                              nis Mar12
  /*--------------------------------------------------------------------------*/
  virtual void SetRange(
    double rangeTol //!< tolerance at which prior function is numerically zero
  ) = 0;

  //! range of basis functions
  double range_;
};

/*============================================================================*/
/*!
 * \brief Maximum-Entropy approximation class
 *///                                                                nis Mar12
/*============================================================================*/
class MaxEntApprox: public MeshfreeApprox
{
private:

  // forward declaration of DualProblem class
  template<int dim> class DualProblem;

public:

  //! ctor
  MaxEntApprox(
    Teuchos::ParameterList const & params, //!< (i) parameter list for maxent basis functions
    int type                               //!< (i) type of basis function: (0) solution or (1) weighting basis functions
  );

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief From base class: evaluating meshfree basis functions and derivatives
   *///                                                           nis Mar12
  /*--------------------------------------------------------------------------*/
  int GetMeshfreeBasisFunction(
    const int                         dim,              //!< (i) dimension
    Teuchos::RCP<LINALG::SerialDenseMatrix>  diffx,     //!< (i) distance vector between node and integration point
    Teuchos::RCP<LINALG::SerialDenseVector>& funct,     //!< (o) basis functions values
    Teuchos::RCP<LINALG::SerialDenseMatrix>& deriv,     //!< (o) spatial derivatives of basis functions
    Teuchos::RCP<const LINALG::SerialDenseVector> cconv = Teuchos::null,//!< (i) convection at cell points for inf-flux weighting basis functions
    Teuchos::RCP<LINALG::SerialDenseVector> sfunct = Teuchos::null,     //!< (i) basis functions values
    Teuchos::RCP<LINALG::SerialDenseMatrix> sderiv = Teuchos::null      //!< (i) spatial derivatives of basis functions
    ) const
  {
    switch (dim)
    {
    case 1: return this->dual1d_->maxent_basisfunction(diffx, funct, deriv, cconv, sfunct, sderiv);
    case 2: return this->dual2d_->maxent_basisfunction(diffx, funct, deriv, cconv, sfunct, sderiv);
    case 3: return this->dual3d_->maxent_basisfunction(diffx, funct, deriv, cconv, sfunct, sderiv);
    default:
    {
      dserror("MaxEnt basis function only built for 1D, 2D, and 3D problems. Code should also work in R^n. Feel free to add!");
      return -1;
    }
    };
  };

private:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief Determines/sets range of prior function for given cut-off tolerance
   *///                                                              nis Mar12
  /*--------------------------------------------------------------------------*/
  void SetRange(
    double rangeTol //!< tolerance at which prior function is numerically zero
  );

  /*==========================================================================*/
  /*!
   * \brief  member class of functions to compute dual problem
   *
   * \author Keijo Nissen (nissen@lnm.mw.tum.de)
   *
   * \date Januar, 2013
   *
   * The rather awkward implementation as a member class with internal
   * classes makes it possible to solve the max-ent problem for different dual
   * functions. This is important for non-standard problems, e.g. for signed
   * basis functions, non-partition-of-unity basis functions, or modification
   * for the weighting basis functions to ensure Newton-stability at high
   * Peclet numbers.
   */
  /*==========================================================================*/
  template<int dim>
  class DualProblem
  {
  public:
    //! ctor
    DualProblem(
      const bool    pu,
      const double  neg,
      MaxEntApprox* maxent
    );

    /*--------------------------------------------------------------------------*/
    /*!
     * \brief  evaluating maxent basis functions and derivatives
     *///                                                           nis Mar12
    /*--------------------------------------------------------------------------*/
    int maxent_basisfunction(
      Teuchos::RCP<LINALG::SerialDenseMatrix>        diffx,  //!< (i) distance vector between node and integration point
      Teuchos::RCP<LINALG::SerialDenseVector>&       funct,  //!< (o) basis functions values
      Teuchos::RCP<LINALG::SerialDenseMatrix>&       deriv,  //!< (o) spatial derivatives of basis functions
      Teuchos::RCP<const LINALG::SerialDenseVector>  cconv,  //!< (i) convection at cell points for inf-flux basis solution functions
      Teuchos::RCP<LINALG::SerialDenseVector>        sfunct, //!< (o) weighting basis functions values
      Teuchos::RCP<LINALG::SerialDenseMatrix>        sderiv  //!< (o) spatial derivatives of basis solution functions
      );

  private:
    /*--------------------------------------------------------------------------*/
    /*!
     * \brief  function to compute the prior function and derivatives
     *///                                                              nis Mar12
    /*--------------------------------------------------------------------------*/
    void SetPriorFunctDeriv(
      LINALG::SerialDenseVector       & q    , //!< (o) prior functions
      LINALG::SerialDenseMatrix       & dxq  , //!< (o) spatial derivatives of prior functions scaled by 1/q
      LINALG::SerialDenseMatrix const & diffx  //!< (i) NOT TO BE CHANGED distance vector between node and integration point
      ) const;

    /*--------------------------------------------------------------------------*/
    /*!
     * \brief  function to compute the compliance condition and derivatives
     *///                                                              nis Mar12
    /*--------------------------------------------------------------------------*/
    void SetComplCondFunctDeriv(
      LINALG::SerialDenseMatrix                   & c    , //!< (o) constraints
      std::vector<LINALG::SerialDenseMatrix>      & dxc  , //!< (o) spatial derivatives of constraints
      const LINALG::SerialDenseMatrix             & diffx, //!< (i) distance vector between node and integration point
      const int                                     type,  //!< (i) type of basis function: 0=solution, 1=weighting
      Teuchos::RCP<const LINALG::SerialDenseVector> conv  //!< (i) convection at evaluation point for inf-flux weighting basis functions
    ) const;

    /*--------------------------------------------------------------------------*/
    /*!
     * \brief  Virtual base class for specific routine for updating parameters
     *         and evaluating spatial derivatives
     *
     * \author Keijo Nissen (nissen@lnm.mw.tum.de)
     *
     * \date Januar, 2013
     */
    /*--------------------------------------------------------------------------*/
    class DualProblemType
    {
    public:
      //! pure virtual function to update dual problem
      virtual int UpdateParams(
        LINALG::SerialDenseVector       & funct, //!< (o) basis functions values
        LINALG::Matrix<dim,1>           & r,     //!< (o) gradient of dual function
        LINALG::Matrix<dim,dim>         & J,     //!< (o) inverted Hessian of dual function
        int                       const   na,    //!< (i) number of nodes
        LINALG::SerialDenseVector const & q,     //!< (i) prior functions
        LINALG::SerialDenseMatrix const & c,     //!< (i) constraints
        LINALG::Matrix<dim,1>     const & lam    //!< (i) argument of dual problem
      ) = 0;

      //! pure virtual function to get spatial derivatives of basis functions
      virtual int GetDerivs(
        LINALG::SerialDenseVector              const & funct, //!< (i) basis functions values
        LINALG::SerialDenseMatrix                    & deriv, //!< (o) spatial derivatives of basis functions
        LINALG::Matrix<dim,1>                  const & lam,   //!< (i) argmax of dual problem
        LINALG::Matrix<dim,dim>                const & Jinv,  //!< (i) inverted Hessian of dual function at maximum
        LINALG::SerialDenseVector              const & q,     //!< (i) prior functions
        LINALG::SerialDenseMatrix              const & dxq,   //!< (i) spatial derivatives of prior functions; = q_{a,i}/q_a
        LINALG::SerialDenseMatrix              const & c,     //!< (i) constraints
        std::vector<LINALG::SerialDenseMatrix> const & dxc    //!< (i) spatial derivatives of constraints
      ) = 0;
    };

    /*--------------------------------------------------------------------------*/
    /*!
     * \brief Class to provide standard routine for updating parameters
     *        and evaluating spatial derivatives (for non-negativity and
     *        partition-of-unity)
     *
     * \author Keijo Nissen (nissen@lnm.mw.tum.de)
     *
     * \date Januar, 2013
     */
    /*--------------------------------------------------------------------------*/
    class DualStandard : public DualProblemType
    {
    public:
      //! standard update of dual problem
      //! (for non-negative basis functions, reduced due to partition-of-unity
      //!  constraint - see Arroyo2006)
      int UpdateParams(
        LINALG::SerialDenseVector       & funct, //!< (o) basis functions values
        LINALG::Matrix<dim,1>           & r    , //!< (o) derivative of dual function
        LINALG::Matrix<dim,dim>         & J    , //!< (o) inverted Hessian of dual function
        int                       const   na   , //!< (i) number of nodes
        LINALG::SerialDenseVector const & q    , //!< (i) prior functions
        LINALG::SerialDenseMatrix const & c    , //!< (i) constraints
        LINALG::Matrix<dim,1>     const & lam    //!< (i) argument of dual problem
      );

      //! computes spatial derivatives for standard problem
      int GetDerivs(
        LINALG::SerialDenseVector              const & funct, //!< (i) basis functions values
        LINALG::SerialDenseMatrix                    & deriv, //!< (o) spatial derivatives of basis functions
        LINALG::Matrix<dim,1>                  const & lam  , //!< (i) argmax of dual problem
        LINALG::Matrix<dim,dim>                const & Jinv , //!< (i) inverted Hessian of dual function at maximum
        LINALG::SerialDenseVector              const & q    , //!< (i) prior functions
        LINALG::SerialDenseMatrix              const & dxq  , //!< (i) spatial derivatives of prior functions
        LINALG::SerialDenseMatrix              const & c    , //!< (i) constraints
        std::vector<LINALG::SerialDenseMatrix> const & dxc    //!< (i) spatial derivatives of constraints
      );
    };

    /*------------------------------------------------------------------------*/

    //! pointer to DualProblemType to specify routine for updating parameters
    //! and evaluating spatial derivatives
    Teuchos::RCP<DualProblemType> dualprob_;

    //! pointer to outer max-ent approximation that owns this DualProblem
    MaxEntApprox* maxent_;

  }; // end member class DualProblem

  /*--------------------------------------------------------------------------*/

  Teuchos::RCP<DualProblem<1> > dual1d_;
  Teuchos::RCP<DualProblem<2> > dual2d_;
  Teuchos::RCP<DualProblem<3> > dual3d_;

  //! type of basis function: (0) solution or (1) weighting basis function
  int type_;
  //! variance of prior
  double var_;
  //! prior type
  int prior_;
  //! compliance type
  int cmpl_;
  //! tolerance at which Newton is considered to be converged
  double newtontol_;
  //! maximum number of Newton steps
  int newtonmaxiter_;
}; // end  class MaxEntApprox

} // end MESHFREE
} // end DRT

#endif // DRT_UTILS_MAXENT_BASISFUNCTION_H
