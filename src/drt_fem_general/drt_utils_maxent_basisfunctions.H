/*!--------------------------------------------------------------------------
\file drt_utils_maxent_basisfunctions.H

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/

#ifndef DRT_UTILS_MAXENT_BASISFUNCTION_H
#define DRT_UTILS_MAXENT_BASISFUNCTION_H

#include <vector>
#include "../linalg/linalg_fixedsizematrix.H"
#include <Teuchos_ParameterList.hpp>
#include "../drt_inpar/inpar_meshfree.H"

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
class MeshfreeApprox{
public:
  virtual ~MeshfreeApprox(){};

  virtual double GetRange(){return range_;};

  virtual int GetMeshfreeBasisFunction(
    LINALG::SerialDenseVector       & funct    , //!< (o) basis functions values
    LINALG::SerialDenseMatrix       & deriv    , //!< (o) spatial derivatives of basis functions
    const LINALG::SerialDenseMatrix & diffx    , //!< (i) distance vector between node and integration point
    const int                         dimension  //!< (i) dimension
  ) const = 0;

protected:

  virtual void SetRange(
    double rangeTol //!< tolerance at which prior function is numerically zero
  ) = 0;

  double range_;
};

/*============================================================================*/
/*!
 * \brief Maximum-Entropy approximation class
 *///                                                                nis Mar12
/*============================================================================*/
class MaxEntApprox: public MeshfreeApprox {
private:

  // forward declaration of DualProblem class
  template<int dim> class DualProblem;

public:
  //! ctor
  MaxEntApprox(
    Teuchos::ParameterList const & params //!< (i) parameter list for maxent basis functions
  );

  int GetMeshfreeBasisFunction(
    LINALG::SerialDenseVector       & funct , //!< (o) basis functions values
    LINALG::SerialDenseMatrix       & deriv , //!< (o) spatial derivatives of basis functions
    const LINALG::SerialDenseMatrix & diffx , //!< (i) distance vector between node and integration point
    const int                         dim     //!< (i) dimension
  ) const
  {
    switch (dim)
    {
    case 1: return this->dual1d_->maxent_basisfunction(funct, deriv, diffx);
    case 2: return this->dual2d_->maxent_basisfunction(funct, deriv, diffx);
    case 3: return this->dual3d_->maxent_basisfunction(funct, deriv, diffx);
    default: {dserror("fucked up!"); return -1;}
    };
  };

private:

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief determines range of prior function for given
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
      bool const     pu,
      double const   neg,
      MaxEntApprox* that
    );

    /*--------------------------------------------------------------------------*/
    /*!
     * \brief  evaluating maxent basis functions and derivatives
     *///                                                           nis Mar12
    /*--------------------------------------------------------------------------*/
    int maxent_basisfunction(
      LINALG::SerialDenseVector       & funct , //!< (o) basis functions values
      LINALG::SerialDenseMatrix       & deriv , //!< (o) spatial derivatives of basis functions
      LINALG::SerialDenseMatrix const & diffx    //!< (i) distance vector between node and integration point
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
      LINALG::SerialDenseMatrix              & c     , //!< (o) constraints
      std::vector<LINALG::SerialDenseMatrix> & dxc   , //!< (o) spatial derivatives of constraints
      LINALG::SerialDenseMatrix        const & diffx   //!< (i) distance vector between node and integration point
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
      virtual void UpdateParams(
        LINALG::SerialDenseVector       & ,
        LINALG::Matrix<dim,1>           & ,
        LINALG::Matrix<dim,dim>         & ,
        int                       const   ,
        LINALG::SerialDenseVector const & ,
        LINALG::SerialDenseMatrix const & ,
        LINALG::Matrix<dim,1>     const &
      ) = 0;

      //! pure virtual function to get spatial derivatives of basis functions
      virtual void GetDerivs(
        LINALG::SerialDenseVector              const & ,
        LINALG::SerialDenseMatrix                    & ,
        LINALG::Matrix<dim,1>                  const & ,
        LINALG::Matrix<dim,dim>                const & ,
        LINALG::SerialDenseVector              const & ,
        LINALG::SerialDenseMatrix              const & , // = q_{a,i}/q_a
        LINALG::SerialDenseMatrix              const & ,
        std::vector<LINALG::SerialDenseMatrix> const &
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
      void UpdateParams(
        LINALG::SerialDenseVector       & funct, //!< (o) basis functions values
        LINALG::Matrix<dim,1>           & r    , //!< (o) derivative of dual function
        LINALG::Matrix<dim,dim>         & J    , //!< (o) inverted Hessian of dual function
        int                       const   na   , //!< (i) number of nodes
        LINALG::SerialDenseVector const & q    , //!< (i) prior functions
        LINALG::SerialDenseMatrix const & c    , //!< (i) constraints
        LINALG::Matrix<dim,1>     const & lam    //!< (i) argmax of dual problem
      );

      //! computes spatial derivatives for standard problem
      void GetDerivs(
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

    //! explicite this_-pointer; this_-pointer of MaxEntApprox is thus also known in DualProblem.
    MaxEntApprox* this_;

    //! pointer to DualProblemType to specify routine for updating parameters
    //! and evaluating spatial derivatives
    Teuchos::RCP<DualProblemType> dualprob_;

  }; // end member class DualProblem

  /*--------------------------------------------------------------------------*/

  Teuchos::RCP<DualProblem<1> > dual1d_;
  Teuchos::RCP<DualProblem<2> > dual2d_;
  Teuchos::RCP<DualProblem<3> > dual3d_;

  //! variance of prior
  double    var_;
  //! prior type
  INPAR::MESHFREE::priortype prior_;
  //! skewness type of prior
  INPAR::MESHFREE::priorskew skew_;
  //! compliance type
  INPAR::MESHFREE::compltype cmpl_;
  //! compliance parameters
  Teuchos::ParameterList* cmpl_param_;
    //! tolerance at which Newton is considered to be converged
  double NewtonTol_;
  //! maximum number of Newton steps
  int NewtonMax_;
}; // end  class MaxEntApprox

} // end MESHFREE
} // end DRT

#endif // DRT_UTILS_MAXENT_BASISFUNCTION_H
