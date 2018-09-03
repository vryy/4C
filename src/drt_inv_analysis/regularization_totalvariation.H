/*----------------------------------------------------------------------*/
/*!
\file regularization_totalvariation.H

\brief Total variation type regularization

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef REGULARIZATION_TOTALVARIATION_H
#define REGULARIZATION_TOTALVARIATION_H

/*----------------------------------------------------------------------*/
/* headers */
#include "regularization_base.H"


namespace INVANA

{
  /*! \class RegularizationTotalVariation
  \brief A class for regularization on a vector 'theta' of parameters based on total variation

  Implementation of the Total Variation Functional
  (see "Nonlinear total variation based noise removal algorithms", Rudin et al. 1992)

  This specific implementation is based on a graph representation of data.
  (see "Nonlinear Multilayered Representation of Graph Signals", Hidane et al, 2013)

  The graph of data is part of the input and not constructed in here.

  The functional is defined as:
  \f[
    f_w(\theta)=\|\nabla\theta\|_1=
    \sum_i\left(\sum_j w_{i,j}(\theta_j-\theta_i)^2+\epsilon\right)^{\frac{1}{2}},
  \f]
  with the adjacency weights \f$ w_{i,j}\f$ of parameter \f$ \theta_i\f$ and
  neighbouring parameters \f$ \theta_j\f$.

  It's gradient is constructed with the contributions of \f$ \theta_i\f$ and \f$ \theta_j\f$.
  These are:
  \f[
    \frac{\partial f_w}{\partial \theta_i}=
    \frac{-\sum_j w_{i,j}(\theta_j-\theta_i)}
    {\sum_j w_{i,j}(\theta_j-\theta_i)^2+\epsilon};
  \f]

  and
  \f[
    \frac{\partial f_w}{\partial \theta_j}=
    \frac{ w_{i,j}(\theta_j-\theta_i)}
    {\sum_j w_{i,j}(\theta_j-\theta_i)^2+\epsilon};
  \f]

  Specific care is taken to importing parameters from different
  procs and to exporting the gradient contributions from all procs
  according to the graph layout given via the adjacency matrix.

  \author kehl \date 10/2014
 */
  class RegularizationTotalVariation : public RegularizationBase
  {
   public:
    //! constructor
    RegularizationTotalVariation();

    //! destructor
    virtual ~RegularizationTotalVariation() {}

    //! setup of the class content
    virtual void Setup(const Teuchos::ParameterList& invp);

    //! @name Evaluation methods
    //@{

    //! Evaluate value of the regularization functional
    virtual void Evaluate(const Epetra_MultiVector& theta,  /// [in] vector to be regularized
        double* value                                       /// [out] functional value
    );

    //! Evaluate gradient of the regularization functional
    virtual void EvaluateGradient(
        const Epetra_MultiVector& theta,           /// [in] vector to be regularized
        Teuchos::RCP<Epetra_MultiVector> gradient  /// [out] directional derivative wrt. theta
    );

    //@}

   private:
    //! differentiation epsilon
    double eps_;

    /*! \brief adjacency matrix
     *
     * graph defining the connectivity of the parameter vector
     * and the weight of every row wrt its neighbour columns.
     * This is usually provided by INVANA::ConnectivityData
     *
     */
    Teuchos::RCP<Epetra_CrsMatrix> adjacency_;
  };

}  // namespace INVANA

/*----------------------------------------------------------------------*/
#endif
