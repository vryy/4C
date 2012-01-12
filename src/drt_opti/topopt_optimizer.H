/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer.H

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef TOPOPT_OPTIMIZER_H_
#define TOPOPT_OPTIMIZER_H_


#include "../drt_inpar/inpar_topopt.H"
#include "Epetra_Vector.h"
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>


using namespace Teuchos;

namespace DRT
{
  class Discretization;
}


namespace TOPOPT
{
  /*
  \brief The optimizer for topology optimization

  The optimizer might uses several methods for optimization. Usually it uses
  the gradient of the objectives every optimization step.

  Possibly optimization algorithms are the MMA (method of moving asymptotes),
  the GCMMA (globally convergent MMA), see Svanberg or the Trust-Region method

  \author winklmaier
  \date 12/11
  */
  class Optimizer
  {
  public:
    /// constructor
    explicit Optimizer(RCP<DRT::Discretization> discret, const ParameterList& params);

    /// desctructor
    virtual ~Optimizer(){};

    RCP<Epetra_Vector> DensityI(){return dens_i_;}; // density field of last optimization step
    RCP<Epetra_Vector> DensityIp(){return dens_ip_;}; // density field of new optimization step
    RCP<DRT::Discretization> Discret(){return discret_;}; // discretization of optimization variable

  protected:

    /// set the initial density field of the optimization variable
    void SetInitialDensityField(const INPAR::TOPOPT::InitialField init,const int startfuncno);


  private:

    RCP<Epetra_Vector> dens_i_; // old density in row map
    RCP<Epetra_Vector> dens_ip_; // new density in row map
    RCP<DRT::Discretization> discret_; // discretization of optimization field
    const ParameterList& params_; // parameters

  }; // end class Optimizer
} // end namespace TOPOPT

#endif // #ifndef TOPOPT_OPTIMIZER_H_
#endif  // #ifdef CCADISCRET
