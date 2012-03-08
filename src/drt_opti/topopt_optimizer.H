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


#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include "../drt_inpar/inpar_topopt.H"


namespace Teuchos
{
  class ParameterList;
}

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
    explicit Optimizer(Teuchos::RCP<DRT::Discretization> discret, const Teuchos::ParameterList& params);

    /// desctructor
    virtual ~Optimizer(){};

    double ComputeObjectiveValue(Teuchos::RCP<Epetra_Vector> porosity);

    void ImportFluidData(Teuchos::RCP<Epetra_Vector> vel,int step);
    void ImportAdjointFluidData(Teuchos::RCP<Epetra_Vector> vel, int step);

    Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > ExportFluidData(){return vel_;}; /// fluid velocity

    Teuchos::RCP<Epetra_Vector> DensityI(){return dens_i_;}; // density field of last optimization step
    Teuchos::RCP<Epetra_Vector> DensityIp(){return dens_ip_;}; // density field of new optimization step

    Teuchos::RCP<DRT::Discretization> Discret(){return discret_;}; // discretization of fluid field
    const Epetra_Map* RowMap(); // row map of optimization variable
    const Epetra_Map* ColMap(); // column map of optimization variable

  protected:

    /// set the initial density field of the optimization variable
    void SetInitialDensityField(const INPAR::TOPOPT::InitialField init,const int startfuncno);

    /// evaluate the objective function
    double EvaluateObjective(
        Teuchos::RCP<Epetra_Vector> porosity,
        const string& action);

  private:

    bool DataComplete() const;

    Teuchos::RCP<Epetra_Vector> dens_i_; // old density in row map
    Teuchos::RCP<Epetra_Vector> dens_ip_; // new density in row map
    Teuchos::RCP<DRT::Discretization> discret_; // discretization of optimization field
    const Teuchos::ParameterList& params_; // parameters

    //! @name fields for optimization

    Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > vel_;
    Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_Vector> > > adjointvel_;

    size_t num_timesteps_;
    //@}

    //! @name parameter of the objective function

    bool dissipation_; //! optimize dissipation
    bool pressure_inlet_; //! optimize pressure at inlet
    bool pressure_drop_; //! optimize pressure drop from inlet to outlet

    double dissipation_fac_; //! dissipation weighting factor
    double pressure_inlet_fac_; //! inlet pressure weighting factor
    double pressure_drop_fac_; //! inlet pressure weighting factor
    //@}

  }; // end class Optimizer
} // end namespace TOPOPT

#endif // #ifndef TOPOPT_OPTIMIZER_H_
#endif  // #ifdef CCADISCRET
