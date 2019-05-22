/*-----------------------------------------------------------*/
/*!

\brief Method to deal with coupling between 3D fluid and 0/1D vascular
problem

\maintainer Martin Kronbichler

\level 3
*/
/*-----------------------------------------------------------*/

//#ifdef D_COUPLED_ARTNET
#ifndef FLUID_COUPLING_RED_MODELS_H_
#define FLUID_COUPLING_RED_MODELS_H_


#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_art_net/art_net_dyn_drt.H"
#include "../drt_red_airways/red_airways_dyn_drt.H"
#include "../drt_io/io.H"


namespace FLD
{
  namespace UTILS
  {
    //--------------------------------------------------------------------
    // Wrapper class (to be called from outside) for coupling 3D/red D  bc
    //--------------------------------------------------------------------

    /*!
    \brief coupling boundary condition wrapper
    this class is meant to do some organisation stuff

    \author ismail (ismail@lnm.mw.tum.de)
    */

    class Fluid_couplingBc;

    class Fluid_couplingWrapperBase
    {
      friend class FluidImplicitTimeInt;

     public:
      /*!
      \brief Standard Constructor
      */
      Fluid_couplingWrapperBase(Teuchos::RCP<DRT::Discretization> dis_3D,
          Teuchos::RCP<DRT::Discretization> dis_redD, IO::DiscretizationWriter& output,
          double dt_3d, double dt_rm);

      /*!
      \brief Destructor
      */
      virtual ~Fluid_couplingWrapperBase() = 0;

      /*!
      \brief Wrapper for Fluid_couplingWrapper::FlowRateCalculation
      */
      void FlowRateCalculation(double time, double dta);

      /*!
      \brief Wrapper for Fluid_couplingWrapper::PressureCalculation
      */
      void PressureCalculation(double time, double dta);

      /*!
      \brief Wrapper for Fluid_couplingWrapper::OutflowBoundary
      */
      void ApplyBoundaryConditions(double time, double dta, double theta);

      /*!
      \brief Wrapper for Fluid_couplingWrapper::UpdateResidual
      */
      void UpdateResidual(Teuchos::RCP<Epetra_Vector> residual);



      void EvaluateDirichlet(
          Teuchos::RCP<Epetra_Vector> velnp, const Epetra_Map& condmap, double time);

      /*!
      \brief Wrapper for Fluid_couplingWrapper::WriteRestart
      */
      void WriteRestart(IO::DiscretizationWriter& output);

      /*!
      \brief Wrapper for Fluid_couplingWrapper::ReadRestart
      */
      void ReadRestart(IO::DiscretizationReader& reader);


      virtual void Integrate(bool flag, Teuchos::RCP<Teuchos::ParameterList>&) = 0;

      virtual void SaveState() = 0;

      virtual void LoadState() = 0;

      /*!
      \brief compute TimeUpdate
      */
      virtual void TimeUpdate() = 0;

     private:
      /*!
      \brief all single coupling conditions
      */
      std::map<const int, Teuchos::RCP<Fluid_couplingBc>> coup_map3D_;

      //! map of coupling variables returned by the reduced-D model at time step n+1
      Teuchos::RCP<std::map<std::string, double>> mapRed_Dnp_;

      //! map of coupling variables returned by the reduced-D model at time step n
      Teuchos::RCP<std::map<std::string, double>> mapRed_Dn_;

      //! map of coupling variables returned by the 3-D model at time step n+1
      Teuchos::RCP<std::map<std::string, double>> map3_Dnp_;

      //! map of coupling variables returned by the 3-D model at time step n
      Teuchos::RCP<std::map<std::string, double>> map3_Dn_;

      //! 3D fluid discretization
      Teuchos::RCP<DRT::Discretization> discret3D_;

      //! Reduced-D artery network discretization
      Teuchos::RCP<DRT::Discretization> discret_redD_;

      //! Reduced-D artery network time integration
      //  Teuchos::RCP<ART::ArtNetExplicitTimeInt>              ArtExpTime_integ_;


      //! the output writer
      IO::DiscretizationWriter& output_;

      //! the fluid 3D time step size
      double dt_f3_;

      //! the reduced model time step size
      double dt_rm_;

    };  // class Fluid_couplingWrapper

    template <class red_D_time_int>

    class Fluid_couplingWrapper : public Fluid_couplingWrapperBase
    {
     public:
      Fluid_couplingWrapper(Teuchos::RCP<DRT::Discretization> dis_3D,
          Teuchos::RCP<DRT::Discretization> dis_redD, Teuchos::RCP<red_D_time_int> time_intg,
          IO::DiscretizationWriter& output, double dt_3d, double dt_rm)
          : Fluid_couplingWrapperBase(dis_3D, dis_redD, output, dt_3d, dt_rm),
            reduced_D_time_integ_(time_intg)
      {
      }

      void Integrate(bool flag, Teuchos::RCP<Teuchos::ParameterList>& params)
      {
        reduced_D_time_integ_->Integrate(true, params);
      }

      void SaveState() { reduced_D_time_integ_->SaveState(); }

      void LoadState() { reduced_D_time_integ_->LoadState(); }

      void TimeUpdate() { reduced_D_time_integ_->TimeUpdate(); }

      ~Fluid_couplingWrapper(){};

     private:
      //! Reduced-D artery network time integration
      Teuchos::RCP<red_D_time_int> reduced_D_time_integ_;
    };

    //--------------------------------------------------------------------
    // Actual coupling bc calculation
    //--------------------------------------------------------------------
    /*!
    \brief coupling boundary condition for vascular outflow boundaries

    \author ismail (ismail@lnm.mw.tum.de)
    */

    class Fluid_couplingBc
    {
      friend class Fluid_couplingWrapperBase;

     public:
      /*!
      \brief Standard Constructor
      */
      Fluid_couplingBc(Teuchos::RCP<DRT::Discretization> dis3D,
          Teuchos::RCP<DRT::Discretization> dis_reD, IO::DiscretizationWriter& output, double dt_3d,
          double dt_rm, int condid, int numcond, int numcond2);

      /*!
      \brief Empty Constructor
      */
      Fluid_couplingBc();

      /*!
      \brief Destructor
      */
      virtual ~Fluid_couplingBc() {}

     protected:
      /*!
      \brief write flowrates_ and flowratespos_ to result files
      */
      void WriteRestart(IO::DiscretizationWriter& output, int condnum);

      /*!
      \brief read flowrates_ and flowratespos_
      */
      void ReadRestart(IO::DiscretizationReader& reader, int condnum);


      /*!
        \brief compute and store flow rate of all previous
        time steps belonging to one cycle
      */
      double FlowRateCalculation(double time, double dta, int condid);

      /*!
        \brief compute and store pressure of all previous
        time steps belonging to one cycle
      */
      double PressureCalculation(double time, double dta, int condid);


      /*!
        \brief compute convolution integral and apply pressure
        to elements
      */
      void OutflowBoundary(double pressure, double time, double dta, double theta, int condid);

      /*!
        \brief compute apply inflow as a Dirichlet BC
      */
      void InflowBoundary(double flowrate, double time, double dta, double theta, int condid);

      void UpdateResidual(Teuchos::RCP<Epetra_Vector> residual);



      void EvaluateDirichlet(
          Teuchos::RCP<Epetra_Vector> velnp, const Epetra_Map& condmap, double time);

      /*!
      \brief compute TimeUpdate
      */
      void TimeUpdate() {}

      void Integrate(bool flag, Teuchos::RCP<Teuchos::ParameterList>& params) {}

      void SaveState(){};

      void LoadState(){};

     private:
      /*!
      \brief calculate area at outflow boundary
      */
      double Area(double& density, double& viscosity, int numcond);



     protected:
      // coupled neumann BC
      Teuchos::RCP<Epetra_Vector> couplingbc_;

     private:
      //! ID of present condition
      int condid_;

      //! 3D fluid time step size
      double dt_f3_;

      //! reduced-D time step size
      double dt_rm_;

      //! coupling error at the boundary
      double max_error_;

      //! number of maximum allowable iterations at the boundary
      double max_itr_;

      //! velocity
      double velocity_;

      //! the processor ID from the communicator
      int myrank_;

      //! 3D fluid discretization
      Teuchos::RCP<DRT::Discretization> discret_3D_;

      //! fluid discretization
      Teuchos::RCP<DRT::Discretization> discret_redD_;

      //! the output writer
      IO::DiscretizationWriter& output_;

      //! flow rate
      double flowrate_;

      //! pressure
      double pressure_;

      //! corrector variable for dirichlet velocity s.t. applied flowrate is correct
      double alfa_;

    };  // class Fluid_couplingBc

  }  // namespace UTILS
}  // namespace FLD

#endif /* FLUID_COUPLING_RED_MODELS_H_ */
