/*-----------------------------------------------------------*/
/*! \file

\brief Method to deal with womersley flow profiles


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FLUIDWOMERSLEYCONDITION_H
#define FLUIDWOMERSLEYCONDITION_H


#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include <Epetra_Map.h>

#include "fluid_utils_mapextractor.H"
#include "fluid_utils.H"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../drt_io/io.H"
#include "../drt_fluid_ele/fluid_ele_action.H"


namespace FLD
{
  namespace UTILS
  {
    //--------------------------------------------------------------------
    // Wrapper class (to be called from outside) for volumetric surface
    // flow
    //--------------------------------------------------------------------

    /*!
    \brief Womersley boundary condition wrapper
    this class is meant to do some organisation stuff
    */
    class FluidVolumetricSurfaceFlowWrapper
    {
      friend class FluidImplicitTimeInt;


     public:
      /*!
      \brief Standard Constructor
      */
      FluidVolumetricSurfaceFlowWrapper(Teuchos::RCP<DRT::Discretization> actdis, double dta);

      /*!
      \brief Destructor
      */
      virtual ~FluidVolumetricSurfaceFlowWrapper();


      /*!
      \brief Wrapper for FluidVolumetricSurfaceFlowBc::EvaluateVelocities
      */
      void EvaluateVelocities(const Teuchos::RCP<Epetra_Vector> velocities, const double time);


      void InsertCondVector(Epetra_Vector& vec1, Epetra_Vector& vec2);

      /*!
      \brief Wrapper for FluidVolumetricSurfaceFlowBc::UpdateResidual
      */
      void UpdateResidual(Teuchos::RCP<Epetra_Vector> residual);


      /*!
      \brief Wrapper for FluidVolumetricSurfaceFlowBc::Output
      */
      void Output(IO::DiscretizationWriter& output);

      /*!
      \brief Wrapper for FluidVolumetricSurfaceFlowBc::ReadRestart
      */
      void ReadRestart(IO::DiscretizationReader& reader);


     private:
      /*!
      \brief all single fluid volumetric surface flow conditions
      */
      std::map<const int, Teuchos::RCP<class FluidVolumetricSurfaceFlowBc>> fvsf_map_;

      //! fluid discretization
      Teuchos::RCP<DRT::Discretization> discret_;

    };  // class FluidWomersleyWrapper

    class TotalTractionCorrector
    {
      friend class FluidImplicitTimeInt;


     public:
      /*!
      \brief Standard Constructor
      */

      TotalTractionCorrector(Teuchos::RCP<DRT::Discretization> actdis, double dta);

      /*!
      \brief Destructor
      */
      virtual ~TotalTractionCorrector();


      /*!
      \brief Wrapper for FluidVolumetricSurfaceFlowBc::EvaluateVelocities
      */
      void EvaluateVelocities(
          Teuchos::RCP<Epetra_Vector> bcdof, double time, double theta, double dta);

      /*!
      \brief export and set boundary values
      */
      void ExportAndSetBoundaryValues(
          Teuchos::RCP<Epetra_Vector> source, Teuchos::RCP<Epetra_Vector> target, std::string name);

      /*!
      \brief Wrapper for FluidVolumetricSurfaceFlowBc::UpdateResidual
      */
      void UpdateResidual(Teuchos::RCP<Epetra_Vector> residual);


      /*!
      \brief Wrapper for FluidVolumetricSurfaceFlowBc::Output
      */
      void Output(IO::DiscretizationWriter& output);

      /*!
      \brief Wrapper for FluidVolumetricSurfaceFlowBc::ReadRestart
      */
      void ReadRestart(IO::DiscretizationReader& reader);


     private:
      /*!
      \brief all single fluid volumetric surface flow conditions
      */
      std::map<const int, Teuchos::RCP<class FluidVolumetricSurfaceFlowBc>> fvsf_map_;

      //! fluid discretization
      Teuchos::RCP<DRT::Discretization> discret_;

    };  // class TotalTractionCorrector



    //--------------------------------------------------------------------
    // Actual Womersley bc calculation stuff
    //--------------------------------------------------------------------
    /*!
    \brief Womersley boundary condition

    */
    class FluidVolumetricSurfaceFlowBc
    {
      friend class FluidVolumetricSurfaceFlowWrapper;
      friend class TotalTractionCorrector;
      //  friend class FluidSurfaceTotalTractionCorrectionWrapper;

     public:
      /*!
      \brief Standard Constructor
      */
      FluidVolumetricSurfaceFlowBc(Teuchos::RCP<DRT::Discretization> actdis, double dta,
          std::string ds_condname, std::string dl_condname, int condid, int surf_numcond,
          int line_numcond);

      /*!
      \brief Empty Constructor
      */
      FluidVolumetricSurfaceFlowBc();

      /*!
      \brief Destructor
      */
      virtual ~FluidVolumetricSurfaceFlowBc() {}


      /*!
      \brief calculates the center of mass
      */
      void CenterOfMassCalculation(Teuchos::RCP<std::vector<double>> coords,
          Teuchos::RCP<std::vector<double>> normal, std::string ds_condname);


      /*!
      \brief calculates the local radii of all nodes
      */
      void EvalLocalNormalizedRadii(std::string ds_condname, std::string dl_condname);

      /*!
      \brief get the node row map of the womersley condition
      */
      void BuildConditionNodeRowMap(Teuchos::RCP<DRT::Discretization> dis, std::string condname,
          int condid, int condnum, Teuchos::RCP<Epetra_Map>& cond_noderowmap);

      /*!
      \brief get the dof row map of the womersley condition
      */
      void BuildConditionDofRowMap(Teuchos::RCP<DRT::Discretization> dis,
          const std::string condname, int condid, int condnum,
          Teuchos::RCP<Epetra_Map>& cond_dofrowmap);

      /*!
      \brief Evaluate velocities
      */
      void EvaluateVelocities(
          const double flowrate, const std::string ds_condname, const double time);


      /*!
      \brief Evaluate flowrate
      */
      double EvaluateFlowrate(const std::string ds_condname, const double time);


      /*!
      \brief UpdateResidual
      */
      void UpdateResidual(Teuchos::RCP<Epetra_Vector> residual);

      /*!
      \brief Evaluate velocities
      */
      void Velocities(Teuchos::RCP<DRT::Discretization> disc, Teuchos::RCP<Epetra_Vector> bcdof,
          Teuchos::RCP<Epetra_Map> cond_noderowmap, Teuchos::RCP<Epetra_Vector> local_radii,
          Teuchos::RCP<Epetra_Vector> border_radii, Teuchos::RCP<std::vector<double>> normal,
          Teuchos::RCP<Teuchos::ParameterList> params);

      /*!
      \brief Polynomail shaped velocity profile
      */
      double PolynomailVelocity(double r, int order);

      /*!
      \brief Womersley shaped velocity profile
      */
      double WomersleyVelocity(double r, double R, double Bn,
          // complex<double> Bn,
          double phi, int n, double t);

      /*!
      \brief Corrects the flow rate
      */
      void CorrectFlowRate(const Teuchos::ParameterList eleparams, const std::string ds_condname,
          const FLD::BoundaryAction action, const double time, const bool force_correction);

      /*!
      \brief Calculate the Flowrate on a boundary
      */
      double FlowRateCalculation(Teuchos::ParameterList eleparams, double time,
          std::string ds_condname, FLD::BoundaryAction action, int condid);

      double PressureCalculation(
          double time, std::string ds_condname, std::string action, int condid);
      /*!
      \brief Calculate the Flowrate on a boundary
      */
      void SetVelocities(const Teuchos::RCP<Epetra_Vector> velocities);

      /*!
      \brief Reset condition velocities
      */
      void ResetVelocities();

      /*!
      \brief evaluate the traction velocity component
      */
      void EvaluateTractionVelocityComp(Teuchos::ParameterList eleparams, std::string ds_condname,
          double flowrate, int condid, double time, double theta, double dta);

      /*!
      \brief export and set boundary values
      */
      void ExportAndSetBoundaryValues(
          Teuchos::RCP<Epetra_Vector> source, Teuchos::RCP<Epetra_Vector> target, std::string name);

      /*!
      \brief reset traction velocity components
      */
      void ResetTractionVelocityComp();

      /*!
      \brief Calculate the Flowrate on a boundary
      */
      void DFT(Teuchos::RCP<std::vector<double>> f,
          Teuchos::RCP<std::vector<std::complex<double>>>& F, int starting_pos);



     protected:
     private:
      /*!
      \brief calculate area at outflow boundary
      */
      double Area(double& density, double& viscosity, std::string ds_condname, int condid);

      /*!
      \brief output
      */
      void Output(IO::DiscretizationWriter& output, std::string ds_condname, int condnum);

      /*!
      \brief Read restart
      */
      void ReadRestart(IO::DiscretizationReader& reader, std::string ds_condname, int condnum);

      /*!
      \brief Bessel function of orders 0 and 1
      */
      std::complex<double> BesselJ01(std::complex<double> z, bool order);

      /*!
      \brief Interpolation function
      */
      void Interpolate(Teuchos::RCP<std::vector<double>> V1, Teuchos::RCP<std::vector<double>> V2,
          int index1, int& index2, double time);

      /*!
      \brief Return prebiasing flag
       */

      std::string PrebiasingFlag() { return prebiasing_flag_; }

     private:
      //! ID of present condition
      int condid_;

      //! Number of present surface condition
      int condnum_s_;

      //! Number of present line condition
      int condnum_l_;

      //! time period of present cyclic problem
      double period_;

      //! fluid viscosity
      double viscosity_;

      //! fluid density
      double density_;

      //! time step size
      double dta_;

      //! the processor ID from the communicator
      int myrank_;

      //! fluid discretization
      Teuchos::RCP<DRT::Discretization> discret_;

      //! Flowrate array for Womersley conditions
      Teuchos::RCP<std::vector<double>> flowrates_;

      //! Position at which the next element should be replaced
      //! initialised to zero as the first element will be replaced first
      int flowratespos_;

      //! center of mass coordinates
      Teuchos::RCP<std::vector<double>> cmass_;

      //! avarage normal of the surface
      Teuchos::RCP<std::vector<double>> normal_;

      //! direction normal of the velocity
      Teuchos::RCP<std::vector<double>> vnormal_;

      //! a Node row map of the nodes that belong to the current condition
      Teuchos::RCP<Epetra_Map> cond_surfnoderowmap_;

      //! a Node row map of the nodes that belong to border of the current condition
      Teuchos::RCP<Epetra_Map> cond_linenoderowmap_;

      //! a Dof row map of the degrees of freedom that belong to the current condition
      Teuchos::RCP<Epetra_Map> cond_dofrowmap_;

      //! A map of the local radii
      Teuchos::RCP<Epetra_Vector> local_radii_;

      //! A map of corresponding border radii
      Teuchos::RCP<Epetra_Vector> border_radii_;

      //! A map of only condition velocites
      Teuchos::RCP<Epetra_Vector> cond_velocities_;

      //! A dof col map of only condition velocites
      Teuchos::RCP<Epetra_Vector> drt_velocities_;

      //! A map of only condition velocites
      Teuchos::RCP<Epetra_Vector> cond_traction_vel_;

      //! initial area of the codition surface
      double area_;

      //! Number of modes
      int n_harmonics_;

      //! order of a polynomial velocity profile
      int order_;

      //! Type of the flow profile
      std::string flowprofile_type_;

      //! Prebiasing flag
      std::string prebiasing_flag_;

      //! is +1 if inflow, else -1
      double flow_dir_;

      //! flag to correct the flowprofile
      bool correct_flow_;


    };  // FluidVolumetricSurfaceFlowBc

  }  // namespace UTILS
}  // namespace FLD

#endif /* FLUIDWOMERSLEYCONDITION_H */
