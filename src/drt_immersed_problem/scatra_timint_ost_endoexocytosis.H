/*!----------------------------------------------------------------------
\file scatra_timint_ost_endoexocytosis.H

\brief One-Step-Theta time-integration scheme for endo-/exocytosis model

\level 2

\maintainer Andreas Rauch

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_OST_ENDOEXOCYTOSIS_H
#define SCATRA_TIMINT_OST_ENDOEXOCYTOSIS_H

#include "../drt_scatra/scatra_timint_ost.H"

namespace SCATRA
{

/*!
\brief A class to model internalization and externalization of biomolecules, which may associate with the membrane in living cells

Internalization, also called endocytosis, is the intake of membrane-bound biomolecules mediated by the formation of vesicles.
This is modeled by an reaction rate equation that describes the sink at the membrane.
The transport of internalized biomolecules to endosomal and lysosomal compartments and the process of recycling is very complex.
The precise real location of these compartments is not known and my vary with time. It makes no sense to try to model this in detail.
Furthermore, we describe the cell as continuum, i.e. the cell organelles are smeared anyway.
Thus, we model all processes between internalization and the start of externalization, also called exocytosis, i.e. the transport
the cell cytoplasm back to the membrane, by a time delay \ref endoexo_delay_.
We model the spatial distribution of newly synthesized biomolecules as a reaction rate equation acting as source at the nucleus surface.
From there the transport back to the membrane starts.


<h3>Implementation</h3>

The vector \ref internalization_vec_ stores the amount of internalized molecules for each time step within the time interval between
intake and the beginning of exocytosis. Thus, its length is \ref endoexo_delay_ divided bay the time step size. We make sure by an
dserror, that this division yields an integer.

The amount of biomolecules is evaluated by integration of the concentration of a dummy species. The intake process is modeled by a
reaction rate equation in the .dat file, e.g.:

MAT 41 MAT_scatra_reaction NUMSCAL 3 STOICH -1 1 0 REACCOEFF 0.2 COUPLING simple_multiplicative // internalization at surface

Here, the first scalar is the externalized (membrane-bound) species. It is internalized with a reaction rate of 0.2 and the total internalized amount is stored in scalar 2.
So the rate of change of scalar 2 equals the rate of internalization. Since, by the boundary integral we determine amount of biomolecules at time
level n+1, also rates of change inferred from these values belong to time level n+1.


<h3>Usage</h3>

For an example on how to use this model see test case 'scatra_cylindrical_cell_endoexocytosis_quad4_hex8.dat'


<h3>Result Filter</h3>

If you filter your results using post_drt_ensight, the maximum dof number of all nodes in the Transport discretization is determined. This gives
the amount of phis which are written. If your membrane wears e.g. 2 dofs and the cell volume only 1, fields phi_1 and phi_2 are written.
For the cell volume phi_2 then has no meaning, i.e. it should just be zero everywhere in the domain.

You can split cell membrane discretization and volume discretization by aplying the Connectivity Filter to the scatra.case. If you apply the Threshold
filter subsequently, you get either the volume or the surface depending on your threshold value (o or 1).

\date 08/2016

\author rauch (rauch@lnm.mw.tum.de)
*/
class TimIntOneStepThetaEndoExocytosis : public TimIntOneStepTheta
{

public:

  //! Constructor
  TimIntOneStepThetaEndoExocytosis(
      Teuchos::RCP<DRT::Discretization>        actdis,        //!< discretization
      Teuchos::RCP<LINALG::Solver>             solver,        //!< linear solver
      Teuchos::RCP<Teuchos::ParameterList>     params,        //!< parameter list
      Teuchos::RCP<Teuchos::ParameterList>     extraparams,   //!< supplementary parameter list
      Teuchos::RCP<IO::DiscretizationWriter>   output,        //!< output writer
      const int                                probnum = 0    //!< global problem number
      );

  //! Destructor
  virtual ~TimIntOneStepThetaEndoExocytosis();

protected:

  //! initialize time integration
  virtual void Init();

  //! Pre Solve Operator
  virtual void PreSolve();

  //! Post Solve Operator
  virtual void PostSolve();

private:

  //! @name various scalar members
  //@{
  double endoexo_delay_;      //!< Time delay between internalization and begin of exocytotic transport
  int internalization_steps_; //!< length of vector internalization_vec_ = Number of steps stored in internalization_vec_
  double exo_surface_area_;   //!< Area of source surface where exocytotic transport starts
  double delta_phi_;          //!< concentration difference of internalized scalar
  //@}

  //! @name various 'std::vector' members
  //@{
  std::vector<double> Delta_phi_;    //!< Internalized scalars
  std::vector<double> mean_scalars_; //!< Vector containing mean scalar values and the domain volume at the last entry
  std::vector<double> internalization_vec_; //<! Vector containing the endoexo_delay_/dta_ step values of internalized molecules
  //@}

  //! @name various global vector
  //@{
  Teuchos::RCP<Epetra_Vector> source_; //!< Source at externalization surface
  //@}

};//class TimIntOneStepThetaEndoExocytosis

}// namespace SCATRA


#endif
