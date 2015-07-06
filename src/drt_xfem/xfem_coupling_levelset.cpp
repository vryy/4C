/*!
\file xfem_coupling_levelset.cpp

\brief manages the different types of level-set based coupling conditions and thereby builds the bridge between the
xfluid class and the cut-library

<pre>
Maintainer: Benedikt Schott & Magnus Winter
            {schott, winter}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/

#include <Teuchos_TimeMonitor.hpp>

#include "xfem_coupling_levelset.H"
#include "xfem_utils.H"

#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_fluid.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

XFEM::LevelSetCoupling::LevelSetCoupling(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
    const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
    const double                        time,      ///< time
    const int                           step       ///< time step
) : CouplingBase(bg_dis, cond_name, bg_dis,time,step)
{
  /// level-set field is given w.r.t background mesh
  /// NOTE: more generally it would be possible cutterdis != bg_dis for the single LevelSetCoupling,
  /// however, the unique bg_phinp vector stored in the ConditionManager has to be given w.r.t bgdis
  cutter_dis_ = bg_dis;

  SetConditionsToCopy();

  SetElementConditions();

  // set the averaging strategy
  SetAveragingStrategy();

  // set coupling discretization
  SetCouplingDiscretization();

  // create node-based vector with level-set values
  phinp_ = Teuchos::rcp(new Epetra_Vector(*cutter_dis_->NodeRowMap()));

  // read initial level-set field
  SetLevelSetField(time_);

  //For output:
  ls_output_ = cutter_dis_->Writer();

}


void XFEM::LevelSetCoupling::SetConditionsToCopy()
{
  // set only the unique given condition name
  conditions_to_copy_.push_back(cond_name_);

  // additional conditions required for the levelset field based on the cutter (background) mesh
  conditions_to_copy_.push_back("XFEMSurfDisplacement");
}


void XFEM::LevelSetCoupling::Output(
    const int step,
    const double time,
    const bool write_restart_data
)
{
  // output for interface
  //ls_output_->NewStep(step,time);

//  ls_output_->WriteVector("phinp_", phinp_);  //This one is not set at the beginning (Output initial field)
//
//  ls_output_->WriteElementData(firstoutputofrun_);  //NEEDED?!?!
//  firstoutputofrun_ = false;

  // write restart

  if (write_restart_data)
  {
    ls_output_->WriteVector("phinp_res", phinp_);
  }
}

// -------------------------------------------------------------------
// Read Restart data for cutter discretization
// -------------------------------------------------------------------
void XFEM::LevelSetCoupling::ReadRestart(
    const int step
)
{

  dserror("Not tested Level Set restart from file. Should most likely work though if this dserror is removed.");

  //-------- boundary discretization
  IO::DiscretizationReader boundaryreader(cutter_dis_, step);

  const double time = boundaryreader.ReadDouble("time");

  if(myrank_ == 0)
  {
    IO::cout << "            RESTART IS PERFORMED FROM FUNCTION IN INPUT FILE!                  " << IO::endl;
    IO::cout << "ReadRestart for Level Set Cut in Xfluid (time="<< time <<" ; step="<< step <<")" << IO::endl;
  }

  SetLevelSetField(time);

}

/*----------------------------------------------------------------------*
 | ... |
 *----------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::SetLevelSetField(const double time)
{
  // initializations
  int err(0);
  double value(0.0);
  std::vector<int> nodedofs;

  // get the function from the first element
  const int lid=0;
  DRT::Condition* cond = cutterele_conds_[lid].second;
  const int func_no = cond->GetInt("levelsetfieldno");

  // check for potential time curve
  const int curvenum  = cond->GetInt("levelsetcurve");

  // initialization of time-curve factor
  double curvefac = 0.0;

  // compute potential time curve or set time-curve factor to one
  if (curvenum >= 0)
  {
    // time factor (negative time indicating error)
    if (time >= 0.0)
      curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
    else dserror("Negative time in function evaluation: time = %f", time);
  }
  else curvefac = 1.0;

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<cutter_dis_->NumMyRowNodes();lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lnode = cutter_dis_->lRowNode(lnodeid);

    // get value
    value=DRT::Problem::Instance()->Funct(func_no-1).Evaluate(0,lnode->X(),time,NULL);

    double final_val = curvefac*value;

    // now copy the values
    err = phinp_->ReplaceMyValue(lnodeid,0,final_val);
    if (err != 0) dserror("error while inserting value into phinp_");
  }

  return;
}

/*----------------------------------------------------------------------*
 | set interface level set field at current time           schott 02/15 |
 *----------------------------------------------------------------------*/
void XFEM::LevelSetCouplingBC::PrepareSolve()
{

  if(myrank_ == 0) IO::cout << "\t set level-set field, time " << time_ << IO::endl;

  SetLevelSetField(time_);
  return;
}


void XFEM::LevelSetCouplingWeakDirichlet::EvaluateCouplingConditions(
    LINALG::Matrix<3,1>& ivel,
    LINALG::Matrix<3,1>& itraction,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond
)
{
  // evaluate interface velocity (given by weak Dirichlet condition)
  EvaluateDirichletFunction(ivel, x, cond, time_);

  // no interface traction to be evaluated
  itraction.Clear();
}

void XFEM::LevelSetCouplingWeakDirichlet::EvaluateCouplingConditionsOldState(
    LINALG::Matrix<3,1>& ivel,
    LINALG::Matrix<3,1>& itraction,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond
)
{
  // evaluate interface velocity (given by weak Dirichlet condition)
  EvaluateDirichletFunction(ivel, x, cond, time_-dt_);

  // no interface traction to be evaluated
  itraction.Clear();
}

void XFEM::LevelSetCouplingNeumann::EvaluateCouplingConditions(
    LINALG::Matrix<3,1>& ivel,
    LINALG::Matrix<3,1>& itraction,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond
)
{
  // no interface velocity to be evaluated
  ivel.Clear();

  // evaluate interface traction (given by Neumann condition)
  EvaluateNeumannFunction(itraction, x, cond, time_);
}

void XFEM::LevelSetCouplingNeumann::EvaluateCouplingConditionsOldState(
    LINALG::Matrix<3,1>& ivel,
    LINALG::Matrix<3,1>& itraction,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond
)
{
  // no interface velocity to be evaluated
  ivel.Clear();

  // evaluate interface traction (given by Neumann condition)
  EvaluateNeumannFunction(itraction, x, cond, time_-dt_);
}

/*----------------------------------------------------------------------*
 | Set the LevelSet Field from a two phase algorithm.                   |
 *----------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::SetLevelSetField(
   Teuchos::RCP<const Epetra_Vector> scalaraf,
   Teuchos::RCP<const Epetra_Vector> curvatureaf,
   Teuchos::RCP<Epetra_MultiVector>  smoothed_gradphiaf,
   Teuchos::RCP<DRT::Discretization> scatradis
   )
{

  //Has the settings for surface tension been set from the Algorithm.
  // WARNING!
  //   This is not nice programming practice. In the future this info might be beneficial to put in
  //   DESIGN XFEM LEVELSET TWOPHASE VOL CONDITIONS in the input file.
  if(not surfacetension_init_)
    dserror("You can't set a LevelSetField without specifying the surface tension specifications.");

  // initializations
  int err(0);
  double value(0.0);
  std::vector<int> nodedofs;

// CUT INFORMATION FROM LEVEL SET
  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<cutter_dis_->NumMyRowNodes();lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lscatranode = scatradis->lRowNode(lnodeid);

    // find out the global dof id of the last(!) dof at the scatra node
    const int numscatradof = scatradis->NumDof(0,lscatranode);
    const int globalscatradofid = scatradis->Dof(0,lscatranode,numscatradof-1);
    const int localscatradofid = scalaraf->Map().LID(globalscatradofid);
    if (localscatradofid < 0)
      dserror("localdofid not found in map for given globaldofid");

    // now copy the values
    value = (*scalaraf)[localscatradofid];
    err = phinp_->ReplaceMyValue(lnodeid,0,value);
    if (err != 0) dserror("error while inserting value into phinp_");
  }

// NODAL CURVATURE!!!!!!
//----------------------------------------------
  //Transfer the vectors onto the NodeColMap.
  if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_nodal_curvature)
  {
    //SAFETY check
    if(curvatureaf==Teuchos::null)
      dserror("Nodal curvature chosen and empty curvatureaf provided.");

      Teuchos::RCP<Epetra_Vector> curvaturenp_rownode = Teuchos::rcp(new Epetra_Vector(*cutter_dis_->NodeRowMap()));

    // loop all column nodes on the processor
    for(int lnodeid=0;lnodeid<cutter_dis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor's local scatra node
      DRT::Node* lscatranode = scatradis->lRowNode(lnodeid);

      // find out the global dof id of the last(!) dof at the scatra node
      const int numscatradof = scatradis->NumDof(0,lscatranode);
      const int globalscatradofid = scatradis->Dof(0,lscatranode,numscatradof-1);

      const int localscatradofid = curvatureaf->Map().LID(globalscatradofid);
      if (localscatradofid < 0)
        dserror("localdofid not found in map for given globaldofid");

      // now copy the values
      value = (*curvatureaf)[localscatradofid];
      err = curvaturenp_rownode->ReplaceMyValue(lnodeid,0,value);
      if (err != 0) dserror("error while inserting value into curvaturenp_rownode");
    }

  curvaturenp_node_ = Teuchos::rcp(new Epetra_Vector(*cutter_dis_->NodeColMap()));
  LINALG::Export(*curvaturenp_rownode,*curvaturenp_node_);
  }
//---------------------------------------------- // NODAL CURVATURE END


// SMOOTHED GRAD PHI!!!!!!
//----------------------------------------------
  // SMoothed gradphi needed for divgrad option or LB with smoothed Projection matrix.
  if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_divgrad_normal or
      (surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_laplacebeltrami and
          (laplacebeltrami_==INPAR::TWOPHASE::matrix_smoothed or laplacebeltrami_==INPAR::TWOPHASE::matrix_mixed_smoothed)))
  {
    //SAFETY check
    if(smoothed_gradphiaf==Teuchos::null)
      dserror("A smoothed grad phi is required, but an empty one is provided!");

    Teuchos::RCP<Epetra_MultiVector> gradphinp_smoothed_rownode = Teuchos::rcp(new Epetra_MultiVector(*cutter_dis_->NodeRowMap(),smoothed_gradphiaf->NumVectors()));
    int numvec = smoothed_gradphiaf->NumVectors();

    // loop all column nodes on the processor
    for(int lnodeid=0;lnodeid<cutter_dis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor's local scatra node
      DRT::Node* lscatranode = scatradis->lRowNode(lnodeid);

      // find out the global dof id of the last(!) dof at the scatra node
      const int numscatradof = scatradis->NumDof(0,lscatranode);
      const int globalscatradofid = scatradis->Dof(0,lscatranode,numscatradof-1);

      const int localscatradofid = smoothed_gradphiaf->Map().LID(globalscatradofid);
      if (localscatradofid < 0)
        dserror("localdofid not found in map for given globaldofid");

      // now copy the values
      for(int i=0; i<numvec; i++)
      {
        value = smoothed_gradphiaf->Pointers()[i][localscatradofid]; //Somehow it is turned around?
        err = gradphinp_smoothed_rownode->ReplaceMyValue(lnodeid,i,value);
        if (err != 0) dserror("error while inserting value into gradphinp_smoothed_rownode");
      }
    }

    gradphinp_smoothed_node_ = Teuchos::rcp(new Epetra_MultiVector(*cutter_dis_->NodeColMap(),smoothed_gradphiaf->NumVectors()));
    LINALG::Export(*gradphinp_smoothed_rownode,*gradphinp_smoothed_node_);
  }
//---------------------------------------------- // SMOOTHED GRAD PHI END

  //SAFETY CHECK
//----------------------------------------------
  //Both empty vectors sent.
  if (curvatureaf!=Teuchos::null and smoothed_gradphiaf!=Teuchos::null)
    if(not (surftensapprox_== INPAR::TWOPHASE::surface_tension_approx_laplacebeltrami
            and laplacebeltrami_ == INPAR::TWOPHASE::matrix_non_smoothed )
            )
    {
      dserror("You can not both have a nodal curvature and a smoothed_gradphinp prescribed at once.");
    }

//---------------------------------------------- // SAFETY CHECK

  //  //Transfer the vectors onto the DofColMap.
//  if(curvatureaf!=Teuchos::null)
//  {
//    curvaturenp_ = Teuchos::rcp(new Epetra_Vector(*scatradis->DofColMap()));
//    LINALG::Export(*curvatureaf,*curvaturenp_);
//  }
//  if(smoothed_gradphiaf!=Teuchos::null)
//  {
//    gradphinp_smoothed_ = Teuchos::rcp(new Epetra_MultiVector(*scatradis->DofColMap(),smoothed_gradphiaf->NumVectors()));
//    LINALG::Export(*smoothed_gradphiaf,*gradphinp_smoothed_);
//  }
//  if (curvaturenp_!=Teuchos::null and gradphinp_smoothed_!=Teuchos::null)
//    dserror("You can not have both a nodal curvature and a gradphinp prescribed at once.");


  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::GetInterfaceSlaveMaterial(
  DRT::Element* actele,
  Teuchos::RCP<MAT::Material> & mat
)
{
  XFEM::UTILS::GetVolumeCellMaterial(actele,mat,GEO::CUT::Point::inside);
}


// -------------------------------------------------------------------
// Read Restart data for ScaTra coupled level set
// -------------------------------------------------------------------
void XFEM::LevelSetCouplingTwoPhase::ReadRestart(
    const int step
)
{

  //-------- boundary discretization
  IO::DiscretizationReader boundaryreader(cutter_dis_, step);

  const double time = boundaryreader.ReadDouble("time");

  if(myrank_ == 0)
  {
    IO::cout << "           RESTART IS PERFORMED FROM STORED VALUES!                            " << IO::endl;
    IO::cout << "ReadRestart for Level Set Cut in Xfluid (time="<< time <<" ; step="<< step <<")" << IO::endl;
  }

  boundaryreader.ReadVector(phinp_,   "phinp_res");

  if (not (cutter_dis_->NodeRowMap())->SameAs(phinp_->Map()))
    dserror("Global node numbering in maps does not match");

}
