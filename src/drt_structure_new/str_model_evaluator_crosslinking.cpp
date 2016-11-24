/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_crosslinking.cpp

\brief model evaluator for crosslinking in biopolymer networks

\maintainer Jonas Eichinger

\date May, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator_crosslinking.H"

#include "str_model_evaluator_data.H"
#include "str_timint_databiopolynetdyn.H"
#include "str_timint_base.H"
#include "str_utils.H"
#include "str_integrator.H"

#include <Epetra_Vector.h>
#include <Epetra_Time.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_fsi/fsi_matrixtransform.H"
#include "../drt_adapter/adapter_coupling.H"

#include "../drt_particle/particle_algorithm.H"
#include "../drt_beam3/beam3_base.H"
#include "../drt_beamcontact/beam3tobeamlinkage.H"
#include "../drt_biopolynet/biopolynet_calc_utils.H"
#include "../drt_biopolynet/crosslinker_node.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Crosslinking::Crosslinking():
    eval_statmech_ptr_(Teuchos::null),
    myrank_(-1),
    numproc_(-1),
    ia_discret_(Teuchos::null),
    coupsia_(Teuchos::null),
    siatransform_(Teuchos::null),
    force_crosslink_(Teuchos::null),
    stiff_crosslink_(Teuchos::null),
    ia_force_crosslink_(Teuchos::null),
    ia_stiff_crosslink_(Teuchos::null),
    ia_disnp_(Teuchos::null),
    binning_(Teuchos::null),
    bindis_(Teuchos::null),
    rowbins_(Teuchos::null),
    bin_beamcontent_(INPAR::BINSTRATEGY::Beam)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Setup()
{
  CheckInit();

  // -------------------------------------------------------------------------
  // setup variables
  // -------------------------------------------------------------------------
  // data
  // todo: this container belongs to browndyn
  eval_statmech_ptr_ = EvalData().StatMechPtr();
  // stiff
  stiff_crosslink_ = Teuchos::rcp(new
      LINALG::SparseMatrix(*GState().DofRowMapView(), 81, true, true));
  // force
  force_crosslink_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(),true));
  // get myrank
  myrank_ = DiscretPtr()->Comm().MyPID();
  // get number of procs
  numproc_ = DiscretPtr()->Comm().NumProc();
  // print logo
  Logo();

  // -------------------------------------------------------------------------
  // clone problem discretization, the idea is simple: we redistribute only
  // the new discretiztion to enable all interactions (including the required
  // search), calculate the resulting force and stiffness contributions, export
  // them no our initial discretization where all evaluation, assembly and
  // solving is done. Therefore the maps of our initial discretization don't
  // change, i.e. there is no need to rebuild the global state.
  // -------------------------------------------------------------------------
  Teuchos::RCP<DRT::UTILS::DiscretizationCreatorBase>  discloner =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreatorBase());
  ia_discret_ = discloner->CreateMatchingDiscretization(DiscretPtr(),"intacdis");

  // get initial (shifted) displacement vector, from now on this vector is based
  // on the maps of the interaction discretization
  ia_disnp_ = Teuchos::rcp(new Epetra_Vector(*GStatePtr()->GetMutableDisNp()));

  // initialize coupling adapter to transform vectors and matrix transform to
  // transform matrices between the two discrets (with distinct parallel
  // distribution)
  coupsia_ = Teuchos::rcp(new ADAPTER::Coupling());
  siatransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);

  // -------------------------------------------------------------------------
  // initialize crosslinker, i.e. add nodes (according to number of crosslinker
  // you want) to bin discretization and set their random reference position
  // -------------------------------------------------------------------------
  InitializeBinDiscret();
  // -------------------------------------------------------------------------
  // initialize particle algorithm, although we don't need/use any of the
  // actual particle "algorithm" or integrator, we are just using some nice
  // methods here.
  // -------------------------------------------------------------------------
  const Teuchos::ParameterList& params = DRT::Problem::Instance()->ParticleParams();
  /// algorithm is created here
  binning_ = Teuchos::rcp(new PARTICLE::Algorithm(DiscretPtr()->Comm(),params));

  // -------------------------------------------------------------------------
  // build periodic boundary conditions in binning strategy
  // -------------------------------------------------------------------------
  if(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()->at(0)>0.0)
    binning_->BuildPeriodicBC(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength());

  // -------------------------------------------------------------------------
  // a complete partitioning including the creation of bins, their weighted
  // distribution to procs as well a new distribution of the beam discret is
  // done here
  // -------------------------------------------------------------------------
  UpdateBinStrategy(false,true);

  // init force vector and stiffness matrix
  ia_stiff_crosslink_ = Teuchos::rcp(new
      LINALG::SparseMatrix(*ia_discret_->DofRowMap(), 81, true, true));
  // force
  ia_force_crosslink_ = Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofRowMap(),true));

  // gather data for all column crosslinker initially
  const int numcolcl = bindis_->NumMyColNodes();
  crosslinker_data_.resize(numcolcl);
  PreComputeCrosslinkerData(numcolcl);

  // gather data for all column beams
  const int numcolbeams = ia_discret_->NumMyColElements();
  beam_data_.resize(numcolbeams);
  PreComputeBeamData(numcolbeams);

  // set flag
  issetup_ = true;

  // that's it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Reset(const Epetra_Vector& x)
{
  CheckInitSetup();

  // get current displacement state and export to interaction discretization dofmap
  UpdateDofMapOfVector(ia_discret_, ia_disnp_, GState().GetMutableDisNp());
  ResetStateOfElementPairs();

  // Zero out force and stiffness contributions
  force_crosslink_->PutScalar(0.0);
  ia_force_crosslink_->PutScalar(0.0);
  stiff_crosslink_->Zero();
  ia_stiff_crosslink_->Zero();

  // that's it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Crosslinking::EvaluateForce()
{
  CheckInitSetup();

  // in case we need to communicate crosslinker
  std::map<int, std::vector<CommForceStiff> > sendforcestiff;
  std::vector<CommForceStiff> recvforcestiff;

  // force and moment exerted on the two connection sites due to the mechanical connection
  LINALG::TMatrix<double,6,1> bspotforce1(true);
  LINALG::TMatrix<double,6,1> bspotforce2(true);

  // resulting discrete element force vectors of the two parent elements
  Epetra_SerialDenseVector ele1force(12);
  Epetra_SerialDenseVector ele2force(12);
  Epetra_SerialDenseMatrix dummystiff(0,0);

  std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> >::const_iterator iter;
  for (iter=doublebondcl_.begin(); iter!=doublebondcl_.end(); ++iter)
  {

    Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> elepairptr = iter->second;

    bspotforce1.Clear();
    bspotforce2.Clear();

    // evaluate beam linkage object to get forces and moments on binding spots
    elepairptr->EvaluateForce(bspotforce1,bspotforce2);


    // ********************** Interpolation *******************************
    DRT::Element* ele1 = ia_discret_->gElement(elepairptr->GetEleGid(0));
    DRT::Element* ele2 = ia_discret_->gElement(elepairptr->GetEleGid(1));

    DRT::ELEMENTS::Beam3Base* beamele1 =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele1);
    DRT::ELEMENTS::Beam3Base* beamele2 =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele2);

    LINALG::TMatrix<double,6,12> trafomat(true);

    beamele1->GetGeneralizedInterpolationMatrixVariationsAtXi(
        trafomat,
        beamele1->GetBindingSpotXi(elepairptr->GetLocBSpotNum(0)));

    LINALG::TMatrix<double,12,1> eleforce;

    eleforce.MultiplyTN(trafomat,bspotforce1);

    for (unsigned int i=0; i<12; ++i)
      ele1force(i) = eleforce(i);


    trafomat.Clear();
    eleforce.Clear();

    beamele2->GetGeneralizedInterpolationMatrixVariationsAtXi(
        trafomat,
        beamele2->GetBindingSpotXi(elepairptr->GetLocBSpotNum(1)));

    eleforce.MultiplyTN(trafomat,bspotforce2);

    for (unsigned int i=0; i<12; ++i)
      ele2force(i) = eleforce(i);
    // ********************** end: Interpolation *******************************

    // assemble the contributions into force vector class variable
    // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
    AssembleEleForceStiffIntoSystemVectorMatrix(
        *ia_discret_,
        elepairptr->GetEleGid(0),
        elepairptr->GetEleGid(1),
        ele1force,
        ele2force,
        dummystiff,
        dummystiff,
        dummystiff,
        dummystiff,
        ia_force_crosslink_,
        Teuchos::null);

    // if needed, communicate force and stiff contributions to involved procs
    const std::set<int>& invprocs = elepairptr->GetInvolvedProcs();
    if(static_cast<int>(invprocs.size())>0)
    {
      // init force and stiff data struct that will be communicated
      CommForceStiff forcestiff;
      forcestiff.elegid1 = elepairptr->GetEleGid(0);
      forcestiff.elegid2 = elepairptr->GetEleGid(1);
      forcestiff.ele1force = ele1force;
      forcestiff.ele2force = ele2force;
      forcestiff.ele11stiff = dummystiff;
      forcestiff.ele12stiff = dummystiff;
      forcestiff.ele21stiff = dummystiff;
      forcestiff.ele22stiff = dummystiff;

      std::set<int>::const_iterator p;
      for(p=invprocs.begin(); p!=invprocs.end(); ++p)
        sendforcestiff[*p].push_back(forcestiff);
    }

  }

  CommunicateForceStiff(sendforcestiff, recvforcestiff);

  AssembleRecvEleForceStiffIntoSystemVectorMatrix(recvforcestiff,ia_force_crosslink_,Teuchos::null);

  // transformation from ia_discret to problem discret
  TransformForceAndStiff(true, false);

  // that is it
return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Crosslinking::EvaluateStiff()
{
  CheckInitSetup();

  // in case we need to communicate crosslinker
  std::map<int, std::vector<CommForceStiff> > sendforcestiff;
  std::vector<CommForceStiff> recvforcestiff;

  /* linearizations, i.e. stiffness contributions due to forces on the two
   * connection sites due to the mechanical connection */
  LINALG::TMatrix<double,6,6> bspotstiff11(true);
  LINALG::TMatrix<double,6,6> bspotstiff12(true);
  LINALG::TMatrix<double,6,6> bspotstiff21(true);
  LINALG::TMatrix<double,6,6> bspotstiff22(true);

  // linearizations, i.e. discrete stiffness contributions to the two parent elements
  // we can't handle this separately for both elements because there are entries which couple the two element stiffness blocks
  Epetra_SerialDenseVector dummyforce(0);
  Epetra_SerialDenseMatrix ele11stiff(12,12);
  Epetra_SerialDenseMatrix ele12stiff(12,12);
  Epetra_SerialDenseMatrix ele21stiff(12,12);
  Epetra_SerialDenseMatrix ele22stiff(12,12);

  ia_stiff_crosslink_->UnComplete();     // Todo check if needed or can be avoided


  // Todo needed?
  // transformation from row to column (after extended ghosting)
  Teuchos::RCP<Epetra_Vector> ia_discolnp =
      Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  LINALG::Export(*ia_disnp_, *ia_discolnp);


  std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> >::const_iterator iter;
  for (iter=doublebondcl_.begin(); iter!=doublebondcl_.end(); ++iter)
  {

    Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> elepairptr = iter->second;

    bspotstiff11.Clear();
    bspotstiff12.Clear();
    bspotstiff21.Clear();
    bspotstiff22.Clear();

    // evaluate beam linkage object to get forces and moments on binding spots
    elepairptr->EvaluateStiff(
        bspotstiff11,
        bspotstiff12,
        bspotstiff21,
        bspotstiff22);

    // ********************** Interpolation *******************************
    DRT::Element* ele1 = ia_discret_->gElement(elepairptr->GetEleGid(0));
    DRT::Element* ele2 = ia_discret_->gElement(elepairptr->GetEleGid(1));

    DRT::ELEMENTS::Beam3Base* beamele1 =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele1);
    DRT::ELEMENTS::Beam3Base* beamele2 =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele2);


    std::vector<double> eledisp;
    GetCurrentElementDis(ele1,ia_discolnp,eledisp);

    LINALG::TMatrix<double,6,12> trafomat(true);

    LINALG::TMatrix<double,12,6> ele11stifftmp;
    LINALG::TMatrix<double,12,6> ele12stifftmp;
    LINALG::TMatrix<double,6,12> ele21stifftmp;
    LINALG::TMatrix<double,6,12> ele22stifftmp;
    LINALG::TMatrix<double,12,12> ele11stifftmp2;
    LINALG::TMatrix<double,12,12> ele12stifftmp2;
    LINALG::TMatrix<double,12,12> ele21stifftmp2;
    LINALG::TMatrix<double,12,12> ele22stifftmp2;



    beamele1->GetGeneralizedInterpolationMatrixVariationsAtXi(
        trafomat,
        beamele1->GetBindingSpotXi(elepairptr->GetLocBSpotNum(0)));

    ele11stifftmp.MultiplyTN(trafomat,bspotstiff11);
    ele12stifftmp.MultiplyTN(trafomat,bspotstiff12);

    trafomat.Clear();


    beamele1->GetGeneralizedInterpolationMatrixIncrementsAtXi(
        trafomat,
        beamele1->GetBindingSpotXi(elepairptr->GetLocBSpotNum(0)),
        eledisp);

    ele11stifftmp2.Multiply(ele11stifftmp,trafomat);

    ele21stifftmp.Multiply(bspotstiff21,trafomat);

    GetCurrentElementDis(ele2,ia_discolnp,eledisp);

    trafomat.Clear();

    beamele2->GetGeneralizedInterpolationMatrixIncrementsAtXi(
        trafomat,
        beamele2->GetBindingSpotXi(elepairptr->GetLocBSpotNum(1)),
        eledisp);

    ele12stifftmp2.Multiply(ele12stifftmp,trafomat);
    ele22stifftmp.Multiply(bspotstiff22,trafomat);

    trafomat.Clear();

    beamele2->GetGeneralizedInterpolationMatrixVariationsAtXi(
        trafomat,
        beamele2->GetBindingSpotXi(elepairptr->GetLocBSpotNum(1)));

    ele21stifftmp2.MultiplyTN(trafomat,ele21stifftmp);
    ele22stifftmp2.MultiplyTN(trafomat,ele22stifftmp);



    // Apply point force to parent elements
    for (unsigned int i=0; i<12; ++i)
      for (unsigned int j=0; j<12; ++j)
      {
        ele11stiff(i,j) = ele11stifftmp2(i,j);
        ele12stiff(i,j) = ele12stifftmp2(i,j);

        ele21stiff(i,j) = ele21stifftmp2(i,j);
        ele22stiff(i,j) = ele22stifftmp2(i,j);
      }

    // ********************** End: Interpolation *******************************


    // assemble the contributions into stiffness matrix class variable
    // stiff_crosslink_ptr_, i.e. in the DOFs of the connected nodes
    AssembleEleForceStiffIntoSystemVectorMatrix(
        *ia_discret_,
        elepairptr->GetEleGid(0),
        elepairptr->GetEleGid(1),
        dummyforce,
        dummyforce,
        ele11stiff,
        ele12stiff,
        ele21stiff,
        ele22stiff,
        Teuchos::null,
        ia_stiff_crosslink_);

    // if needed, communicate force and stiff contributions to involved procs
    const std::set<int>& invprocs = elepairptr->GetInvolvedProcs();
    if(static_cast<int>(invprocs.size())>0)
    {
      // init force and stiff data struct that will be communicated
      CommForceStiff forcestiff;
      forcestiff.elegid1 = elepairptr->GetEleGid(0);
      forcestiff.elegid2 = elepairptr->GetEleGid(1);
      forcestiff.ele1force = dummyforce;
      forcestiff.ele2force = dummyforce;
      forcestiff.ele11stiff = ele11stiff;
      forcestiff.ele12stiff = ele12stiff;
      forcestiff.ele21stiff = ele21stiff;
      forcestiff.ele22stiff = ele22stiff;

      std::set<int>::const_iterator p;
      for(p=invprocs.begin(); p!=invprocs.end(); ++p)
        sendforcestiff[*p].push_back(forcestiff);
    }

  }

  CommunicateForceStiff(sendforcestiff, recvforcestiff);

  AssembleRecvEleForceStiffIntoSystemVectorMatrix(recvforcestiff,Teuchos::null,ia_stiff_crosslink_);

  if (not ia_stiff_crosslink_->Filled())
    ia_stiff_crosslink_->Complete();

  TransformForceAndStiff(false, true);

  if (not stiff_crosslink_->Filled())
    stiff_crosslink_->Complete();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Crosslinking::EvaluateForceStiff()
{
  CheckInitSetup();

  // in case we need to communicate crosslinker
  std::map<int, std::vector<CommForceStiff> > sendforcestiff;
  std::vector<CommForceStiff> recvforcestiff;

  // force and moment exerted on the two connection sites due to the mechanical connection
  LINALG::TMatrix<double,6,1> bspotforce1(true);
  LINALG::TMatrix<double,6,1> bspotforce2(true);
  // corresponding linearizations, i.e. stiffness contributions
  LINALG::TMatrix<double,6,6> bspotstiff11(true);
  LINALG::TMatrix<double,6,6> bspotstiff12(true);
  LINALG::TMatrix<double,6,6> bspotstiff21(true);
  LINALG::TMatrix<double,6,6> bspotstiff22(true);

  // resulting discrete element force vectors of the two parent elements
  Epetra_SerialDenseVector ele1force(12);
  Epetra_SerialDenseVector ele2force(12);
  // linearizations, i.e. discrete stiffness contributions to the two parent elements
  // we can't handle this separately for both elements because there are entries which couple the two element stiffness blocks
  Epetra_SerialDenseMatrix ele11stiff(12,12);
  Epetra_SerialDenseMatrix ele12stiff(12,12);
  Epetra_SerialDenseMatrix ele21stiff(12,12);
  Epetra_SerialDenseMatrix ele22stiff(12,12);

  ia_stiff_crosslink_->UnComplete();     // Todo check if needed or can be avoided



  // Todo needed?
  // transformation from row to column (after extended ghosting)
  Teuchos::RCP<Epetra_Vector> ia_discolnp =
      Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  LINALG::Export(*ia_disnp_, *ia_discolnp);



  std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> >::const_iterator iter;
  for (iter=doublebondcl_.begin(); iter!=doublebondcl_.end(); ++iter)
  {

    Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> elepairptr = iter->second;

    bspotforce1.Clear();
    bspotforce2.Clear();
    bspotstiff11.Clear();
    bspotstiff12.Clear();
    bspotstiff21.Clear();
    bspotstiff22.Clear();

    // evaluate beam linkage object to get forces and moments on binding spots
    elepairptr->EvaluateForceStiff(bspotforce1,
                                   bspotforce2,
                                   bspotstiff11,
                                   bspotstiff12,
                                   bspotstiff21,
                                   bspotstiff22);


    // ********************** Interpolation *******************************
    DRT::Element* ele1 = ia_discret_->gElement(elepairptr->GetEleGid(0));
    DRT::Element* ele2 = ia_discret_->gElement(elepairptr->GetEleGid(1));


    DRT::ELEMENTS::Beam3Base* beamele1 =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele1);
    DRT::ELEMENTS::Beam3Base* beamele2 =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele2);


    std::vector<double> eledisp;
    GetCurrentElementDis(ele1,ia_discolnp,eledisp);

    LINALG::TMatrix<double,6,12> trafomat(true);

    LINALG::TMatrix<double,12,1> ele1forcetmp;
    LINALG::TMatrix<double,12,1> ele2forcetmp;

    LINALG::TMatrix<double,12,6> ele11stifftmp;
    LINALG::TMatrix<double,12,6> ele12stifftmp;
    LINALG::TMatrix<double,6,12> ele21stifftmp;
    LINALG::TMatrix<double,6,12> ele22stifftmp;
    LINALG::TMatrix<double,12,12> ele11stifftmp2;
    LINALG::TMatrix<double,12,12> ele12stifftmp2;
    LINALG::TMatrix<double,12,12> ele21stifftmp2;
    LINALG::TMatrix<double,12,12> ele22stifftmp2;



    beamele1->GetGeneralizedInterpolationMatrixVariationsAtXi(
        trafomat,
        beamele1->GetBindingSpotXi(elepairptr->GetLocBSpotNum(0)));

    ele1forcetmp.MultiplyTN(trafomat,bspotforce1);

    ele11stifftmp.MultiplyTN(trafomat,bspotstiff11);
    ele12stifftmp.MultiplyTN(trafomat,bspotstiff12);

    trafomat.Clear();


    // ************* DEBUG **************************************
//    std::cout << "\neledisp:";
//    for (unsigned int i=0; i<eledisp.size(); ++i)
//      std::cout << " " << eledisp[i];
//    std::cout << "\n";
    // **********************************************************


    beamele1->GetGeneralizedInterpolationMatrixIncrementsAtXi(
        trafomat,
        beamele1->GetBindingSpotXi(elepairptr->GetLocBSpotNum(0)),
        eledisp);

    ele11stifftmp2.Multiply(ele11stifftmp,trafomat);

    ele21stifftmp.Multiply(bspotstiff21,trafomat);


    // ************** DEBUG *************************************
//    std::cout << "\nxi= "<< beamele1->GetBindingSpotXi(locbspotnum1) << ", trafomat: ";
//    trafomat.Print(std::cout);
    // **********************************************************

    GetCurrentElementDis(ele2,ia_discolnp,eledisp);

    trafomat.Clear();

    beamele2->GetGeneralizedInterpolationMatrixIncrementsAtXi(
        trafomat,
        beamele2->GetBindingSpotXi(elepairptr->GetLocBSpotNum(1)),
        eledisp);

    ele12stifftmp2.Multiply(ele12stifftmp,trafomat);
    ele22stifftmp.Multiply(bspotstiff22,trafomat);

    trafomat.Clear();

    beamele2->GetGeneralizedInterpolationMatrixVariationsAtXi(
        trafomat,
        beamele2->GetBindingSpotXi(elepairptr->GetLocBSpotNum(1)));

    ele2forcetmp.MultiplyTN(trafomat,bspotforce2);

    ele21stifftmp2.MultiplyTN(trafomat,ele21stifftmp);
    ele22stifftmp2.MultiplyTN(trafomat,ele22stifftmp);



    // Apply point force to parent elements
    for (unsigned int i=0; i<12; ++i)
    {
      ele1force(i) = ele1forcetmp(i);
      ele2force(i) = ele2forcetmp(i);

      for (unsigned int j=0; j<12; ++j)
      {
        ele11stiff(i,j) = ele11stifftmp2(i,j);
        ele12stiff(i,j) = ele12stifftmp2(i,j);

        ele21stiff(i,j) = ele21stifftmp2(i,j);
        ele22stiff(i,j) = ele22stifftmp2(i,j);
      }
    }

    // ********************** End: Interpolation *******************************

    // assemble the contributions into force and stiffness class variables
    // f_crosslink_np_ptr_, stiff_crosslink_ptr_, i.e. in the DOFs of the connected nodes
    AssembleEleForceStiffIntoSystemVectorMatrix(
        *ia_discret_,
        elepairptr->GetEleGid(0),
        elepairptr->GetEleGid(1),
        ele1force,
        ele2force,
        ele11stiff,
        ele12stiff,
        ele21stiff,
        ele22stiff,
        ia_force_crosslink_,
        ia_stiff_crosslink_);

    // if needed, communicate force and stiff contributions to involved procs
    const std::set<int>& invprocs = elepairptr->GetInvolvedProcs();
    if(static_cast<int>(invprocs.size())>0)
    {
      // init force and stiff data struct that will be communicated
      CommForceStiff forcestiff;
      forcestiff.elegid1 = elepairptr->GetEleGid(0);
      forcestiff.elegid2 = elepairptr->GetEleGid(1);
      forcestiff.ele1force = ele1force;
      forcestiff.ele2force = ele2force;
      forcestiff.ele11stiff = ele11stiff;
      forcestiff.ele12stiff = ele12stiff;
      forcestiff.ele21stiff = ele21stiff;
      forcestiff.ele22stiff = ele22stiff;

      std::set<int>::const_iterator p;
      for(p=invprocs.begin(); p!=invprocs.end(); ++p)
        sendforcestiff[*p].push_back(forcestiff);
    }

  }

  CommunicateForceStiff(sendforcestiff, recvforcestiff);

  AssembleRecvEleForceStiffIntoSystemVectorMatrix(recvforcestiff,ia_force_crosslink_,ia_stiff_crosslink_);

  if (not ia_stiff_crosslink_->Filled())
    ia_stiff_crosslink_->Complete();

  TransformForceAndStiff();

  if (not stiff_crosslink_->Filled())
    stiff_crosslink_->Complete();

  // that's it
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Crosslinking::AssembleForce(Epetra_Vector& f,
    const double & timefac_np) const
{
  CheckInitSetup();

  STR::AssembleVector(1.0,f,timefac_np,*force_crosslink_);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Crosslinking::AssembleJacobian(
    LINALG::SparseOperator& jac,
    const double & timefac_np) const
{
  CheckInitSetup();

  Teuchos::RCP<LINALG::SparseMatrix> jac_dd_ptr =
      GState().ExtractDisplBlock(jac);
  jac_dd_ptr->Add(*stiff_crosslink_,false,timefac_np,1.0);

  // no need to keep it
  stiff_crosslink_->Zero();
  ia_stiff_crosslink_->Zero();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::WriteRestart(
    IO::DiscretizationWriter& iowriter,
    const bool& forced_writerestart) const
{
  CheckInitSetup();

  return;
} // WriteRestart()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::ReadRestart(
    IO::DiscretizationReader& ioreader)
{
  CheckInitSetup();

  return;
} // ReadRestart()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::RecoverState(
    const Epetra_Vector& xold,
    const Epetra_Vector& dir,
    const Epetra_Vector& xnew)
{
 // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateStepState(
    const double& timefac_n)
{
  CheckInitSetup();

  // add the old time factor scaled contributions to the residual
  Teuchos::RCP<Epetra_Vector>& fstructold_ptr =
      GState().GetMutableFstructureOld();

  fstructold_ptr->Update(-timefac_n,*force_crosslink_,1.0);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateStepElement()
{
  // -------------------------------------------------------------------------
  // update all binding states and redistribute bindis and intactdis
  // -------------------------------------------------------------------------
  UpdateCrosslinking();

  return;

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::DetermineStressStrain()
{

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::DetermineEnergy()
{
  CheckInitSetup();
  dserror("Not yet implemented!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();

  // mesh is not written to disc, only maximum node id is important for output
  bindis_->Writer()->ParticleOutput(GState().GetStepN(), GState().GetTimeN(), false);
  bindis_->Writer()->NewStep(GState().GetStepN(), GState().GetTimeN());
  Teuchos::RCP<Epetra_Vector> dis = LINALG::CreateVector(*bindis_->DofRowMap(),true);
  Teuchos::RCP<Epetra_Vector> numbond = LINALG::CreateVector(*bindis_->NodeRowMap(),true);

  //todo: this is of course not nice, this needs to be done somewhere else
  for(int i=0;i<bindis_->NumMyRowNodes();++i)
  {
    CROSSLINKING::CrosslinkerNode *crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(bindis_->lRowNode(i));
    // std::vector holding gids of dofs
    std::vector<int> dofnode  = bindis_->Dof(crosslinker_i);

    // loop over all dofs
    for(int dim=0;dim<3;++dim)
    {
      int doflid = dis->Map().LID(dofnode[dim]);
      (*dis)[doflid] = crosslinker_i->X()[dim];
    }

    (*numbond)[i] = crosslinker_i->ClData()->GetNumberOfBonds();
  }
  bindis_->Writer()->WriteVector("displacement", dis);
  bindis_->Writer()->WriteVector("numbond", numbond, bindis_->Writer()->nodevector);
  // as we know that our maps have changed every time we write output, we can empty
  // the map cache as we can't get any advantage saving the maps anyway
  bindis_->Writer()->ClearMapCache();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::Crosslinking::
    GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Crosslinking::
    GetCurrentSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::Crosslinking::
    GetLastTimeStepSolutionPtr() const
{
  // there are no model specific solution entries
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::PostOutput()
{
  CheckInitSetup();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::ResetStepState()
{
  CheckInitSetup();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::ResetStateOfElementPairs()
{
  CheckInit();

  // get current displacement state and export to interaction discretization dofmap
  UpdateDofMapOfVector(ia_discret_, ia_disnp_, GState().GetMutableDisNp());

  // transformation from row to column
  Teuchos::RCP<Epetra_Vector> ia_discolnp =
      Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  LINALG::Export(*ia_disnp_, *ia_discolnp);

  std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> >::const_iterator iter;
  for (iter=doublebondcl_.begin(); iter!=doublebondcl_.end(); ++iter)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> elepairptr = iter->second;

    // init positions and triads
    std::vector<LINALG::Matrix<3,1> > pos(2);
    std::vector<LINALG::Matrix<3,3> > triad(2);

    for(int i=0; i<2; ++i)
    {
      int elegid = elepairptr->GetEleGid(i);
      int locbspotnum = elepairptr->GetLocBSpotNum(i);
      DRT::Element* ele = ia_discret_->gElement(elegid);

      GetPosAndTriadOfBindingSpot(ele, ia_discolnp,locbspotnum,pos[i],triad[i]);
    }
    // finally reset state
    elepairptr->ResetState(pos,triad);
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateBinStrategy(bool transfer,
                                                          bool partition,
                                                          bool repartition,
                                                          bool createxaabb,
                                                          bool setcutoff,
                                                          bool cltoclinteraction)
{
  CheckInit();

  // fixme:
  cltoclinteraction=true;

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::UpdateBinStrategy");

  if(repartition)
    dserror("repartioning not yet functional, will follow soon");

  // some safety checks
  if((transfer && partition) || (transfer && repartition) || (partition && repartition))
    dserror("You can either do a transfer, partitioning or repartitioning.");
  if((repartition || transfer) && (createxaabb || setcutoff))
    dserror("If you want to reset XAABB and/or cutoff you need to do a full partitioning");

  // store structure discretization in vector
  std::vector<Teuchos::RCP<DRT::Discretization> > discret_vec(1);
  discret_vec[0] = ia_discret_;
  // displacement vector according to periodic boundary conditions
  std::vector<Teuchos::RCP<Epetra_Vector> > disnp(1);
  disnp[0] = ia_disnp_;
  // todo: unshifted current positions are needed here
  if(setcutoff)
    dserror("unshifted configuration is needed (not yet here) for calculation of cutoff.");

  // create XAABB and optionally set cutoff radius
  if(createxaabb)
    binning_->CreateXAABB(discret_vec,disnp,setcutoff);
  // just set cutoff radius
  else if (setcutoff)
    binning_->ComputeMaxCutoff(discret_vec,disnp);

  // -------------------------------------------------------------------------
  // Create bins according to XAABB and cutoff set in constructor (read
  // of input file)
  // -------------------------------------------------------------------------
  if(partition)
    binning_->CreateBins(Teuchos::null);

  // -------------------------------------------------------------------------
  // assign bins to procs according to a weighted partitioning, i.e. bins
  // are weighted with respect to the number of nodes they contain and then
  // distributed to the procs. Then an optimal distribution of bins to
  // procs can be obtained
  // -------------------------------------------------------------------------
  // nodes, that are owned by a proc, are distributed to the bins of this proc
  std::vector<std::map<int, std::vector<int> > > nodesinbin(1);
  if(partition || repartition)
  {
    // get optimal row distribution of bins to procs
    rowbins_ = binning_->WeightedDistributionOfBinsToProcs(discret_vec,
                                                           disnp,
                                                           nodesinbin,
                                                           repartition);
  }

  // -------------------------------------------------------------------------
  // build element rowmap of binning discretization according to rowbins
  // -------------------------------------------------------------------------
  // completely new bins are created
  if(partition)
  {
    // extract noderowmap because it will be called Reset() after adding elements
    Teuchos::RCP<Epetra_Map> noderowmap =
        Teuchos::rcp(new Epetra_Map(*bindis_->NodeRowMap()));

    // delete old bins
    bindis_->DeleteElements();
    // loop over all bins on this proc, bindis than contains bins as elements
    for(int i=0; i<rowbins_->NumMyElements(); i++)
    {
      const int gid = rowbins_->GID(i);
      Teuchos::RCP<DRT::Element> bin =
          DRT::UTILS::Factory("MESHFREEMULTIBIN","dummy",gid,myrank_);
      bindis_->AddElement(bin);
    }

    // -----------------------------------------------------------------------
    // now node (=crosslinker) to bin (=element) relation needs to be
    // established in binning discretization. Therefore some nodes need to
    // change their owner according to the bins owner they are in
    // -----------------------------------------------------------------------
    // i) create a list of homeless particles (owned by this proc) that do not
    //    reside in a bin owned by this proc
    std::list<Teuchos::RCP<DRT::Node> > homelesscrosslinker;
    for (int lid = 0; lid < noderowmap->NumMyElements(); ++lid)
    {
      DRT::Node* node = bindis_->gNode(noderowmap->GID(lid));
      const double* currpos = node->X();
      binning_->PlaceNodeCorrectly(
          Teuchos::rcp(node,false), currpos, homelesscrosslinker);
    }
    // ii) round robin loop to assign homeless crosslinker to correct bin and owner
    binning_->FillParticlesIntoBinsRoundRobin(homelesscrosslinker);

  }
  // bin gids are the same
  else if(repartition)
  {
    // export row elements/bins to new layout
    bindis_->ExportRowElements(*rowbins_);

    // export row nodes to new layout
    // create a set of row crosslinker IDs for each proc
    std::set<int> crosslinker;
    for (int lid=0; lid<rowbins_->NumMyElements(); ++lid)
    {
      DRT::Element* bin = bindis_->gElement(rowbins_->GID(lid));
      const int* crosslinkerids = bin->NodeIds();
      for(int icl=0; icl<bin->NumNode(); ++icl)
        crosslinker.insert(crosslinkerids[icl]);
    }
    // copy crosslinkergids to a vector and create crosslinkerrowmap
    std::vector<int> rowcrosslinker(crosslinker.begin(),crosslinker.end());
    Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(
        -1,(int)rowcrosslinker.size(),&rowcrosslinker[0],0,bindis_->Comm()));

    // place all nodes on the correct processor
    bindis_->ExportRowNodes(*noderowmap);
  }
  else if (transfer)
  {
    // transfer crosslinker to their new bins
    binning_->TransferParticles();
  }
  else
  {
    dserror("You should not be here");
  }

  // call fill complete to build new node row map
  if(cltoclinteraction)
    bindis_->FillComplete(false,false,false);

  // -------------------------------------------------------------------------
  // each owner of a bin gets owner of the nodes (of the structure discret) this
  // bin contains. All other nodes of elements, of which proc is owner of at
  // least one node, are ghosted.
  // -------------------------------------------------------------------------

  // as extended ghosting is applied to discret_vec, colmaps of standard ghosting
  // are given back separately if needed
  Teuchos::RCP<Epetra_Map> stdelecolmap;
  Teuchos::RCP<Epetra_Map> stdnodecolmap;
  binning_->StandardGhosting(ia_discret_,
                             rowbins_,
                             ia_disnp_,
                             stdelecolmap,
                             stdnodecolmap,
                             nodesinbin[0]);

#ifdef DEBUG
  // print distribution after standard ghosting
  if(myrank_ == 0 && (partition || repartition))
  {
    IO::cout<<"\n+--------------------------------------------------+"<<IO::endl;
    IO::cout<<"   parallel distribution with standard ghosting   " << IO::endl;
    IO::cout<<"+--------------------------------------------------+"<<IO::endl;
  }
  DRT::UTILS::PrintParallelDistribution(*ia_discret_);
#endif

  // ----------------------------------------------------------------------
  // extended ghosting means the following here: Each proc ghosts
  // all elements whose XAABB cuts a bin that is next to a bin that is
  // owned by a proc an not empty. All associated nodes are ghosted as well
  // ----------------------------------------------------------------------

  // to distribute eles to bins correctly, we need a column map for disnp here
  // as the owner of a element does not need to be the owner of all nodes,
  // but we need the information of all nodes for correct distribution
  // (note this vector is col format before extended ghosting, later we
  // a col vector based on maps of extended ghosting)
  Teuchos::RCP<Epetra_Vector> iadiscolnp =
      Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  LINALG::Export(*ia_disnp_, *iadiscolnp);

  // extend ghosting
  binning_->ExtendBinGhosting(ia_discret_,
                              rowbins_,
                              iadiscolnp,
                              extbintoelemap_,
                              cltoclinteraction,
                              true);

  // build element to bin map according to extended ghosting
  // extbintoelemap_ than contains for each col element the bins that are
  // somehow touched by an axis aligned bounding box around the element
  BuildEleToBinMap();

#ifdef DEBUG
  // print distribution after extended ghosting
  if(myrank_ == 0 && (partition || repartition))
  {
    IO::cout<<"\n+--------------------------------------------------+"<<IO::endl;
    IO::cout<<"   parallel distribution with extended ghosting   " << IO::endl;
    IO::cout<<"+--------------------------------------------------+"<<IO::endl;
  }
  DRT::UTILS::PrintParallelDistribution(*ia_discret_);
#endif

  // -------------------------------------------------------------------------
  // one layer ghosting of bins and particles, i.e. each proc ghosts the
  // bins (multibin elements, elemap) and particles (nodes, nodemap) of
  // its 26 neighbored bins -> final FillComplete() for maps included
  // -------------------------------------------------------------------------
  binning_->BuildBinColMap(rowbins_,extbintoelemap_);

#ifdef DEBUG
  // print distribution after extended ghosting
  if (myrank_ == 0 && (partition || repartition))
  {
    IO::cout<<"\n+--------------------------------------------------+"<<IO::endl;
    IO::cout<<"           particles after ghosting     " << IO::endl;
    IO::cout<<"+--------------------------------------------------+"<<IO::endl;
  }
  DRT::UTILS::PrintParallelDistribution(*bindis_);
#endif

  //--------------------------------------------------------------------------
  // update vectors and matrices as maps have changed during redistribution
  //--------------------------------------------------------------------------
  UpdateMaps();

  //--------------------------------------------------------------------------
  // assign beam elements to bins to enable crosslinker to beam interaction
  //--------------------------------------------------------------------------
  // in case have not build new bins we need to clear the content
  if(!partition)
    binning_->RemoveElesFromBins(bin_beamcontent_);
  // note: extbintoelemap contains all necessary information for this step as
  // there is no assigning to do for empty bins (not in extbintoelemap) on a proc
  // bins that just have crosslinker are contained as well but do no damage
  // loop over all bins and remove assigned wall elements
  binning_->AssignElesToBins(ia_discret_, extbintoelemap_, bin_beamcontent_);

#ifdef DEBUG
  // safety check if a crosslinker got lost
  if(bindis_->NumGlobalNodes()!=eval_statmech_ptr_->GetDataSMDynPtr()->NumCrosslink())
    dserror("A crosslinker got lost, something went wrong.");
#endif

  //--------------------------------------------------------------------------
  // reset transformation member variables (eg. exporter) by rebuilding
  // and provide new maps for coupling adapter
  //--------------------------------------------------------------------------
  siatransform_ = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  coupsia_->SetupCoupling(*ia_discret_, Discret());

  // -------------------------------------------------------------------------
  // in case we have double bonded crosslinker on myrank we have to check if
  // myrank is still owner of all its crosslinker (if not, set up double bond on
  // other proc that is now responsible)
  UpdateMyDoubleBondsAfterRedistribution();

  // that's it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::BuildEleToBinMap()
{
  CheckInit();

  // delete old map
  exteletobinmap_.clear();
  // loop over bins
  std::map<int, std::set<int> >::const_iterator biniter;
  for(biniter = extbintoelemap_.begin(); biniter!=extbintoelemap_.end(); ++biniter)
  {
    // loop over ele content of this bin
    std::set<int>::const_iterator eleiter;
    for(eleiter = biniter->second.begin(); eleiter!=biniter->second.end(); ++eleiter)
    {
      int elegid = *eleiter;
      int bingid = biniter->first;
      // assign bins to elements
      exteletobinmap_[elegid].insert(bingid);
    }
  }

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateMyDoubleBondsAfterRedistribution()
{
  CheckInit();

  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> > > dbcltosend;
  std::set<int> dbtoerase;

  // loop over all double bonds on myrank
  std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> >::const_iterator iter;;
  for(iter=doublebondcl_.begin(); iter!=doublebondcl_.end(); ++iter)
  {
    const int clgid = iter->first;
    DRT::Node* doublebondedcl_i = bindis_->gNode(clgid);

#ifdef DEBUG
    // safety check
    if(doublebondedcl_i==NULL)
      dserror("Crosslinker moved further than the bin length in one time step, "
              "this is not allowed (maybe increase cutoff radius). ");
#endif

    // check ownership
    int owner = doublebondedcl_i->Owner();
    if(owner != myrank_)
    {
#ifdef DEBUG
      if(not doublebondcl_.count(clgid))
        dserror("willing to delete not existing entry, something went wrong");
#endif
      dbtoerase.insert(clgid);
      dbcltosend[owner].push_back(iter->second);
    }
  }

  std::set<int>::const_iterator i;
  for(i=dbtoerase.begin(); i!=dbtoerase.end(); ++i)
    doublebondcl_.erase(*i);

  // add new double bonds
  CommunicateBeamToBeamLinkageAfterRedistribution(dbcltosend);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::InitializeBinDiscret()
{
  CheckInit();
  // -------------------------------------------------------------------------
  // add bin discretization manually (as it is not part of the input file)
  // to the global problem to enable output and other stuff (e.g. that particle
  // methods work here)
  // -------------------------------------------------------------------------
  Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(DiscretPtr()->Comm().Clone());
  bindis_ = Teuchos::rcp(new DRT::Discretization("particle" ,com));
  // create discretization writer - in constructor set into and owned by
  // corresponding discret
  bindis_->SetWriter(
      Teuchos::rcp(new IO::DiscretizationWriter(bindis_)));
  // this is necessary for particle algorithm constructor
  DRT::Problem::Instance()->AddDis("particle", bindis_);
  // -------------------------------------------------------------------------
  // set range for uniform random number generator
  // -------------------------------------------------------------------------
  DRT::Problem::Instance()->Random()->SetRandRange(0.0,1.0);
  // -------------------------------------------------------------------------
  // add nodes to bin discretization
  // only proc 0 is doing this (as the number of crosslinker is manageable)
  // -------------------------------------------------------------------------
  if(myrank_==0)
  {
    for (int i=0; i<eval_statmech_ptr_->GetDataSMDynPtr()->NumCrosslink(); i++)
    {
      // random reference position of crosslinker in bounding box
      std::vector<double> X(3);
      for (int dim=0; dim<3; dim++)
        X[dim] = eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()->at(dim)
                 * DRT::Problem::Instance()->Random()->Uni();

      Teuchos::RCP<DRT::Node> newcrosslinker =
          Teuchos::rcp(new CROSSLINKING::CrosslinkerNode(i,&X[0],myrank_));

      // todo: put next two calls in constructor of CrosslinkerNode?
      // init crosslinker data container
      Teuchos::RCP<CROSSLINKING::CrosslinkerNode> clnode =
          Teuchos::rcp_dynamic_cast<CROSSLINKING::CrosslinkerNode>(newcrosslinker);
      clnode->InitializeDataContainer();
      // set material
      // fixme: assign matnum to crosslinker type in crosslinker section in input file
      //       for now, only one linker type with matnum 2
      clnode->SetMaterial(2);

      // add crosslinker to bin discretization
      bindis_->AddNode(newcrosslinker);
    }
  }

  // set row map of newly created particle discretization
  bindis_->FillComplete(false,false,false);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::DiffuseCrosslinker()
{
  CheckInit();

  // get standard deviation and mean value for crosslinker that are free to
  // diffuse
  double standarddev = sqrt(eval_statmech_ptr_->GetDataSMDynPtr()->KT() /
                       (2*M_PI * eval_statmech_ptr_->GetDataSMDynPtr()->Eta()
                       * eval_statmech_ptr_->GetDataSMDynPtr()->RLink())
                       * (*GState().GetDeltaTime())[0]);
  double meanvalue = 0.0;
  // Set mean value and standard deviation of normal distribution
  DRT::Problem::Instance()->Random()->SetMeanVariance(meanvalue,standarddev);

  // transformation from row to column (after extended ghosting)
  Teuchos::RCP<Epetra_Vector> ia_discolnp =
      Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  LINALG::Export(*ia_disnp_, *ia_discolnp);

  // loop over all row crosslinker (beam binding status not touched here)
  const int numrowcl = bindis_->NumMyRowNodes();
  for(int rowcli=0; rowcli<numrowcl; ++rowcli)
  {
    // get current linker
    CROSSLINKING::CrosslinkerNode* crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(bindis_->lRowNode(rowcli));
    const int clcollid = crosslinker_i->LID();
    CrosslinkerData& cldata_i = crosslinker_data_[clcollid];

    // different treatment according to number of bonds a crosslinker has
    switch(cldata_i.clnumbond)
    {
      case 0:
      {
        // crosslinker has zero bonds, i.e. is free to diffuse according to
        // brownian dynamics
        DiffuseUnboundCrosslinker(crosslinker_i);
        break;
      }
      case 1:
      {
        // get clbspot that is currently bonded
        int occbspotid = 0;
        GetOccupiedClBspot(occbspotid, cldata_i.clbspots);

        // get current position of binding spot of filament partner
        // note: we can not use our beam data container, as bspot position is not current position (as this
        // is the result of a sum, you can not have a reference to that)
        const int elegid = cldata_i.clbspots[occbspotid].first;

#ifdef DEBUG
        // safety check
        const int colelelid = ia_discret_->ElementColMap()->LID(elegid);
        if(colelelid<0)
          dserror("Crosslinker has %i bonds but his binding partner with gid %i "
                  "is \nnot ghosted/owned on proc %i (owner of crosslinker)",cldata_i.clnumbond,elegid,myrank_);
#endif

        DRT::ELEMENTS::Beam3Base* ele =
            dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ia_discret_->gElement(elegid));

        // get current position of filament binding spot
        LINALG::Matrix<3,1> bbspotpos;
        std::vector<double> eledisp;
        GetCurrentElementDis(ele,ia_discolnp,eledisp);
        ele->GetPosOfBindingSpot(bbspotpos,eledisp,cldata_i.clbspots[occbspotid].second,
            *(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()));

        SetCrosslinkerPosition(crosslinker_i, bbspotpos);

        break;
      }

      // crosslinker has two bonds (cl gets current mid position between the filament
      // binding spot it is attached to)
      case 2:
      {
        // -----------------------------------------------------------------
        // partner one
        // -----------------------------------------------------------------
        int elegid = cldata_i.clbspots[0].first;

#ifdef DEBUG
        // safety check
        int colelelid = ia_discret_->ElementColMap()->LID(elegid);
        if(colelelid<0)
          dserror("Crosslinker has %i bonds but his binding partner with gid %i "
                  "is not \nghosted/owned on proc %i (owner of crosslinker)",cldata_i.clnumbond,elegid,myrank_);
#endif

        DRT::ELEMENTS::Beam3Base* ele =
            dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ia_discret_->gElement(elegid));

        // get current position of filament binding spot
        LINALG::Matrix<3,1> bbspotposone;
        std::vector<double> eledisp;
        GetCurrentElementDis(ele,ia_discolnp,eledisp);
        ele->GetPosOfBindingSpot(bbspotposone,eledisp,cldata_i.clbspots[0].second,
            *(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()));

        // -----------------------------------------------------------------
        // partner two
        // -----------------------------------------------------------------
        elegid = cldata_i.clbspots[1].first;

#ifdef DEBUG
        // safety check
        colelelid = ia_discret_->ElementColMap()->LID(elegid);
        if(colelelid<0)
          dserror("Crosslinker has %i bonds but his binding partner with gid %i "
                  "is \nnot ghosted/owned on proc %i (owner of crosslinker)",cldata_i.clnumbond,elegid,myrank_);
#endif

        ele = dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ia_discret_->gElement(elegid));

        // get current position of filament binding spot
        LINALG::Matrix<3,1> bbspotpostwo;
        GetCurrentElementDis(ele,ia_discolnp,eledisp);
        ele->GetPosOfBindingSpot(bbspotpostwo,eledisp,cldata_i.clbspots[1].second,
            *(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()));

        // fixme: this is wrong if on binding partner leaves the periodic bounding box
        // cl gets current mid position between the filament binding spots it is attached to
        cldata_i.clpos.Update(0.5,bbspotposone,0.5,bbspotpostwo);
        SetCrosslinkerPosition(crosslinker_i, cldata_i.clpos);
        //fixme:
        static LINALG::Matrix<3,1> dist_vec;
        dist_vec.Update(1.0, bbspotposone, -1.0, cldata_i.clpos);
        const double distance = dist_vec.Norm2();
        if(abs(distance)>0.3*(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()->at(0)))
          dserror("You should have fixed this a long time ago (one binding partner left the domain");

        break;
      }
      default:
      {
        dserror("Unrealistic number %i of bonds for a crosslinker.", cldata_i.clnumbond);
        break;
      }
    }
  }

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::DiffuseUnboundCrosslinker(
    CROSSLINKING::CrosslinkerNode* crosslinker) const
{
  CheckInit();

  // diffuse crosslinker according to brownian dynamics
  std::vector<double> randvec;
  int count = 3;
  DRT::Problem::Instance()->Random()->Normal(randvec,count);
  // note: check for compliance with periodic boundary conditions is
  // done during crosslinker transfer in UpdateBinStrategy()
  crosslinker->ChangePos(randvec);

  // that's it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateCrosslinking()
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::UpdateCrosslinking");

  // -------------------------------------------------------------------------
  // 1) crosslinker diffusion:
  //    - according to browninan dyn for free cl
  //    - according to beams for single and double bonded
  // -------------------------------------------------------------------------
  // note: it is possible that a crosslinker leaves the computational domain
  // at this point, it gets shifted back in in UpdateBinStrategy which is called
  // right after this
  DiffuseCrosslinker();

  // -------------------------------------------------------------------------
  // 2) update parallel distribution for bindis and intacdis as crosslinker
  //    and beams have moved during the timestep. As we call export in the
  //    following as well, ghosted crosslinker and beams get current information
  //    regarding their binding status, position and so on.
  // -------------------------------------------------------------------------
  // get current beam positions based on maps of intactdis
  UpdateDofMapOfVector(ia_discret_, ia_disnp_, GState().GetMutableDisNp());

  // todo: do we want this ? (how to choose between (re-) and partitioning?)
  // do a complete new weighted partitioning every 100th step
  int dlbevery = eval_statmech_ptr_->GetDataSMDynPtr()->DynLoadBalanceEvery();
  if(dlbevery)
  {
    if(GState().GetStepN()%dlbevery == 0 and numproc_ != 1)
      UpdateBinStrategy(false,true);
  }
  else
  {
    // transfer beams and crosslinker
    UpdateBinStrategy();
  }

  // erase temporary data container as both distributions as well the col data
  // might have changed
  crosslinker_data_.clear();
  beam_data_.clear();

  // -------------------------------------------------------------------------
  // 3) now we manage binding events, this includes:
  //    - find potential binding events on each myrank
  //    - make a decision by asking other procs
  //    - set bonds and adapt states accordingly
  // -------------------------------------------------------------------------
  // intended bonds row cl to row ele
  std::map<int, BindEventData > mybonds; // key is clgid
  // intended bond col crosslinker to row element
  std::map<int, std::vector<BindEventData> > undecidedbonds; // key is owner!=myrank

  // fill binding event maps
  LookForBindingEvents(mybonds,undecidedbonds);

  // bind events where myrank only owns the elements, cl are taken care of by their owner
  std::map<int, BindEventData > myelebonds; // key is clgid

  // now each row owner of a linker gets requests, makes a random decision and
  // informs back its requesters
  ManageBindingInParallel(mybonds,undecidedbonds,myelebonds);

  // actual update of binding states is done here
  BindMyCrosslinkerAndElements(mybonds,myelebonds);

  // -------------------------------------------------------------------------
  // 4) unbinding events if probability check is passed
  // -------------------------------------------------------------------------
  UnBindCrosslinker();

  // -------------------------------------------------------------------------
  // 5) check which procs need to be informed about stiff and force contributions
  //    of double bonds on myrank
  // -------------------------------------------------------------------------
  GetProcsInvolvedInMyDoubleBonds();

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::LookForBindingEvents(
    std::map<int, BindEventData >&              mybonds,
    std::map<int, std::vector<BindEventData> >& undecidedbonds)
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::LookForBindingEvents");

  // gather data for all column crosslinker initially
  const int numcolcl = bindis_->NumMyColNodes();
  crosslinker_data_.resize(numcolcl);
  PreComputeCrosslinkerData(numcolcl);

  // gather data for all column beams
  const int numcolbeams = ia_discret_->NumMyColElements();
  beam_data_.resize(numcolbeams);
  PreComputeBeamData(numcolbeams);

  // store bins, which have already been examined
  std::set<int> examinedbins;
  // loop over all column crosslinker in random order
  // create random order of indices
  std::vector<int> rordercolcl = BIOPOLYNET::UTILS::Permutation(numcolcl);
  std::vector<int>::const_iterator icl;
  for(icl=rordercolcl.begin(); icl!=rordercolcl.end(); ++icl)
  {
    DRT::Node *currcrosslinker = bindis_->lColNode(*icl);

    // get bin that contains this crosslinker (can only be one)
    DRT::Element* CurrentBin = currcrosslinker->Elements()[0];
    const int currbinId = CurrentBin->Id();

    // if a bin has already been examined --> continue with next crosslinker
    if(examinedbins.find(currbinId) != examinedbins.end())
      continue;
    //else: bin is examined for the first time --> new entry in examinedbins_
    else
      examinedbins.insert(currbinId);

    // get neighboring bins
    // note: interaction distance cl to beam needs to be smaller than bin size
    std::vector<int> neighboring_binIds;
    neighboring_binIds.reserve(27);
    // do not check on existence here -> shifted to GetBinContent
    binning_->GetNeighborAndOwnBinIds(currbinId,neighboring_binIds);

    // get set of neighbouring beam elements (i.e. elements that somehow touch nb bins)
    // as explained above, we only need row elements
    std::set<DRT::Element*> neighboring_beams;
    binning_->GetBinContent(neighboring_beams,bin_beamcontent_,neighboring_binIds,true);

    // get all crosslinker in current bin
    DRT::Node **ClInCurrentBin = CurrentBin->Nodes();
    const int numcrosslinker = CurrentBin->NumNode();

    // obtain random order in which crosslinker are addressed
    std::vector<int> randorder = BIOPOLYNET::UTILS::Permutation(numcrosslinker);

    // loop over all crosslinker in CurrentBin in random order
    std::vector<int>::const_iterator randcliter;
    for(randcliter=randorder.begin(); randcliter!=randorder.end(); ++randcliter)
    {
      // get random crosslinker in current bin
      DRT::Node *crosslinker_i = ClInCurrentBin[*randcliter];
      // get all potential binding events on myrank
      PrepareBinding(crosslinker_i,neighboring_beams,mybonds,undecidedbonds);
    }
  }

  // that is it
  return;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::PrepareBinding(
    DRT::Node*                                  crosslinker_i,
    const std::set<DRT::Element*>&              neighboring_beams,
    std::map<int, BindEventData >&              mybonds,
    std::map<int, std::vector<BindEventData> >& undecidedbonds)
{
  CheckInit();

  // -------------------------------------------------------------------------
  // get precomputed data of crosslinker i
  // -------------------------------------------------------------------------
  const int clcollid = crosslinker_i->LID();
  const CrosslinkerData& cldata_i = crosslinker_data_[clcollid];

  // -------------------------------------------------------------------------
  // 1. criterion: in case crosslinker is double bonded, we can leave here
  // -------------------------------------------------------------------------
  if(cldata_i.clnumbond==2) return;

  // spherical shell around crosslinker center where binding spots have to lie
  // for potential binding event
  // todo: this needs to go crosslinker material
  double rmin = (eval_statmech_ptr_->GetDataSMDynPtr()->RLink()
               - eval_statmech_ptr_->GetDataSMDynPtr()->DeltaRLink()) / 2.0;
  double rmax = (eval_statmech_ptr_->GetDataSMDynPtr()->RLink()
               + eval_statmech_ptr_->GetDataSMDynPtr()->DeltaRLink()) / 2.0;
  // probability with which a crosslinker is established between crosslink
  // molecule and neighbor binding spot
  const double konstart = eval_statmech_ptr_->GetDataSMDynPtr()->KOnStart();
  double plink = 1.0 - exp( -(*GState().GetDeltaTime())[0]* konstart);

  // -------------------------------------------------------------------------
  // look for potential interaction of crosslinker i and a binding spot of an
  // element, i.e. distance \Delta = R +/- \Delta R
  // -------------------------------------------------------------------------
  // loop over all neighboring beam elements in random order (keep in mind
  // we are only looping over row elements)
  std::vector<DRT::Element*> beamvec(neighboring_beams.begin(),neighboring_beams.end());
  const int numbeams = beamvec.size();
  std::vector<int> randorder = BIOPOLYNET::UTILS::Permutation(numbeams);
  std::vector<int> ::const_iterator randiter;
  for(randiter=randorder.begin(); randiter!=randorder.end();  ++randiter)
  {
    // get neighboring (nb) beam element
    DRT::ELEMENTS::Beam3Base* nbbeam =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(beamvec.at(*randiter));

#ifdef DEBUG
      if(nbbeam == NULL)
        dserror("Dynamic cast to beam3base failed");
#endif

    // -----------------------------------------------------------------------
    // get pre computed data of current nbbeam
    // -----------------------------------------------------------------------
    const BeamData& beamdata_i = beam_data_[nbbeam->LID()];

    // 2. criterion:
    // exclude binding of a single bonded crosslinker in close proximity on the
    // same filament (i.e. element cloud of old element binding partner is excluded)
    if(cldata_i.clnumbond==1)
    {
      if(CheckCrosslinkOfAdjacentElements(nbbeam, cldata_i.bnodegids_))
        continue;
    }

    // loop over all binding spots of current element in random order
    std::vector<int> randbspot = BIOPOLYNET::UTILS::Permutation(beamdata_i.bbspotstatus.size());
    std::vector<int> ::const_iterator rbspotiter;
    for(rbspotiter=randbspot.begin(); rbspotiter!=randbspot.end(); ++rbspotiter)
    {
      // get local number of binding spot in element
      const int locnbspot = *rbspotiter;

      // -----------------------------------------------------------------------
      // we are now doing some checks if a binding event could happen
      // -----------------------------------------------------------------------
      {
        // 3. criterion:
        // first check if binding spot is free, if not, check next bspot on curr ele
        // note: bspotstatus in bonded case holds cl gid, otherwise -1 (meaning free)
        if(beamdata_i.bbspotstatus.at(locnbspot) != -1)
          continue;

        // 4. criterion:
        // if binding spot not in binding range, continue with next binding spot
        // compute distance between current binding spot and center (if free) or one
        // end (if single bonded) of crosslinker i
        const LINALG::Matrix<3,1> currbbspos = beamdata_i.bbspotpos.at(locnbspot);
        static LINALG::Matrix<3,1> dist_vec;
        dist_vec.Update(1.0, currbbspos, -1.0, cldata_i.clpos);
        const double distance = dist_vec.Norm2();

        if(distance>rmax || distance<rmin)
          continue;

        // 5. criterion:
        // a crosslink is set if and only if it passes a probability check
        if(DRT::Problem::Instance()->Random()->Uni() > plink)
          continue;
      }

      // ---------------------------------------------------------------------
      // if we came this far, we can add this potential binding event to its
      // corresponding map
      // ---------------------------------------------------------------------
      BindEventData bindeventdata;
      bindeventdata.clgid = crosslinker_i->Id();
      bindeventdata.elegid = nbbeam->Id();
      bindeventdata.bspotlocn = locnbspot;
      bindeventdata.requestproc = myrank_;
      // this is default, is changed if owner of cl has something against it
      bindeventdata.permission = true;

      // in case myrank is owner, we add it to the mybonds map
      if(cldata_i.clowner == myrank_)
      {
        mybonds[bindeventdata.clgid] = bindeventdata;
      }
      // myrank is not owner, we add it to the map of events that need to be
      // communicated to make a decision
      else
      {
        undecidedbonds[cldata_i.clowner].push_back(bindeventdata);
      }

      // as we allow only one binding event for each cl in one time step,
      // we are done here, if we made it so far (i.e met 1., 2., 3., 4. and 5. criterion)
      return;
    }
  }

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Crosslinking::CheckCrosslinkOfAdjacentElements(
    DRT::Element* nbbeam,
    const std::vector<std::pair<int, int> >& bnodegids) const
{
  CheckInit();

  int occbspotid = 0;
  GetOccupiedClBspot(occbspotid, bnodegids);

  // check if two considered eles share nodes
  for (int i=0; i<2; i++)
    if(nbbeam->NodeIds()[i]==bnodegids[occbspotid].first ||
       nbbeam->NodeIds()[i]==bnodegids[occbspotid].second)
      return true;

  // that is it
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::BindMyCrosslinkerAndElements(
    std::map<int, BindEventData >& mybonds,
    std::map<int, BindEventData >& myelebonds)
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::BindCrosslinker");

  // map key is crosslinker gid to be able to uniquely address one entry over all procs
  std::map<int, NewDoubleBonds> mynewdbondcl;

  // myrank owner of crosslinker and most elements
  BindMyCrosslinker(mybonds, mynewdbondcl);

  // myrank only owner of current binding partner ele
  BindMyElements(myelebonds);

  // setup new double bonds and insert them in doublebondcl_
  SetupNewDoubleBonds(mynewdbondcl, true);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::BindMyCrosslinker(
    const std::map<int, BindEventData >& mybonds,
    std::map<int, NewDoubleBonds>&       mynewdbondcl)
{
  CheckInit();

  std::map<int, BindEventData>::const_iterator cliter;
  for(cliter=mybonds.begin(); cliter!=mybonds.end(); ++cliter)
  {
    // get binding event data
    BindEventData binevdata = cliter->second;

    // get current linker
    const int clcollid = bindis_->NodeColMap()->LID(cliter->first);
#ifdef DEBUG
    // safety checks
    if(clcollid<0)
      dserror("Crosslinker not even ghosted, should be owned here.");
#endif
    CROSSLINKING::CrosslinkerNode *crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(bindis_->lColNode(clcollid));
    // get crosslinker data
    CrosslinkerData& cldata_i = crosslinker_data_[clcollid];

    // get beam data
    const int colelelid = ia_discret_->ElementColMap()->LID(binevdata.elegid);
#ifdef DEBUG
    // safety checks
    if(colelelid<0)
      dserror("Element not ghosted.");
#endif
    DRT::ELEMENTS::Beam3Base* beamele_i =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ia_discret_->lColElement(colelelid));
    BeamData& beamdata_i = beam_data_[colelelid];

#ifdef DEBUG
    // safety checks
    if(cliter->first!=binevdata.clgid)
      dserror("Map key does not match crosslinker gid of current binding event.");

    const int& clowner_i = cldata_i.clowner;
    if(clowner_i!=myrank_)
      dserror("Only row owner of crosslinker is changing its status");

    if(colelelid<0)
      dserror("Binding element partner of current row crosslinker is not ghosted, "
              "this must be the case though.");
#endif

    // -------------------------------------------------------------------------
    // different treatment according to number of bonds crosslinker had before
    // this binding event
    // -------------------------------------------------------------------------
    switch(cldata_i.clnumbond)
    {
      // crosslinker with zero bonds before this binding event
      case 0:
      {
        // -----------------------------------------------------------------
        // update crosslinker status
        // -----------------------------------------------------------------
        // store gid and bspot local number of this element, first binding spot
        // always bonded first
        cldata_i.clbspots[0].first  = binevdata.elegid;
        cldata_i.clbspots[0].second = binevdata.bspotlocn;
        crosslinker_i->ClData()->SetClBSpotStatus(cldata_i.clbspots);

        // store gid of first and second node of new binding partner
        cldata_i.bnodegids_[0].first  = beamele_i->NodeIds()[0];
        cldata_i.bnodegids_[0].second = beamele_i->NodeIds()[1];
        crosslinker_i->ClData()->SetBeamNodeGids(cldata_i.bnodegids_);

        // update number of bonds
        cldata_i.clnumbond = 1;
        crosslinker_i->ClData()->SetNumberOfBonds(cldata_i.clnumbond);

        // update position
        cldata_i.clpos = beamdata_i.bbspotpos.at(binevdata.bspotlocn);
        SetCrosslinkerPosition(crosslinker_i, cldata_i.clpos);

        // -----------------------------------------------------------------
        // update beam status
        // -----------------------------------------------------------------
        // store crosslinker gid in status of beam binding spot if myrank
        // is owner of beam
        if(beamdata_i.bowner == myrank_)
        {
          beamdata_i.bbspotstatus.at(binevdata.bspotlocn) = binevdata.clgid;
          beamele_i->SetBindingSpotStatus(beamdata_i.bbspotstatus);
        }

#ifdef DEBUG
        // safety check
        if(not (cldata_i.clbspots[1].first < 0))
          dserror("Numbond does not fit to clbspot vector.");
#endif

        break;
      }
      // crosslinker with one bond before this binding event
      case 1:
      {
        // get clbspot that is currently bonded
        int occbspotid = 0;
        GetOccupiedClBspot(occbspotid, cldata_i.clbspots);
        int freebspotid = 1;
        if(occbspotid==1)
          freebspotid = 0;

        // -----------------------------------------------------------------
        // update crosslinker status
        // -----------------------------------------------------------------
        // store gid and bspot local number of this element
        cldata_i.clbspots[freebspotid].first = binevdata.elegid;
        cldata_i.clbspots[freebspotid].second = binevdata.bspotlocn;
        crosslinker_i->ClData()->SetClBSpotStatus(cldata_i.clbspots);

        // store gid of first and second node of this element
        cldata_i.bnodegids_[freebspotid].first = beamele_i->NodeIds()[0];
        cldata_i.bnodegids_[freebspotid].second = beamele_i->NodeIds()[1];
        crosslinker_i->ClData()->SetBeamNodeGids(cldata_i.bnodegids_);

        // update number of bonds
        cldata_i.clnumbond = 2;
        crosslinker_i->ClData()->SetNumberOfBonds(cldata_i.clnumbond);

        // update position
        // fixme
        cldata_i.clpos.Update(0.5,beamdata_i.bbspotpos.at(cldata_i.clbspots[freebspotid].second),0.5);
        SetCrosslinkerPosition(crosslinker_i,cldata_i.clpos);

        static LINALG::Matrix<3,1> dist_vec;
        dist_vec.Update(1.0, beamdata_i.bbspotpos.at(cldata_i.clbspots[freebspotid].second), -1.0, cldata_i.clpos);
        const double distance = dist_vec.Norm2();
        if(abs(distance)>0.3*(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()->at(0)))
          dserror("You should have fixed this a long time ago (one binding partner left the domain");

        // create double bond cl data
        NewDoubleBonds dbondcl;
        dbondcl.id = binevdata.clgid;
        dbondcl.eleids.push_back(cldata_i.clbspots[freebspotid]);
        dbondcl.eleids.push_back(cldata_i.clbspots[occbspotid]);

        // insert pair in mypairs
        mynewdbondcl[dbondcl.id] = dbondcl;

        // first check if myrank is owner of element of current binding event
        // (additionally to being owner of cl)
        if(beamdata_i.bowner == myrank_)
        {
          // update beam data
          // store crosslinker gid in status of beam binding spot
          beamdata_i.bbspotstatus.at(binevdata.bspotlocn) = binevdata.clgid;
          beamele_i->SetBindingSpotStatus(beamdata_i.bbspotstatus);
        }

        break;
      }
      default:
      {
        dserror("You should not be here, crosslinker has unrealistic number "
                "%i of bonds.", cldata_i.clnumbond);
        break;
      }
    }
  }

  //that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::BindMyElements(
    const std::map<int, BindEventData >& myelebonds)
{
  CheckInit();

  /*
   * 1| 2__|2  or  1| 3__|2  or  1| 2__|1
   * 1|    |2      1|    |2      1|    |1
   * legend: | = beam; __= cl; 2,3 = owner; 1=myrank
   */
  // loop through all binding events
  std::map<int, BindEventData>::const_iterator cliter;
  for(cliter=myelebonds.begin(); cliter!=myelebonds.end(); ++cliter)
  {
    // get binding event data
    BindEventData binevdata = cliter->second;

    // ---------------------------------------------------------------------
    // get linker data
    // ---------------------------------------------------------------------
    const int clcollid = bindis_->NodeColMap()->LID(cliter->first);
    #ifdef DEBUG
    // safety check
    if(clcollid<0)
     dserror("Crosslinker needs to be ghosted, but this isn't the case.");
    #endif
    CrosslinkerData& cldata_i = crosslinker_data_[clcollid];

    // ---------------------------------------------------------------------
    // get beam data
    // ---------------------------------------------------------------------
    const int colelelid = ia_discret_->ElementColMap()->LID(binevdata.elegid);
    #ifdef DEBUG
    // safety check
    if(colelelid<0)
     dserror("element with gid %i not ghosted on proc %i",binevdata.elegid,myrank_);
    #endif
    // get beam element of current binding event
    DRT::ELEMENTS::Beam3Base* ele_i =
       dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ia_discret_->lColElement(colelelid));
    BeamData& beamdata_i = beam_data_[colelelid];

    #ifdef DEBUG
    // safety checks
    if(beamdata_i.bowner!=myrank_)
     dserror("Only row owner of element is allowed to change its status");
    const int& clowner_i = cldata_i.clowner;
    if(clowner_i==myrank_)
     dserror("myrank should not be owner of this crosslinker");
    #endif

    // different treatment according to number of bonds crosslinker has before
    // this binding event
    switch(cldata_i.clnumbond)
    {
      case 0:
      {
        // update beam data
        // store crosslinker gid in status of beam binding spot
        beamdata_i.bbspotstatus.at(binevdata.bspotlocn) = binevdata.clgid;
        ele_i->SetBindingSpotStatus(beamdata_i.bbspotstatus);
        break;
      }
      case 1:
      {
        // update beam data
        // store crosslinker gid in status of beam binding spot
        beamdata_i.bbspotstatus.at(binevdata.bspotlocn) = binevdata.clgid;
        ele_i->SetBindingSpotStatus(beamdata_i.bbspotstatus);

        break;
      }
      default:
      {
        dserror("You should not be here, crosslinker has to many bonds.");
        break;
      }
    }
  }

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::SetupNewDoubleBonds(
    const std::map<int, NewDoubleBonds>& mynewdbondcl,
    bool precomputed)
{
  CheckInit();

  // transformation from row to column
  Teuchos::RCP<Epetra_Vector> ia_discolnp =
      Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  LINALG::Export(*ia_disnp_, *ia_discolnp);

  std::map<int, NewDoubleBonds>::const_iterator iter;
  for (iter=mynewdbondcl.begin(); iter!=mynewdbondcl.end(); ++iter)
  {
    // init positions and triads
    std::vector<LINALG::Matrix<3,1> > pos(2);
    std::vector<LINALG::Matrix<3,3> > triad(2);

    const NewDoubleBonds& newdoublebond_i = iter->second;

    for(int i=0; i<2; ++i)
    {
      int elegid = newdoublebond_i.eleids[i].first;
      int locbspotnum = newdoublebond_i.eleids[i].second;
      DRT::Element* ele = ia_discret_->gElement(elegid);


#ifdef DEBUG
      // safety checks
      if(ele==NULL)
        dserror("Element not there on rank %i", myrank_);
#endif

      pos[i] = beam_data_[ele->LID()].bbspotpos.at(locbspotnum);
      triad[i] = beam_data_[ele->LID()].bbspottriad.at(locbspotnum);
    }

    // ToDo specify and pass material parameters for crosslinker element

    // create and initialize objects of beam-to-beam connections
    // Todo introduce enum for type of linkage (only linear Beam3r element possible so far)
    //      and introduce corresponding input parameter or even condition for mechanical
    //      links between beams in general
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> elepairptr =
      BEAMINTERACTION::BeamToBeamLinkage::Create();

    // finally initialize and setup object
    elepairptr->Init(iter->first,newdoublebond_i.eleids,pos,triad);
    elepairptr->Setup();

    // add to my double bonds
    doublebondcl_[elepairptr->Id()] = elepairptr;
  }

  // print some information
  if(mynewdbondcl.size()>0 or doublebondcl_.size()>0)
  {
    IO::cout(IO::standard) <<"\n************************************************"<<IO::endl;
    IO::cout(IO::standard) << "PID " << myrank_ << ": added " << mynewdbondcl.size()
        << " new db crosslinkers. Now have " << doublebondcl_.size() <<IO::endl;
  }

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::SetCrosslinkerPosition(
    DRT::Node* crosslinker,
    const LINALG::Matrix<3,1>& newclpos) const
{
  CheckInit();

  std::vector<double> newpos(3,0.0);
  for(int dim=0; dim<3; ++dim)
    newpos[dim] = newclpos(dim);
  crosslinker->SetPos(newpos);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::SetPositionOfNewlyFreeCrosslinker(
    DRT::Node* crosslinker,
    LINALG::Matrix<3,1>& clpos) const
{
  CheckInit();

  //generate vector in random direction of length R_LINK to "reset" crosslink
  // molecule position: it may now reenter or leave the bonding proximity
  // todo: does this make sense?
  LINALG::Matrix<3,1> cldeltapos_i;
  std::vector<double> randunivec(3);
  int count = 3;
  DRT::Problem::Instance()->Random()->Uni(randunivec, count);
  for (int dim=0; dim<3; ++dim)
    cldeltapos_i(dim) = randunivec[dim];
  cldeltapos_i.Scale(eval_statmech_ptr_->GetDataSMDynPtr()->RLink() / cldeltapos_i.Norm2());

  clpos.Update(1.0,cldeltapos_i,1.0);
  SetCrosslinkerPosition(crosslinker, clpos);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UnBindCrosslinker()
{
  CheckInit();

  // todo: this needs to go somewhere else
  //------------------------------------------------------------------------------
  // get current off-rate for crosslinkers
  double koff = 0;
  double starttime = eval_statmech_ptr_->GetDataSMDynPtr()->ActionTime()->at(1);
  const double timenp = GState().GetTimeNp();
  const double dt = (*GState().GetDeltaTime())[0];

  if (timenp <= starttime || (timenp>starttime && fabs(starttime)<dt/1e4))
    koff = eval_statmech_ptr_->GetDataSMDynPtr()->KOffStart();
  else
    koff = eval_statmech_ptr_->GetDataSMDynPtr()->KOffEnd();

  // probability with which a crosslink breaks up in the current time step
  double p_unlink = 1.0 - exp(-dt * koff);

  //------------------------------------------------------------------------------

  // data containing information about elements that need to be updated on
  // procs != myrank
  std::map<int, std::vector<UnBindEventData> > sendunbindevents;
  // elements that need to be updated on myrank
  std::vector<UnBindEventData> myrankunbindevents;
  // double bonded crosslinker ("elements") that needs to be deleted on myrank
  std::vector<int> doublebondstodelete;

  // loop over all row linker (in random order) and dissolve bond if probability
  // criterion is met
  /* note: we loop over all row crosslinker, i.e. myrank needs to update all
   * crosslinker information. As it possible that a row crosslinker is linked
   * to col element, we potentially need to communicate if such an element
   * needs to be updated*/
  const int numrowcl = bindis_->NumMyRowNodes();
  std::vector<int> rorderrowcl = BIOPOLYNET::UTILS::Permutation(numrowcl);
  std::vector<int>::const_iterator rowcli;
  for(rowcli=rorderrowcl.begin(); rowcli!=rorderrowcl.end(); ++rowcli)
  {
    // get current linker
    CROSSLINKING::CrosslinkerNode *crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(bindis_->lRowNode(*rowcli));

    // get linker data
    const int clcollid = crosslinker_i->LID();
    CrosslinkerData& cldata_i = crosslinker_data_[clcollid];

    // different treatment according to number of bonds of a crosslinker
    switch(cldata_i.clnumbond)
    {
      case 0:
      {
        // nothing to do here
        break;
      }
      case 1:
      {
        // if probability criterion is not met, we are done here
        if (DRT::Problem::Instance()->Random()->Uni() > p_unlink)
          break;

        // dissolve bond
        int occbspotid = 0;
        GetOccupiedClBspot(occbspotid, cldata_i.clbspots);

        // -----------------------------------------------------------------
        // update beam status (first check which rank is responsible)
        // -----------------------------------------------------------------

        // store unbinding event data
        UnBindEventData unbindevent;
        unbindevent.eletoupdate = cldata_i.clbspots[occbspotid];

        // owner of beam
        const int beamowner =
            ia_discret_->gElement(unbindevent.eletoupdate.first)->Owner();

        // check who needs to update the element status
        if(beamowner == myrank_)
          myrankunbindevents.push_back(unbindevent);
        else
          sendunbindevents[beamowner].push_back(unbindevent);

        // -----------------------------------------------------------------
        // update crosslinker status
        // -----------------------------------------------------------------
        // update binding status of linker
        cldata_i.clbspots[occbspotid].first = -1;
        cldata_i.clbspots[occbspotid].second = -1;
        crosslinker_i->ClData()->SetClBSpotStatus(cldata_i.clbspots);

        // store gid of first and second node of new binding partner
        cldata_i.bnodegids_[occbspotid].first = -1;
        cldata_i.bnodegids_[occbspotid].second = -1;
        crosslinker_i->ClData()->SetBeamNodeGids(cldata_i.bnodegids_);

        // update number of bonds
        cldata_i.clnumbond = 0;
        crosslinker_i->ClData()->SetNumberOfBonds(cldata_i.clnumbond);

        // update position of crosslinker
        SetPositionOfNewlyFreeCrosslinker(crosslinker_i,cldata_i.clpos);

        break;
      }
      case 2:
      {
        // get id of freed and still occupied bspot
        int freedbspotid = -1;
        int stayocc = -1;
        // loop through crosslinker bonds in random order
        std::vector<int> ro = BIOPOLYNET::UTILS::Permutation(cldata_i.clnumbond);
        std::vector<int>::const_iterator clbspotiter;
        for(clbspotiter=ro.begin(); clbspotiter!=ro.end(); ++clbspotiter)
        {
          // if probability criterion isn't met, go to next spot
          if (DRT::Problem::Instance()->Random()->Uni() > p_unlink)
            continue;

          // get id of freed and still occupied bspot
          freedbspotid = *clbspotiter;
          stayocc = 0;
          if(freedbspotid == 0)
            stayocc = 1;

          // add crosslinker to list of deleted double bonds
          doublebondstodelete.push_back(crosslinker_i->Id());

          // -----------------------------------------------------------------
          // update beam status (first check which rank is responsible)
          // -----------------------------------------------------------------
          // initialize ubindevent
          UnBindEventData unbindevent;
          unbindevent.eletoupdate = cldata_i.clbspots[freedbspotid];

          // owner of beam element bond that gets dissolved
          const int freedbeamowner =
              ia_discret_->gElement(cldata_i.clbspots[freedbspotid].first)->Owner();

          if(freedbeamowner == myrank_)
            myrankunbindevents.push_back(unbindevent);
          else
            sendunbindevents[freedbeamowner].push_back(unbindevent);

          // -----------------------------------------------------------------
          // update crosslinker status
          // -----------------------------------------------------------------
          // reset binding status of freed crosslinker binding spot
          cldata_i.clbspots[freedbspotid].first = -1;
          cldata_i.clbspots[freedbspotid].second = -1;
          crosslinker_i->ClData()->SetClBSpotStatus(cldata_i.clbspots);

          // store gid of first and second node of new binding partner
          cldata_i.bnodegids_[freedbspotid].first = -1;
          cldata_i.bnodegids_[freedbspotid].second = -1;
          crosslinker_i->ClData()->SetBeamNodeGids(cldata_i.bnodegids_);

          // update number of bonds
          cldata_i.clnumbond = 1;
          crosslinker_i->ClData()->SetNumberOfBonds(cldata_i.clnumbond);

          // update postion
          const int collidoccbeam =
              ia_discret_->ElementColMap()->LID(cldata_i.clbspots[stayocc].first);
#ifdef DEBUG
          // safety check
          if(collidoccbeam<0)
            dserror("element with gid %i not ghosted on proc %i",cldata_i.clbspots[stayocc].first,myrank_);
#endif
          BeamData& beamdata_i = beam_data_[collidoccbeam];
          cldata_i.clpos = beamdata_i.bbspotpos.at(cldata_i.clbspots[stayocc].second);
          SetCrosslinkerPosition(crosslinker_i, cldata_i.clpos);

          // we only want to dissolve one bond per timestep, therefore we go to
          // next crosslinker if we made it so far (i.e. a bond got dissolved)
          break;
        }

        break;
      }
      default:
      {
        dserror("Unrealistic number %i of bonds for a crosslinker.", cldata_i.clnumbond);
        break;
      }
    }
  }

  // communicate which elements that need to be updated on rank!=myrank
  CommunicateCrosslinkerUnbinding(sendunbindevents,myrankunbindevents);

  // update binding status of beam binding partners on myrank
  UpdateBeamBindingStatusAfterUnbinding(myrankunbindevents);

  // update list of double bonds on myrnak
  UpdateDoubleBondsAfterUbinding(doublebondstodelete);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::GetProcsInvolvedInMyDoubleBonds()
{
  CheckInit();

  // loop over all double bonds on myrank
  std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> >::const_iterator iter;
  for(iter=doublebondcl_.begin(); iter!=doublebondcl_.end(); ++iter)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> elepairptr = iter->second;

    // get elements
    DRT::Element* ele1 = ia_discret_->gElement(elepairptr->GetEleGid(0));
    DRT::Element* ele2 = ia_discret_->gElement(elepairptr->GetEleGid(1));

    // those procs need to be informed
    std::set<int> procs;
    procs.clear();

    // get owner of element dofs
    std::vector<int> lm, lmowner, lmstride;
    ele1->LocationVector(*ia_discret_,lm,lmowner,lmstride);

    // insert in proc list
    procs.insert(lmowner.begin(),lmowner.end());

    // same for second element
    ele2->LocationVector(*ia_discret_,lm,lmowner,lmstride);
    procs.insert(lmowner.begin(),lmowner.end());

    // delete myrank id of set
    procs.erase(myrank_);

    elepairptr->SetInvolvedProcs(procs);
  }

  // that's it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateBeamBindingStatusAfterUnbinding(
    const std::vector<UnBindEventData>& unbindevent )
{
  CheckInit();

  // loop through all unbinding events on myrank
  std::vector<UnBindEventData>::const_iterator iter;
  for(iter=unbindevent.begin(); iter!=unbindevent.end(); ++ iter)
  {
    // get data
    const int elegidtoupdate = iter->eletoupdate.first;
    const int bspotlocn = iter->eletoupdate.second;
    const int colelelid = ia_discret_->ElementColMap()->LID(elegidtoupdate);

#ifdef DEBUG
    // safety check
    if(colelelid<0)
      dserror("element with gid %i not owned by proc %i",elegidtoupdate,myrank_);
#endif

    // get beam element of current binding event
    DRT::ELEMENTS::Beam3Base* ele_i =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ia_discret_->lColElement(colelelid));

    BeamData& beamdata_i = beam_data_[colelelid];
    std::map<int, int>& bbspotstatus_i = beamdata_i.bbspotstatus;

    // update beam data
    bbspotstatus_i.at(bspotlocn) = -1;
    ele_i->SetBindingSpotStatus(bbspotstatus_i);
    }

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateDoubleBondsAfterUbinding(
    const std::vector<int>& doublebondstodelete)
{
  CheckInit();

  // loop through all unbinding events on myrank
  std::vector<int>::const_iterator iter;
  for(iter=doublebondstodelete.begin(); iter!=doublebondstodelete.end(); ++ iter)
  {
    // get data
    const int doublebondtoerase = *iter;

#ifdef DEBUG
    // safety check
    if(not doublebondcl_.count(doublebondtoerase))
      dserror("willing to delete not existing entry, something went wrong");
#endif

    // delete double bonded crosslinker
    doublebondcl_.erase(doublebondtoerase);

    }

  if(doublebondstodelete.size()>0)
  {
    IO::cout(IO::standard) << "PID " << myrank_ << ": deleted " << doublebondstodelete.size()
        << " db crosslinkers. Now have " << doublebondcl_.size() <<IO::endl;
    IO::cout(IO::standard) <<"************************************************\n"<<IO::endl;
  }

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::PreComputeCrosslinkerData(
    const int numcolcl)
{
  CheckInit();

#ifdef DEBUG
  // sanity check
  if(static_cast<int>(crosslinker_data_.size())!= bindis_->NumMyColNodes())
    dserror("temporary crosslinker data container has wrong size");
#endif

  for (int i=0; i<numcolcl; ++i)
  {
    // crosslinker i for which data will be collected
    CROSSLINKING::CrosslinkerNode *crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(bindis_->lColNode(i));

#ifdef DEBUG
      if(crosslinker_i == NULL)
        dserror("Dynamic cast to CrosslinkerNode failed");
#endif

    // col lid of this crosslinker
    const int collid = crosslinker_i->LID();
    // store date of crosslinker i according to column lid
    CrosslinkerData& cldata = crosslinker_data_[collid];

    // store positions
    for(int dim=0; dim<3; ++dim)
      cldata.clpos(dim) = crosslinker_i->X()[dim];
    // get current binding spot status of crosslinker
    cldata.clbspots = crosslinker_i->ClData()->GetClBSpotStatus();
    // get current binding spot status of crosslinker
    cldata.bnodegids_ = crosslinker_i->ClData()->GetBeamNodeGids();
    // get number of bonds
    cldata.clnumbond = crosslinker_i->ClData()->GetNumberOfBonds();
    // get type of crosslinker (i.e. its material)
    cldata.clmat = crosslinker_i->GetMaterial();
    // get owner
    cldata.clowner = crosslinker_i->Owner();

#ifdef DEBUG
      if(crosslinker_i->NumElement() != 1)
        dserror("More than one element for this crosslinker");
#endif
  }

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::PreComputeBeamData(
    const int numcolbeams)
{
  CheckInit();

#ifdef DEBUG
  // sanity check
  if(static_cast<int>(beam_data_.size())!= ia_discret_->NumMyColElements())
    dserror("temporary beam data container has wrong size");
#endif

  // transformation from row to column (after extended ghosting)
  Teuchos::RCP<Epetra_Vector> ia_discolnp =
      Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofColMap()));
  LINALG::Export(*ia_disnp_, *ia_discolnp);

  // loop over all column beams elements
  for (int i=0; i<numcolbeams; ++i)
  {
    // beam element i for which data will be collected
    DRT::ELEMENTS::Beam3Base* ele_i =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ia_discret_->lColElement(i));

#ifdef DEBUG
      if(ele_i == NULL)
        dserror("Dynamic cast to Beam3Base failed");
#endif

    // get elelid
    const int elelid = ele_i->LID();
    // store data
    BeamData& bdata = beam_data_[elelid];

    std::vector<double> eledisp;
    GetCurrentElementDis(ele_i,ia_discolnp,eledisp);

    // loop over all binding spots of current element
    const int numbbspot = static_cast<int>(ele_i->GetBindingSpotStatus().size());
    for(int j=0; j<numbbspot; ++j)
      GetPosAndTriadOfBindingSpot(ele_i,
                                  ia_discolnp,
                                  j,
                                  bdata.bbspotpos[j],
                                  bdata.bbspottriad[j],
                                  eledisp);

    // get status of beam binding spots
    bdata.bbspotstatus = ele_i->GetBindingSpotStatus();
    // get owner
    bdata.bowner = ele_i->Owner();
  }

  // that is it
  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::GetCurrentElementDis(
    const DRT::Element* ele,
    const Teuchos::RCP<Epetra_Vector> ia_discolnp,
    std::vector<double>& eledisp) const
{
  CheckInit();

  // clear
  eledisp.clear();

  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;

  ele->LocationVector(*ia_discret_,lm,lmowner,lmstride);
  DRT::UTILS::ExtractMyValues(*ia_discolnp,eledisp,lm);

  // that's it
  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::GetPosAndTriadOfBindingSpot(
    DRT::Element* ele,
    const Teuchos::RCP<Epetra_Vector> ia_discolnp,
    const int& locbspotnum,
    LINALG::Matrix<3,1>& bspotpos,
    LINALG::Matrix<3,3>& bspottriad) const
{
  CheckInit();

  std::vector<double> eledisp;
  GetCurrentElementDis(ele,ia_discolnp,eledisp);

  GetPosAndTriadOfBindingSpot(ele,ia_discolnp,locbspotnum,
      bspotpos,bspottriad,eledisp);

  // that's it
  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::GetPosAndTriadOfBindingSpot(
    DRT::Element* ele,
    const Teuchos::RCP<Epetra_Vector> ia_discolnp,
    const int& locbspotnum,
    LINALG::Matrix<3,1>& bspotpos,
    LINALG::Matrix<3,3>& bspottriad,
    std::vector<double>& eledisp) const
{
  CheckInit();

  // cast to beambase element
  DRT::ELEMENTS::Beam3Base* beamele =
      dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele);

#ifdef DEBUG
      if(beamele == NULL)
        dserror("Dynamic cast to beam3base failed");
#endif

  // get current position at binding spot xi
  beamele->GetPosOfBindingSpot(bspotpos,eledisp,locbspotnum,
      *(eval_statmech_ptr_->GetDataSMDynPtr()->PeriodLength()));

  // get current triad at binding spot xi
  beamele->GetTriadOfBindingSpot(bspottriad,eledisp,locbspotnum);

  // that's it
  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::GetOccupiedClBspot(
    int& occbspotid,
    const std::vector<std::pair<int, int> >& clbspots) const
{
  CheckInit();

  occbspotid = 0;
  // check, which clbspot is free
  if(clbspots[0].first < 0)
    occbspotid  = 1;

#ifdef DEBUG
    // safety check
  if(clbspots[occbspotid].first < 0)
    dserror("determined binding spot for singly bonded crosslinker is not occupied");
#endif

  // that's it
  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::ManageBindingInParallel(
    std::map<int, BindEventData >&              mybonds,
    std::map<int, std::vector<BindEventData> >& undecidedbonds,
    std::map<int, BindEventData >&              myelebonds) const
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::ManageBindingInParallel");

  // variable for safety check
  int numrecrequest;
  // exporter
  DRT::Exporter exporter(bindis_->Comm());

  // -------------------------------------------------------------------------
  // 1) each procs makes his requests and receives the request of other procs
  // -------------------------------------------------------------------------
  // store requested cl and its data
  std::map<int, std::vector<BindEventData> > requestedcl;
  CommunicateUndecidedBonds(exporter,undecidedbonds,numrecrequest,requestedcl);

  // -------------------------------------------------------------------------
  // 2) now myrank needs to decide which proc is allowed to set the requested
  //    link
  // -------------------------------------------------------------------------
  std::map<int, std::vector<BindEventData> > decidedbonds;
  DecideBindingInParallel(requestedcl,mybonds,decidedbonds);

  // -------------------------------------------------------------------------
  // 3) communicate the binding decisions made on myrank, receive decisions
  //    made for its own requests and create colbondmap accordingly
  // -------------------------------------------------------------------------
  int answersize = static_cast<int>(undecidedbonds.size());
  CommunicateDecidedBonds(exporter,decidedbonds,myelebonds,numrecrequest,answersize);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::CommunicateUndecidedBonds(
    DRT::Exporter& exporter,
    std::map<int, std::vector<BindEventData> >& undecidedbonds,
    int& numrecrequest,
    std::map<int, std::vector<BindEventData> >& requestedcl) const
{
  CheckInit();

  // -----------------------------------------------------------------------
  // unblocking send
  // -----------------------------------------------------------------------
  std::vector<MPI_Request> request;
  ISend(exporter,request,undecidedbonds);

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  std::vector<int> summedtargets;
  PrepareReceivingProcs(undecidedbonds,summedtargets);

  numrecrequest = summedtargets[myrank_];
  for(int rec=0; rec<numrecrequest; ++rec)
  {
      std::vector<char> rdata;
      int length = 0;
      int tag = -1;
      int from = -1;
      exporter.ReceiveAny(from,tag,rdata,length);
      if (tag != 1234)
       dserror("Received on proc %i data with wrong tag from proc %i", myrank_, from);

      // store received data
      std::vector<char>::size_type position = 0;
      while (position < rdata.size())
      {
        // ---- extract received data -----
        BindEventData reccldata;
        UnPack(position,rdata,reccldata);

        // create map holding all requests
        requestedcl[reccldata.clgid].push_back(reccldata);
      }

      if (position != rdata.size())
        dserror("Mismatch in size of data %d <-> %d",static_cast<int>(rdata.size()),position);
  }

  // wait for all communication to finish
  Wait(exporter,request,static_cast<int>(undecidedbonds.size()));

  // that's it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::CommunicateDecidedBonds(
    DRT::Exporter& exporter,
    std::map<int, std::vector<BindEventData> >& decidedbonds,
    std::map<int, BindEventData >&              myelebonds,
    const int& numrecrequest,
    const int& answersize) const
{
  CheckInit();

  // -----------------------------------------------------------------------
  // send back decisions for all requests that were made
  // -----------------------------------------------------------------------
  // store requested cl and its data
  std::vector<MPI_Request> request;

  // safety check
  // todo: if passed in release, do it only in debug mode
  if(static_cast<int>(decidedbonds.size()) != numrecrequest)
    dserror("Number of received requests %i unequal to number of answers %i",
        numrecrequest,decidedbonds.size());

  // unblocking send
  ISend(exporter,request,decidedbonds);

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  // store requested cl and its data
  for(int rec=0; rec<answersize; ++rec)
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from,tag,rdata,length);
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", myrank_, from);

    // store received data
    std::vector<char>::size_type position = 0;
    while (position < rdata.size())
    {
      // ---- extract received data -----
      BindEventData reccldata;
      UnPack(position,rdata,reccldata);

      // add binding events to new colbond map
      if(reccldata.permission)
        myelebonds[reccldata.clgid] = reccldata;
    }

    if (position != rdata.size())
      dserror("Mismatch in size of data %d <-> %d",static_cast<int>(rdata.size()),position);
  }

  // wait for all communication to finish
  Wait(exporter,request,static_cast<int>(decidedbonds.size()));

  // that's it
  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::DecideBindingInParallel(
    std::map<int, std::vector<BindEventData> >& requestedcl,
    std::map<int, BindEventData >&              mybonds,
    std::map<int, std::vector<BindEventData> >& decidedbonds) const
{
  CheckInit();

   std::map<int, std::vector<BindEventData> >::iterator cliter;
   // loop over all requested cl (note myrank is owner of these)
   for(cliter=requestedcl.begin(); cliter!=requestedcl.end(); ++cliter)
   {
     // check if myrank wants to bind this crosslinker
     const bool myrankbond = (mybonds.find(cliter->first) != mybonds.end());

     // ---------------------------------------------------------------------
     // if only one request and myrank does not want to bind this cl,
     // requesting proc gets the permission to do so
     // ---------------------------------------------------------------------
     if(static_cast<int>(cliter->second.size()) == 1 and not myrankbond)
     {
       // we send back the permission to the relevant proc, because myrank as row
       // owner of bspot needs to set the respective stuff for the element of this
       // binding event
       // note: permission = true was send as default, so this can be sent back
       // without changes
       decidedbonds[cliter->second[0].requestproc].push_back(cliter->second[0]);

       // insert this new binding event in map of myrank, because as row owner of
       // this cl he is responsible to set the respective stuff for the crosslinker
       // of this binding event
       mybonds[cliter->first] = cliter->second[0];

       // go to next crosslinker
       continue;
     }

     // ---------------------------------------------------------------------
     // in case number of requesting procs >1 for this cl or myrank wants to
     // set it itself
     // ---------------------------------------------------------------------
     int numrequprocs = cliter->second.size();
     if(myrankbond)
       numrequprocs += 1;

     // get random proc out of affected ones
     DRT::Problem::Instance()->Random()->SetRandRange(0.0,1.0);
     // todo: what if random number exactly = 1?
     int randowner = static_cast<int>(floor(numrequprocs* DRT::Problem::Instance()->Random()->Uni()));

     // myrank is allowed to set link
     if(myrankbond and randowner==numrequprocs-1)
     {
       // note: this means link is set between row cl and row ele on myrank,
       // all relevant information for myrank is stored in mybonds
       // loop over all requesters and store their veto
       std::vector<BindEventData>::iterator iter;
       for(iter=cliter->second.begin(); iter!=cliter->second.end(); ++iter)
       {
         iter->permission = false;
         decidedbonds[iter->requestproc].push_back(*iter);
       }
     }
     // certain requester is allowed to set the link
     else
     {
       // loop over all requesters and store veto for all requester except for one
       std::vector<BindEventData>::iterator iter;

       int counter = 0;
       for(iter=cliter->second.begin(); iter!=cliter->second.end(); ++iter)
       {
         if(randowner==counter)
         {
           // permission for this random proc
           decidedbonds[iter->requestproc].push_back(*iter);

           // erase old binding event
           if(myrankbond)
             mybonds.erase(cliter->first);

           // insert new binding event
           mybonds[cliter->first] = *iter;
         }
         else
         {
           iter->permission = false;
           decidedbonds[iter->requestproc].push_back(*iter);
         }
         counter++;
       }
     }
   }

  // that is it
  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::CommunicateBeamToBeamLinkageAfterRedistribution(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> > >& dbondcltosend)
{
  CheckInit();

  // build exporter
  DRT::Exporter exporter(bindis_->Comm());

  // -----------------------------------------------------------------------
  // send
  // -----------------------------------------------------------------------
  // number of messages
  const int length = dbondcltosend.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> > >::const_iterator p;
  for(p=dbondcltosend.begin(); p!=dbondcltosend.end(); ++p)
  {
    // ---- pack data for sending -----
    std::vector<char> sdata;
    std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> >::const_iterator iter;
    DRT::PackBuffer data;
    for(iter=p->second.begin(); iter!=p->second.end(); ++iter)
    {
     (*iter)->Pack(data);
    }
    data.StartPacking();
    for(iter=p->second.begin(); iter!=p->second.end(); ++iter)
    {
      (*iter)->Pack(data);
    }
    swap(sdata,data());
     // unblocking send
    exporter.ISend(myrank_, p->first, &(sdata[0]), static_cast<int>(sdata.size()), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  std::vector<int> summedtargets;
  PrepareReceivingProcs(dbondcltosend, summedtargets);

  // myrank receive all packs that are sent to him
  for(int rec=0; rec<summedtargets[myrank_]; ++rec)
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from,tag,rdata,length);
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", myrank_, from);

    // store received data
    std::vector<char>::size_type position = 0;
    while (position < rdata.size())
    {
      std::vector<char> data;
      DRT::ParObject::ExtractfromPack(position,rdata,data);
      // this Teuchos::rcp holds the memory of the node
      Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data),true);
      Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> beamtobeamlink =
          Teuchos::rcp_dynamic_cast<BEAMINTERACTION::BeamToBeamLinkage>(object);
      if (beamtobeamlink == Teuchos::null) dserror("Received object is not a beam to beam linkage");

      // insert new double bonds in my list
      doublebondcl_[beamtobeamlink->Id()] = beamtobeamlink;
    }

    if (position != rdata.size())
      dserror("Mismatch in size of data %d <-> %d",static_cast<int>(rdata.size()),position);
  }

  // wait for all communication to finish
  Wait(exporter,request,static_cast<int>(dbondcltosend.size()));

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::CommunicateForceStiff(
    std::map<int, std::vector<CommForceStiff> >& sendforcestiff,
    std::vector<CommForceStiff>&                 recvforcestiff) const
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::CommunicateForceStiff");

  ISendRecvAny(sendforcestiff,recvforcestiff);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::CommunicateCrosslinkerUnbinding(
    std::map<int, std::vector<UnBindEventData> >& sendunbindevent,
    std::vector<UnBindEventData>&                 myrankunbindevent) const
{
  CheckInit();

  ISendRecvAny(sendunbindevent,myrankunbindevent);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::GetNeighboringEles()
{
  CheckInit();

  // measure time for evaluating this function
  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::GetNeighboringEles");

#ifdef DEBUG
  if(static_cast<int>(exteletobinmap_.size()) != ia_discret_->ElementColMap()->NumMyElements())
    dserror("std::map does not equal elecolmap (check e.g. if extended ghosting contains "
            " standard ghosting). Therefore not every contact can be detected");
#endif

  // loop over all column elements
  std::map<int, std::set<int> >::const_iterator coleleiter;
  for(coleleiter=exteletobinmap_.begin(); coleleiter!=exteletobinmap_.end(); ++coleleiter)
  {
    const int elegid = coleleiter->first;
    DRT::Element* currele = ia_discret_->gElement(elegid);

    // (unique) set of neighboring bins for all col bins assigned to current element
    std::set<int> neighboring_binIds;

    // loop over all row bins
    std::set<int>::const_iterator colbiniter;
    for(colbiniter=exteletobinmap_[elegid].begin(); colbiniter!=exteletobinmap_[elegid].end(); ++colbiniter)
    {
      std::vector<int> loc_neighboring_binIds;
      loc_neighboring_binIds.reserve(27);

      // do not check on existence here -> shifted to GetBinContent
      binning_->GetNeighborAndOwnBinIds(*colbiniter,loc_neighboring_binIds);

      // build up comprehensive unique set of neighboring bins
      neighboring_binIds.insert(loc_neighboring_binIds.begin(), loc_neighboring_binIds.end());
    }
    // get unique vector of comprehensive neighboring bins
    std::vector<int> glob_neighboring_binIds(neighboring_binIds.begin(), neighboring_binIds.end());

    // set of elements that lie in neighboring bins
    std::set<DRT::Element*> neighboring_elements;
    binning_->GetBinContent(neighboring_elements,bin_beamcontent_,glob_neighboring_binIds);

    // sort out elements that should not be considered in contact evaluation
    VerifyNeighbors(currele, neighboring_elements);

  } // loop over all row elements

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::VerifyNeighbors(
    DRT::Element*            currele,
    std::set<DRT::Element*>& neighbors) const
{
  CheckInit();

  // sort out elements that should not be considered in contact evaluation
  std::set<DRT::Element*>::const_iterator eiter;
  for(eiter=neighbors.begin(); eiter!=neighbors.end(); ++eiter)
  {
    // 1) ensure each contact only evalutated once on myrank
    if(not(currele->Id() < (*eiter)->Id()))
    {
      neighbors.erase(*eiter);
      continue;
    }

    // 2) ensure that two elements sharing the same node do not get into contact
    for (int i=0; i<2; i++)
      for (int j=0; j<2; j++)
        if((*eiter)->NodeIds()[i]==currele->NodeIds()[j])
          neighbors.erase(*eiter);
  }

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateMaps()
{
  CheckInit();

  // todo: performance improvement by using the same exporter object every time
  // and not doing the safety checks in Linalg::Export. See in particle_timint
  // how this can be done.

  //todo: check if update is necessary (->SameAs())

  // beam displacement
  UpdateDofMapOfVector(ia_discret_, ia_disnp_);

  ia_stiff_crosslink_ = Teuchos::rcp(new
      LINALG::SparseMatrix(*ia_discret_->DofRowMap(), 81, true, true));
  // force
  ia_force_crosslink_ = Teuchos::rcp(new Epetra_Vector(*ia_discret_->DofRowMap(),true));

  // that is it
  return;

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UpdateDofMapOfVector(
    Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<Epetra_Vector>&      dofmapvec,
    Teuchos::RCP<Epetra_Vector>       old)
{
  CheckInit();

  // todo: performance improvement by using the same exporter object every time
  // and not doing the safety checks in Linalg::Export. See in particle_timint
  // how this can be done.

  if (dofmapvec != Teuchos::null)
  {
    if(old==Teuchos::null)
      old = dofmapvec;
    dofmapvec = LINALG::CreateVector(*discret->DofRowMap(),true);
    LINALG::Export(*old, *dofmapvec);
  }

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::AssembleEleForceStiffIntoSystemVectorMatrix(
  const DRT::Discretization&         discret,
  const int                          elegid1,
  const int                          elegid2,
  const Epetra_SerialDenseVector&    elevec1,
  const Epetra_SerialDenseVector&    elevec2,
  const Epetra_SerialDenseMatrix&    elemat11,
  const Epetra_SerialDenseMatrix&    elemat12,
  const Epetra_SerialDenseMatrix&    elemat21,
  const Epetra_SerialDenseMatrix&    elemat22,
  Teuchos::RCP<Epetra_Vector>        sysvec,
  Teuchos::RCP<LINALG::SparseMatrix> sysmat
  ) const
{
  // the entries of elevec1 belong to the Dofs of this element
  DRT::Element* ele1 = discret.gElement(elegid1);
  // the entries of elevec2 belong to the Dofs of this element
  DRT::Element* ele2 = discret.gElement(elegid2);

  // get element location vector and ownerships
  std::vector<int> lmrow1;
  std::vector<int> lmrow2;
  std::vector<int> lmrowowner1;
  std::vector<int> lmrowowner2;
  std::vector<int> lmstride;

  ele1->LocationVector(discret,lmrow1,lmrowowner1,lmstride);
  ele2->LocationVector(discret,lmrow2,lmrowowner2,lmstride);

  // assemble both element vectors into global system vector
  if(sysvec!=Teuchos::null)
  {
    LINALG::Assemble(*sysvec,elevec1,lmrow1,lmrowowner1);
    LINALG::Assemble(*sysvec,elevec2,lmrow2,lmrowowner2);
  }

  // and finally also assemble stiffness contributions
  if(sysmat!=Teuchos::null)
  {
    sysmat->Assemble(0,elemat11,lmrow1,lmrowowner1,lmrow1);
    sysmat->Assemble(0,elemat12,lmrow1,lmrowowner1,lmrow2);
    sysmat->Assemble(0,elemat21,lmrow2,lmrowowner2,lmrow1);
    sysmat->Assemble(0,elemat22,lmrow2,lmrowowner2,lmrow2);
  }

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::AssembleRecvEleForceStiffIntoSystemVectorMatrix(
  const std::vector<CommForceStiff>& recvforcestiff,
  Teuchos::RCP<Epetra_Vector>        sysvec,
  Teuchos::RCP<LINALG::SparseMatrix> sysmat
  ) const
{
  // assemble received stiffness and force contributions
  std::vector<CommForceStiff>::const_iterator fsiter;
  for(fsiter=recvforcestiff.begin(); fsiter!=recvforcestiff.end(); ++fsiter)
  {
    // safety check (todo: put this in debug mode as soon as extensively tested in release)
    if(ia_discret_->ElementColMap()->LID(fsiter->elegid1)<0)
      dserror("Element with gid %i for which proc %i received forcestiff not ghosted",fsiter->elegid1,myrank_);
    if(ia_discret_->ElementColMap()->LID(fsiter->elegid2)<0)
      dserror("Element with gid %i for which proc %i received forcestiff not ghosted",fsiter->elegid2,myrank_);

    AssembleEleForceStiffIntoSystemVectorMatrix(*ia_discret_,
                                                fsiter->elegid1,
                                                fsiter->elegid2,
                                                fsiter->ele1force,
                                                fsiter->ele2force,
                                                fsiter->ele11stiff,
                                                fsiter->ele12stiff,
                                                fsiter->ele21stiff,
                                                fsiter->ele22stiff,
                                                sysvec,
                                                sysmat);
  }

  // that is it
  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::TransformForceAndStiff(
    bool force,
    bool stiff)
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR("STR::MODELEVALUATOR::Crosslinking::TransformForceAndStiff");

  if(force)
  {
    // transform force vector to problem discret layout/distribution
    force_crosslink_ = coupsia_->MasterToSlave(ia_force_crosslink_);
  }
  // transform stiffness matrix to problem discret layout/distribution
  if(stiff)
  {
    stiff_crosslink_->UnComplete();
    // transform stiffness matrix to problem discret layout/distribution
    (*siatransform_)(*ia_stiff_crosslink_,
                     1.0,
                     ADAPTER::CouplingMasterConverter(*coupsia_),
                     *stiff_crosslink_,
                     false);
  }

  // that is it
  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template<typename T>
void STR::MODELEVALUATOR::Crosslinking::ISend(
    DRT::Exporter& exporter,
    std::vector<MPI_Request>& request,
    const std::map<int, std::vector<T> >& send) const
{
  CheckInit();

  // number of messages
  const int length = send.size();
  request.resize(length);
  int tag = 0;
  typename std::map<int, std::vector<T> >::const_iterator p;
  for(p=send.begin(); p!=send.end(); ++p)
  {
    // ---- pack data for sending -----
    std::vector<char> sdata;
    typename std::vector<T>::const_iterator iter;
    DRT::PackBuffer data;
    for(iter=p->second.begin(); iter!=p->second.end(); ++iter)
    {
     Pack(data,*iter);
    }
    data.StartPacking();
    for(iter=p->second.begin(); iter!=p->second.end(); ++iter)
    {
     Pack(data,*iter);
    }
    swap(sdata,data());

    // unblocking send
    exporter.ISend(myrank_, p->first, &(sdata[0]), static_cast<int>(sdata.size()), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // that is it
  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template<typename T>
void STR::MODELEVALUATOR::Crosslinking::PrepareReceivingProcs(
    const std::map<int, std::vector<T> >& datasenttorank,
    std::vector<int>& summedtargets) const
{
  CheckInit();

  // get number of procs from which myrank receives data
  std::vector<int> targetprocs(numproc_,0);
  typename std::map<int, std::vector<T> >::const_iterator prociter;
  for(prociter=datasenttorank.begin(); prociter!=datasenttorank.end(); ++prociter)
    targetprocs[prociter->first] = 1;
  // store number of messages myrank receives
  summedtargets.resize(numproc_,0);
  bindis_->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc_);

  // that is it
  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template<typename T>
void STR::MODELEVALUATOR::Crosslinking::RecvAny(
    DRT::Exporter&  exporter,
    const int& receivesize,
    std::vector<T>& recv) const
{
  CheckInit();

  // myrank receive all packs that are sent to him
  for(int rec=0; rec<receivesize; ++rec)
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from,tag,rdata,length);
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", myrank_, from);

    // store received data
    std::vector<char>::size_type position = 0;
    while (position < rdata.size())
    {
      // ---- extract received data -----
      T recdata;
      UnPack(position,rdata,recdata);

      // add received data to list of unbindevents on myrank
      recv.push_back(recdata);
    }

    if (position != rdata.size())
      dserror("Mismatch in size of data %d <-> %d",static_cast<int>(rdata.size()),position);
  }

  // that is it
  return;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template<typename T>
void STR::MODELEVALUATOR::Crosslinking::ISendRecvAny(
  const std::map<int, std::vector<T> >& send,
  std::vector<T>&                       recv) const
{
  CheckInit();

  // build exporter
  DRT::Exporter exporter(bindis_->Comm());

  // -----------------------------------------------------------------------
  // send
  // -----------------------------------------------------------------------
  // unblocking send
  std::vector<MPI_Request> request;
  ISend(exporter,request,send);

  // -----------------------------------------------------------------------
  // prepare receive
  // -----------------------------------------------------------------------
  std::vector<int> summedtargets;
  PrepareReceivingProcs(send, summedtargets);

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  int receivesize = summedtargets[myrank_];
  RecvAny(exporter,receivesize,recv);

  // wait for all communication to finish
  Wait(exporter,request,static_cast<int>(send.size()));

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Wait(
    DRT::Exporter& exporter,
    std::vector<MPI_Request>& request,
    const int& length) const
{
  CheckInit();

  // wait for all communication to finish
  for (int i=0; i<length; ++i)
    exporter.Wait(request[i]);

  // note: if we have done everything correct, this should be a no time operation
  bindis_->Comm().Barrier(); // I feel better this way ;-)

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Pack(
  DRT::PackBuffer&     data,
  const BindEventData& bindeventdata) const
{
  CheckInit();

  // pack data that is communicated
  DRT::ParObject::AddtoPack(data,bindeventdata.clgid);
  DRT::ParObject::AddtoPack(data,bindeventdata.elegid);
  DRT::ParObject::AddtoPack(data,bindeventdata.bspotlocn);
  DRT::ParObject::AddtoPack(data,bindeventdata.requestproc);
  DRT::ParObject::AddtoPack(data,bindeventdata.permission);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Pack(
  DRT::PackBuffer&       data,
  const UnBindEventData& unbindeventdata) const
{
  CheckInit();

  // pack data that is communicated
  DRT::ParObject::AddtoPack(data,unbindeventdata.eletoupdate);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Pack(
  DRT::PackBuffer&      data,
  const CommForceStiff& forcestiffdata) const
{
  CheckInit();

  // pack data that is communicated
  DRT::ParObject::AddtoPack(data,forcestiffdata.elegid1);
  DRT::ParObject::AddtoPack(data,forcestiffdata.elegid2);
  DRT::ParObject::AddtoPack(data,forcestiffdata.ele1force);
  DRT::ParObject::AddtoPack(data,forcestiffdata.ele2force);
  DRT::ParObject::AddtoPack(data,forcestiffdata.ele11stiff);
  DRT::ParObject::AddtoPack(data,forcestiffdata.ele12stiff);
  DRT::ParObject::AddtoPack(data,forcestiffdata.ele21stiff);
  DRT::ParObject::AddtoPack(data,forcestiffdata.ele22stiff);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UnPack(
  std::vector<char>::size_type& position,
  std::vector<char>             data,
  BindEventData&                bindeventdata) const
{
  CheckInit();

  // extract data
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.clgid);
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.elegid);
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.bspotlocn);
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.requestproc);
  DRT::ParObject::ExtractfromPack(position,data,bindeventdata.permission);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UnPack(
  std::vector<char>::size_type& position,
  std::vector<char>             data,
  UnBindEventData&              unbindeventdata) const
{
  CheckInit();

  // extract data
  DRT::ParObject::ExtractfromPack(position,data,unbindeventdata.eletoupdate);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::UnPack(
  std::vector<char>::size_type& position,
  std::vector<char>             data,
  CommForceStiff&               forcestiffdata) const
{
  CheckInit();

  // extract data
  DRT::ParObject::ExtractfromPack(position,data,forcestiffdata.elegid1);
  DRT::ParObject::ExtractfromPack(position,data,forcestiffdata.elegid2);
  DRT::ParObject::ExtractfromPack(position,data,forcestiffdata.ele1force);
  DRT::ParObject::ExtractfromPack(position,data,forcestiffdata.ele2force);
  DRT::ParObject::ExtractfromPack(position,data,forcestiffdata.ele11stiff);
  DRT::ParObject::ExtractfromPack(position,data,forcestiffdata.ele12stiff);
  DRT::ParObject::ExtractfromPack(position,data,forcestiffdata.ele21stiff);
  DRT::ParObject::ExtractfromPack(position,data,forcestiffdata.ele22stiff);

  // that is it
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Crosslinking::Logo() const
{
  CheckInit();

  if(myrank_==0)
  {
    IO::cout << "\n****************************************************************" << IO::endl;
    IO::cout << "*                                                              *" << IO::endl;
    IO::cout << "*          Welcome to a Biopolymer Network Simulation          *" << IO::endl;
    IO::cout << "*                                                              *" << IO::endl;
    IO::cout << "****************************************************************" << IO::endl;
    IO::cout << "                                                                  " << IO::endl;
    IO::cout << "                                                                  " << IO::endl;
    IO::cout << "                      0=========================0                 " << IO::endl;
    IO::cout << "                    //|   \\            /       /||                " << IO::endl;
    IO::cout << "                   // |    \\ |       |/       //||                " << IO::endl;
    IO::cout << "                  //  |  /  \\|       /       // ||                " << IO::endl;
    IO::cout << "                 //   |  \\   \\   /  /|\\     //  ||                " << IO::endl;
    IO::cout << "                //    |  /   |\\ /  / | \\   //   ||                " << IO::endl;
    IO::cout << "               //     |  \\   | \\     |  \\ //  / ||                " << IO::endl;
    IO::cout << "              //  \\  /|  /   |/      |   //  /  ||                " << IO::endl;
    IO::cout << "              0=========================0 \\ /   ||                " << IO::endl;
    IO::cout << "             ||    /\\ |____          |  || \\    ||                " << IO::endl;
    IO::cout << "             ||   /  \\|    \\   ------   ||/ \\   ||                " << IO::endl;
    IO::cout << "             ||  /    |                 ||      ||                " << IO::endl;
    IO::cout << "             || /     0----------/------||------0-                " << IO::endl;
    IO::cout << "             ||      /   /       \\      ||     //                 " << IO::endl;
    IO::cout << "             ||     /___/  \\     /    / ||    //                  " << IO::endl;
    IO::cout << "             ||    /        \\    \\   /  ||   //                   " << IO::endl;
    IO::cout << "             ||   /  \\/\\/\\/  \\   /  /   ||  //                    " << IO::endl;
    IO::cout << "             ||  /      /     \\  \\ /    || //                     " << IO::endl;
    IO::cout << "             || /      /         /      ||//                      " << IO::endl;
    IO::cout << "             ||/                       /||/                       " << IO::endl;
    IO::cout << "              0=========================0                         " << IO::endl;
    IO::cout << "                                                                     " << IO::endl;
    IO::cout << "                                                                     " << IO::endl;
  }

  // that is it
  return;
}

//-----------------------------------------------------------------------------
// explicit template instantiation (to please every compiler)
//-----------------------------------------------------------------------------
template void STR::MODELEVALUATOR::Crosslinking::ISend(
    DRT::Exporter&,std::vector<MPI_Request>&,const std::map<int, std::vector<BindEventData> >&) const;
template void STR::MODELEVALUATOR::Crosslinking::ISend(
    DRT::Exporter&,std::vector<MPI_Request>&,const std::map<int, std::vector<UnBindEventData> >&) const;
template void STR::MODELEVALUATOR::Crosslinking::ISend(
    DRT::Exporter&,std::vector<MPI_Request>&,const std::map<int, std::vector<CommForceStiff> >&) const;

template void STR::MODELEVALUATOR::Crosslinking::PrepareReceivingProcs(
    const std::map<int, std::vector<BindEventData> >&,std::vector<int>&) const;
template void STR::MODELEVALUATOR::Crosslinking::PrepareReceivingProcs(
    const std::map<int, std::vector<UnBindEventData> >&,std::vector<int>&) const;
template void STR::MODELEVALUATOR::Crosslinking::PrepareReceivingProcs(
    const std::map<int, std::vector<CommForceStiff> >&,std::vector<int>&) const;

template void STR::MODELEVALUATOR::Crosslinking::RecvAny(
    DRT::Exporter&,const int&,std::vector<BindEventData>&) const;
template void STR::MODELEVALUATOR::Crosslinking::RecvAny(
    DRT::Exporter&,const int&,std::vector<UnBindEventData>&) const;
template void STR::MODELEVALUATOR::Crosslinking::RecvAny(
    DRT::Exporter&,const int&,std::vector<CommForceStiff>&) const;

template void STR::MODELEVALUATOR::Crosslinking::ISendRecvAny(
    const std::map<int, std::vector<BindEventData> >&,std::vector<BindEventData>&) const;
template void STR::MODELEVALUATOR::Crosslinking::ISendRecvAny(
    const std::map<int, std::vector<UnBindEventData> >&,std::vector<UnBindEventData>&) const;
template void STR::MODELEVALUATOR::Crosslinking::ISendRecvAny(
    const std::map<int, std::vector<CommForceStiff> >&,std::vector<CommForceStiff>&) const;


///*-----------------------------------------------------------------------------*
// | nodes are checked and transferred if necessary              eichinger 09/16 |
// *-----------------------------------------------------------------------------*/
//void PARTICLE::Algorithm::TransferNodes(
//    Teuchos::RCP<DRT::Discretization> discret,
//    Teuchos::RCP<Epetra_Vector>       disnp)
//{
//  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::TransferNodes");
//
//  // list of homeless nodes
//  std::list<Teuchos::RCP<DRT::Node> > homelessnodes;
//
//  // check in each bin whether nodes have moved out
//  // first run over nodes and then process whole bin in which node is located
//  // until all particles have been checked
//  std::map<int, std::set<int> > rownodetobin;
//  std::map<int, std::set<int> >::const_iterator rownodeiter;
//
//  // loop over all row nodes
//  for(rownodeiter=rownodetobin.begin(); rownodeiter!=rownodetobin.end(); ++rownodeiter)
//  {
//    // get current node
//    DRT::Node *currnode = discret->gNode(rownodeiter->first);
//    // as checked above, there is only one element in currele array
//    const int binId = rownodeiter->second;
//
//    // get current node position
//    double currpos[3] = {0.0,0.0,0.0};
//    GetCurrentNodePos(discret,currnode,disnp,currpos);
//
//    // get current bin of node
//    const int gidofbin = ConvertPosToGid(currpos);
//    if(gidofbin != binId) // node has left current bin
//    {
//      // gather all node Ids that will be removed and remove them afterwards
//      // (looping over nodes and deleting at the same time is detrimental)
//      tobemoved.push_back(currnode->Id());
//      // find new bin for particle
//      /*bool placed = */PlaceNodeCorrectly(Teuchos::rcp(currnode,false), currpos, homelessparticles);
//    }
//
//
//
//    // finally remove nodes from their old bin
//    for(size_t iter=0; iter<tobemoved.size(); iter++)
//    {
//      dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currbin)->DeleteNode(tobemoved[iter]);
//    }
//
//  } // end for ibin
//
//#ifdef DEBUG
//  if(homelessparticles.size())
//    std::cout << "There are " << homelessparticles.size() << " homeless particles on proc" << myrank_ << std::endl;
//#endif
//
//  // homeless particles are sent to their new processors where they are inserted into their correct bin
//  FillParticlesIntoBinsRemoteIdList(homelessparticles);
//
//  // check whether all procs have a filled bindis_,
//  // oldmap in ExportColumnElements must be Reset() on every proc or nowhere
//  bindis_->CheckFilledGlobally();
//
//  // new ghosting if necessary
//  if (ghosting)
//    bindis_->ExtendedGhosting(*bincolmap_,true,false,true,false);
//  else
//    bindis_->FillComplete(true, false, true);
//
//  // reconstruct element -> bin pointers for fixed particle wall elements and fluid elements
//  bool rebuildwallpointer = true;
//  if(moving_walls_)
//    rebuildwallpointer = false;
//  BuildElementToBinPointers(rebuildwallpointer);
//
//  // update state vectors in time integrator to the new layout
//  if(updatestates)
//  {
//    TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::TransferParticles::UpdateStates");
//    particles_->UpdateStatesAfterParticleTransfer();
//    UpdateStates();
//  }
//
//  return;
//}
