/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_multiphase.cpp

 \brief material for multiphase porous fluid

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/




#include <vector>
#include "fluidporo_multiphase.H"
#include "fluidporo_singlephase.H"

#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_calc_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

#include <Epetra_SerialDenseSolver.h>

/*----------------------------------------------------------------------*
 | constructor of paramter class                            vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroMultiPhase::FluidPoroMultiPhase(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: MatList(matdata),
  permeability_(matdata->GetDouble("PERMEABILITY")),
  numfluidphases_(matdata->GetInt("NUMFLUIDPHASES")),
  dof2pres_(Teuchos::null),
  constraintphaseID_(-1),
  isinit_(false)
{
}

/*----------------------------------------------------------------------*
 | create a poro multiphase material                        vuong 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::FluidPoroMultiPhase::CreateMaterial()
{
  return Teuchos::rcp(new MAT::FluidPoroMultiPhase(this));
}

/*----------------------------------------------------------------------*
 | initialize                                               vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::PAR::FluidPoroMultiPhase::Initialize()
{
  //  matrix holding the conversion from pressures and dofs
  dof2pres_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(numfluidphases_,numfluidphases_));

  //  matrix holding the conversion from pressures and dofs
  // reset
  dof2pres_->Scale(0.0);

  for(int iphase=0;iphase<(int)matids_->size();iphase++)
  {
    // get the single phase material by its ID
    const int matid = (*matids_)[iphase];
    Teuchos::RCP< MAT::Material> singlemat = MaterialById(matid);

    if(iphase < numfluidphases_)
    {
      // safety check and cast
      if(singlemat->MaterialType() != INPAR::MAT::m_fluidporo_singlephase)
        dserror("You have chosen %d fluidphases, however your material number %d is no poro singlephase material", numfluidphases_, iphase+1);
      const MAT::FluidPoroSinglePhase& singlephase = static_cast<const MAT::FluidPoroSinglePhase&>(*singlemat);

      if(singlephase.PoroPhaseLawType() == INPAR::MAT::m_fluidporo_phaselaw_constraint)
      {
        if(constraintphaseID_!=-1)
          dserror("More than one constraint phase law defined. Are you sure this makes sense?");
        constraintphaseID_ = iphase;
      }

      // fill the coefficients into matrix
      singlephase.FillDoFMatrix(*dof2pres_,iphase);
    }
    else
    {
      // safety check and cast
      if(singlemat->MaterialType() != INPAR::MAT::m_fluidporo_singlevolfrac)
        dserror("You have chosen %d volume fractions, however your material number %d is no poro volume fraction material");
    }
  }

  // check
  if(constraintphaseID_==-1)
    dserror("No constraint phase law defined. Are you sure this makes sense?");

  // invert dof2pres_ to get conversion from dofs to prbessures
  {
    Epetra_SerialDenseSolver inverse;
    inverse.SetMatrix(*dof2pres_);
    int err = inverse.Invert();
    if (err != 0)
      dserror("Inversion of matrix for DOF transform failed with errorcode %d. Is your system of DOFs linear independent?",err);
  }

  isinit_=true;
  return;
}

/*----------------------------------------------------------------------*
 | global instance of parameter class                       vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroMultiPhaseType MAT::FluidPoroMultiPhaseType::instance_;

/*----------------------------------------------------------------------*
 | create material from data                                vuong 08/16 |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::FluidPoroMultiPhaseType::Create( const std::vector<char> & data )
{
  MAT::FluidPoroMultiPhase* FluidPoroMultiPhase = new MAT::FluidPoroMultiPhase();
  FluidPoroMultiPhase->Unpack(data);
  return FluidPoroMultiPhase;
}


/*----------------------------------------------------------------------*
 | construct empty material object                          vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroMultiPhase::FluidPoroMultiPhase()
  : MatList(),
    paramsporo_(NULL)
{
}

/*----------------------------------------------------------------------*
 | construct the material object given material parameter   vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroMultiPhase::FluidPoroMultiPhase(MAT::PAR::FluidPoroMultiPhase* params)
  : MatList(params),
    paramsporo_(params)
{
}

/*----------------------------------------------------------------------*
 | reset everything                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::Clear()
{
  paramsporo_ = NULL;
  return;
}

/*----------------------------------------------------------------------*
 | initialize                                               vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::Initialize()
{
  std::map<int,Teuchos::RCP<MAT::Material> >* materials;

  if (Parameter() != NULL) // params is null pointer in post-process mode
  {
    if(Parameter()->local_)
      materials = MaterialMapWrite();
    else
      materials = Parameter()->MaterialMapWrite();

    std::map<int,Teuchos::RCP<MAT::Material> >::iterator it;
    for(it=materials->begin();it!=materials->end();it++)
    {
      Teuchos::RCP<MAT::FluidPoroSinglePhaseBase> actphase =
          Teuchos::rcp_dynamic_cast<FluidPoroSinglePhaseBase>(it->second,true);
      actphase->Initialize();
    }

    if(not paramsporo_->isinit_)
      paramsporo_->Initialize();
  }
  return;
}

/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (paramsporo_ != NULL) matid = paramsporo_->Id();  // in case we are in post-process mode

  AddtoPack(data,matid);

  // Pack base class material
  MAT::MatList::Pack(data);
}

/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  Clear();

  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover paramsporo_
  int matid(-1);
  ExtractfromPack(position,data,matid);
  paramsporo_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        //Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are in a diamond inheritance structure
        paramsporo_ = dynamic_cast<MAT::PAR::FluidPoroMultiPhase*>(mat);
      }
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  MAT::MatList::ExtractfromPack(position,data,basedata);
  MAT::MatList::Unpack(basedata);

  // in the postprocessing mode, we do not unpack everything we have packed
  // -> position check cannot be done in this case
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

/*----------------------------------------------------------------------*
 *  Evaluate generalized pressure of all phases            vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::EvaluateGenPressure(
    std::vector<double>& genpressure,
    const std::vector<double>& phinp) const
{
  // evaluate the pressures
  for(int iphase=0; iphase<NumFluidPhases(); iphase++)
  {
    //get the single phase material
    const MAT::FluidPoroSinglePhase& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMultiMaterial(*this,iphase);

    // evaluate generalized pressure (i.e. some kind of linear combination of the true pressures)
    genpressure[iphase] = singlephasemat.EvaluateGenPressure(iphase,phinp);
  }
  return;
}

/*----------------------------------------------------------------------*
 *  Evaluate generalized pressure of all phases            vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::GetVolFracPressure(
    std::vector<double>& volfracpress) const
{
  // evaluate the pressures
  for(int ivolfrac=0; ivolfrac<NumMat()-NumFluidPhases(); ivolfrac++)
  {
    //get the single phase material
    const MAT::FluidPoroSingleVolFrac& singlevolfracmat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSingleVolFracMatFromMaterial(*this,ivolfrac+NumFluidPhases());

    // evaluate volume fraction pressure
    volfracpress[ivolfrac] = singlevolfracmat.Pressure();
  }
  return;
}

/*----------------------------------------------------------------------*
 *   Evaluate saturation of the phase                       vuong 08/16 |
*----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::EvaluateSaturation(
    std::vector<double>& saturation,
    const std::vector<double>& phinp,
    const std::vector<double>& pressure) const
{
  // get the number of the phase, which saturation is calculated by the saturation constraint
  const int constraintsaturationphase = paramsporo_->constraintphaseID_;

  // the constraint saturation is calculated from 1- sum(all other saturations)
  saturation[constraintsaturationphase] = 1.0;
  for(int iphase=0; iphase<NumFluidPhases(); iphase++)
  {
    if(iphase != constraintsaturationphase)
    {
      // get the single phase material
      const MAT::FluidPoroSinglePhase& singlephasemat =
          POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMultiMaterial(*this,iphase);

      saturation[iphase] = singlephasemat.EvaluateSaturation(iphase,phinp,pressure);
      // the saturation of the last phase is 1.0- (sum of all saturations)
      saturation[constraintsaturationphase] -= saturation[iphase];
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | transform generalized pressures to true pressures        vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::TransformGenPresToTruePres(
    const std::vector<double>& phinp,
    std::vector<double>& phi_transformed) const
{
  // get trafo matrix
  const Epetra_SerialDenseMatrix& dof2pres = *paramsporo_->dof2pres_;
  //simple matrix vector product
  phi_transformed.resize(phinp.size());
  for(int i=0;i<NumFluidPhases();i++)
    for(int j=0;j<NumFluidPhases();j++)
      phi_transformed[i] += dof2pres(i,j)*phinp[j];
  return;
}

/*----------------------------------------------------------------------------------------*
 * Evaluate derivative of degree of freedom with respect to pressure          vuong 08/16 |
*----------------------------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::EvaluateDerivOfDofWrtPressure(
    Epetra_SerialDenseMatrix& derivs,
    const std::vector<double>& state) const
{
  for(int iphase=0; iphase<NumFluidPhases(); iphase++)
  {
    // get the single phase material by its ID
    const int matid = MatID(iphase);
    Teuchos::RCP< MAT::Material> singlemat = MaterialById(matid);
    const MAT::FluidPoroSinglePhase& singlephase = static_cast<const MAT::FluidPoroSinglePhase&>(*singlemat);

    for(int jphase=0; jphase<NumFluidPhases(); jphase++)
    {
      derivs(iphase,jphase)   = singlephase.EvaluateDerivOfDofWrtPressure(iphase,jphase,state);
    }
  }
  return;
}

/*--------------------------------------------------------------------------*
 *  Evaluate derivative of saturation w.r.t. pressure           vuong 08/16 |
*---------------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::EvaluateDerivOfSaturationWrtPressure(
    Epetra_SerialDenseMatrix& derivs,
    const std::vector<double>& pressure) const
{
  // get the number of the phase, which saturation is calculated by the saturation constraint
  const int constraintsaturationphase = paramsporo_->constraintphaseID_;

  for(int iphase=0; iphase<NumFluidPhases(); iphase++)
  {
    // skip constraint saturation phase
    if(iphase == constraintsaturationphase)
      continue;

    // get the single phase material by its ID
    const int matid = MatID(iphase);
    Teuchos::RCP< MAT::Material> singlemat = MaterialById(matid);
    const MAT::FluidPoroSinglePhase& singlephase = static_cast<const MAT::FluidPoroSinglePhase&>(*singlemat);

    for(int jphase=0; jphase<NumFluidPhases(); jphase++)
    {
      const double saturationderiv = singlephase.EvaluateDerivOfSaturationWrtPressure(iphase,jphase,pressure);
      derivs(iphase,jphase) = saturationderiv;
      // the saturation of the last phase is 1.0- (sum of all saturations)
      // -> the derivative of this saturation = -1.0 (sum of all saturation derivatives)
      derivs(constraintsaturationphase,jphase) += -1.0*saturationderiv;
    }
  }
  return;
}

/*--------------------------------------------------------------------------*
 *  Evaluate 2nd derivative of saturation w.r.t. pressure  kremheller 05/17 |
*---------------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::EvaluateSecondDerivOfSaturationWrtPressure(
    std::vector<Epetra_SerialDenseMatrix>& derivs,
    const std::vector<double>& pressure) const
{
  // get the number of the phase, which saturation is calculated by the saturation constraint
  const int constraintsaturationphase = paramsporo_->constraintphaseID_;

  for(int iphase=0; iphase<NumFluidPhases(); iphase++)
  {
    // skip constraint saturation phase
    if(iphase == constraintsaturationphase)
      continue;

    // get the single phase material by its ID
    const int matid = MatID(iphase);
    Teuchos::RCP< MAT::Material> singlemat = MaterialById(matid);
    const MAT::FluidPoroSinglePhase& singlephase = static_cast<const MAT::FluidPoroSinglePhase&>(*singlemat);

    for(int jphase=0; jphase<NumFluidPhases(); jphase++)
    {
      for(int kphase=0; kphase<NumFluidPhases(); kphase++)
      {
        const double saturationderivderiv
                    = singlephase.EvaluateSecondDerivOfSaturationWrtPressure(iphase,jphase,kphase,pressure);
        derivs[iphase](jphase,kphase) = saturationderivderiv;
        // the saturation of the last phase is 1.0- (sum of all saturations)
        // -> the derivative of this saturation = -1.0 (sum of all saturation derivatives)
        derivs[constraintsaturationphase](jphase,kphase) += -1.0*saturationderivderiv;
      }
    }
  }
  return;
}

