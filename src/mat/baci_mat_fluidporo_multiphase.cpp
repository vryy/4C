/*----------------------------------------------------------------------*/
/*! \file
 \brief material for multiphase porous fluid

   \level 3

 *----------------------------------------------------------------------*/



#include "baci_mat_fluidporo_multiphase.H"

#include "baci_lib_globalproblem.H"
#include "baci_mat_fluidporo_singlephase.H"
#include "baci_mat_par_bundle.H"
#include "baci_porofluidmultiphase_ele_calc_utils.H"

#include <Teuchos_SerialDenseSolver.hpp>

#include <vector>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor of paramter class                            vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroMultiPhase::FluidPoroMultiPhase(Teuchos::RCP<MAT::PAR::Material> matdata)
    : MatList(matdata),
      permeability_(matdata->GetDouble("PERMEABILITY")),
      numfluidphases_(matdata->GetInt("NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE")),
      numvolfrac_(-1),
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
  dof2pres_ = Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix(numfluidphases_, numfluidphases_));

  //  matrix holding the conversion from pressures and dofs
  // reset
  dof2pres_->putScalar(0.0);

  // get number of volume fractions
  numvolfrac_ = (int)(((int)matids_->size() - numfluidphases_) / 2);

  // safety check
  if ((int)matids_->size() != (int)(numvolfrac_ * 2 + numfluidphases_))
    dserror(
        "You have chosen %i materials, %i fluidphases and %f volume fractions, check your input "
        "definition\n"
        "Your Input should always look like (for example: 4 fluid phases, 2 volume fractions):\n"
        "MAT 0 MAT_FluidPoroMultiPhase LOCAL No PERMEABILITY 1.0 NUMMAT 8 MATIDS    1 2 3 4 5 6 7 "
        "8 NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE 4 END\n"
        "with: 4 fluid phases in multiphase pore space: materials have to be "
        "MAT_FluidPoroSinglePhase \n"
        "      2 volume fractions: materials have to be MAT_FluidPoroSingleVolFrac \n"
        "      2 volume fraction pressures: materials have to be MAT_FluidPoroVolFracPressure ",
        (int)matids_->size(), numfluidphases_,
        (double)(((double)matids_->size() - (double)numfluidphases_) / 2.0));

  for (int iphase = 0; iphase < (int)matids_->size(); iphase++)
  {
    // get the single phase material by its ID
    const int matid = (*matids_)[iphase];
    Teuchos::RCP<MAT::Material> singlemat = MaterialById(matid);

    // fluidphases at [0...numfluidphases-1]
    if (iphase < numfluidphases_)
    {
      // safety check and cast
      if (singlemat->MaterialType() != INPAR::MAT::m_fluidporo_singlephase)
        dserror(
            "You have chosen %i fluidphases, however your material number %i is no poro "
            "singlephase material\n"
            "Your Input should always look like (for example: 4 fluid phases, 2 volume "
            "fractions):\n"
            "MAT 0 MAT_FluidPoroMultiPhase LOCAL No PERMEABILITY 1.0 NUMMAT 8 MATIDS    1 2 3 4 5 "
            "6 7 8 NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE 4 END\n"
            "with: 4 fluid phases in multiphase pore space: materials have to be "
            "MAT_FluidPoroSinglePhase \n"
            "      2 volume fractions: materials have to be MAT_FluidPoroSingleVolFrac \n"
            "      2 volume fraction pressures: materials have to be MAT_FluidPoroVolFracPressure ",
            numfluidphases_, iphase + 1);

      const MAT::FluidPoroSinglePhase& singlephase =
          static_cast<const MAT::FluidPoroSinglePhase&>(*singlemat);

      if (singlephase.PoroPhaseLawType() == INPAR::MAT::m_fluidporo_phaselaw_constraint)
      {
        if (constraintphaseID_ != -1)
          dserror("More than one constraint phase law defined. Are you sure this makes sense?");
        constraintphaseID_ = iphase;
      }

      // fill the coefficients into matrix
      singlephase.FillDoFMatrix(*dof2pres_, iphase);
    }
    // volume fractions at [numfluidphases-1...numfluidphases-1+numvolfrac]
    else if (iphase < numfluidphases_ + numvolfrac_)
    {
      // safety check
      if (singlemat->MaterialType() != INPAR::MAT::m_fluidporo_singlevolfrac)
        dserror(
            "You have chosen %i fluid phases and %i volume fractions, however your material number "
            "%i is no poro volume fraction material\n"
            "Your Input should always look like (for example: 4 fluid phases, 2 volume "
            "fractions):\n"
            "MAT 0 MAT_FluidPoroMultiPhase LOCAL No PERMEABILITY 1.0 NUMMAT 8 MATIDS    1 2 3 4 5 "
            "6 7 8 NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE 4 END\n"
            "with: 4 fluid phases in multiphase pore space: materials have to be "
            "MAT_FluidPoroSinglePhase \n"
            "      2 volume fractions: materials have to be MAT_FluidPoroSingleVolFrac \n"
            "      2 volume fraction pressures: materials have to be MAT_FluidPoroVolFracPressure ",
            numfluidphases_, (int)matids_->size() - numfluidphases_, iphase + 1);
    }
    // volume fraction pressures at [numfluidphases-1+numvolfrac...numfluidphases-1+2*numvolfrac]
    else if (iphase < numfluidphases_ + 2 * numvolfrac_)
    {
      // safety check
      if (singlemat->MaterialType() != INPAR::MAT::m_fluidporo_volfracpressure)
        dserror(
            "You have chosen %i fluid phases and %i volume fractions, however your material number "
            "%i is no poro volume fraction pressure material\n"
            "Your Input should always look like (for example: 4 fluid phases, 2 volume "
            "fractions):\n"
            "MAT 0 MAT_FluidPoroMultiPhase LOCAL No PERMEABILITY 1.0 NUMMAT 8 MATIDS    1 2 3 4 5 "
            "6 7 8 NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE 4 END\n"
            "with: 4 fluid phases in multiphase pore space: materials have to be "
            "MAT_FluidPoroSinglePhase \n"
            "      2 volume fractions: materials have to be MAT_FluidPoroSingleVolFrac \n"
            "      2 volume fraction pressures: materials have to be MAT_FluidPoroVolFracPressure ",
            numfluidphases_, (int)matids_->size() - numfluidphases_, iphase + 1);
    }
    else
      dserror("something went wrong here, why is iphase = %i", iphase);
  }

  // check
  if (constraintphaseID_ == -1 && numfluidphases_ > 0)
    dserror(
        "No constraint phase law defined but NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE > 0. Are you "
        "sure this makes sense?");

  // invert dof2pres_ to get conversion from dofs to pressures for the fluid phases
  if (numfluidphases_ > 0)
  {
    using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
    using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
    Teuchos::SerialDenseSolver<ordinalType, scalarType> inverse;
    inverse.setMatrix(dof2pres_);
    int err = inverse.invert();
    if (err != 0)
      dserror(
          "Inversion of matrix for DOF transform failed with errorcode %d. Is your system of DOFs "
          "linear independent?",
          err);
  }

  isinit_ = true;
  return;
}

/*----------------------------------------------------------------------*
 | global instance of parameter class                       vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroMultiPhaseType MAT::FluidPoroMultiPhaseType::instance_;

/*----------------------------------------------------------------------*
 | create material from data                                vuong 08/16 |
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* MAT::FluidPoroMultiPhaseType::Create(const std::vector<char>& data)
{
  MAT::FluidPoroMultiPhase* FluidPoroMultiPhase = new MAT::FluidPoroMultiPhase();
  FluidPoroMultiPhase->Unpack(data);
  return FluidPoroMultiPhase;
}


/*----------------------------------------------------------------------*
 | construct empty material object                          vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroMultiPhase::FluidPoroMultiPhase() : MatList(), paramsporo_(nullptr) {}

/*----------------------------------------------------------------------*
 | construct the material object given material parameter   vuong 08/16 |
 *----------------------------------------------------------------------*/
MAT::FluidPoroMultiPhase::FluidPoroMultiPhase(MAT::PAR::FluidPoroMultiPhase* params)
    : MatList(params), paramsporo_(params)
{
}

/*----------------------------------------------------------------------*
 | reset everything                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::Clear()
{
  paramsporo_ = nullptr;
  return;
}

/*----------------------------------------------------------------------*
 | initialize                                               vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::Initialize()
{
  std::map<int, Teuchos::RCP<MAT::Material>>* materials;

  if (Parameter() != nullptr)  // params is null pointer in post-process mode
  {
    if (Parameter()->local_)
      materials = MaterialMapWrite();
    else
      materials = Parameter()->MaterialMapWrite();

    std::map<int, Teuchos::RCP<MAT::Material>>::iterator it;
    for (it = materials->begin(); it != materials->end(); it++)
    {
      Teuchos::RCP<MAT::FluidPoroSinglePhaseBase> actphase =
          Teuchos::rcp_dynamic_cast<FluidPoroSinglePhaseBase>(it->second, true);
      actphase->Initialize();
    }

    if (not paramsporo_->isinit_) paramsporo_->Initialize();
  }
  return;
}

/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (paramsporo_ != nullptr) matid = paramsporo_->Id();  // in case we are in post-process mode

  AddtoPack(data, matid);

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

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover paramsporo_
  int matid(-1);
  ExtractfromPack(position, data, matid);
  paramsporo_ = nullptr;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction
        // are in a diamond inheritance structure
        paramsporo_ = dynamic_cast<MAT::PAR::FluidPoroMultiPhase*>(mat);
      }
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  MAT::MatList::ExtractfromPack(position, data, basedata);
  MAT::MatList::Unpack(basedata);

  // in the postprocessing mode, we do not unpack everything we have packed
  // -> position check cannot be done in this case
  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
 *  Evaluate generalized pressure of all phases            vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::EvaluateGenPressure(
    std::vector<double>& genpressure, const std::vector<double>& phinp) const
{
  // evaluate the pressures
  for (int iphase = 0; iphase < NumFluidPhases(); iphase++)
  {
    // get the single phase material
    const MAT::FluidPoroSinglePhase& singlephasemat =
        POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMultiMaterial(*this, iphase);

    // evaluate generalized pressure (i.e. some kind of linear combination of the true pressures)
    genpressure[iphase] = singlephasemat.EvaluateGenPressure(iphase, phinp);
  }
  return;
}

/*----------------------------------------------------------------------*
 *   Evaluate saturation of the phase                       vuong 08/16 |
 *----------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::EvaluateSaturation(std::vector<double>& saturation,
    const std::vector<double>& phinp, const std::vector<double>& pressure) const
{
  // get the number of the phase, which saturation is calculated by the saturation constraint
  const int constraintsaturationphase = paramsporo_->constraintphaseID_;

  // the constraint saturation is calculated from 1- sum(all other saturations)
  saturation[constraintsaturationphase] = 1.0;
  for (int iphase = 0; iphase < NumFluidPhases(); iphase++)
  {
    if (iphase != constraintsaturationphase)
    {
      // get the single phase material
      const MAT::FluidPoroSinglePhase& singlephasemat =
          POROFLUIDMULTIPHASE::ELEUTILS::GetSinglePhaseMatFromMultiMaterial(*this, iphase);

      saturation[iphase] = singlephasemat.EvaluateSaturation(iphase, phinp, pressure);
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
    const std::vector<double>& phinp, std::vector<double>& phi_transformed) const
{
  // get trafo matrix
  const CORE::LINALG::SerialDenseMatrix& dof2pres = *paramsporo_->dof2pres_;
  // simple matrix vector product
  phi_transformed.resize(phinp.size());
  for (int i = 0; i < NumFluidPhases(); i++)
    for (int j = 0; j < NumFluidPhases(); j++) phi_transformed[i] += dof2pres(i, j) * phinp[j];
  return;
}

/*----------------------------------------------------------------------------------------*
 * Evaluate derivative of degree of freedom with respect to pressure          vuong 08/16 |
 *----------------------------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::EvaluateDerivOfDofWrtPressure(
    CORE::LINALG::SerialDenseMatrix& derivs, const std::vector<double>& state) const
{
  for (int iphase = 0; iphase < NumFluidPhases(); iphase++)
  {
    // get the single phase material by its ID
    const int matid = MatID(iphase);
    Teuchos::RCP<MAT::Material> singlemat = MaterialById(matid);
    const MAT::FluidPoroSinglePhase& singlephase =
        static_cast<const MAT::FluidPoroSinglePhase&>(*singlemat);

    for (int jphase = 0; jphase < NumFluidPhases(); jphase++)
    {
      derivs(iphase, jphase) = singlephase.EvaluateDerivOfDofWrtPressure(iphase, jphase, state);
    }
  }
  return;
}

/*--------------------------------------------------------------------------*
 *  Evaluate derivative of saturation w.r.t. pressure           vuong 08/16 |
 *---------------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::EvaluateDerivOfSaturationWrtPressure(
    CORE::LINALG::SerialDenseMatrix& derivs, const std::vector<double>& pressure) const
{
  // get the number of the phase, which saturation is calculated by the saturation constraint
  const int constraintsaturationphase = paramsporo_->constraintphaseID_;

  for (int iphase = 0; iphase < NumFluidPhases(); iphase++)
  {
    // skip constraint saturation phase
    if (iphase == constraintsaturationphase) continue;

    // get the single phase material by its ID
    const int matid = MatID(iphase);
    Teuchos::RCP<MAT::Material> singlemat = MaterialById(matid);
    const MAT::FluidPoroSinglePhase& singlephase =
        static_cast<const MAT::FluidPoroSinglePhase&>(*singlemat);

    for (int jphase = 0; jphase < NumFluidPhases(); jphase++)
    {
      const double saturationderiv =
          singlephase.EvaluateDerivOfSaturationWrtPressure(iphase, jphase, pressure);
      derivs(iphase, jphase) = saturationderiv;
      // the saturation of the last phase is 1.0- (sum of all saturations)
      // -> the derivative of this saturation = -1.0 (sum of all saturation derivatives)
      derivs(constraintsaturationphase, jphase) += -1.0 * saturationderiv;
    }
  }
  return;
}

/*--------------------------------------------------------------------------*
 *  Evaluate 2nd derivative of saturation w.r.t. pressure  kremheller 05/17 |
 *---------------------------------------------------------------------------*/
void MAT::FluidPoroMultiPhase::EvaluateSecondDerivOfSaturationWrtPressure(
    std::vector<CORE::LINALG::SerialDenseMatrix>& derivs, const std::vector<double>& pressure) const
{
  // get the number of the phase, which saturation is calculated by the saturation constraint
  const int constraintsaturationphase = paramsporo_->constraintphaseID_;

  for (int iphase = 0; iphase < NumFluidPhases(); iphase++)
  {
    // skip constraint saturation phase
    if (iphase == constraintsaturationphase) continue;

    // get the single phase material by its ID
    const int matid = MatID(iphase);
    Teuchos::RCP<MAT::Material> singlemat = MaterialById(matid);
    const MAT::FluidPoroSinglePhase& singlephase =
        static_cast<const MAT::FluidPoroSinglePhase&>(*singlemat);

    for (int jphase = 0; jphase < NumFluidPhases(); jphase++)
    {
      for (int kphase = 0; kphase < NumFluidPhases(); kphase++)
      {
        const double saturationderivderiv = singlephase.EvaluateSecondDerivOfSaturationWrtPressure(
            iphase, jphase, kphase, pressure);
        derivs[iphase](jphase, kphase) = saturationderivderiv;
        // the saturation of the last phase is 1.0- (sum of all saturations)
        // -> the derivative of this saturation = -1.0 (sum of all saturation derivatives)
        derivs[constraintsaturationphase](jphase, kphase) += -1.0 * saturationderivderiv;
      }
    }
  }
  return;
}

BACI_NAMESPACE_CLOSE
