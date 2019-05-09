/*-----------------------------------------------------------------------------------------------*/
/*!
\brief constitutive relations for beam cross-section resultants (hyperelastic stored energy
function)

\maintainer Maximilian Grill

\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_elasthyper.H"
#include "beam_elasthyper_parameter.H"

#include "../drt_mat/matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Sacado.hpp>


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::BeamElastHyperMaterialType MAT::BeamElastHyperMaterialType::instance_;



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
DRT::ParObject* MAT::BeamElastHyperMaterialType::Create(const std::vector<char>& data)
{
  MAT::BeamElastHyperMaterial* matobject = new MAT::BeamElastHyperMaterial();
  matobject->Unpack(data);
  return matobject;
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::BeamElastHyperMaterial::BeamElastHyperMaterial() : params_(NULL)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::BeamElastHyperMaterial::BeamElastHyperMaterial(
    MAT::PAR::BeamElastHyperMaterialParameterGeneric* params)
    : params_(params)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void MAT::BeamElastHyperMaterial::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void MAT::BeamElastHyperMaterial::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;

  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);

      /* the idea is that we have a generic type of material (this class), but various
       * possible sets of material parameters to 'feed' these very general constitutive relations */
      if (mat->Type() == INPAR::MAT::m_beam_reissner_elast_hyper or
          mat->Type() == INPAR::MAT::m_beam_reissner_elast_hyper_bymodes or
          mat->Type() == INPAR::MAT::m_beam_kirchhoff_elast_hyper or
          mat->Type() == INPAR::MAT::m_beam_kirchhoff_elast_hyper_bymodes or
          mat->Type() == INPAR::MAT::m_beam_kirchhoff_torsionfree_elast_hyper or
          mat->Type() == INPAR::MAT::m_beam_kirchhoff_torsionfree_elast_hyper_bymodes)
        params_ = static_cast<MAT::PAR::BeamElastHyperMaterialParameterGeneric*>(mat);
      else
        dserror("Type of material parameter %d does not fit to type of material law %d",
            mat->Type(), MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::Parameter* MAT::BeamElastHyperMaterial::Parameter() const { return params_; }

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
const MAT::PAR::BeamElastHyperMaterialParameterGeneric& MAT::BeamElastHyperMaterial::Params() const
{
  ThrowErrorIfParamsPointerIsNull();

  return *params_;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void MAT::BeamElastHyperMaterial::ThrowErrorIfParamsPointerIsNull() const
{
  if (params_ == NULL) dserror("pointer to parameter class is not set!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial::GetConstitutiveMatrixOfForcesMaterialFrame(
    LINALG::TMatrix<T, 3, 3>& C_N) const
{
  // defining material constitutive matrix CN between Gamma and N
  // according to Jelenic 1999, section 2.4
  C_N.Clear();

  C_N(0, 0) = Params().GetAxialRigidity();
  C_N(1, 1) = Params().GetShearRigidity2();
  C_N(2, 2) = Params().GetShearRigidity3();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial::GetConstitutiveMatrixOfMomentsMaterialFrame(
    LINALG::TMatrix<T, 3, 3>& C_M) const
{
  // defining material constitutive matrix CM between curvature and moment
  // according to Jelenic 1999, section 2.4
  C_M.Clear();

  C_M(0, 0) = Params().GetTorsionalRigidity();
  C_M(1, 1) = Params().GetBendingRigidity2();
  C_M(2, 2) = Params().GetBendingRigidity3();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
double MAT::BeamElastHyperMaterial::GetTranslationalMassInertiaFactor() const
{
  return Params().GetTranslationalMassInertia();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial::GetMassMomentOfInertiaTensorMaterialFrame(
    LINALG::TMatrix<T, 3, 3>& J) const
{
  J.Clear();

  J(0, 0) = Params().GetPolarMassMomentOfInertia();
  J(1, 1) = Params().GetMassMomentOfInertia2();
  J(2, 2) = Params().GetMassMomentOfInertia3();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
double MAT::BeamElastHyperMaterial::GetInteractionRadius() const
{
  return Params().GetInteractionRadius();
}

// explicit template instantiations
template void MAT::BeamElastHyperMaterial::GetConstitutiveMatrixOfForcesMaterialFrame<double>(
    LINALG::TMatrix<double, 3, 3>&) const;
template void
MAT::BeamElastHyperMaterial::GetConstitutiveMatrixOfForcesMaterialFrame<Sacado::Fad::DFad<double>>(
    LINALG::TMatrix<Sacado::Fad::DFad<double>, 3, 3>&) const;

template void MAT::BeamElastHyperMaterial::GetConstitutiveMatrixOfMomentsMaterialFrame<double>(
    LINALG::TMatrix<double, 3, 3>&) const;
template void
MAT::BeamElastHyperMaterial::GetConstitutiveMatrixOfMomentsMaterialFrame<Sacado::Fad::DFad<double>>(
    LINALG::TMatrix<Sacado::Fad::DFad<double>, 3, 3>&) const;

template void MAT::BeamElastHyperMaterial::GetMassMomentOfInertiaTensorMaterialFrame<double>(
    LINALG::TMatrix<double, 3, 3>&) const;
template void
MAT::BeamElastHyperMaterial::GetMassMomentOfInertiaTensorMaterialFrame<Sacado::Fad::DFad<double>>(
    LINALG::TMatrix<Sacado::Fad::DFad<double>, 3, 3>&) const;
