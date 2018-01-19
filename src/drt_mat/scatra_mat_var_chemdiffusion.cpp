/*----------------------------------------------------------------------*/
/*!
 \file scatra_mat_var_chemdiffusion.cpp

 \brief scatra material for chemical diffusion under a variational framework

   \level 3

   \maintainer  Jorge De Anda Salazar
                DeAnda@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/




#include <vector>
#include "scatra_mat_var_chemdiffusion.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraMatVarChemDiffusion::ScatraMatVarChemDiffusion(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: ScatraMat(matdata),
  chemdiff_model_(StringToModel(*matdata->Get<std::string>("MODEL"))),
  refMu_(matdata->GetDouble("REFMU")),
  refC_(matdata->GetDouble("REFC")),
  refTemp_(matdata->GetDouble("REFTEMP")),
  gasConstant_(matdata->GetDouble("GASCON"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ScatraMatVarChemDiffusion::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraMatVarChemDiffusion(this));
}


MAT::ScatraMatVarChemDiffusionType MAT::ScatraMatVarChemDiffusionType::instance_;

DRT::ParObject* MAT::ScatraMatVarChemDiffusionType::Create( const std::vector<char> & data )
{
  MAT::ScatraMatVarChemDiffusion* scatra_mat = new MAT::ScatraMatVarChemDiffusion();
  scatra_mat->Unpack(data);
  return scatra_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatVarChemDiffusion::ScatraMatVarChemDiffusion()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatVarChemDiffusion::ScatraMatVarChemDiffusion(MAT::PAR::ScatraMatVarChemDiffusion* params)
  :   ScatraMat(params),
      params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatVarChemDiffusion::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL)
    matid = params_->Id(); // in case we are in post-process mode
  AddtoPack(data, matid);

  // add base class material
  ScatraMat::Pack(data);

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatVarChemDiffusion::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ScatraMatVarChemDiffusion*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  ScatraMat::Unpack(basedata);

}

/*---------------------------------------------------------------------------*
 | convert string to model for chemical diffusion          deanda 08/17 |
 *---------------------------------------------------------------------------*/
MAT::PAR::Models MAT::PAR::ScatraMatVarChemDiffusion::StringToModel(const std::string& modelstring) const
{
  Models modelenum(chemdiff_undefined);

  // linear mass action law + Ficks law
  if(modelstring == "Linear")
    modelenum = chemdiff_linear;

  // logarithmic mass action law + Ficks law
  else if(modelstring == "Fickean")
    modelenum = chemdiff_fickean;

  // unknown model
  if ( modelenum == chemdiff_undefined )
    dserror("The variational chemical diffusion Model does not exist");

  return modelenum;
}

/*----------------------------------------------------------------------*
 | compute internal energy for chemical diffusion           deanda 08/17 |
 *----------------------------------------------------------------------*/
double MAT::ScatraMatVarChemDiffusion::ComputeInternalEnergy(
  const double concentration,     //!< Concentration
  const double refMu,             //!< Reference chemical potential
  const double refC,              //!< Reference concentration
  const double rt,                //!< factor RT
  const int orderderivative       //!< Order of the derivative to use
  ) const
{
  double IntEnergy(0.);

  switch(params_->chemdiff_model_)
  {
    // Computes linearized version of the Fickean Model
    case MAT::PAR::chemdiff_linear:
    {
      // Selects the corresponding derivative of the function to compute
        switch (orderderivative)
        {
        // Computes internal energy function
        case 0:
        {
          IntEnergy = refMu*concentration + (1/2)*(rt/refC)*(concentration - refC)*(concentration - refC);
          break;
        }
        // Computes internal energy first derivative
        case 1:
        {
          IntEnergy = refMu + (rt/refC)*(concentration - refC);
          break;
        }
        // Computes internal energy second derivative
        case 2:
        {
          IntEnergy = rt/refC;
          break;
        }
        // Computes internal energy higher derivatives
        default:
        {
          IntEnergy = 0.0;
          break;
        }
        }//switch (orderderivative)
        break;
      }//case MAT::PAR::chemdiff_linear:
    // Computes the fickean model
    case MAT::PAR::chemdiff_fickean:
      {
        // Selects the corresponding derivative of the function to compute
        switch (orderderivative)
        {
        // Computes internal energy function
        case 0:
        {
          const double arg_log = concentration/refC ;
          if (arg_log < 0)
          {
            std::cout<< "Concentration value is "<< concentration << std::endl;
            dserror("Negative argument for logarithm while computing Internal energy");
          }
          IntEnergy = refMu*concentration + (rt)*(concentration*log(arg_log) -concentration + refC);
          break;
        }
        // Computes internal energy first derivative
        case 1:
        {
          const double arg_log = concentration/refC;
          if (arg_log < 0)
          {
            std::cout<< "Concentration value is "<< concentration << std::endl;
            dserror("Negative argument for logarithm while computing Internal energy first derivative");
          }
          IntEnergy = refMu + (rt)*log(arg_log) ;
          break;
        }
        // Computes internal energy second derivative
        case 2:
        {
          IntEnergy = rt/concentration;
        break;
        }
        // Computes internal energy higher derivatives
        default:
        {
          dserror("Higher derivatives of the Internal energy have not being defined!");
            break;
        }
        }//switch (orderderivative)
        break;
    }//case MAT::PAR::chemdiff_fickean:
    default:
    {
      dserror("Constitutive Model for chemical diffusion not recognized!");
      break;
    }
  }// switch(params_->model_)

  return IntEnergy;
} // MAT::ScatraMatVarChemDiffusion::ComputeInternalEnergy


/*----------------------------------------------------------------------*
 | compute dual internal energy for chemical diffusion      deanda 08/17 |
 *----------------------------------------------------------------------*/
double MAT::ScatraMatVarChemDiffusion::ComputeDualInternalEnergy(
  const double ChemPot,       //!< Chemical Potential
  const double refMu,         //!< Reference chemical potential
  const double refC,          //!< Reference concentration
  const double rt,            //!< factor RT
  const int orderderivative   //!< Order of the derivative to use
  ) const
{
  double DualIntEnergy(0.);

  switch(params_->chemdiff_model_)
  {
    // Computes linearized version of the Fickean Model
    case MAT::PAR::chemdiff_linear:
    {
    // Selects the corresponding derivative of the function to compute
      switch (orderderivative)
      {
      // Computes internal energy function
      case 0:
      {
        DualIntEnergy = (ChemPot - refMu)*refC + (1/2)*(refC/rt)*(ChemPot - refMu)*(ChemPot - refMu);
        break;
      }
      // Computes internal energy first derivative
      case 1:
      {
        DualIntEnergy = refC + (refC/rt)*(ChemPot - refMu);
        break;
      }
      // Computes internal energy second derivative
      case 2:
      {
        DualIntEnergy = (refC/rt);
        break;
      }
      // Computes internal energy higher derivatives
      default:
      {
        DualIntEnergy = 0.0;
        break;
      }
      }//switch (orderderivative)
      break;
    }//case MAT::PAR::chemdiff_linear:

    // Computes the fickean model
    case MAT::PAR::chemdiff_fickean:
    {
        // Selects the corresponding derivative of the function to compute
          switch (orderderivative)
          {
          // Computes internal energy function
          case 0:
          {
          DualIntEnergy = (rt*refC)*exp((ChemPot - refMu)/rt) - 1.0;
            break;
          }
          // Computes internal energy first derivative
          case 1:
          {
          DualIntEnergy = refC*exp((ChemPot - refMu)/rt);
          break;
          }
          // Computes internal energy second derivative
          case 2:
          {
          DualIntEnergy = (1/rt)*refC*exp((ChemPot - refMu)/rt);
          break;
          }
          // Computes internal energy higher derivatives
          default:
          {
            dserror("Higher derivatives of the Dual Internal energy have not being defined!");
              break;
          }
          }//switch (orderderivative)
      break;
    }//case MAT::PAR::chemdiff_fickean:
    default:
    {
      dserror("Constitutive Model for chemical diffusion not recognized!");
      break;
    }
  }// switch(params_->model_)

  return DualIntEnergy;
} // MAT::ScatraMatVarChemDiffusion::ComputeDualInternalEnergy

///*----------------------------------------------------------------------*
// | compute dissipation function for chemical diffusion      deanda 01/18 |
// *----------------------------------------------------------------------*/

//! Set Dissipation Potential
template<int NSD>
void MAT::ScatraMatVarChemDiffusion::ComputeDissipationPot(
   double concentration                        //!< Concentration at t_(n+1)
  ,const double& refC                          //!< Reference concentration
  ,const LINALG::Matrix<NSD,1>& chemicalfield  //!< Chemical field at t_(n+1)
  ,const double& mobility                      //!< mobility
  ,double& DissipationPot                      //!< Output
  )const
{
  if ( params_->chemdiff_model_ == MAT::PAR::chemdiff_linear)
    concentration = refC;
  DissipationPot = (0.5)*(concentration*mobility)* chemicalfield.Dot(chemicalfield);
  return;
}

//! Set Dissipation Potential 1st derivative
template<int NSD>
void MAT::ScatraMatVarChemDiffusion::ComputeDissipationPot_D1(
   double concentration                         //!< Concentration at t_(n+1)
  ,const double& refC                           //!< Reference concentration
  ,const LINALG::Matrix<NSD,1>& chemicalfield   //!< Chemical field at t_(n+1)
  ,const double& mobility                       //!< mobility
  ,LINALG::Matrix<NSD,1>& DissipationPot1deriv  //!< Output
)const
{
  DissipationPot1deriv.PutScalar(0.);
  if ( params_->chemdiff_model_ == MAT::PAR::chemdiff_linear)
    concentration = refC;
  DissipationPot1deriv.Update(concentration*mobility,chemicalfield,0.);
  return;
}

//! Set Dissipation Potential 2nd derivative
template<int NSD>
void MAT::ScatraMatVarChemDiffusion::ComputeDissipationPot_D2(
   double concentration                             //!< Concentration at t_(n+1)
  ,const double& refC                               //!< Reference concentration
  ,const LINALG::Matrix<NSD,1>& chemicalfield       //!< Chemical field at t_(n+1)
  ,const double&  mobility                          //!< mobility
 ,LINALG::Matrix<NSD,NSD>& DissipationPot2deriv     //!< Output
)const
{
  DissipationPot2deriv.PutScalar(0.);
  if ( params_->chemdiff_model_ == MAT::PAR::chemdiff_linear)
    concentration = refC;
  for (int d=0; d<NSD; ++d)
    DissipationPot2deriv(d,d) = (concentration*mobility);
  return;
}

//ComputeDissipationPot
template void MAT::ScatraMatVarChemDiffusion::ComputeDissipationPot<1>(double, const double&, const LINALG::Matrix<1,1>&, const double&, double& )const;//
template void MAT::ScatraMatVarChemDiffusion::ComputeDissipationPot<2>(double, const double&, const LINALG::Matrix<2,1>&, const double&, double& )const;
template void MAT::ScatraMatVarChemDiffusion::ComputeDissipationPot<3>(double, const double&, const LINALG::Matrix<3,1>&, const double&, double& )const;
//ComputeDissipationPot_D1
template void MAT::ScatraMatVarChemDiffusion::ComputeDissipationPot_D1<1>(double, const double&, const LINALG::Matrix<1,1>&, const double&, LINALG::Matrix<1,1>&)const;
template void MAT::ScatraMatVarChemDiffusion::ComputeDissipationPot_D1<2>(double, const double&, const LINALG::Matrix<2,1>&, const double&, LINALG::Matrix<2,1>&)const;
template void MAT::ScatraMatVarChemDiffusion::ComputeDissipationPot_D1<3>(double, const double&, const LINALG::Matrix<3,1>&, const double&, LINALG::Matrix<3,1>&)const;
//ComputeDissipationPot_D2
template void MAT::ScatraMatVarChemDiffusion::ComputeDissipationPot_D2<1>(double, const double&, const LINALG::Matrix<1,1>&, const double& , LINALG::Matrix<1,1>&)const;
template void MAT::ScatraMatVarChemDiffusion::ComputeDissipationPot_D2<2>(double, const double&, const LINALG::Matrix<2,1>&, const double& , LINALG::Matrix<2,2>&)const;
template void MAT::ScatraMatVarChemDiffusion::ComputeDissipationPot_D2<3>(double, const double&, const LINALG::Matrix<3,1>&, const double& , LINALG::Matrix<3,3>&)const;
