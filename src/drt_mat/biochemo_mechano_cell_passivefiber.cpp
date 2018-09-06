/*!----------------------------------------------------------------------
\file biochemo_mechano_cell_passivefiber.cpp

\brief Implementation of Biochemo-Mechano Coupled passive, viscoelastic material model for the cell.

\level 3

\maintainer Andreas Rauch

*----------------------------------------------------------------------*/

#include "biochemo_mechano_cell_passivefiber.H"

#include "material_service.H"
#include "matpar_bundle.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::BioChemoMechanoCellPassiveFiber::BioChemoMechanoCellPassiveFiber(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      idmatelast_(matdata->GetInt("IDMATELAST")),
      mu_(matdata->GetDouble("VISC")),
      analyticalmaterialtangent_(true)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  if (problem->ProblemType() == prb_immersed_cell)
  {
    if (DRT::INPUT::IntegralValue<int>(
            problem->CellMigrationParams().sublist("STRUCTURAL DYNAMIC"), "MATERIALTANGENT"))
      analyticalmaterialtangent_ = false;
  }
  else
  {
    if (DRT::INPUT::IntegralValue<int>(problem->StructuralDynamicParams(), "MATERIALTANGENT"))
      analyticalmaterialtangent_ = false;
  }
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::BioChemoMechanoCellPassiveFiber::CreateMaterial()
{
  return Teuchos::rcp(new MAT::BioChemoMechanoCellPassiveFiber(this));
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::BioChemoMechanoCellPassiveFiberType MAT::BioChemoMechanoCellPassiveFiberType::instance_;


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::BioChemoMechanoCellPassiveFiberType::Create(const std::vector<char>& data)
{
  MAT::BioChemoMechanoCellPassiveFiber* material = new MAT::BioChemoMechanoCellPassiveFiber();
  material->Unpack(data);
  return material;
}


/*----------------------------------------------------------------------*
 |  Constructor                                            rauch  01/17 |
 *----------------------------------------------------------------------*/
MAT::BioChemoMechanoCellPassiveFiber::BioChemoMechanoCellPassiveFiber()
    : params_(NULL), isinit_(false)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                                       rauch  01/17 |
 *----------------------------------------------------------------------*/
MAT::BioChemoMechanoCellPassiveFiber::BioChemoMechanoCellPassiveFiber(
    MAT::PAR::BioChemoMechanoCellPassiveFiber* params)
    : params_(params), isinit_(false)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                                   rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::Pack(DRT::PackBuffer& data) const
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

  // pack history data
  int histsize = -1;
  // if material is not initialized, i.e. start simulation, nothing to pack
  if (!Initialized())
  {
    histsize = 0;
  }
  else
  {
    // if material is initialized (restart): size equates number of gausspoints
    histsize = histdefgrdlast_->size();
  }

  AddtoPack(data, histsize);  // length of history vector(s)

  for (int var = 0; var < histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data, histdefgrdlast_->at(var));
  }

  // Pack data of elastic material
  if (matelast_ != Teuchos::null)
  {
    matelast_->Pack(data);
  }

  return;

}  // Pack()


/*----------------------------------------------------------------------*
 |  Unpack                                                 rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::Unpack(const std::vector<char>& data)
{
  // construct current deformation gradient history variable
  histdefgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  isinit_ = true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid = -1;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::BioChemoMechanoCellPassiveFiber*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // history data
  int histsize = -1;
  ExtractfromPack(position, data, histsize);

  // if system is not yet initialized, the history vectors have to be initialized
  if (histsize == 0)
  {
    isinit_ = false;
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
    return;
  }

  // construct vector of old deformation gradient matrices
  histdefgrdlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  for (int var = 0; var < histsize; ++var)
  {
    // initialize
    LINALG::Matrix<3, 3> tmp_matrix3x3(true);

    // matrices of last converged state are unpacked
    ExtractfromPack(position, data, tmp_matrix3x3);
    histdefgrdlast_->push_back(tmp_matrix3x3);

    // initialize current deformation gradient
    histdefgrdcurr_->push_back(tmp_matrix3x3);
  }

  // Unpack data of passive elastic material
  std::vector<char> dataelastic;
  ExtractfromPack(position, data, dataelastic);
  if (dataelastic.size() > 0)
  {
    DRT::ParObject* o = DRT::UTILS::Factory(dataelastic);  // Unpack is done here
    MAT::So3Material* matel = dynamic_cast<MAT::So3Material*>(o);
    if (matel == NULL) dserror("failed to unpack passive material");
    matelast_ = Teuchos::rcp(matel);
  }
  else
    matelast_ = Teuchos::null;

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;

}  // Unpack()


/*----------------------------------------------------------------------*
 | initialize / allocate internal variables (public)       rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // construct history variables
  histdefgrdlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  histdefgrdlast_->resize(numgp);

  histdefgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  histdefgrdcurr_->resize(numgp);

  LINALG::Matrix<3, 3> emptymat3x3(true);
  for (int i = 0; i < 3; i++) emptymat3x3(i, i) = 1.0;

  for (int i = 0; i < numgp; i++)
  {
    histdefgrdlast_->at(i) = emptymat3x3;
    histdefgrdcurr_->at(i) = emptymat3x3;
  }

  // Setup of passive material
  matelast_ =
      Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Material::Factory(params_->idmatelast_));
  matelast_->Setup(numgp, linedef);

  isinit_ = true;
  return;
}  // Setup()


/*----------------------------------------------------------------------*
 |  ResetAll                                               rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::ResetAll(const int numgp)
{
  // construct history variables
  histdefgrdlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  histdefgrdlast_->resize(numgp);
  histdefgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  histdefgrdcurr_->resize(numgp);

  LINALG::Matrix<3, 3> emptymat3x3(true);
  for (int i = 0; i < 3; i++) emptymat3x3(i, i) = 1.0;

  for (int i = 0; i < numgp; i++)
  {
    histdefgrdlast_->at(i) = emptymat3x3;
    histdefgrdcurr_->at(i) = emptymat3x3;
  }

  matelast_->ResetAll(numgp);
  isinit_ = false;
  return;
}  // ResetAll()


/*----------------------------------------------------------------------*
 |  Update internal variables                              rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::Update()
{
  // make current values at time step t^{n+1} to old values at time  t^n
  histdefgrdlast_ = histdefgrdcurr_;
  histdefgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  histdefgrdcurr_->resize(histdefgrdlast_->size());

  matelast_->Update();
  return;
}  // Update()


/*----------------------------------------------------------------------*
 |  Reset internal variables                               rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::ResetStep()
{
  matelast_->ResetStep();
  return;
}  // ResetStep()


/*----------------------------------------------------------------------*
 |  Evaluate Material                                      rauch  01/17 |
 *--------------------------------------------------------------------- */
void MAT::BioChemoMechanoCellPassiveFiber::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int eleGID)
{
  // get gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");
  // get time algorithmic parameters
  double dt = params.get<double>("delta time", -1.0);
  if (dt == -1.0) dserror("no time step size provided in material");


  // copy deformation gradient in histroy variable
  histdefgrdcurr_->at(gp) = *defgrd;


  // calc the viscosity tensor
  const double viscosity = params_->mu_;


  // setup inverse of deformation gradient
  LINALG::Matrix<3, 3> invdefgrd(*defgrd);
  invdefgrd.Invert();
  // setup deformation gradient rate, rotation tensor, strain rate and rotation rate
  // \dot{F} = \frac {F^n - F^{n-1}} {\Delta t}
  LINALG::Matrix<3, 3> defgrdrate(true);
  // R = F * U^{-1}
  LINALG::Matrix<3, 3> R;
  // \dot{\epsilon} = d = 0.5 * (\dot{F}F^{-1} + (\dot{F}F^{-1})^{T}
  LINALG::Matrix<6, 1> strainrate;
  // calc the rates
  SetupRates(*defgrd, invdefgrd, params, defgrdrate, R, strainrate, gp, dt);



  // temporary viscous stress tensor
  LINALG::Matrix<6, 1> visc_stress_cauchy;
  // calc viscous stress contribution
  visc_stress_cauchy(0) = strainrate(0);
  visc_stress_cauchy(1) = strainrate(1);
  visc_stress_cauchy(2) = strainrate(2);
  visc_stress_cauchy(3) = strainrate(3);
  visc_stress_cauchy(4) = strainrate(4);
  visc_stress_cauchy(5) = strainrate(5);
  visc_stress_cauchy.Scale(2.0 * viscosity);

  // Transform Cauchy stress to PK2 stress
  // S = J * F^{-1} \sigma F^{-T}
  LINALG::Matrix<NUM_STRESS_3D, 1> visc_stress_PK2(true);  // 6x1
  LINALG::Matrix<3, 3> visc_stress_cauchy_mat(true);       // 3x3
  CauchytoPK2(visc_stress_PK2, visc_stress_cauchy_mat, *defgrd, invdefgrd, visc_stress_cauchy);


  /////////////////////////////////////////////////////////////////
  // Calculate constitutive tensor
  //
  //  Cmat = 2*dS/dC =
  //
  //   d ( J * F^{-1} * \sigma * F^{-T} )
  //  ------------------------------------   =
  //   d               C
  //
  // dJ/dC * F^{-1} * \sigma * F^{-T}  +
  //
  // J * (d F^{-1}/dC) * \sigma * F^{-T}  +
  //
  // J * F^{-1} * (d \sigma / dC) * F^{-T}  +
  //
  // J * F^{-1} * \sigma * (d F^{-T} / dC)
  //
  //
  /////////////////////////////////////////////////////////////////

  // viscous part of material tangent
  LINALG::Matrix<6, 6> visc_cmat;

  if (params_->analyticalmaterialtangent_)
  {
    // Jacobi Determinant
    const double detF = defgrd->Determinant();

    // theta is hard coded to be 1.0 (backward-euler)
    const double theta = 1.0;

    // Right-Cauchy-Green tensor(3x3): C = F^{T} * F
    LINALG::Matrix<3, 3> C;
    C.MultiplyTN(*defgrd, *defgrd);
    // Inverse of C: C^{-1}
    LINALG::Matrix<3, 3> Cinv(C);
    Cinv.Invert();
    // Root of C: \sqrt{C}
    LINALG::Matrix<3, 3> RootC(C);
    MatrixRoot3x3(RootC);
    // Inverse of sqrt(C): \sqrt(C)^{-1}
    LINALG::Matrix<3, 3> RootCInv(RootC);
    RootCInv.Invert();

    // Derivative of \sqrt{C} with respect to C: DerviC =  d sqrt(C) / d C
    LINALG::Matrix<6, 6> dsqrtCdC(true);  // 6x6 Voigt matrix
    double dsqrtCdC_Tensor[3][3][3][3] = {{{{0.}}}};
    double dsqrtCinvdC_Tensor[3][3][3][3] = {{{{0.}}}};
    MatrixRootDerivativeSym3x3(C, dsqrtCdC);
    Setup4Tensor(dsqrtCdC_Tensor, dsqrtCdC);  // 3x3x3x3 Tensor
    // -sqrt(C)^{-1}*(d sqrt(C)/dC)*sqrt(C)^{-1}
    MatrixInverseOfRootDerivative(dsqrtCdC_Tensor,  // d sqrt(C)/dC
        RootCInv,                                   // sqrt(C)^{-1}
        dsqrtCinvdC_Tensor                          // result
    );

    //    ///////////////////////////////////
    //    // FD for d sqrt{C}^{-1}/ d C
    //    ///////////////////////////////////
    //    double delta    = 1e-08;
    //    double invdelta = 1e+08;
    //    LINALG::Matrix<3,3> Ccopy(C);
    //    LINALG::Matrix<3,3> RootCcopy(true);
    //    LINALG::Matrix<3,3> RootCcopy2(true);
    //    LINALG::Matrix<3,3> InvRootCcopy(true);
    //    LINALG::Matrix<3,3> InvRootCcopy2(true);
    //
    //    double dsqrtCinvdC_Tensor_fd[3][3][3][3]  = {{{{0.}}}};
    //    double dsqrtCdC_Tensor_fd[3][3][3][3]     = {{{{0.}}}};
    //
    //    for(int k=0;k<3;++k)
    //    {
    //      for(int l=0;l<3;++l)
    //      {
    //        // toggle C_kl with positive direction
    //        Ccopy(k,l)+=delta/2.0;
    //        if(l != k)
    //          Ccopy(l,k)+=delta/2.0;
    //
    //        // calc root of toggled copy of C
    //        RootCcopy.Update(Ccopy);
    //        MatrixRoot3x3(RootCcopy);
    //        // calc inverse
    //        InvRootCcopy.Update(RootCcopy);
    //        InvRootCcopy.Invert();
    //
    //        // toggle C_kl with negative direction
    //        Ccopy(k,l)-=delta;
    //        if(l != k)
    //          Ccopy(l,k)-=delta;
    //
    //        // calc root of toggled copy of C
    //        RootCcopy2.Update(Ccopy);
    //        MatrixRoot3x3(RootCcopy2);
    //        // calc inverse
    //        InvRootCcopy2.Update(RootCcopy2);
    //        InvRootCcopy2.Invert();
    //
    //        // finite difference
    //        InvRootCcopy.Update(-1.0,InvRootCcopy2,1.0);
    //        InvRootCcopy.Scale(invdelta);
    //
    //        RootCcopy.Update(-1.0,RootCcopy2,1.0);
    //        RootCcopy.Scale(invdelta);
    //
    //        // fill 4-Tensor
    //        for(int i=0;i<3;++i)
    //        {
    //          for (int j=0;j<3;++j)
    //          {
    //            dsqrtCinvdC_Tensor_fd [i][j][k][l] = InvRootCcopy(i,j);
    //            dsqrtCdC_Tensor_fd [i][j][k][l]    = RootCcopy(i,j);
    //          } // j loop
    //        } // i loop
    //
    //        // reset Ccopy
    //        Ccopy.Update(C);
    //
    //      } // l loop
    //    } // k loop

    //    //////////////
    //    // FD CHECK //
    //    //////////////
    //    for (int i=0; i<3; i++)
    //      for (int j=0; j<3; j++)
    //        for (int k=0; k<3; k++)
    //          for (int l=0; l<3; l++)
    //            if(abs(dsqrtCinvdC_Tensor[i][j][k][l])-abs(dsqrtCinvdC_Tensor_fd[i][j][k][l])>(delta*5.0))
    //            {
    //              std::cout<<"dsqrt(C)^{-1}/dC = "<<dsqrtCinvdC_Tensor[i][j][k][l]<<" Approx =
    //              "<<dsqrtCinvdC_Tensor_fd[i][j][k][l]<<" at ijkl="<<i<<j<<k<<l<<" err
    //              ="<<abs(dsqrtCinvdC_Tensor[i][j][k][l])-abs(dsqrtCinvdC_Tensor_fd[i][j][k][l])<<std::endl;
    //            }

    //    //////////////
    //    // FD CHECK //
    //    //////////////
    //    for (int i=0; i<3; i++)
    //      for (int j=0; j<3; j++)
    //        for (int k=0; k<3; k++)
    //          for (int l=0; l<3; l++)
    //            if(abs(dsqrtCdC_Tensor[i][j][k][l])-abs(dsqrtCdC_Tensor_fd[i][j][k][l])>(delta*5.0))
    //            {
    //              std::cout<<"dsqrt(C)/dC = "<<dsqrtCdC_Tensor[i][j][k][l]<<" Approx =
    //              "<<dsqrtCdC_Tensor_fd[i][j][k][l]<<" at ijkl="<<i<<j<<k<<l<<" err
    //              ="<<abs(dsqrtCdC_Tensor[i][j][k][l])-abs(dsqrtCdC_Tensor_fd[i][j][k][l])<<std::endl;
    //            }


    // Setup transposed matrices
    LINALG::Matrix<3, 3> Rtrans;
    Rtrans.UpdateT(R);
    LINALG::Matrix<3, 3> defgrdratetrans;
    defgrdratetrans.UpdateT(defgrdrate);
    LINALG::Matrix<3, 3> invdefgrdtrans;
    invdefgrdtrans.UpdateT(invdefgrd);

    // 3x3x3x3 Tensor auxiliary variables
    double tens1[3][3][3][3] = {{{{0.}}}};
    double tens2[3][3][3][3] = {{{{0.}}}};
    double tens3[3][3][3][3] = {{{{0.}}}};
    double tens4[3][3][3][3] = {{{{0.}}}};
    double tens5[3][3][3][3] = {{{{0.}}}};
    double tens6[3][3][3][3] = {{{{0.}}}};
    double tens7[3][3][3][3] = {{{{0.}}}};
    double auxtens[3][3][3][3] = {{{{0.}}}};
    // 3x3 matrix auxiliary variables
    LINALG::Matrix<3, 3> tempmat1;
    LINALG::Matrix<3, 3> tempmat2;
    // 6x6 matrix auxiliary variables
    LINALG::Matrix<6, 6> auxvoigt;

    // F^{-1} * \sigma
    tempmat1.MultiplyNN(invdefgrd, visc_stress_cauchy_mat);

    // F^{-1} * \sigma * R
    tempmat2.MultiplyNN(tempmat1, R);

    // F^{-1} * \sigma * (dF^{-T}/ dC) =
    // F^{-1} * \sigma * R * (d sqrt(C)^{-1} / dC) =
    // F^{-1} * \sigma * R * (-sqrt(C)^{-1}*(d sqrt(C)/dC)*sqrt(C)^{-1})
    MultMatrixFourTensor(tens1, tempmat2, dsqrtCinvdC_Tensor, false);

    // (F^{-1}*\sigma*F^{-T}) * (dJ/dC) :
    // (F^{-1} \sigma F^{-T}) dyad 0.5*C^{-1} (dyadic product to obtain 4-Tensor)
    tempmat2.MultiplyNT(tempmat1, invdefgrd);
    MAT::ElastSymTensorMultiply(auxvoigt, 0.5, tempmat2, Cinv, 0.0);
    Setup4Tensor(tens2, auxvoigt);

    // (d F^{-1}/dC) * \sigma * F^{-T} =
    // (d sqrt(C)^{-1} / dC) * R^{T} * \sigma * F^{-T}
    tempmat1.MultiplyNT(visc_stress_cauchy_mat, invdefgrd);
    tempmat2.MultiplyTN(R, tempmat1);
    MultFourTensorMatrix(tens3, tempmat2, dsqrtCinvdC_Tensor, false);

    /////////////////////////////////////////////////////////////////////
    // velocity gradient with respect to right cauchy-green
    // d d / d C = d 0.5*(l+l^T) / dC
    //
    // l   = \dot F F^{-1}
    // l^T = F^{-T} \dot F^T
    //
    // dl / dC = (d\dotF / dC)F^{-1}  +  \dot F (dF^{-1} / dC)
    //
    // dl^T / dC = (dF^{-T} / dC) \dot F^T  +  F^{-T} (d \dot F^T / dC)
    //
    /////////////////////////////////////////////////////////////////////

    // F^{-1} [\dot F  (d F^{-1} / dC)] F^{-T}
    // F^{-1} {[F^dot * [d sqrt(C)^{-1} / dC] * R^{T}]} F^{-T}
    tempmat1.MultiplyNN(invdefgrd, defgrdrate);
    tempmat2.MultiplyTN(R, invdefgrdtrans);
    MultMatrixFourTensor(auxtens, tempmat1, dsqrtCinvdC_Tensor, false);
    MultFourTensorMatrix(tens4, tempmat2, auxtens, false);

    // F^{-1} [(d\dotF / dC) F^{-1}] F^{-T} =
    // 1/dt F^{-1} R (d sqrt(C) / dC) F^{-1} F^{-T}
    tempmat1.MultiplyNN(invdefgrd, R);
    tempmat1.Scale(1.0 / (theta * dt));
    tempmat2.MultiplyNN(invdefgrd, invdefgrdtrans);
    MultMatrixFourTensor(auxtens, tempmat1, dsqrtCdC_Tensor, true);
    MultFourTensorMatrix(tens5, tempmat2, auxtens, false);

    // F^{-1} [ (d F^{-T} / dC) \dot F^T ] F^{-T}
    // F^{-1} {R [d sqrt(C)^{-1} / dC] \dot F^T} F^{-T}
    tempmat1.MultiplyNN(invdefgrd, R);
    tempmat2.MultiplyTN(defgrdrate, invdefgrdtrans);
    MultMatrixFourTensor(auxtens, tempmat1, dsqrtCinvdC_Tensor, true);
    MultFourTensorMatrix(tens6, tempmat2, auxtens, false);

    // F^{-1} [ F^{-T} (d\dotF^T / dC) ] F^{-T} =
    // 1/dt F^{-1} [ F^{-T} (d sqrt(C) / dC) R^T] F^{-T}
    tempmat1.MultiplyNN(invdefgrd, invdefgrdtrans);
    tempmat1.Scale(1.0 / (theta * dt));
    tempmat2.MultiplyTN(R, invdefgrdtrans);
    MultMatrixFourTensor(auxtens, tempmat1, dsqrtCdC_Tensor, true);
    MultFourTensorMatrix(tens7, tempmat2, auxtens, false);

    // FINALLY SUM IT ALL UP

    // put together viscous constitutive tensor
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
            auxtens[i][j][k][l] =
                tens1[i][j][k][l] + tens2[i][j][k][l] + tens3[i][j][k][l] +
                (tens4[i][j][k][l] + tens5[i][j][k][l] + tens6[i][j][k][l] + tens7[i][j][k][l]) *
                    viscosity;

    // convert tensor to voigt notation
    Setup6x6VoigtMatrix(visc_cmat, auxtens);

    // every part of the linearization needs to be scaled with 2*J
    // this factor has been neglected in the single terms.
    visc_cmat.Scale(2.0 * detF);

  }  // if analyticalmaterialtangent

  // Evaluate elastic stress and material tangent
  matelast_->Evaluate(defgrd, glstrain, params, stress, cmat, eleGID);

  // sum up total stress
  stress->Update(1.0, visc_stress_PK2, 1.0);
  if (params_->analyticalmaterialtangent_) cmat->Update(1.0, visc_cmat, 1.0);

  return;
}  // MAT::BioChemoMechanoCellPassiveFiber::Evaluate


/*----------------------------------------------------------------------------------*
 | Calculation of deformation gradient rate, rotation tensor, strain rate and       |
 | rotation rate (finite difference scheme)                            rauch  01/17 |
 *----------------------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::SetupRates(LINALG::Matrix<3, 3> defgrd,
    LINALG::Matrix<3, 3>& invdefgrd, Teuchos::ParameterList& params,
    LINALG::Matrix<3, 3>& defgrdrate, LINALG::Matrix<3, 3>& R, LINALG::Matrix<6, 1>& strainrate,
    const int& gp, const double& dt)
{
  // Read history
  LINALG::Matrix<3, 3> defgrdlast;
  defgrdlast = histdefgrdlast_->at(gp);

  // Rate of deformation gradient: \dot{F} = \frac {F^{n+1} - F^{n}} {\Delta t}
  defgrdrate.Update(1.0, defgrd, 0.0);
  defgrdrate.Update(-1.0, defgrdlast, 1.0);
  defgrdrate.Scale(1.0 / dt);

  // Calculate velocity gradient l = \dot{F}F^{-1}
  LINALG::Matrix<3, 3> velgradient;
  velgradient.MultiplyNN(defgrdrate, invdefgrd);

  // Rate of strain/symmetric part of velocity gradient
  // d = 0.5 * (l + l^{T}) = 0.5 * (\dot{F}F^{-1} + (\dot{F}F^{-1})^{T})
  // Remark: strain-like 6-Voigt vector
  strainrate(0) = velgradient(0, 0) + velgradient(0, 0);
  strainrate(1) = velgradient(1, 1) + velgradient(1, 1);
  strainrate(2) = velgradient(2, 2) + velgradient(2, 2);
  strainrate(3) = velgradient(0, 1) + velgradient(1, 0);
  strainrate(4) = velgradient(1, 2) + velgradient(2, 1);
  strainrate(5) = velgradient(0, 2) + velgradient(2, 0);
  strainrate.Scale(0.5);

  // this means a lot of computational cost.
  // we do this only if analytical material tangent is requested.
  if (params_->analyticalmaterialtangent_)
  {
    // Rate of rotation tensor (!= skew part w of velocity gradient l, see Holzapfel S.99)
    // Determine rotation tensor R from F (F=R*U) -> polar decomposition of displacement based F
    LINALG::Matrix<3, 3> Q(true);
    LINALG::Matrix<3, 3> S(true);
    LINALG::Matrix<3, 3> VT(true);

    // Calculate rotcurr from defgrd
    // Singular Value Decomposition,analogously to micromaterial_evaluate.cpp lines 81ff
    LINALG::SVD<3, 3>(defgrd, Q, S, VT);
    R.MultiplyNN(Q, VT);
  }

}  // SetupRates


/*----------------------------------------------------------------------*
 | pull back of spatial stresses                           rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::CauchytoPK2(LINALG::Matrix<6, 1>& Sactive,
    LINALG::Matrix<3, 3>& cauchystress, LINALG::Matrix<3, 3> defgrd,
    LINALG::Matrix<3, 3>& invdefgrd, LINALG::Matrix<6, 1>& sigma)
{
  // calculate the Jacobi-determinant
  const double detF = defgrd.Determinant();

  // Convert stress like 6x1-Voigt vector to 3x3 matrix
  cauchystress(0, 0) = sigma(0);
  cauchystress(0, 1) = sigma(3);
  cauchystress(0, 2) = sigma(5);
  cauchystress(1, 0) = cauchystress(0, 1);
  cauchystress(1, 1) = sigma(1);
  cauchystress(1, 2) = sigma(4);
  cauchystress(2, 0) = cauchystress(0, 2);
  cauchystress(2, 1) = cauchystress(1, 2);
  cauchystress(2, 2) = sigma(2);

  // S = J * F^{-1} * sigma * F^{-T}
  LINALG::Matrix<3, 3> temp(true);
  LINALG::Matrix<3, 3> S(true);
  temp.MultiplyNN(invdefgrd, cauchystress);
  S.MultiplyNT(temp, invdefgrd);
  S.Scale(detF);

  // Sactive is stress like 6x1-Voigt vector
  Sactive(0) = S(0, 0);
  Sactive(1) = S(1, 1);
  Sactive(2) = S(2, 2);
  Sactive(3) = S(0, 1);
  Sactive(4) = S(1, 2);
  Sactive(5) = S(0, 2);

}  // CauchytoPK2()


/*----------------------------------------------------------------------*
 |  Names of gp data to be visualized                      rauch  01/17 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::VisNames(std::map<std::string, int>& names)
{
  matelast_->VisNames(names);
  return;
}  // VisNames()


/*----------------------------------------------------------------------*
 |  gp data to be visualized                               rauch  01/17 |
 *----------------------------------------------------------------------*/
bool MAT::BioChemoMechanoCellPassiveFiber::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  return matelast_->VisData(name, data, numgp, eleID);
}


/*-------------------------------------------------------------------------------------*
 |  Setup 4-Tensor from 6x6 Voigt notation                                rauch  07/14 |
 *-------------------------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::Setup4Tensor(
    double FourTensor[3][3][3][3], LINALG::Matrix<6, 6>& VoigtMatrix)
{
  // Setup 4-Tensor from 6x6 Voigt matrix
  // Voigt matrix has to be the representative of a 4 tensor with
  // at least minor symmetries.
  FourTensor[0][0][0][0] = VoigtMatrix(0, 0);  // C1111
  FourTensor[0][0][1][1] = VoigtMatrix(0, 1);  // C1122
  FourTensor[0][0][2][2] = VoigtMatrix(0, 2);  // C1133
  FourTensor[0][0][0][1] = VoigtMatrix(0, 3);  // C1112
  FourTensor[0][0][1][0] = VoigtMatrix(0, 3);  // C1121
  FourTensor[0][0][1][2] = VoigtMatrix(0, 4);  // C1123
  FourTensor[0][0][2][1] = VoigtMatrix(0, 4);  // C1132
  FourTensor[0][0][0][2] = VoigtMatrix(0, 5);  // C1113
  FourTensor[0][0][2][0] = VoigtMatrix(0, 5);  // C1131

  FourTensor[1][1][0][0] = VoigtMatrix(1, 0);  // C2211
  FourTensor[1][1][1][1] = VoigtMatrix(1, 1);  // C2222
  FourTensor[1][1][2][2] = VoigtMatrix(1, 2);  // C2233
  FourTensor[1][1][0][1] = VoigtMatrix(1, 3);  // C2212
  FourTensor[1][1][1][0] = VoigtMatrix(1, 3);  // C2221
  FourTensor[1][1][1][2] = VoigtMatrix(1, 4);  // C2223
  FourTensor[1][1][2][1] = VoigtMatrix(1, 4);  // C2232
  FourTensor[1][1][0][2] = VoigtMatrix(1, 5);  // C2213
  FourTensor[1][1][2][0] = VoigtMatrix(1, 5);  // C2231

  FourTensor[2][2][0][0] = VoigtMatrix(2, 0);  // C3311
  FourTensor[2][2][1][1] = VoigtMatrix(2, 1);  // C3322
  FourTensor[2][2][2][2] = VoigtMatrix(2, 2);  // C3333
  FourTensor[2][2][0][1] = VoigtMatrix(2, 3);  // C3312
  FourTensor[2][2][1][0] = VoigtMatrix(2, 3);  // C3321
  FourTensor[2][2][1][2] = VoigtMatrix(2, 4);  // C3323
  FourTensor[2][2][2][1] = VoigtMatrix(2, 4);  // C3332
  FourTensor[2][2][0][2] = VoigtMatrix(2, 5);  // C3313
  FourTensor[2][2][2][0] = VoigtMatrix(2, 5);  // C3331

  FourTensor[0][1][0][0] = VoigtMatrix(3, 0);
  FourTensor[1][0][0][0] = VoigtMatrix(3, 0);  // C1211 = C2111
  FourTensor[0][1][1][1] = VoigtMatrix(3, 1);
  FourTensor[1][0][1][1] = VoigtMatrix(3, 1);  // C1222 = C2122
  FourTensor[0][1][2][2] = VoigtMatrix(3, 2);
  FourTensor[1][0][2][2] = VoigtMatrix(3, 2);  // C1233 = C2133
  FourTensor[0][1][0][1] = VoigtMatrix(3, 3);
  FourTensor[1][0][0][1] = VoigtMatrix(3, 3);  // C1212 = C2112
  FourTensor[0][1][1][0] = VoigtMatrix(3, 3);
  FourTensor[1][0][1][0] = VoigtMatrix(3, 3);  // C1221 = C2121
  FourTensor[0][1][1][2] = VoigtMatrix(3, 4);
  FourTensor[1][0][1][2] = VoigtMatrix(3, 4);  // C1223 = C2123
  FourTensor[0][1][2][1] = VoigtMatrix(3, 4);
  FourTensor[1][0][2][1] = VoigtMatrix(3, 4);  // C1232 = C2132
  FourTensor[0][1][0][2] = VoigtMatrix(3, 5);
  FourTensor[1][0][0][2] = VoigtMatrix(3, 5);  // C1213 = C2113
  FourTensor[0][1][2][0] = VoigtMatrix(3, 5);
  FourTensor[1][0][2][0] = VoigtMatrix(3, 5);  // C1231 = C2131

  FourTensor[1][2][0][0] = VoigtMatrix(4, 0);
  FourTensor[2][1][0][0] = VoigtMatrix(4, 0);  // C2311 = C3211
  FourTensor[1][2][1][1] = VoigtMatrix(4, 1);
  FourTensor[2][1][1][1] = VoigtMatrix(4, 1);  // C2322 = C3222
  FourTensor[1][2][2][2] = VoigtMatrix(4, 2);
  FourTensor[2][1][2][2] = VoigtMatrix(4, 2);  // C2333 = C3233
  FourTensor[1][2][0][1] = VoigtMatrix(4, 3);
  FourTensor[2][1][0][1] = VoigtMatrix(4, 3);  // C2312 = C3212
  FourTensor[1][2][1][0] = VoigtMatrix(4, 3);
  FourTensor[2][1][1][0] = VoigtMatrix(4, 3);  // C2321 = C3221
  FourTensor[1][2][1][2] = VoigtMatrix(4, 4);
  FourTensor[2][1][1][2] = VoigtMatrix(4, 4);  // C2323 = C3223
  FourTensor[1][2][2][1] = VoigtMatrix(4, 4);
  FourTensor[2][1][2][1] = VoigtMatrix(4, 4);  // C2332 = C3232
  FourTensor[1][2][0][2] = VoigtMatrix(4, 5);
  FourTensor[2][1][0][2] = VoigtMatrix(4, 5);  // C2313 = C3213
  FourTensor[1][2][2][0] = VoigtMatrix(4, 5);
  FourTensor[2][1][2][0] = VoigtMatrix(4, 5);  // C2331 = C3231

  FourTensor[0][2][0][0] = VoigtMatrix(5, 0);
  FourTensor[2][0][0][0] = VoigtMatrix(5, 0);  // C1311 = C3111
  FourTensor[0][2][1][1] = VoigtMatrix(5, 1);
  FourTensor[2][0][1][1] = VoigtMatrix(5, 1);  // C1322 = C3122
  FourTensor[0][2][2][2] = VoigtMatrix(5, 2);
  FourTensor[2][0][2][2] = VoigtMatrix(5, 2);  // C1333 = C3133
  FourTensor[0][2][0][1] = VoigtMatrix(5, 3);
  FourTensor[2][0][0][1] = VoigtMatrix(5, 3);  // C1312 = C3112
  FourTensor[0][2][1][0] = VoigtMatrix(5, 3);
  FourTensor[2][0][1][0] = VoigtMatrix(5, 3);  // C1321 = C3121
  FourTensor[0][2][1][2] = VoigtMatrix(5, 4);
  FourTensor[2][0][1][2] = VoigtMatrix(5, 4);  // C1323 = C3123
  FourTensor[0][2][2][1] = VoigtMatrix(5, 4);
  FourTensor[2][0][2][1] = VoigtMatrix(5, 4);  // C1332 = C3132
  FourTensor[0][2][0][2] = VoigtMatrix(5, 5);
  FourTensor[2][0][0][2] = VoigtMatrix(5, 5);  // C1313 = C3113
  FourTensor[0][2][2][0] = VoigtMatrix(5, 5);
  FourTensor[2][0][2][0] = VoigtMatrix(5, 5);  // C1331 = C3131

}  // Setup4Tensor()


/*------------------------------------------------------------------------------------*
 |  Set every value of 4-Tensor(3x3x3x3) to zero                         rauch  07/14 |
 *------------------------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::Clear4Tensor(double FourTensor[3][3][3][3])
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++) FourTensor[i][j][k][l] = 0.0;

}  // Clear4Tensor()


/*------------------------------------------------------------------------------------*
 |  Multiply: Matrix(3x3) * 4-Tensor(3x3x3x3)                            rauch  07/14 |
 *------------------------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::MultMatrixFourTensor(double FourTensorResult[3][3][3][3],
    const LINALG::Matrix<3, 3>& Matrix,   // B^i_m
    const double FourTensor[3][3][3][3],  // A^mjkl
    const bool clearresulttensor)
{
  if (clearresulttensor) Clear4Tensor(FourTensorResult);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          for (int m = 0; m < 3; m++)
            FourTensorResult[i][j][k][l] +=
                Matrix(i, m) * FourTensor[m][j][k][l];  // C^ijkl = B^i_m * A^mjkl

}  // MultMatrixFourTensor()


/*------------------------------------------------------------------------------------*
 |  Multiply: 4-Tensor(3x3x3x3) * Matrix(3,3)                            rauch  07/14 |
 *------------------------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::MultFourTensorMatrix(double FourTensorResult[3][3][3][3],
    LINALG::Matrix<3, 3>& Matrix,   // B^m_l
    double FourTensor[3][3][3][3],  // A^ijkm
    bool clearresulttensor)
{
  if (clearresulttensor) Clear4Tensor(FourTensorResult);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          for (int m = 0; m < 3; m++)
            FourTensorResult[i][j][k][l] +=
                FourTensor[i][j][k][m] * Matrix(m, l);  // C^ijkl = A^ijkm B_ml

}  // MultMatrixFourTensor()


/*------------------------------------------------------------------------------------------*
 |  Setup 6x6 matrix in Voigt notation from 4-Tensor                           rauch  07/14 |
 *------------------------------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::Setup6x6VoigtMatrix(
    LINALG::Matrix<6, 6>& VoigtMatrix, double FourTensor[3][3][3][3])
{
  ///*  [      C1111                 C1122                C1133                0.5*(C1112+C1121)
  /// 0.5*(C1123+C1132)                0.5*(C1113+C1131)      ]
  //    [      C2211                 C2222                C2233                0.5*(C2212+C2221)
  //    0.5*(C2223+C2232)                0.5*(C2213+C2231)      ] [      C3311                 C3322
  //    C3333                0.5*(C3312+C3321)               0.5*(C3323+C3332) 0.5*(C3313+C3331) ]
  //    [0.5*(C1211+C2111)    0.5*(C1222+C2122)    0.5*(C1233+C2133) 0.5*(C1212+C2112+C1221+C2121)
  //    0.5*(C1223+C2123+C1232+C2132)    0.5*(C1213+C2113+C1231+C2131)] [0.5*(C2311+C3211)
  //    0.5*(C2322+C3222)    0.5*(C2333+C3233)    0.5*(C2312+C3212+C2321+C3221)
  //    0.5*(C2323+C3223+C2332+C3232)    0.5*(C2313+C3213+C2331+C3231)] [0.5*(C1322+C3122)
  //    0.5*(C1322+C3122)    0.5*(C1333+C3133)    0.5*(C1312+C3112+C1321+C3121)
  //    0.5*(C1323+C3123+C1332+C3132)    0.5*(C1313+C3113+C1331+C3131)] */

  // Setup 4-Tensor from 6x6 Voigt matrix
  VoigtMatrix(0, 0) = FourTensor[0][0][0][0];                                   // C1111
  VoigtMatrix(0, 1) = FourTensor[0][0][1][1];                                   // C1122
  VoigtMatrix(0, 2) = FourTensor[0][0][2][2];                                   // C1133
  VoigtMatrix(0, 3) = 0.5 * (FourTensor[0][0][0][1] + FourTensor[0][0][1][0]);  // 0.5*(C1112+C1121)
  VoigtMatrix(0, 4) = 0.5 * (FourTensor[0][0][1][2] + FourTensor[0][0][2][1]);  // 0.5*(C1123+C1132)
  VoigtMatrix(0, 5) = 0.5 * (FourTensor[0][0][0][2] + FourTensor[0][0][2][0]);  // 0.5*(C1113+C1131)

  VoigtMatrix(1, 0) = FourTensor[1][1][0][0];                                   // C2211
  VoigtMatrix(1, 1) = FourTensor[1][1][1][1];                                   // C2222
  VoigtMatrix(1, 2) = FourTensor[1][1][2][2];                                   // C2233
  VoigtMatrix(1, 3) = 0.5 * (FourTensor[1][1][0][1] + FourTensor[1][1][1][0]);  // 0.5*(C2212+C2221)
  VoigtMatrix(1, 4) = 0.5 * (FourTensor[1][1][1][2] + FourTensor[1][1][2][1]);  // 0.5*(C2223+C2232)
  VoigtMatrix(1, 5) = 0.5 * (FourTensor[1][1][0][2] + FourTensor[1][1][2][0]);  // 0.5*(C2213+C2231)

  VoigtMatrix(2, 0) = FourTensor[2][2][0][0];                                   // C3311
  VoigtMatrix(2, 1) = FourTensor[2][2][1][1];                                   // C3322
  VoigtMatrix(2, 2) = FourTensor[2][2][2][2];                                   // C3333
  VoigtMatrix(2, 3) = 0.5 * (FourTensor[2][2][0][1] + FourTensor[2][2][1][0]);  // 0.5*(C3312+C3321)
  VoigtMatrix(2, 4) = 0.5 * (FourTensor[2][2][1][2] + FourTensor[2][2][2][1]);  // 0.5*(C3323+C3332)
  VoigtMatrix(2, 5) = 0.5 * (FourTensor[2][2][0][2] + FourTensor[2][2][2][0]);  // 0.5*(C3313+C3331)

  VoigtMatrix(3, 0) = 0.5 * (FourTensor[0][1][0][0] + FourTensor[1][0][0][0]);  // 0.5*(C1211+C2111)
  VoigtMatrix(3, 1) = 0.5 * (FourTensor[0][1][1][1] + FourTensor[1][0][1][1]);  // 0.5*(C1222+C2122)
  VoigtMatrix(3, 2) = 0.5 * (FourTensor[0][1][2][2] + FourTensor[1][0][2][2]);  // 0.5*(C1233+C2133)
  VoigtMatrix(3, 3) =
      0.25 * (FourTensor[0][1][0][1] + FourTensor[1][0][0][1] + FourTensor[0][1][1][0] +
                 FourTensor[1][0][1][0]);  // 0.5*(C1212+C2112+C1221+C2121)
  VoigtMatrix(3, 4) =
      0.25 * (FourTensor[0][1][1][2] + FourTensor[1][0][1][2] + FourTensor[0][1][2][1] +
                 FourTensor[1][0][2][1]);  // 0.5*(C1223+C2123+C1232+C2132)
  VoigtMatrix(3, 5) =
      0.25 * (FourTensor[0][1][0][2] + FourTensor[1][0][0][2] + FourTensor[0][1][2][0] +
                 FourTensor[1][0][2][0]);  // 0.5*(C1213+C2113+C1231+C2131)

  VoigtMatrix(4, 0) = 0.5 * (FourTensor[1][2][0][0] + FourTensor[2][1][0][0]);  // 0.5*(C2311+C3211)
  VoigtMatrix(4, 1) = 0.5 * (FourTensor[1][2][1][1] + FourTensor[2][1][1][1]);  // 0.5*(C2322+C3222)
  VoigtMatrix(4, 2) = 0.5 * (FourTensor[1][2][2][2] + FourTensor[2][1][2][2]);  // 0.5*(C2333+C3233)
  VoigtMatrix(4, 3) =
      0.25 * (FourTensor[1][2][0][1] + FourTensor[2][1][0][1] + FourTensor[1][2][1][0] +
                 FourTensor[2][1][1][0]);  // 0.5*(C2312+C3212+C2321+C3221)
  VoigtMatrix(4, 4) =
      0.25 * (FourTensor[1][2][1][2] + FourTensor[2][1][1][2] + FourTensor[1][2][2][1] +
                 FourTensor[2][1][2][1]);  // 0.5*(C2323+C3223+C2332+C3232)
  VoigtMatrix(4, 5) =
      0.25 * (FourTensor[1][2][0][2] + FourTensor[2][1][0][2] + FourTensor[1][2][2][0] +
                 FourTensor[2][1][2][0]);  // 0.5*(C2313+C3213+C2331+C3231)

  VoigtMatrix(5, 0) = 0.5 * (FourTensor[0][2][0][0] + FourTensor[2][0][0][0]);  // 0.5*(C1311+C3111)
  VoigtMatrix(5, 1) = 0.5 * (FourTensor[0][2][1][1] + FourTensor[2][0][1][1]);  // 0.5*(C1322+C3122)
  VoigtMatrix(5, 2) = 0.5 * (FourTensor[0][2][2][2] + FourTensor[2][0][2][2]);  // 0.5*(C1333+C3133)
  VoigtMatrix(5, 3) =
      0.25 * (FourTensor[0][2][0][1] + FourTensor[2][0][0][1] + FourTensor[0][2][1][0] +
                 FourTensor[2][0][1][0]);  // 0.5*(C1312+C3112+C1321+C3121)
  VoigtMatrix(5, 4) =
      0.25 * (FourTensor[0][2][1][2] + FourTensor[2][0][1][2] + FourTensor[0][2][2][1] +
                 FourTensor[2][0][2][1]);  // 0.5*(C1323+C3123+C1332+C3132)
  VoigtMatrix(5, 5) =
      0.25 * (FourTensor[0][2][0][2] + FourTensor[2][0][0][2] + FourTensor[0][2][2][0] +
                 FourTensor[2][0][2][0]);  // 0.5*(C1313+C3113+C1331+C3131)

}  // Setup6x6VoigtMatrix()


/*---------------------------------------------------------------*
 |  matrix root                                     rauch  07/14 |
 *---------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::MatrixRoot3x3(LINALG::Matrix<3, 3>& MatrixInOut)
{
  double Norm = MatrixInOut.Norm2();
  // direct calculation for zero-matrix
  if (Norm == 0.)
  {
    MatrixInOut.Clear();
    return;
  }
  else
  {
    LINALG::Matrix<3, 3> EV(MatrixInOut);
    LINALG::Matrix<3, 3> EW;

    // MatixInOut = EV * EW * EV^{-1}
    LINALG::SYEV(EV, EW, EV);

    // sqrt(MatrixInOut) = EV * sqrt(EW) * EVT
    // loop over all eigenvalues
    for (int a = 0; a < 3; a++) EW(a, a) = sqrt(EW(a, a));

    MatrixInOut.Clear();

    // temp = sqrt(EW) * EVT
    LINALG::Matrix<3, 3> temp;
    temp.MultiplyNT(EW, EV);

    // sqrt(MatrixInOut) = EV * sqrt(EW) * EV^{-1} = EV * temp
    MatrixInOut.MultiplyNN(EV, temp);
  }

  return;

}  // MatrixRoot3x3()


/*-------------------------------------------------------------------------------------*
 |  matrix root derivative of a symmetric 3x3 matrix                      rauch  07/14 |
 *-------------------------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::MatrixRootDerivativeSym3x3(
    const LINALG::Matrix<3, 3>& MatrixIn, LINALG::Matrix<6, 6>& MatrixRootDeriv)
{
  double Norm = MatrixIn.Norm2();

  LINALG::Matrix<6, 6> id4sharp(true);  // souza S.31 eq. (2.110)???
  for (int i = 0; i < 3; i++) id4sharp(i, i) = 1.0;
  for (int i = 3; i < 6; i++) id4sharp(i, i) = 0.5;

  // direct calculation for zero-matrix
  if (Norm == 0.)
  {
    MatrixRootDeriv = id4sharp;
    dserror("d sqrt(C)/ d C not defined for C==0");
    return;
  }

  else
  {
    double EWtolerance = 1.e-12;  // see souza S.737 Remark A.2

    LINALG::Matrix<3, 3> EV(MatrixIn);
    LINALG::Matrix<3, 3> EW;
    LINALG::SYEV(EV, EW, EV);

    MatrixRootDeriv.Clear();
    // souza eq. (A.52)
    // note: EW stored in ascending order

    //  d X^2 / d X  =  1/2 * (  delta_jk X_lj + delta_il X_kj
    //                         + delta_jl X_ik + delta_kj X_il )    souza eq. (A.46)
    //
    // y_i = sqrt(x_i)
    // dy_i / dx_j = delta_ij 1/(2*sqrt(x_i))

    LINALG::Matrix<3, 3> id2(true);
    for (int i = 0; i < 3; i++) id2(i, i) = 1.0;
    //  // --------------------------------- switch by number of equal eigenvalues

    if (abs(EW(0, 0) - EW(1, 1)) < EWtolerance &&
        abs(EW(1, 1) - EW(2, 2)) < EWtolerance)  // ------------------ x_a == x_b == x_c
    {
      // calculate derivative
      MatrixRootDeriv = id4sharp;
      MatrixRootDeriv.Scale(1.0 / (2.0 * sqrt(EW(0, 0))));
    }

    else if ((abs(EW(0, 0) - EW(1, 1)) < EWtolerance && abs(EW(1, 1) - EW(2, 2)) > EWtolerance) ||
             (abs(EW(0, 0) - EW(1, 1)) > EWtolerance &&
                 abs(EW(1, 1) - EW(2, 2)) <
                     EWtolerance))  // ---- x_a != x_b == x_c or x_a == x_b != x_c
    {
      // scalar factors
      double s1 = 0.0;
      double s2 = 0.0;
      double s3 = 0.0;
      double s4 = 0.0;
      double s5 = 0.0;
      double s6 = 0.0;

      int a = 0;
      int c = 0;

      // switch which two EW are equal
      if (abs(EW(0, 0) - EW(1, 1)) < EWtolerance &&
          abs(EW(1, 1) - EW(2, 2)) > EWtolerance)  // ----------------------- x_a == x_b != x_c
      {
        a = 2;
        c = 0;
      }
      else if (abs(EW(0, 0) - EW(1, 1)) > EWtolerance &&
               abs(EW(1, 1) - EW(2, 2)) < EWtolerance)  // ------------------ x_a != x_b == x_c
      {
        a = 0;
        c = 2;
      }
      else
        dserror("you should not end up here");

      // in souza eq. (A.53):
      s1 = (sqrt(EW(a, a)) - sqrt(EW(c, c))) / (pow(EW(a, a) - EW(c, c), 2.0)) -
           (1.0 / (2.0 * sqrt(EW(c, c)))) / (EW(a, a) - EW(c, c));
      s2 = 2.0 * EW(c, c) * (sqrt(EW(a, a)) - sqrt(EW(c, c))) / (pow(EW(a, a) - EW(c, c), 2.0)) -
           (EW(a, a) + EW(c, c)) / (EW(a, a) - EW(c, c)) * (1.0 / (2.0 * sqrt(EW(c, c))));
      s3 = 2.0 * (sqrt(EW(a, a)) - sqrt(EW(c, c))) / (pow(EW(a, a) - EW(c, c), 3.0)) -
           ((1.0 / (2.0 * sqrt(EW(a, a)))) + (1.0 / (2.0 * sqrt(EW(c, c))))) /
               (pow(EW(a, a) - EW(c, c), 2.0));
      s4 = EW(c, c) * s3;
      s5 = s4;
      s6 = EW(c, c) * EW(c, c) * s3;

      // calculate derivative
      // + s_1 (d X^2 / d X)
      MAT::AddToCmatDerivTensorSquare(MatrixRootDeriv, s1, MatrixIn, 1.);
      // - s_2 I_s
      MatrixRootDeriv.Update(-s2, id4sharp, 1.);
      // - s_3 (X \dyad X)
      MAT::ElastSymTensorMultiply(MatrixRootDeriv, -1. * s3, MatrixIn, MatrixIn, 1.);
      // + s_4 (X \dyad I)
      MAT::ElastSymTensorMultiply(MatrixRootDeriv, s4, MatrixIn, id2, 1.);
      // + s_5 (I \dyad X)
      MAT::ElastSymTensorMultiply(MatrixRootDeriv, s5, id2, MatrixIn, 1.);
      // - s_6 (I \dyad I)
      MAT::ElastSymTensorMultiply(MatrixRootDeriv, -1. * s6, id2, id2, 1.);
    }

    else if (abs(EW(0, 0) - EW(1, 1)) > EWtolerance &&
             abs(EW(1, 1) - EW(2, 2)) > EWtolerance)  // ----------------- x_a != x_b != x_c
    {
      for (int a = 0; a < 3; a++)  // loop over all eigenvalues
      {
        // a=0 || a=1 || a=2
        int b = (a + 1) % 3;  // b=1 || b=2 || b=0     even (cyclic) permutations of (a,b,c)
        int c = (a + 2) % 3;  // c=2 || c=0 || c=1

        LINALG::Matrix<3, 1> ea;
        LINALG::Matrix<3, 1> eb;
        LINALG::Matrix<3, 1> ec;
        for (int i = 0; i < 3; i++)
        {
          ea(i) = EV(i, a);
          eb(i) = EV(i, b);
          ec(i) = EV(i, c);
        }
        LINALG::Matrix<3, 3> Ea;
        Ea.MultiplyNT(ea, ea);  // souza S.26 eq. (2.63)
        LINALG::Matrix<3, 3> Eb;
        Eb.MultiplyNT(eb, eb);
        LINALG::Matrix<3, 3> Ec;
        Ec.MultiplyNT(ec, ec);

        double fac = sqrt(EW(a, a)) / ((EW(a, a) - EW(b, b)) * (EW(a, a) - EW(c, c)));

        // calculate derivative
        // + d X^2 / d X
        MAT::AddToCmatDerivTensorSquare(MatrixRootDeriv, fac, MatrixIn, 1.);
        // - (x_b + x_c) I_s
        MatrixRootDeriv.Update(-1. * (EW(b, b) + EW(c, c)) * fac, id4sharp, 1.);
        // - [(x_a - x_b) + (x_a - x_c)] (E_a \dyad E_a)
        MAT::ElastSymTensorMultiply(MatrixRootDeriv,
            -1. * fac * ((EW(a, a) - EW(b, b)) + (EW(a, a) - EW(c, c))), Ea, Ea, 1.);
        // - (x_b - x_c) (E_b \dyad E_b)
        MAT::ElastSymTensorMultiply(MatrixRootDeriv, -1. * fac * (EW(b, b) - EW(c, c)), Eb, Eb, 1.);
        // + (x_b - x_c) (E_c \dyad E_c)
        MAT::ElastSymTensorMultiply(MatrixRootDeriv, fac * (EW(b, b) - EW(c, c)), Ec, Ec, 1.);
        // dy / dx_a (E_a \dyad E_a)
        MAT::ElastSymTensorMultiply(MatrixRootDeriv, 1.0 / (2.0 * sqrt(EW(a, a))), Ea, Ea, 1.);
      }  // end loop over all eigenvalues
    }

    else
      dserror("you should not end up here.");
  }

  return;
}  // MatrixRootDerivativeSym3x3()


/*------------------------------------------------------------------------------------------*
 |  Get the real derivative of the inverse of the square root.                 rauch  07/14 |
 *------------------------------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::MatrixInverseOfRootDerivative(
    double MatrixRootDerivative[3][3][3][3], const LINALG::Matrix<3, 3>& MatrixRootInverse,
    double result[3][3][3][3])
{
  double temp[3][3][3][3] = {{{{0.}}}};

  // -sqrt(C)^{-1}_ij (d sqrt(C)/dC)_klmn sqrt(C)^{-1}_qp =
  // B_ilmp = -sqrt(C)^{-1}_ij (d sqrt(C)/dC)_jlmq sqrt(C)^{-1}_qp
  for (int i = 0; i < 3; ++i)
    for (int l = 0; l < 3; ++l)
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
          for (int j = 0; j < 3; ++j)
            temp[i][l][m][n] += -MatrixRootInverse(i, j) * MatrixRootDerivative[j][l][m][n];

  for (int i = 0; i < 3; ++i)
    for (int l = 0; l < 3; ++l)
      for (int m = 0; m < 3; ++m)
        for (int p = 0; p < 3; ++p)
          for (int q = 0; q < 3; ++q)
            result[i][l][m][p] += temp[i][l][m][q] * MatrixRootInverse(q, p);
  return;
}


/*------------------------------------------------------------------------------------------*
 |  Print Four Tensor                                                          rauch  07/14 |
 *------------------------------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellPassiveFiber::PrintFourTensor(double FourTensor[3][3][3][3])
{
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          std::cout << "ELEMENT " << i << j << k << l << " : " << FourTensor[i][j][k][l]
                    << std::endl;
  return;
}
