/*----------------------------------------------------------------------*/
/*!
\file elast_remodelfiber.cpp
\brief


the input line should read
  MAT 1 ELAST_RemodelFiber NUMMAT 1 MATIDS 100 TDECAY 1.0 SIGMAPRE 1.0

<pre>
Maintainer: Fabian Bräu
            braeu@lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_remodelfiber.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::RemodelFiber::RemodelFiber(
    Teuchos::RCP<MAT::PAR::Material> matdata
)
: Parameter(matdata),
  nummat_(matdata->GetInt("NUMMAT")),
  matids_(matdata->Get<std::vector<int> >("MATIDS")),
  tdecay_(matdata->GetDouble("TDECAY")),
  sigmapre_(matdata->GetDouble("SIGMAPRE")),
  k_growth_(matdata->GetDouble("GROWTHFAC")),
  cur_w_collagen_(matdata->GetMutable<std::vector<double> >("COLMASSFRAC"))
{
  // check if sizes fit
  if (nummat_ != (int)matids_->size())
    dserror("number of materials %d does not fit to size of material vector %d", nummat_, matids_->size());

  // check decay time validity
  if (tdecay_<=0.)
    dserror("decay time must be positive");
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   fb         09/15 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::RemodelFiber::RemodelFiber(MAT::ELASTIC::PAR::RemodelFiber* params)
  : params_(params),
    potsumfiber_(0)
{
  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;
  for (m=params_->matids_->begin(); m!=params_->matids_->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to allocate");
    potsumfiber_.push_back(sum);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::PackSummand(DRT::PackBuffer& data) const
{
  int num_fiber = 0;
  num_fiber = last_ilambda_r_.size();

  AddtoPack(data,num_fiber);

  for(int i=0;i<num_fiber;++i)
  {
    AddtoPack(data,last_ilambda_r_[i]);
    AddtoPack(data,current_ilambda_r_[i]);
    AddtoPack(data,current_w_collagen_[i]);
    AddtoPack(data,stress_[i]);
  }

  if (params_ != NULL) // summands are not accessible in postprocessing mode
    // loop map of associated potential summands
    for (unsigned int p=0; p<potsumfiber_.size(); ++p)
     potsumfiber_[p]->PackSummand(data);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::UnpackSummand(const std::vector<char>& data,
                                               std::vector<char>::size_type& position)
{
  int num_fiber=0;
  ExtractfromPack(position,data,num_fiber);

  last_ilambda_r_.resize(num_fiber);
  current_ilambda_r_.resize(num_fiber);
  current_w_collagen_.resize(num_fiber);
  stress_.resize(num_fiber);

  for(int i=0;i<num_fiber;++i)
  {
    ExtractfromPack(position,data,last_ilambda_r_[i]);
    ExtractfromPack(position,data,current_ilambda_r_[i]);
    ExtractfromPack(position,data,current_w_collagen_[i]);
    ExtractfromPack(position,data,stress_[i]);
  }

  // loop map of associated potential summands
  for (unsigned int p=0; p<potsumfiber_.size(); ++p)
    potsumfiber_[p]->UnpackSummand(data,position);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::Setup(int numgp,DRT::INPUT::LineDefinition* linedef)
{
  // setup fiber and inelastic history variable
  last_ilambda_r_.resize(potsumfiber_.size());
  current_ilambda_r_.resize(potsumfiber_.size());
  current_w_collagen_.resize(potsumfiber_.size());
  stress_.resize(potsumfiber_.size());

  for(unsigned p=0;p<potsumfiber_.size();++p)
  {
    last_ilambda_r_[p].resize(numgp,1.0);
    current_ilambda_r_[p].resize(numgp,1.0);
    current_w_collagen_[p].resize(numgp,1.0);
    stress_[p].resize(numgp,1.0);

    potsumfiber_[p]->Setup(linedef);
  }


  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::Update()
{
  // update history variable
  for(unsigned p=0;p<potsumfiber_.size();++p)
    for(unsigned gp=0;gp<current_ilambda_r_[p].size();++gp)
      last_ilambda_r_[p][gp] = current_ilambda_r_[p][gp];

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::AddStressCmatRemodel(
  const LINALG::Matrix<3,3>* defgrd,
  Teuchos::ParameterList& params,
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress,
  const int eleGID
  )
{
  // Clear variables
  stress.Clear();
  cmat.Clear();

  // current Gauß-Point
  int gp = params.get<int>("gp");

  // time step
  double dt = params.get<double>("delta time");

  // Calculate the inelastic stretch in fiberdirection
  // by using Eq. A13 in "A homogenized constrained mixture (and mechanical analog) model
  // for grwoth and remodeling of soft tissue
  //
  // Implicit time integration scheme is used
  //
  // REMARK: 1D formulation is used for fiber material

  // fiberdirection
  std::vector<LINALG::Matrix<3,1> > fibervecs;

  // Structural Tensor in "stress-like" Voigt notation
  LINALG::Matrix<6,1> A(true);

  // right Cauchy-Green
  LINALG::Matrix<3,3> C(true);

  // some kinematic quantities
  // stretch in fiberdirection
  double lambda = 0.0;

  // elastic stretch in fiberdirection
  double lambda_e = 0.0;

  // last inelastic stretch in fiberdirection
  double lambda_r_last = 1.0;

  // current inverse inelastic stretch in fiberdirection
  double ilambda_r_cur = 1.0;

  // current mass fraction of collagen
  double cur_w_col = 0.0;

  // principal invariants (I4,I5)
  LINALG::Matrix<2,1> prinv(true);

  // first derivatives of strain energy function w.r.t. I4
  LINALG::Matrix<2,1> dPI(true);

  // second derivatives of strain energy function w.r.t. I4
  LINALG::Matrix<3,1> ddPII(true);

  // third derivatives of strain energy function w.r.t. I4
  LINALG::Matrix<4,1> dddPIII(true);

  // Derivative of the evolution eq. for the inelastic stretch w.r.t. the inelastic stretch in fiberdirection
  double dRdlambda_r = 0.0;

  // Derivative of the inelastic stretch rate w.r.t. the quadratic stretch in fiberdirection
  LINALG::Matrix<1,6> dlambda_rdC(true);

  // Derivative of the 2nd Piola Kirchhoff stress w.r.t. the inelastic stretch in fiberdirection
  double dSdlambda_r = 0.0;

  // Residual
  double R = 0.0;

  // Cauchy stress
  double sig = 0.0;

  // determinant of deformation gradient
  double J = defgrd->Determinant();

  // lambda_r increment
  double delta_lambda_r = 0.0;

  // temporary variables
  LINALG::Matrix<3,1> tmp_vec(true);
  LINALG::Matrix<1,1> tmp_scal(true);

  for(unsigned p=0;p<potsumfiber_.size();++p)
  {
    // last inelastic stretch in fiberdirection
    lambda_r_last = 1.0/last_ilambda_r_[p][gp];

    // number of loops
    int nr_loop = 0;

    // start newton with initialized parameter
    ilambda_r_cur = current_ilambda_r_[p][gp];

    // get fiberdirection
    potsumfiber_[p]->GetFiberVecs(fibervecs);

    // stretch in fiberdirection
    C.MultiplyTN(1.0,*defgrd,*defgrd,0.0);
    tmp_vec.MultiplyNN(1.0,C,fibervecs[p],0.0);
    tmp_scal.MultiplyTN(1.0,fibervecs[p],tmp_vec,0.0);
    lambda = sqrt(tmp_scal(0));

    // initialize mass fraction collagen
    cur_w_col = current_w_collagen_[p][gp];

    do
    {
      if(nr_loop != 0)
      {
        // lambda_r increment
        delta_lambda_r = -R/dRdlambda_r;

        // Update lambda_r_cur
        ilambda_r_cur = 1.0/(delta_lambda_r + (1.0/ilambda_r_cur));
      }
      // Clear some variables
      prinv.Clear();
      dPI.Clear();
      ddPII.Clear();
      dddPIII.Clear();

      // elastic stretch in fiberdirection
      lambda_e = lambda * ilambda_r_cur;

      // principal invariant (I4)
      prinv(0) = lambda_e * lambda_e;

      // get first,second and third derivative of starin energy function
      potsumfiber_[p]->GetDerivativesAniso(dPI,ddPII,dddPIII,prinv,eleGID);

      // Cauchy stress
      sig = (2.0/J)*prinv(0)*cur_w_col*dPI(0);


      /*-----------------d((sig-sig_pre)/T)/dlambda_r-----------------*/
      dRdlambda_r = -(4.0/J)*cur_w_col*ilambda_r_cur*prinv(0)*(dPI(0)+prinv(0)*ddPII(0))*(1.0/params_->tdecay_);

      /*------------------(-2.0)*(dsig/dlambda_e^2)/dlambda_r*(lambda_e^2*(1-last_lambda_r*ilambda_r)/dt)-----------------*/
      dRdlambda_r += (8.0/J)*cur_w_col*ilambda_r_cur*prinv(0)*(2.0*ddPII(0)+prinv(0)*dddPIII(0))*(prinv(0)/dt)*(1.0-lambda_r_last*ilambda_r_cur);

      /*-----------------(-2.0)*(dsig/dlambda_e^2)*d(lambda_e^2*(1-last_lambda_r*ilambda_r)/dt)/dlambda_r-----------------*/
      dRdlambda_r += -(4.0/J)*cur_w_col*(dPI(0)+prinv(0)*ddPII(0))*ilambda_r_cur*prinv(0)*(1.0/dt)*(3.0*lambda_r_last*ilambda_r_cur-2.0);

      /*-----------------Residual-----------------*/
      R = (sig-params_->sigmapre_)/params_->tdecay_ - (4.0/J)*cur_w_col*(dPI(0)+prinv(0)*ddPII(0))*prinv(0)*(1.0/dt)*(1.0-lambda_r_last*ilambda_r_cur);

      nr_loop++;
      if(nr_loop > 50)
        dserror("no convergence!!!");
    }
    while(fabs(R) > 1.0e-8);

    // update inelastic stretch in fiber direction
    current_ilambda_r_[p][gp] = ilambda_r_cur;

    // update current Cauchy stress
    stress_[p][gp] = sig;

    /*-----------------------Update 2nd Piola-Kirchhoff Stress and Elasticity tensor-----------------------*/
    // Update stress and (one part!) of elasticity tensor
    // structural tensor (stress-like Voigt notation)
    for (int i = 0; i < 3; ++i)
      A(i) = fibervecs[p](i)*fibervecs[p](i);
    A(3) = fibervecs[p](0)*fibervecs[p](1);
    A(4) = fibervecs[p](1)*fibervecs[p](2);
    A(5) = fibervecs[p](0)*fibervecs[p](2);


    stress.Update(2.0*cur_w_col*ilambda_r_cur*ilambda_r_cur*dPI(0),A,1.0);

    cmat.MultiplyNT(4.0*cur_w_col*ilambda_r_cur*ilambda_r_cur*ilambda_r_cur*ilambda_r_cur*ddPII(0),A,A,1.0);


    /*-----------------------Update (extra part!) of Elasticity tensor-----------------------*/
    // derivative of the residual R w.r.t. C
    LINALG::Matrix<1,6> dRdC(true);

    //d((sig-sig_pre)/T)/dC
    LINALG::Matrix<6,1> iC_voigt(true);
    LINALG::Matrix<3,3> tmp3x3(true);

    tmp3x3.Invert(C);
    for (int i = 0; i < 3; i++)
      iC_voigt(i) = tmp3x3(i,i);
    iC_voigt(3) = 0.5*(tmp3x3(0,1)+tmp3x3(1,0));
    iC_voigt(4) = 0.5*(tmp3x3(1,2)+tmp3x3(2,1));
    iC_voigt(5) = 0.5*(tmp3x3(0,2)+tmp3x3(2,0));



    // update dRdC
    // d((sig-sig_pre)/T)dC
    dRdC.UpdateT(-cur_w_col*(1.0/J)*prinv(0)*dPI(0)*(1.0/params_->tdecay_),iC_voigt,0.0);
    dRdC.UpdateT(2.0*cur_w_col*(1.0/J)*(1.0/params_->tdecay_)*ilambda_r_cur*ilambda_r_cur*(dPI(0)+prinv(0)*ddPII(0)),A,1.0);


    // (-2.0)*d(dsig/dlambda_e^2)/dC*(lambda_e^2*(1-last_lambda_r*ilambda_r)/dt)
    dRdC.UpdateT(2.0*cur_w_col*(1.0/J)*(dPI(0)+prinv(0)*ddPII(0))*prinv(0)*(1.0/dt)*(1.0-lambda_r_last*ilambda_r_cur),iC_voigt,1.0);
    dRdC.UpdateT(-4.0*cur_w_col*(1.0/J)*ilambda_r_cur*ilambda_r_cur*(2.0*ddPII(0)+prinv(0)*dddPIII(0))*(1.0/dt)*(1.0-lambda_r_last*ilambda_r_cur)*ilambda_r_cur*ilambda_r_cur,A,1.0);


    // (-2.0)*(dsig/dlambda_e^2)*d(lambda_e^2*(1-last_lambda_r*ilambda_r)/dt)/dC
    dRdC.UpdateT(-(4.0/J)*cur_w_col*(dPI(0)+prinv(0)*ddPII(0))*(1.0/dt)*(1.0-lambda_r_last*ilambda_r_cur)*ilambda_r_cur*ilambda_r_cur,A,1.0);



    // evaluate derivative of inelastic stretch w.r.t. right Cauchy Green tensor
    dlambda_rdC.Update(-1.0/dRdlambda_r,dRdC,0.0);

    // evaluate derivative of 2nd Piola Kirchhoff stress w.r.t. inelastic stretch
    dSdlambda_r = -4.0*cur_w_col*ilambda_r_cur*ilambda_r_cur*ilambda_r_cur*(dPI(0)+prinv(0)*ddPII(0));


    // Update elasticity tensor
    cmat.MultiplyNN(2.0*dSdlambda_r,A,dlambda_rdC,1.0);
  }

  return;
}



/*---------------------------------------------------------------------*
 | return names of visualization data (public)                         |
 *---------------------------------------------------------------------*/
void MAT::ELASTIC::RemodelFiber::VisNames(std::map<std::string,int>& names)
{
  std::string inelastic_defgrd = "inelastic_defgrd_fiber";
  std::string result_inelastic_defgrad;

  for (unsigned int p=0; p<potsumfiber_.size(); ++p)
  {
    std::stringstream sstm;
    sstm << inelastic_defgrd <<"_" << p;
    result_inelastic_defgrad = sstm.str();

    names[result_inelastic_defgrad] = 1;
  }


  std::string fiber_cauchy_stress = "fiber_cauchy_stress";
  std::string result_fiber_cauchy_stress;

  for (unsigned int p=0; p<potsumfiber_.size(); ++p)
  {
    std::stringstream sstm;
    sstm << fiber_cauchy_stress <<"_" << p;
    result_fiber_cauchy_stress = sstm.str();

    names[result_fiber_cauchy_stress] = 1;
  }


  std::string fiberdirection = "fiberdirection";
  std::string result_fiber;

  for (unsigned int p=0; p<potsumfiber_.size(); ++p)
  {
    std::stringstream sstm;
    sstm << fiberdirection <<"_" << p;
    result_fiber = sstm.str();

    names[result_fiber] = 3;
  }


  std::string mass_fraction_col = "mass_fraction_col";
  std::string result_mass_fraction_col;

  for (unsigned int p=0; p<potsumfiber_.size(); ++p)
  {
    std::stringstream sstm;
    sstm << mass_fraction_col <<"_" << p;
    result_mass_fraction_col = sstm.str();

    names[result_mass_fraction_col] = 1;
  }
}  // VisNames()


/*---------------------------------------------------------------------*
 | return visualization data (public)                                  |
 *---------------------------------------------------------------------*/
bool MAT::ELASTIC::RemodelFiber::VisData(
  const std::string& name,
  std::vector<double>& data,
  int numgp,
  int eleID
)
{
  if (name == "inelastic_defgrd_fiber_0")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<last_ilambda_r_[0].size(); gp++)
    {
      data[0] += 1.0/last_ilambda_r_[0][gp];
    }
    data[0] = data[0]/last_ilambda_r_[0].size();

    return true;
  }
  if (name == "inelastic_defgrd_fiber_1")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<last_ilambda_r_[1].size(); gp++)
    {
      data[0] += 1.0/last_ilambda_r_[1][gp];
    }
    data[0] = data[0]/last_ilambda_r_[1].size();

    return true;
  }
  if (name == "inelastic_defgrd_fiber_2")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<last_ilambda_r_[2].size(); gp++)
    {
      data[0] += 1.0/last_ilambda_r_[2][gp];
    }
    data[0] = data[0]/last_ilambda_r_[2].size();

    return true;
  }



  if (name == "fiber_cauchy_stress_0")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<stress_[0].size(); gp++)
    {
      data[0] += stress_[0][gp];
    }
    data[0] = data[0]/stress_[0].size();

    return true;
  }
  if (name == "fiber_cauchy_stress_1")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<stress_[1].size(); gp++)
    {
      data[0] += stress_[1][gp];
    }
    data[0] = data[0]/stress_[1].size();

    return true;
  }
  if (name == "fiber_cauchy_stress_2")
  {
    if (data.size()!= 1) dserror("size mismatch");
    for(unsigned gp=0; gp<stress_[2].size(); gp++)
    {
      data[0] += stress_[2][gp];
    }
    data[0] = data[0]/stress_[2].size();

    return true;
  }



std::vector<LINALG::Matrix<3,1> > fiberdirection;
if(name == "fiberdirection_0")
{
  if (data.size()!= 3) dserror("size mismatch");

  potsumfiber_[0]->GetFiberVecs(fiberdirection);

  for(unsigned i=0;i<data.size();i++)
    data[i] = fiberdirection[0](i);

  return true;
}
if(name == "fiberdirection_1")
{
  if (data.size()!= 3) dserror("size mismatch");

  potsumfiber_[1]->GetFiberVecs(fiberdirection);

  for(unsigned i=0;i<data.size();i++)
    data[i] = fiberdirection[0](i);

  return true;
}
if(name == "fiberdirection_2")
{
  if (data.size()!= 3) dserror("size mismatch");

  potsumfiber_[2]->GetFiberVecs(fiberdirection);

  for(unsigned i=0;i<data.size();i++)
    data[i] = fiberdirection[0](i);

  return true;
}


if(name == "mass_fraction_col_0")
{
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<current_w_collagen_[0].size();++i)
  {
    data[0] += current_w_collagen_[0][i];
  }
  data[0] = data[0]/current_w_collagen_[0].size();


  return true;
}
if(name == "mass_fraction_col_1")
{
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<current_w_collagen_[0].size();++i)
  {
    data[0] += current_w_collagen_[1][i];
  }
  data[0] = data[0]/current_w_collagen_[1].size();


  return true;
}
if(name == "mass_fraction_col_2")
{
  if (data.size()!= 1) dserror("size mismatch");
  for(unsigned i=0;i<current_w_collagen_[0].size();++i)
  {
    data[0] += current_w_collagen_[2][i];
  }
  data[0] = data[0]/current_w_collagen_[2].size();


  return true;
}

dserror("The output is only implemented for three different fiber directions!!!");
return false;
}  // VisData()



