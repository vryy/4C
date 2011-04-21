/*----------------------------------------------------------------------*/
/*!
\file plasticlinelast.cpp
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material for a 3D hex element
       following perfectly von Mises plasticity and a linear elastic material law
       (St.Venant Kirchhoff).

       perfect plasticity:
        - no hardening allowed
        - independent yield stress level of degree of plastification
        - constant uniaxial yield stress \f \sigma_y \,=\, const.\f

       geometric linear, for small strains

       example input line:
       MAT 1 MAT_Struct_PlasticLinElast YOUNG 206.9 NUE 0.29 DENS 0.0
         YIELD 0.45

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*
 | Definitions                                                    04/11 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 | Headers                                                        04/11 |
 *----------------------------------------------------------------------*/
#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "plasticlinelast.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*
 | constructor (public)                                           04/11 |
 *----------------------------------------------------------------------*/
MAT::PAR::PlasticLinElast::PlasticLinElast(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  poissonratio_(matdata->GetDouble("NUE")),
  density_(matdata->GetDouble("DENS")),
  yield_(matdata->GetDouble("YIELD"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::PlasticLinElast::CreateMaterial()
{
  return Teuchos::rcp(new MAT::PlasticLinElast(this));
}


MAT::PlasticLinElastType MAT::PlasticLinElastType::instance_;


DRT::ParObject* MAT::PlasticLinElastType::Create( const std::vector<char> & data )
{
  MAT::PlasticLinElast* plastic = new MAT::PlasticLinElast();
  plastic->Unpack(data);
  return plastic;
}


/*----------------------------------------------------------------------*
 | constructor (public)                                          04/11 |
 *----------------------------------------------------------------------*/
MAT::PlasticLinElast::PlasticLinElast()
  : params_(NULL)
{
  isinit_=false;

  /// plastic history deformation
  // old strain \f \varepsilon_n \f
  Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > strain_n;
  // old plastic strain \f \varepsilon^p_n \f
  Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > strainp_n;
  // old elastic strain \f \varepsilon^e_n\f
  Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > straine_n;
}


/*----------------------------------------------------------------------*
 | copy-constructor (public)                                      04/11 |
 *----------------------------------------------------------------------*/
MAT::PlasticLinElast::PlasticLinElast(MAT::PAR::PlasticLinElast* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 | pack (public)                                                 04/11 |
 *----------------------------------------------------------------------*/
void MAT::PlasticLinElast::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  // in case we are in post-process mode
  if (params_ != NULL) matid = params_->Id();
  AddtoPack(data,matid);

  // pack history data
  int histsize;
  if (!Initialized())
  {
    histsize=0;
  }
  else
  {
    // TODO 21.04 check
    histsize = strain_n->size();
  }
  // TODO 18.04.11 why to use "2*"??
  AddtoPack(data,2*histsize); // Length of history vector(s)
  for (int var = 0; var < histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data,strain_n->at(var));
    AddtoPack(data,strainp_n->at(var));
    AddtoPack(data,straine_n->at(var));
  }

  return;
} // Pack()


/*----------------------------------------------------------------------*
 | unpack (public)                                               04/11 |
 *----------------------------------------------------------------------*/
void MAT::PlasticLinElast::Unpack(const vector<char>& data)
{
  isinit_=true;
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::PlasticLinElast*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // history data
  int twicehistsize;
  ExtractfromPack(position,data,twicehistsize);

  if (twicehistsize == 0) isinit_=false;

  // unpack strain vectors
  strain_n  = rcp( new vector<LINALG::Matrix<NUM_STRESS_3D,1> > );
  strainp_n = rcp( new vector<LINALG::Matrix<NUM_STRESS_3D,1> > );
  straine_n = rcp( new vector<LINALG::Matrix<NUM_STRESS_3D,1> > );
  strain_np  = rcp( new vector<LINALG::Matrix<NUM_STRESS_3D,1> > );
  strainp_np = rcp( new vector<LINALG::Matrix<NUM_STRESS_3D,1> > );
  straine_np = rcp( new vector<LINALG::Matrix<NUM_STRESS_3D,1> > );

  for (int var=0; var<twicehistsize; var+=2)
  {
    LINALG::Matrix<NUM_STRESS_3D,1> tmp1(true);
    // current vectors have to be initialized
    strain_np->push_back(tmp1);
    strainp_np->push_back(tmp1);
    straine_np->push_back(tmp1);

    // last vectors are unpacked
    ExtractfromPack(position,data,tmp1);
    strain_n->push_back(tmp1);
    ExtractfromPack(position,data,tmp1);
    strainp_n->push_back(tmp1);
    ExtractfromPack(position,data,tmp1);
    straine_n->push_back(tmp1);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;

} // Unpack()


/*---------------------------------------------------------------------*
 | initialise / allocate internal stress variables (public)      04/11 |
 *---------------------------------------------------------------------*/
void MAT::PlasticLinElast::Setup(const int numgp)
{
  // TODO 21.04.11 check
  // initialize hist variables
  strain_n   = rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  strainp_n  = rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  straine_n  = rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  strain_np  = rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  strainp_np = rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  straine_np = rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);


  LINALG::Matrix<NUM_STRESS_3D,1> emptymat(true);
  strain_n->resize(numgp);
  strainp_n->resize(numgp);
  straine_n->resize(numgp);
  strain_np->resize(numgp);
  strainp_np->resize(numgp);
  straine_np->resize(numgp);
  // 21.04.11 check
  for (int i=0; i<numgp; i++)
  {
    strain_n->at(i)  = emptymat;
    strainp_n->at(i) = emptymat;
    straine_n->at(i) = emptymat;
    strain_np->at(i)  = emptymat;
    strainp_np->at(i) = emptymat;
    straine_np->at(i) = emptymat;
  }

  isinit_=true;
  return;

}  // Setup()


/*---------------------------------------------------------------------*
 | update internal stress variables (public)                     04/11 |
 *---------------------------------------------------------------------*/
void MAT::PlasticLinElast::Update()
{
  // make current values at time step t_n+1 to values of last step t_n
  strain_n  = strain_np;
  strainp_n = strainp_np;
  straine_n = straine_np;

  return;
}  // Update()


/*---------------------------------------------------------------------*
 | reset internal stress variables (public)                      04/11 |
 *---------------------------------------------------------------------*/
void MAT::PlasticLinElast::Reset()
{
  // do nothing,
  // because #histplasticrcgcurr_ and #histeplasticscurr_ are recomputed
  // anyway at every iteration based upon #histplasticrcglast_ and
  // #histeplasticslast_ untouched within time step

  return;
}  // Reset()


/*----------------------------------------------------------------------*
 | evaluate material (public)                                     04/11 |
 *----------------------------------------------------------------------*/
void MAT::PlasticLinElast::Evaluate(
  const LINALG::Matrix<6,1>& strain,  //!< linear strain vector
  const int gp, //!< current Gauss point
  Teuchos::ParameterList& params,  //!< parameter list for communication & HISTORY
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat, //!< material stiffness matrix
  LINALG::Matrix<NUM_STRESS_3D,1>& stress //!< 2nd PK-stress
  )
{
  // get material parameters
  // Young's modulus
  double young = params_->youngs_;
  // Poisson's ratio
  double nu = params_->poissonratio_;
  // yield stress
  double sigma_y = params_->yield_;

  // initialize scalars
  // lame constant
  double lambda = 0.0;
  lambda = nu * young / ( (1.0 + nu) * (1.0 - 2.0*nu) );
  // lame constant
  // shear modulus parameter mu == G
  double G = 0.0;
  G = young / ( 2.0*(1.0 + nu) );
  // bulk modulus K:= kappa = lambda + 2 * mu
  double K = 0.0;
  K = lambda + ( 2.0 / 3.0 * G );

  // yield function
  double Phi = 0.0;
  // incremental plastic multiplier Delta gamma
  double Dgamma = 0.0;

  // strain (in): independent variable
  //  strain_p: evolution is determined by the flow rule
  //  strain_e: definition of additive decomposition: strain_e = strain - strain_p

  // calculate strain increment
  // deltastrain = strain_n+1 - strain_n
  // copy total strain on deltastrain
  LINALG::Matrix<NUM_STRESS_3D,1> deltastrain(strain);
  // deltastrain = deltstrain - epsilon_n
  LINALG::Matrix<NUM_STRESS_3D,1> epsilon(true);
  for (int i=0; i<6; i++)
    epsilon(i,0) = strain_n->at(gp)(i,0);
  deltastrain.Update( (-1.0), epsilon, 1.0 );

  //-------------------------------------------------------------------
  // elastic predictor (trial values)
  //-------------------------------------------------------------------

  // ---------------------------------------------- elastic trial strain

  // epsilon^{e,trial}_{n+1} = epsilon^e_n + Delta epsilon
  LINALG::Matrix<NUM_STRESS_3D,1> epsilon_e(true);
  for (int i=0; i<6; i++)
    epsilon_e(i,0) = straine_n->at(gp)(i,0);
  epsilon_e.Update( 1.0, deltastrain, 1.0 );

  // split the elastic strain into deviatoric and volumetric part

  // volumetrically elastic strain
  // epsilon_v^{e,trial}_{n+1} = 1/3 . (epsilon^{e,trial}_{n+1} : I) . I
  //                           = tr( epsilon^{e,trial}_{n+1} ) . I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i=0; i<3; i++) Id(i) = 1.0;
  LINALG::Matrix<NUM_STRESS_3D,1> epsilon_ev(true);
  for (int i=0; i<6; i++)
    epsilon_ev(i,0) = epsilon_e(i,0) * Id(i,0);
  epsilon_ev.Scale(1.0/3.0);

  // deviatorically elastic strain
  // epsilon_d^{e,trial}_{n+1}
  //   = epsilon^{e,trial}_{n+1} - epsilon_v^{e,trial}_{n+1}
  LINALG::Matrix<NUM_STRESS_3D,1> epsilon_ed(true);
  for (int i=0; i<6; i++)
    epsilon_ed(i,0) = epsilon_e(i,0) - epsilon_ev(i,0);

  // plastic trial strain
  // epsilon^{p,trial}_{n+1} = epsilon^p_n
  // equivalent plastic strain
  LINALG::Matrix<NUM_STRESS_3D,1> epsilon_p(true);
  for (int i=0; i<6; i++)
    epsilon_p(i,0) = strainp_n->at(gp)(i,0);

  // ----------------------------------------------------- trial stress
  // sigma^{trial}_{n+1} = cmat^e . epsilon^{e,trial}_{n+1}

  // consistent tangent modulus
  // for linear elastic material: cmat == cmat^e == cmat_stvenantkirchhoff
  SetupCmat(cmat);

  // evaluate stresses
  // sigma^{trial}_{n+1} = C . epsilon^{e,trial}_{n+1}
  LINALG::Matrix<NUM_STRESS_3D,1> sigmatr(true);
  sigmatr.MultiplyNN(cmat,epsilon_e);

  // split the trial stress into distortional and volumetric stress
  // deviatoric/distortional stress
  // G == mu
  // sigma^D^{trial}_{n+1} := s
  // s = 2 . G . epsilon_d^{e,trial}_{n+1}
  LINALG::Matrix<NUM_STRESS_3D,1> s(true);
  const double coeff = 2.0 * G;
  for (int i=0; i<6; i++)
    s(i) = coeff * epsilon_ed(i);

  // volumetric stress - pressure
  // p^{trial}_{n+1} = K . epsilon_v^{e,trial}_{n+1}
  LINALG::Matrix<NUM_STRESS_3D,1> press(true);
  for(int i=0; i<6; i++)
    press(i) = K * epsilon_ev(i);

  // q^(trial)_{n+1} := q(s^(trial)_{n+1}) = \sqrt{ 3/2 . s : s }
  double q = 0.0;
  for (int i=0; i<6; i++)
    q = sqrt( (3.0/2.0) * s(i) * s(i) );

  //-------------------------------------------------------------------
  // check plastic admissibility
  //-------------------------------------------------------------------

  // calculate the yield function
  // Phi = \sqrt{ 3.0 . J2 } - sigma_y

  // J2 = 0.5 . s : s
  double J2 = 0.0;
  for (int i=0; i<6; i++)
    J2 += s(i)*s(i);
    J2 *= 0.5;

  // Phi = sqrt(3.0 * J2) - sigma_y;
  Phi = sqrt(3.0 * J2) - sigma_y;

  //-------------------------------------------------------------------
  // IF elastic step (Phi < 0.0, Dgamma = 0.0)
  //-------------------------------------------------------------------
  if (Phi < 0.0)
  {
    // trial state vectors = result vectors of time step n+1
    // sigma^e_n+1 = sigma^(e,trial)_n+1
    stress.Update(1.0,sigmatr,0.0);

    // --------------------------------------------------- update history

    // epsilon^e_n+1 = epsilon^(e,trial)_n+1
    straine_np->at(gp) = epsilon_e;

    // epsilon^p_n+1 = epsilon^(p,trial)_n+1 = epsilon^p_n
    strainp_np->at(gp) = epsilon_p;

  }  // elastic step

  //-------------------------------------------------------------------
  // ELSE IF plastic step ( Phi = 0.0, Dgamma >= 0.0 )
  //-------------------------------------------------------------------
  else
  {
    // --------------------------------------------------- return mapping

    // Phi(stress) = Phi(s_n+1^(trial)) = sqrt{ 3 J2(s_{n+1}) } - sigma_y = 0
    // Phi = 0 = q^(trial)_{n+1} - 3 . G . Delta gamma - sigma_y

    // calculate incremental plastic multiplier
    // Delta gamma = ( 1 / (3G) ) . ( q_n+1^{trial} - sigma_y )
    Dgamma = (1.0/ (3*G) ) * ( q - sigma_y );

    // deviatoric stress
    // s = ( 1 - (3 . G . Delta gamma) / ( q_{n+1}^{trial} ) ) s_{n+1}^{trial}
    const double devfac = 1 - ( 3 * G * Dgamma ) / q;
    s.Scale(devfac);

    // flow vector N (Prandtl-Reuss)
    // (using the updated deviatoric stress s_n+1)
    // N = sqrt{3/2} . ( s_{n+1} / || s_{n+1} || )
    LINALG::Matrix<NUM_STRESS_3D,1> N(true);
    N.Update(1.0,s,0.0);
    const double fac = -sqrt(3.0/2.0) * ( 1.0 / ( s.Norm2() ) );
    N.Scale(fac);

    // epsilon^e_{n+1} = epsilon^(e,trial)_{n+1} - Dgamma . N
    epsilon_e.Update((-Dgamma),N,1.0);

    // epsilon^p_{n+1} = epsilon^p_n + Dgamma . N
    epsilon_p.Update((-Dgamma),N,1.0);

    // ------------------------------------------------ stress correction

    // total stress
    // sigma_{n+1} = s_{n+1} + p_{n+1} . I
    // pressure/volumetric stress no influence due to plasticity
    stress.Update(1.0,s,0.0);
    stress.Update(1.0,press,1.0);

    // --------------------------------------------------- update history

    // epsilon^e_n+1 = epsilon^(e,trial)_n+1
    straine_np->at(gp) = epsilon_e;

    // epsilon^p_n+1 = epsilon^(p,trial)_n+1 = epsilon^p_n
    strainp_np->at(gp) = epsilon_p;

  }  // plastic corrector

  return;

}  // Evaluate()


/*----------------------------------------------------------------------*
 | computes isotropic elasticity tensor in matrix notion     dano 04/11 |
 | for 3d                                                               |
 *----------------------------------------------------------------------*/
void MAT::PlasticLinElast::SetupCmat(LINALG::Matrix<6,6>& cmat)
{
  // get material parameters
  // Young's modulus (modulus of elasticity)
  double young = params_->youngs_;
  // Poisson's ratio
  double nu = params_->poissonratio_;

  // isotropic elasticity tensor C in Voigt matrix notation
  //                       [ 1-nu     nu     nu |          0    0    0 ]
  //                       [        1-nu     nu |          0    0    0 ]
  //           E           [               1-nu |          0    0    0 ]
  //   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~  ~~~ ]
  //       (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0    0 ]
  //                       [                    |      (1-2*nu)/2    0 ]
  //                       [ symmetric          |           (1-2*nu)/2 ]
  //
  const double mfac = young/((1.0+nu)*(1.0-2.0*nu));  // factor

  // clear the material tangent
  cmat.Clear();
  // write non-zero components
  cmat(0,0) = mfac*(1.0-nu);
  cmat(0,1) = mfac*nu;
  cmat(0,2) = mfac*nu;
  cmat(1,0) = mfac*nu;
  cmat(1,1) = mfac*(1.0-nu);
  cmat(1,2) = mfac*nu;
  cmat(2,0) = mfac*nu;
  cmat(2,1) = mfac*nu;
  cmat(2,2) = mfac*(1.0-nu);
  // ~~~
  cmat(3,3) = mfac*0.5*(1.0-2.0*nu);
  cmat(4,4) = mfac*0.5*(1.0-2.0*nu);
  cmat(5,5) = mfac*0.5*(1.0-2.0*nu);

}  // SetupCmat()


/*----------------------------------------------------------------------*/
#endif // CCADISCRET


