/*----------------------------------------------------------------------*/
/*!
\file viscogenmax.cpp

<pre>
Maintainer: Aline Bel
            brunon@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/



/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

// unnecessary
//#include <vector>
#include "viscogenmax.H"
#include "elasthyper.H"
#include "../drt_matelast/elast_summand.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ViscoGenMax::ViscoGenMax(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: ElastHyper(matdata),
//  nummat_(matdata->GetInt("NUMMAT")),
//  matids_(matdata->Get<std::vector<int> >("MATIDS")),
//  density_(matdata->GetDouble("DENS")),
//  gamma_(matdata->GetDouble("GAMMA")),
  relax_(matdata->GetDouble("RELAX")),
  beta_(matdata->GetDouble("BETA"))
//  init_mode_(matdata->GetInt("INIT_MODE"))
{
//  // check if sizes fit
//  if (nummat_ != (int)matids_->size())
//    dserror("number of materials %d does not fit to size of material vector %d", nummat_, matids_->size());
//
//  // make sure the referenced materials in material list have quick access parameters
//  std::vector<int>::const_iterator m;
//  for (m=matids_->begin(); m!=matids_->end(); ++m)
//  {
//    const int matid = *m;
//    Teuchos::RCP<MAT::ELASTIC::Summand> potsum = MAT::ELASTIC::Summand::Factory(matid);
//    if (potsum == Teuchos::null) dserror("Failed to allocate");
//    potsum_.insert(std::pair<int,Teuchos::RCP<MAT::ELASTIC::Summand> >(matid,potsum));
//  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ViscoGenMax::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ViscoGenMax(this));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const MAT::ELASTIC::Summand> MAT::PAR::ViscoGenMax::MaterialById(const int id) const
{
  std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::const_iterator m = potsum_.find(id);
  if (m == potsum_.end())
  {
    dserror("Material %d could not be found", id);
    return Teuchos::null;
  }
  else
  {
    return m->second;
  }
}


MAT::ViscoGenMaxType MAT::ViscoGenMaxType::instance_;


DRT::ParObject* MAT::ViscoGenMaxType::Create( const std::vector<char> & data )
{
  MAT::ViscoGenMax* viscogenmax = new MAT::ViscoGenMax();
  viscogenmax->Unpack(data);
  return viscogenmax;
}

///*----------------------------------------------------------------------*
// |  initialise static arrays                                 bborn 08/09|
// *----------------------------------------------------------------------*/
//// 6-Voigt C-index                              0 1 2  3 4 5
//const int MAT::ElastHyper::VOIGT6ROW_[6] = {0,1,2, 0,1,2};
//const int MAT::ElastHyper::VOIGT6COL_[6] = {0,1,2, 1,2,0};
//
//// tensor indices ij = 11, 12, 13, 21, 22, 23, 31, 32, 33
//// C indices           00, 01, 02, 10, 11, 12, 20, 21, 22
//// Access : 3*i+j
//// 6-Voigt C-indices    0   3   5   3   1   4   5   4   2
//const int MAT::ElastHyper::VOIGT3X3SYM_[9] = {0,3,5, 3,1,4, 5,4,2};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ViscoGenMax::ViscoGenMax()
  : params_(NULL)
{
  isinit_=false;
  histstressisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisoprinclast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinclast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressisomodisocurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodisolast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisocurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisolast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressisomodvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodvollast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvollast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressanisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressanisoprinclast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinclast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ViscoGenMax::ViscoGenMax(MAT::PAR::ViscoGenMax* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoGenMax::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  MAT::ElastHyper::Pack(data);

//  // matid
//  int matid = -1;
//  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
//  AddtoPack(data,matid);
//  AddtoPack(data,anisotropic_);
//  AddtoPack(data,a1_);
//  AddtoPack(data,a2_);
//  AddtoPack(data,A1_);
//  AddtoPack(data,A2_);
//  AddtoPack(data,A1A2_);
//  AddtoPack(data,HU_);
//  AddtoPack(data,HUlumen_);
//  AddtoPack(data,normdist_);

  //  pack history data
  int histsize;
  if (!Initialized())
  {
    histsize=0;
  }
  else
  {
    histsize = histstressisoprinclast_->size();
  }
  AddtoPack(data,histsize);  // Length of history vector(s)
  for (int var = 0; var < histsize; ++var)
  {
    AddtoPack(data,histstressisoprinclast_->at(var));
    AddtoPack(data,artstressisoprinclast_->at(var));
    AddtoPack(data,histstressisomodisolast_->at(var));
    AddtoPack(data,artstressisomodisolast_->at(var));
    AddtoPack(data,histstressisomodvollast_->at(var));
    AddtoPack(data,artstressisomodvollast_->at(var));
    AddtoPack(data,histstressanisoprinclast_->at(var));
    AddtoPack(data,artstressanisoprinclast_->at(var));
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoGenMax::Unpack(const std::vector<char>& data)
{
  isinit_=true;
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  MAT::ElastHyper::Unpack(basedata);

//  // matid and recover params_
//  int matid;
//  ExtractfromPack(position,data,matid);
//  params_ = NULL;
//  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
//    if (DRT::Problem::Instance()->Materials()->Num() != 0)
//    {
//      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
//      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
//      if (mat->Type() == MaterialType())
//        params_ = static_cast<MAT::PAR::ElastHyper*>(mat);
//      else
//        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
//    }
//
//  int anisotropic;
//
//  ExtractfromPack(position,data,anisotropic);
//  anisotropic_ = anisotropic != 0;
//
//  ExtractfromPack(position,data,a1_);
//  ExtractfromPack(position,data,a2_);
//  ExtractfromPack(position,data,A1_);
//  ExtractfromPack(position,data,A2_);
//  ExtractfromPack(position,data,A1A2_);
//  ExtractfromPack(position,data,HU_);
//  ExtractfromPack(position,data,HUlumen_);
//  ExtractfromPack(position,data,normdist_);

  // history data
  int histsize;
  ExtractfromPack(position,data,histsize);

  if (histsize == 0) isinit_=false;

  histstressisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisoprinclast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinclast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressisomodisocurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodisolast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisocurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisolast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressisomodvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodvollast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvollast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressanisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressanisoprinclast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinclast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  for (int var=0; var<histsize; var+=1)
  {
    LINALG::Matrix<NUM_STRESS_3D,1> tmp(true);
    histstressisoprinccurr_->push_back(tmp);
    artstressisoprinccurr_->push_back(tmp);
    histstressisomodisocurr_->push_back(tmp);
    artstressisomodisocurr_->push_back(tmp);
    histstressisomodvolcurr_->push_back(tmp);
    artstressisomodvolcurr_->push_back(tmp);
    histstressanisoprinccurr_->push_back(tmp);
    artstressanisoprinccurr_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    histstressisoprinclast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    artstressisoprinclast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    histstressisomodisolast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    artstressisomodisolast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    histstressisomodvollast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    artstressisomodvollast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    histstressanisoprinclast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    artstressanisoprinclast_->push_back(tmp);

  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
//int MAT::ViscoGenMax::MatID(
//  const unsigned index
//) const
//{
//  if ((int)index < params_->nummat_)
//    return params_->matids_->at(index);
//  else
//  {
//    dserror("Index too large");
//    return -1;
//  }
//}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
//        double MAT::ViscoGenMax::ShearMod() const
//        {
//          // principal coefficients
//          bool haveshearmod = false;
//          double shearmod = 0.0;
//          {
//            // loop map of associated potential summands
//            std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
//            std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
//            for (p=pot.begin(); p!=pot.end(); ++p)
//            {
//              p->second->AddShearMod(haveshearmod,shearmod);
//            }
//          }
//
//          if (haveshearmod)
//          {
//            return shearmod;
//          }
//          else
//          {
//            dserror("Cannot provide shear modulus equivalent");
//            return -1.0;
//          }
//        }

//      /*----------------------------------------------------------------------*/
//      /*----------------------------------------------------------------------*/
//      void MAT::ViscoGenMax::SetupAAA(Teuchos::ParameterList& params)
//      {
//        normdist_ = params.get("iltthick meanvalue",-999.0);
//        const ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
//        HUlumen_ = pslist.get<int>("MAXHULUMEN");
//        return;
//      }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoGenMax::Setup(const int numgp,DRT::INPUT::LineDefinition* linedef)
{
  MAT::ElastHyper::Setup(linedef);

//  if (linedef->HaveNamed("HU"))
//  {
//    linedef->ExtractDouble("HU",HU_);
//  }
//  else
//  {
//    HU_ = -999.0;
//  }
//  anisotropic_ = true;
//
//  // fibers aligned in local element cosy with gamma_i around circumferential direction
//  vector<double> rad;
//  vector<double> axi;
//  vector<double> cir;
//  if (not linedef->HaveNamed("RAD") or
//      not linedef->HaveNamed("AXI") or
//      not linedef->HaveNamed("CIR"))
//  {
//    //dserror("Reading of element local cosy failed");
//    anisotropic_=false;
//  }
//
//  else
//  {
//    // read local (cylindrical) cosy-directions at current element
//    // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
//    LINALG::Matrix<3,3> locsys(true);
//    linedef->ExtractDoubleVector("RAD",rad);
//    linedef->ExtractDoubleVector("AXI",axi);
//    linedef->ExtractDoubleVector("CIR",cir);
//    double radnorm=0.; double axinorm=0.; double cirnorm=0.;
//
//    for (int i = 0; i < 3; ++i)
//    {
//      radnorm += rad[i]*rad[i]; axinorm += axi[i]*axi[i]; cirnorm += cir[i]*cir[i];
//    }
//    radnorm = sqrt(radnorm); axinorm = sqrt(axinorm); cirnorm = sqrt(cirnorm);
//
//    for (int i=0; i<3; ++i)
//    {
//      locsys(i,0) = rad[i]/radnorm;
//      locsys(i,1) = axi[i]/axinorm;
//      locsys(i,2) = cir[i]/cirnorm;
//    }
//    // INIT_MODE = 0 : Fiber direction derived from local cosy
//    if( 0 == params_->init_mode_)
//    {
//      // alignment angles gamma_i are read from first entry of then unnecessary vectors a1 and a2
//      if ((params_->gamma_<0) || (params_->gamma_ >90)) dserror("Fiber angle not in [0,90]");
//      //convert
//      const double gamma = (params_->gamma_*PI)/180.;
//
//      for (int i = 0; i < 3; ++i)
//      {
//        // a1 = cos gamma e3 + sin gamma e2
//        a1_(i) = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
//        // a2 = cos gamma e3 - sin gamma e2
//        a2_(i) = cos(gamma)*locsys(i,2) - sin(gamma)*locsys(i,1);
//      }
//    }
//    // INIT_MODE = 1 : Fiber direction aligned to local cosy
//    else if (1 == params_->init_mode_)
//    {
//      for (int i = 0; i < 3; ++i)
//      {
//        a1_(i) = locsys(i,0);
//        a2_(i) = locsys(i,1);
//      }
//    }
//    // INIT_MODE = -1 = default value; usage of fiber direction without initialization mode
//    else if (-1 == params_->init_mode_)
//    {
//      dserror("Forgotten to give INIT_MODE in .dat-file");
//    }
//    else
//    {
//      dserror("Problem with fiber initialization");
//    }
//    for (int i = 0; i < 3; ++i) {
//      A1_(i) = a1_(i)*a1_(i);
//      A2_(i) = a2_(i)*a2_(i);
//      for (int j=0; j<3; j++)
//      {
//        A1A2_(j,i) = a1_(j)*a2_(i);
//      }
//    }
//    A1_(3) = a1_(0)*a1_(1); A1_(4) = a1_(1)*a1_(2); A1_(5) = a1_(0)*a1_(2);
//    A2_(3) = a2_(0)*a2_(1); A2_(4) = a2_(1)*a2_(2); A2_(5) = a2_(0)*a2_(2);
//  }
//  //return;

  // Initialise/allocate internal stress variables

  histstressisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisoprinclast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinclast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressisomodisocurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodisolast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisocurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisolast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressisomodvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodvollast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvollast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressanisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressanisoprinclast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinclast_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  const LINALG::Matrix<NUM_STRESS_3D,1> emptyvec(true);
  histstressisoprinccurr_->resize(numgp);
  histstressisoprinclast_->resize(numgp);
  artstressisoprinccurr_->resize(numgp);
  artstressisoprinclast_->resize(numgp);

  histstressisomodisocurr_->resize(numgp);
  histstressisomodisolast_->resize(numgp);
  artstressisomodisocurr_->resize(numgp);
  artstressisomodisolast_->resize(numgp);

  histstressisomodvolcurr_->resize(numgp);
  histstressisomodvollast_->resize(numgp);
  artstressisomodvolcurr_->resize(numgp);
  artstressisomodvollast_->resize(numgp);

  histstressanisoprinccurr_->resize(numgp);
  histstressanisoprinclast_->resize(numgp);
  artstressanisoprinccurr_->resize(numgp);
  artstressanisoprinclast_->resize(numgp);

  for (int j=0; j<numgp; ++j)
  {
    histstressisoprinccurr_->at(j) = emptyvec;
    histstressisoprinclast_->at(j) = emptyvec;
    artstressisoprinccurr_->at(j) = emptyvec;
    artstressisoprinclast_->at(j) = emptyvec;
    histstressisomodisocurr_->at(j) = emptyvec;
    histstressisomodisolast_->at(j) = emptyvec;
    artstressisomodisocurr_->at(j) = emptyvec;
    artstressisomodisolast_->at(j) = emptyvec;
    histstressisomodvolcurr_->at(j) = emptyvec;
    histstressisomodvollast_->at(j) = emptyvec;
    artstressisomodvolcurr_->at(j) = emptyvec;
    artstressisomodvollast_->at(j) = emptyvec;
    histstressanisoprinccurr_->at(j) = emptyvec;
    histstressanisoprinclast_->at(j) = emptyvec;
    artstressanisoprinccurr_->at(j) = emptyvec;
    artstressanisoprinclast_->at(j) = emptyvec;
  }
  isinit_ = true;
  return ;
}

/*----------------------------------------------------------------------*
 |  Update internal stress variables              (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoGenMax::Update()
{
  histstressisoprinclast_=histstressisoprinccurr_;
  artstressisoprinclast_=artstressisoprinccurr_;
  histstressisomodisolast_=histstressisomodisocurr_;
  artstressisomodisolast_=artstressisomodisocurr_;
  histstressisomodvollast_=histstressisomodvolcurr_;
  artstressisomodvollast_=artstressisomodvolcurr_;
  histstressanisoprinclast_=histstressanisoprinccurr_;
  artstressanisoprinclast_=artstressanisoprinccurr_;
  const LINALG::Matrix<NUM_STRESS_3D,1> emptyvec(true);
  histstressisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodisocurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisocurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvolcurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressanisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinccurr_=rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  const int numgp=histstressisoprinclast_->size();
  histstressisoprinccurr_->resize(numgp);
  artstressisoprinccurr_->resize(numgp);
  histstressisomodisocurr_->resize(numgp);
  artstressisomodisocurr_->resize(numgp);
  histstressisomodvolcurr_->resize(numgp);
  artstressisomodvolcurr_->resize(numgp);
  histstressanisoprinccurr_->resize(numgp);
  artstressanisoprinccurr_->resize(numgp);
  for (int j=0; j<numgp; ++j)
  {
    histstressisoprinccurr_->at(j) = emptyvec;
    artstressisoprinccurr_->at(j) = emptyvec;
    histstressisomodisocurr_->at(j) = emptyvec;
    artstressisomodisocurr_->at(j) = emptyvec;
    histstressisomodvolcurr_->at(j) = emptyvec;
    artstressisomodvolcurr_->at(j) = emptyvec;
    histstressanisoprinccurr_->at(j) = emptyvec;
    artstressanisoprinccurr_->at(j) = emptyvec;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Reset internal stress variables               (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoGenMax::Reset()
{
  // do nothing,
  // because #histstresscurr_ and #artstresscurr_ are recomputed anyway at every iteration
  // based upon #histstresslast_ and #artstresslast_ untouched within time step
  return;
}

//            /*----------------------------------------------------------------------*/
//            /*----------------------------------------------------------------------*/
//            void MAT::ViscoGenMax::InvariantsPrincipal(
//              LINALG::Matrix<3,1>& prinv,
//              const LINALG::Matrix<6,1>& rcg
//              )
//            {
//              // 1st invariant, trace
//              prinv(0) = rcg(0) + rcg(1) + rcg(2);
//              // 2nd invariant
//              prinv(1) = 0.5*( prinv(0)*prinv(0)
//                               - rcg(0)*rcg(0) - rcg(1)*rcg(1) - rcg(2)*rcg(2)
//                               - .5*rcg(3)*rcg(3) - .5*rcg(4)*rcg(4) - .5*rcg(5)*rcg(5) );
//              // 3rd invariant, determinant
//              prinv(2) = rcg(0)*rcg(1)*rcg(2)
//                + 0.25 * rcg(3)*rcg(4)*rcg(5)
//                - 0.25 * rcg(1)*rcg(5)*rcg(5)
//                - 0.25 * rcg(2)*rcg(3)*rcg(3)
//                - 0.25 * rcg(0)*rcg(4)*rcg(4);
//
//              return;
//            }
//
//            /*----------------------------------------------------------------------*/
//            /*----------------------------------------------------------------------*/
//            void MAT::ElastHyper::InvariantsModified(
//              LINALG::Matrix<3,1>& modinv,
//              const LINALG::Matrix<3,1>& prinv
//              )
//            {
//              // 1st invariant, trace
//              modinv(0) = prinv(0)*std::pow(prinv(2),-1./3.);
//              // 2nd invariant
//              modinv(1) = prinv(1)*std::pow(prinv(2),-2./3.);
//              // J
//              modinv(2) = std::pow(prinv(2),1./2.);
//
//              return;
//            }
//
//            /*----------------------------------------------------------------------*/
//            /*----------------------------------------------------------------------*/
//            void MAT::ElastHyper::InvariantsPrincipalAniso(
//              LINALG::Matrix<6,1>& prinv,
//              const LINALG::Matrix<6,1>& rcg
//              )
//            {
//              // 1st invariant, trace
//              prinv(0) = rcg(0) + rcg(1) + rcg(2);
//              // 2nd invariant
//              prinv(1) = 0.5*( prinv(0)*prinv(0)
//                               - rcg(0)*rcg(0) - rcg(1)*rcg(1) - rcg(2)*rcg(2)
//                               - .5*rcg(3)*rcg(3) - .5*rcg(4)*rcg(4) - .5*rcg(5)*rcg(5) );
//              // 3rd invariant, determinant
//              prinv(2) = rcg(0)*rcg(1)*rcg(2)
//                + 0.25 * rcg(3)*rcg(4)*rcg(5)
//                - 0.25 * rcg(1)*rcg(5)*rcg(5)
//                - 0.25 * rcg(2)*rcg(3)*rcg(3)
//                - 0.25 * rcg(0)*rcg(4)*rcg(4);
//
//              prinv(3) =  A1_(0)*rcg(0) + A1_(1)*rcg(1) + A1_(2)*rcg(2)
//                        + A1_(3)*rcg(3) + A1_(4)*rcg(4) + A1_(5)*rcg(5);
//
//              prinv(4) =  A2_(0)*rcg(0) + A2_(1)*rcg(1) + A2_(2)*rcg(2)
//                        + A2_(3)*rcg(3) + A2_(4)*rcg(4) + A2_(5)*rcg(5);
//
//              prinv(5) =  A1A2_(0,0)*rcg(0) + A1A2_(1,1)*rcg(1) + A1A2_(2,2)*rcg(2)
//                        + 0.5*(A1A2_(0,1)*rcg(3) + A1A2_(1,2)*rcg(4) + A1A2_(0,2)*rcg(5))
//                        + 0.5*(A1A2_(1,0)*rcg(3) + A1A2_(2,1)*rcg(4) + A1A2_(2,0)*rcg(5));
//
//              return;
//            }
//
//            /*----------------------------------------------------------------------*/
//            /*----------------------------------------------------------------------*/
//            void MAT::ElastHyper::StretchesPrincipal(
//              LINALG::Matrix<3,1>& prstr,
//              LINALG::Matrix<3,3>& prdir,
//              const LINALG::Matrix<6,1>& rcg
//              )
//            {
//              // create right Cauchy-Green 2-tensor
//              LINALG::Matrix<3,3> rcgt(false);
//              rcgt(0,0) = rcg(0);
//              rcgt(1,1) = rcg(1);
//              rcgt(2,2) = rcg(2);
//              rcgt(0,1) = rcgt(1,0) = 0.5*rcg(3);
//              rcgt(1,2) = rcgt(2,1) = 0.5*rcg(4);
//              rcgt(2,0) = rcgt(0,2) = 0.5*rcg(5);
//
//              // eigenvalue decomposition
//              LINALG::Matrix<3,3> prstr2;  // squared principal stretches
//              LINALG::SYEV(rcgt,prstr2,prdir);
//
//              // THE principal stretches
//              for (int al=0; al<3; ++al) prstr(al) = std::sqrt(prstr2(al,al));
//
//              // bye
//              return;
//            }
//
//            /*----------------------------------------------------------------------*/
//            /*----------------------------------------------------------------------*/
//            void MAT::ElastHyper::StretchesModified(
//              LINALG::Matrix<3,1>& modstr,
//              const LINALG::Matrix<3,1>& prstr
//              )
//            {
//              // determinant of deformation gradient
//              const double detdefgrad = prstr(0)*prstr(1)*prstr(2);
//
//              // determine modified principal stretches
//              modstr.Update(std::pow(detdefgrad,-1.0/3.0),prstr);
//
//              return;
//            }
//
//            /*----------------------------------------------------------------------*/
//            /*----------------------------------------------------------------------*/
//            bool MAT::ElastHyper::HaveCoefficientsStretchesPrincipal()
//            {
//              // set default
//              bool havecoeff = false;
//
//              // loop map of associated potential summands and see
//              {
//                std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
//                std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
//                for (p=pot.begin(); p!=pot.end(); ++p)
//                {
//                  havecoeff = havecoeff or p->second->HaveCoefficientsStretchesPrincipal();
//                }
//              }
//
//              // deliver
//              return havecoeff;
//            }
//
//            /*----------------------------------------------------------------------*/
//            /*----------------------------------------------------------------------*/
//            bool MAT::ElastHyper::HaveCoefficientsStretchesModified()
//            {
//              // set default
//              bool havecoeff = false;
//
//              // loop map of associated potential summands and see
//              {
//                std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >& pot = params_->potsum_;
//                std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> >::iterator p;
//                for (p=pot.begin(); p!=pot.end(); ++p)
//                {
//                  havecoeff = havecoeff or p->second->HaveCoefficientsStretchesModified();
//                }
//              }
//
//              // deliver
//              return havecoeff;
//            }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoGenMax::Evaluate(
  const LINALG::Matrix<6,1>& glstrain,
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress,
  const int gp,
  Teuchos::ParameterList& params
  )
{
  //initialize scalars
  double artscalar1;
  double artscalar2;
  double scalarvisco;

  //  Get relaxation time and time integration constant
  double tau  = params_->relax_;
  double beta = params_->beta_;
  const double theta= 0.5 ;//params_->theta_;

  // get time algorithmic parameters
  // NOTE: dt can be zero (in restart of STI) for Generalized Maxwell model
  // there is no special treatment required. Adaptation for Kelvin-Voigt were necessary.
  double dt = params.get<double>("delta time"); // TIMESTEP in the .dat file

//#define GEN_MAXWELL
//#ifdef GEN_MAXWELL
//  double tau1=tau;
  // evaluate "alpha" factors which distribute stress or stiffness between parallel springs
  // sum_0^i alpha_j = 1
//  alpha0 = E_s / E_f;
//  alpha1 = 1.0 - alpha0;
//
//  // evaluate Lame constants, bulk modulus
//  lambda = nue*E_f / ((1.0+nue)*(1.0-2.0*nue));
//  mue = E_f / (2.0*(1.0+nue));
//  kappa = lambda + 2.0/3.0 * mue;

  // WE HAVE TO FIND A WAY TO CALCULATE alpha0, alpha1, tau and beta
//  alpha0 = 1.;
//  alpha1 = 1.;

  // evaluate scalars to compute
  // Q^(n+1) = tau/(tau+theta*dt) [(tau-dt+theta*dt)/tau Q + beta(S^(n+1) - S^n)]
  artscalar1=(tau - dt + theta*dt)/tau;
  artscalar2=tau/(tau + theta*dt);

  // factor to calculate visco stiffness matrix from elastic stiffness matrix
  scalarvisco = 1+beta*exp(-dt/(2*tau));//+alpha1*tau/(tau+theta*dt);

//#else
//  //in this case we have a parallel layout of a spring and a dashpot,
//  //so no stress distribution between parallel springs
////  alpha0 = 1.;
////  alpha1 = 1.;
//
//  // do we have to propagate in time?
//  if (dt >0.0)
//  {
//    // evaluate scalars to compute
//    // Q^(n+1) = tau/(theta*dt) [(-dt+theta*dt)/tau Q + S^(n+1) - S^n]
//    artscalar1=(tau-dt+theta*dt)/tau;
//    artscalar2=tau/(tau+theta*dt);
//
//    // factor to calculate visco stiffness matrix from elastic stiffness matrix
//    scalarvisco = 1.0+tau/(tau+theta*dt);
//  }
//  else
//  {
//    // in case we do not want to propagate in time, Q^{n+1} = Q^{n}
//    artscalar1 = 1.0;
//    artscalar2 = 1.0;
//    // factor to calculate visco stiffness matrix from elastic stiffness matrix
//    scalarvisco = 2.0;
//  }
//
//#endif

//  bool havecoeffprinc_ = false;
//  bool havecoeffmodi_ = false;
//  bool havecoeffpraniso_ = false;

  LINALG::Matrix<6,1> id2(true) ;
  LINALG::Matrix<6,1> rcg(true) ;
  LINALG::Matrix<6,1> scg(true) ;
  LINALG::Matrix<6,1> icg(true) ;
  LINALG::Matrix<6,6> id4(true) ;
  LINALG::Matrix<6,6> id4sharp(true) ;

  LINALG::Matrix<3,1> prinv(true);
  LINALG::Matrix<3,1> modinv(true);
  LINALG::Matrix<6,1> pranisoinv(true);

  LINALG::Matrix<3,1> gamma(true);
  LINALG::Matrix<8,1> delta(true);
  LINALG::Matrix<3,1> modgamma(true);
  LINALG::Matrix<5,1> moddelta(true);
  LINALG::Matrix<3,1> anisogamma(true);
  LINALG::Matrix<15,1> anisodelta(true);

  EvaluateKinQuant(glstrain,id2,scg,rcg,icg,id4,id4sharp,prinv,modinv,pranisoinv);
  EvaluateGammaDelta(rcg,prinv,modinv,pranisoinv,gamma,delta,modgamma,moddelta,anisogamma,anisodelta);

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress.Clear();
  cmat.Clear();

  if (isoprinc_)
    {
    LINALG::Matrix<NUM_STRESS_3D,1> stressisoprinc(true) ;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisoprinc(true) ;
    EvaluateIsotropicPrinc(stressisoprinc,cmatisoprinc,scg,id2,icg,id4sharp,gamma,delta);
    stress.Update(1.0, stressisoprinc, 1.0);
    cmat.Update(1.0,cmatisoprinc,1.0);

    // viscous part of the stress
    // read history
    LINALG::Matrix<NUM_STRESS_3D,1> S_n (histstressisoprinclast_->at(gp));
    S_n.Scale(-beta);
    LINALG::Matrix<NUM_STRESS_3D,1> Q_n (artstressisoprinclast_->at(gp));

    // artificial visco stresses
    LINALG::Matrix<NUM_STRESS_3D,1> Q(Q_n);
    Q.Scale(artscalar1);
    stressisoprinc.Scale(beta);
    Q += stressisoprinc;
    Q += S_n;
    Q.Scale(artscalar2);  // Q^(n+1) = artscalar2* [artscalar1* Q + beta*(S^(n+1) - S^n)]

    // update history
    histstressisoprinccurr_->at(gp) = stressisoprinc;
    artstressisoprinccurr_->at(gp) = Q;

    // add the viscous to the elastic contribution
    stress.Update(1.0,Q,1.0);

    // constitutive tensor
    // viscous contribution : Cmat = Cmatx(1+delta) with (1+delta)=scalarvisco
    cmat.Scale(scalarvisco);
    }


  if (isomod_)
  {
    LINALG::Matrix<NUM_STRESS_3D,1> stressisomodiso(true) ;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisomodiso(true) ;
    LINALG::Matrix<NUM_STRESS_3D,1> stressisomodvol(true) ;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisomodvol(true) ;
    EvaluateIsotropicMod(stressisomodiso,stressisomodvol,cmatisomodiso,cmatisomodvol,rcg,id2,icg,id4,id4sharp,modinv,prinv,modgamma,moddelta);
    stress.Update(1.0, stressisomodiso, 1.0);
    stress.Update(1.0, stressisomodvol, 1.0);
    cmat.Update(1.0,cmatisomodiso,1.0);
    cmat.Update(1.0,cmatisomodvol,1.0);

    // viscous part of the isochoric contribution
    //-----------------------------------------------------------------------
    // read history
    LINALG::Matrix<NUM_STRESS_3D,1> Siso_n (histstressisomodisolast_->at(gp));
    Siso_n.Scale(-beta);
    LINALG::Matrix<NUM_STRESS_3D,1> Qiso_n (artstressisomodisolast_->at(gp));

    // artificial visco stresses
    LINALG::Matrix<NUM_STRESS_3D,1> Qiso(Qiso_n);
    Qiso.Scale(artscalar1);
    stressisomodiso.Scale(beta);
    Qiso += stressisomodiso;
    Qiso += Siso_n;
    Qiso.Scale(artscalar2);  // Q^(n+1) = artscalar2* [artscalar1* Q + beta*(S^(n+1) - S^n)]

    // update history
    histstressisomodisocurr_->at(gp) = stressisomodiso;
    artstressisomodisocurr_->at(gp) = Qiso;

    // add the viscous to the elastic contribution
    stress.Update(1.0,Qiso,1.0); //--> at this point, stress includes the isochoric viscoelastic contributions
    //----------------------------------------------------------------------


    // viscous part of the volumetric contribution
    //-----------------------------------------------------------------------
    // read history
    LINALG::Matrix<NUM_STRESS_3D,1> Svol_n (histstressisomodvollast_->at(gp));
    Svol_n.Scale(-beta);
    LINALG::Matrix<NUM_STRESS_3D,1> Qvol_n (artstressisomodvollast_->at(gp));

    // artificial visco stresses
    LINALG::Matrix<NUM_STRESS_3D,1> Qvol(Qvol_n);
    Qvol.Scale(artscalar1);
    stressisomodvol.Scale(beta);
    Qvol += stressisomodvol;
    Qvol += Svol_n;
    Qvol.Scale(artscalar2);  // Q^(n+1) = artscalar2* [artscalar1* Q + S^(n+1) - S^n]

    // update history
    histstressisomodvolcurr_->at(gp) = stressisomodvol;
    artstressisomodvolcurr_->at(gp) = Qvol;

    // add the viscous contribution
    stress.Update(1.0,Qvol,1.0);//--> at this point, stress includes the isochoric and volumetric viscoelastic contributions
    //----------------------------------------------------------------------

    // constitutive tensor

    // Cmat = Cmat_iso_inf(1+delta_iso) + Cmat_vol_inf(1+delta_vol)
    // 1+delta = scalarvisco ; could be two different scalars for iso and vol contributions
    // overall contributions by taking into account the viscous contribution
    cmat.Scale(scalarvisco);
  }

  /*----------------------------------------------------------------------*/
  // coefficients in principal stretches
  const bool havecoeffstrpr = HaveCoefficientsStretchesPrincipal();
  const bool havecoeffstrmod = HaveCoefficientsStretchesModified();
  if (havecoeffstrpr or havecoeffstrmod) {
    ResponseStretches(cmat,stress,rcg,havecoeffstrpr,havecoeffstrmod);
  }

  /*----------------------------------------------------------------------*/
  //Do all the anisotropic stuff!
  if (anisoprinc_)
  {
      LINALG::Matrix<NUM_STRESS_3D,1> stressanisoprinc(true) ;
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatanisoprinc(true) ;
      EvaluateAnisotropicPrinc(stressanisoprinc,cmatanisoprinc,scg,id2,icg,anisogamma,anisodelta);
      stress.Update(1.0, stressanisoprinc, 1.0);
      cmat.Update(1.0, cmatanisoprinc, 1.0);

      // viscous part of the stress
      // read history
      LINALG::Matrix<NUM_STRESS_3D,1> S_n (histstressanisoprinclast_->at(gp));
      S_n.Scale(-beta);
      LINALG::Matrix<NUM_STRESS_3D,1> Q_n (artstressanisoprinclast_->at(gp));

      // artificial visco stresses
      LINALG::Matrix<NUM_STRESS_3D,1> Q(Q_n);
      Q.Scale(artscalar1);
      stressanisoprinc.Scale(beta);
      Q += stressanisoprinc;
      Q += S_n;
      Q.Scale(artscalar2);  // Q^(n+1) = artscalar2* [artscalar1* Q + S^(n+1) - S^n]

      // update history
      histstressanisoprinccurr_->at(gp) = stressanisoprinc;
      artstressanisoprinccurr_->at(gp) = Q;

      // add the viscous to the elastic contribution
      stress.Update(1.0,Q,1.0);

      // constitutive tensor
      // viscous contribution : Cmat = Cmatx(1+delta) with (1+delta)=scalarvisco
      cmat.Scale(scalarvisco);
  }

}
///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//void MAT::ElastHyper::Evaluate(const Epetra_SerialDenseVector* glstrain_e,
//                               Epetra_SerialDenseMatrix* cmat_e,
//                               Epetra_SerialDenseVector* stress_e,
//                               const int gp,   ///< current Gauss point
//                               Teuchos::ParameterList& params)  ///< parameter list for communication
//{
//  // this is temporary as long as the material does not have a
//  // Matrix-type interface
//  const LINALG::Matrix<6,1> glstrain(glstrain_e->A(),true);
//  LINALG::Matrix<6,6> cmat(cmat_e->A(),true);
//  LINALG::Matrix<6,1> stress(stress_e->A(),true);
//
//  Evaluate(glstrain,cmat,stress,gp,params);
//}

#endif
