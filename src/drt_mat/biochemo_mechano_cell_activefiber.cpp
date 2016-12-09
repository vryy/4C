/*!----------------------------------------------------------------------
\file biochemo_mechano_cell_activefiber.cpp

\brief Implementation of Biochemo-Mechano Coupled Stress Fiber Formation and Dissociation.

\level 3

\maintainer Andreas Rauch

*----------------------------------------------------------------------*/

#include "biochemo_mechano_cell_activefiber.H"

#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"

#include "../drt_io/io_control.H"

#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_so3/so3_defines.H"
#include "../drt_so3/so_hex8.H"

#include "../drt_immersed_problem/immersed_field_exchange_manager.H"

#include <Epetra_SerialDenseSolver.h>


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::BioChemoMechanoCellActiveFiber::BioChemoMechanoCellActiveFiber(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(    matdata                        ),
  density_(     matdata->GetDouble("DENS")     ),
  idmatpassive_(matdata->GetInt("IDMATPASSIVE")),
  kforwards_(   matdata->GetDouble("KFOR")     ),
  kbackwards_(  matdata->GetDouble("KBACK")    ),
  kRockEta_(  matdata->GetDouble("KROCKETA")    ),
  kActin_(  matdata->GetDouble("KACTIN")    ),
  ratemax_(  matdata->GetDouble("RATEMAX")    ),
  nmax_(   matdata->GetDouble("NMAX")     ),
  kstress_(    matdata->GetDouble("KSTRESS")   ),
  sourceconst_(    matdata->GetDouble("SOURCE")   ),
  myintmethod_(*(matdata->Get<std::string>("METHOD")))
{
  if(DRT::Problem::Instance()->ProblemType()!=prb_immersed_cell)
    dserror("Material 'BioChemoMechanoCellActiveFiber' is only available in problems of type Immersed_CellMigration.\n"
            "If you want to use it in a different context, have a look at the necessary setup in immersed_problem_dyn.cpp.");

  if(DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->CellMigrationParams().sublist("STRUCTURAL DYNAMIC"),"MATERIALTANGENT"))
    analyticalmaterialtangent_ = false;
  else
    analyticalmaterialtangent_ = true;
}


Teuchos::RCP<MAT::Material> MAT::PAR::BioChemoMechanoCellActiveFiber::CreateMaterial()
{
  return Teuchos::rcp(new MAT::BioChemoMechanoCellActiveFiber(this));
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::BioChemoMechanoCellActiveFiberType MAT::BioChemoMechanoCellActiveFiberType::instance_;


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::BioChemoMechanoCellActiveFiberType::Create( const std::vector<char> & data )
{
  MAT::BioChemoMechanoCellActiveFiber* cellcontraction = new MAT::BioChemoMechanoCellActiveFiber();
  cellcontraction->Unpack(data);
  return cellcontraction;
}


/*----------------------------------------------------------------------*
 |  Constructor                                             rauch  01/16|
 *----------------------------------------------------------------------*/
MAT::BioChemoMechanoCellActiveFiber::BioChemoMechanoCellActiveFiber()
  : params_(NULL),
    isinit_(false)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                                        rauch  01/16|
 *----------------------------------------------------------------------*/
MAT::BioChemoMechanoCellActiveFiber::BioChemoMechanoCellActiveFiber(MAT::PAR::BioChemoMechanoCellActiveFiber* params)
  : params_(params),
    isinit_(false)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                                    rauch  01/16|
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
    AddtoPack(data,matid);

  // pack history data
  int histsize;
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

  AddtoPack(data,histsize);  // length of history vector(s)
  for (int var=0; var<histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data,histdefgrdlast_->at(var));
    AddtoPack(data,histdefgrdlastsurf_->at(var));
    AddtoPack(data,strainratelast_->at(var));
    AddtoPack(data,strainratelastsurf_->at(var));
    AddtoPack(data,etalast_->at(var));

    AddtoPack(data,etalastphi_->at(var));
    AddtoPack(data,etalastphisurf_->at(var));

    AddtoPack(data,Nfil_->at(var));
    AddtoPack(data,rate_->at(var));
    AddtoPack(data,D_->at(var));
    AddtoPack(data,clast_->at(var));
    AddtoPack(data,nsfhor_->at(var));
    AddtoPack(data,nsfver_->at(var));
    AddtoPack(data,nsfdiagup_->at(var));
    AddtoPack(data,nsfdiagdown_->at(var));
  }

  // Pack data of passive elastic material
  if (matpassive_!=Teuchos::null) {
    matpassive_->Pack(data);
  }

  return;

}   // Pack()


/*----------------------------------------------------------------------*
 |  Unpack                                                  rauch  01/16|
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::Unpack(const std::vector<char>& data)
{
  isinit_=true;
  std::vector<char>::size_type position = 0;
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
        params_ = static_cast<MAT::PAR::BioChemoMechanoCellActiveFiber*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // history data
  int histsize;
  ExtractfromPack(position,data,histsize);

  // if system is not yet initialized, the history vectors have to be initialized
  if (histsize == 0){
    isinit_=false;
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
    return;
  }

  // unpack deformation gradient matrices
  histdefgrdlast_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );
  histdefgrdlastsurf_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );
  histdefgrdcurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );
  histdefgrdcurrsurf_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );
  strainratelast_= Teuchos::rcp( new std::vector<LINALG::Matrix<6,1> > );
  strainratelastsurf_= Teuchos::rcp( new std::vector<LINALG::Matrix<6,1> > );
  strainratecurr_= Teuchos::rcp( new std::vector<LINALG::Matrix<6,1> > );
  strainratecurrsurf_= Teuchos::rcp( new std::vector<LINALG::Matrix<6,1> > );

  // unpack formation level eta
  etalast_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );
  etacurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );
//#ifdef option3
  etalastphi_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );
  etalastphisurf_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );
  etacurrphi_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );
  etacurrphisurf_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );
//#endif

  Nfil_  = Teuchos::rcp( new std::vector<double> );
  rate_  = Teuchos::rcp( new std::vector<double> );
  D_  = Teuchos::rcp( new std::vector<double> );
  clast_  = Teuchos::rcp( new std::vector<double> );
  ccurr_  = Teuchos::rcp( new std::vector<double> );
  nsfhor_= Teuchos::rcp( new std::vector<double> );
  nsfver_= Teuchos::rcp( new std::vector<double> );
  nsfdiagup_= Teuchos::rcp( new std::vector<double> );
  nsfdiagdown_= Teuchos::rcp( new std::vector<double> );
  for (int var=0; var<histsize; ++var)
  {
    // initialize
    LINALG::Matrix<3,3> tmp_matrix3x3(true);
    LINALG::Matrix<6,1> tmp_matrix6x1(true);
    LINALG::Matrix<nphi,1> tmp_matrix(true);

    //LINALG::Matrix<numbgp,twice> tmp_matrix(true);
    double tmp_scalar = 0.0;

    // matrices of last converged state are unpacked
    ExtractfromPack(position,data,tmp_matrix3x3);
    histdefgrdlast_->push_back(tmp_matrix3x3);
    ExtractfromPack(position,data,tmp_matrix3x3);
    histdefgrdlastsurf_->push_back(tmp_matrix3x3);
    ExtractfromPack(position,data,tmp_matrix6x1);
    strainratelast_->push_back(tmp_matrix6x1);
    ExtractfromPack(position,data,tmp_matrix6x1);
    strainratelastsurf_->push_back(tmp_matrix6x1);
    //  of last converged state are unpacked
    ExtractfromPack(position,data,tmp_matrix3x3);
    etalast_->push_back(tmp_matrix3x3);


//#ifdef option3
    ExtractfromPack(position,data,tmp_matrix);
    etalastphi_->push_back(tmp_matrix);
    ExtractfromPack(position,data,tmp_matrix);
    etalastphisurf_->push_back(tmp_matrix);
//#endif
//    ExtractfromPack(position, data, tmp_matrix);
//    sigmaomegaphilast_->push_back(tmp_matrix);
    ExtractfromPack(position, data, tmp_scalar);
    Nfil_->push_back(tmp_scalar);
    ExtractfromPack(position, data, tmp_scalar);
    rate_->push_back(tmp_scalar);
    ExtractfromPack(position, data, tmp_scalar);
    D_->push_back(tmp_scalar);
    ExtractfromPack(position, data, tmp_scalar);
    clast_->push_back(tmp_scalar);
    ExtractfromPack(position, data, tmp_scalar);
    nsfhor_->push_back(tmp_scalar);
    ExtractfromPack(position, data, tmp_scalar);
    nsfver_->push_back(tmp_scalar);
    ExtractfromPack(position, data, tmp_scalar);
    nsfdiagup_->push_back(tmp_scalar);
    ExtractfromPack(position, data, tmp_scalar);
    nsfdiagdown_->push_back(tmp_scalar);
    // current vectors have to be initialized
    histdefgrdcurr_->push_back(tmp_matrix3x3);
    histdefgrdcurrsurf_->push_back(tmp_matrix3x3);
    strainratecurr_->push_back(tmp_matrix6x1);
    strainratecurrsurf_->push_back(tmp_matrix6x1);
    etacurr_->push_back(tmp_matrix3x3);
//#ifdef option3
//    LINALG::Matrix<nphi,1> tmp_matrix(true);
    etacurrphi_->push_back(tmp_matrix);
    etacurrphisurf_->push_back(tmp_matrix);
    ccurr_->push_back(tmp_scalar);
//#endif
    //sigmaomegaphicurr_->push_back(tmp_matrix);
  }

  // Unpack data of passive elastic material (these lines are copied from drt_element.cpp)
  std::vector<char> dataelastic;
  ExtractfromPack(position,data,dataelastic);
  if (dataelastic.size()>0) {
    DRT::ParObject* o = DRT::UTILS::Factory(dataelastic);  // Unpack is done here
    MAT::So3Material* matel = dynamic_cast<MAT::So3Material*>(o);
    if (matel==NULL)
      dserror("failed to unpack passive material");
    matpassive_ = Teuchos::rcp(matel);
  } else matpassive_ = Teuchos::null;


  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;

}   // Unpack()


/*----------------------------------------------------------------------*
 | initialize / allocate internal variables (public)        rauch  01/16|
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // initialize history variables
  histdefgrdlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  histdefgrdlastsurf_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  histdefgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  histdefgrdcurrsurf_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  strainratelast_= Teuchos::rcp( new std::vector<LINALG::Matrix<6,1> > );
  strainratelastsurf_= Teuchos::rcp( new std::vector<LINALG::Matrix<6,1> > );
  strainratecurr_= Teuchos::rcp( new std::vector<LINALG::Matrix<6,1> > );
  strainratecurrsurf_= Teuchos::rcp( new std::vector<LINALG::Matrix<6,1> > );

  etalast_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );
  etacurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );

  etalastphi_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );
  etalastphisurf_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );
  etacurrphi_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );
  etacurrphisurf_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );

  Nfil_ = Teuchos::rcp( new std::vector<double> );
  rate_ = Teuchos::rcp( new std::vector<double> );
  D_ = Teuchos::rcp( new std::vector<double> );
  clast_=Teuchos::rcp( new std::vector<double> );
  ccurr_=Teuchos::rcp( new std::vector<double> );
  nsfhor_=Teuchos::rcp( new std::vector<double> );
  nsfver_=Teuchos::rcp( new std::vector<double> );
  nsfdiagup_=Teuchos::rcp( new std::vector<double> );
  nsfdiagdown_=Teuchos::rcp( new std::vector<double> );
  // set all history variables to zero
  histdefgrdlast_->resize(numgp);
  histdefgrdlastsurf_->resize(numgp);
  histdefgrdcurr_->resize(numgp);
  histdefgrdcurrsurf_->resize(numgp);
  strainratelast_->resize(numgp);
  strainratelastsurf_->resize(numgp);
  strainratecurr_->resize(numgp);
  strainratecurrsurf_->resize(numgp);

  etalast_->resize(numgp);
  etacurr_->resize(numgp);

  etalastphi_->resize(numgp);
  etalastphisurf_->resize(numgp);
  etacurrphi_->resize(numgp);
  etacurrphisurf_->resize(numgp);
  Nfil_->resize(numgp);
  rate_->resize(numgp);
  D_->resize(numgp);

  clast_->resize(numgp);
  ccurr_->resize(numgp);
  nsfhor_->resize(numgp);
  nsfver_->resize(numgp);
  nsfdiagup_->resize(numgp);
  nsfdiagdown_->resize(numgp);
  LINALG::Matrix<6,1> emptymat6x1(true);
  LINALG::Matrix<3,3> emptymat3x3(true);
  LINALG::Matrix<3,3> emptymatempty(true);
  LINALG::Matrix<nphi,1> emptymat(true);
  for (int i=0; i<3; i++)
    emptymat3x3(i,i) = 1.0;

  for (int i=0; i<numgp; i++)
  {
    histdefgrdlast_->at(i) = emptymat3x3;
    histdefgrdlastsurf_->at(i) = emptymat3x3;
    histdefgrdcurr_->at(i) = emptymat3x3;
    histdefgrdcurrsurf_->at(i) = emptymat3x3;
    strainratelast_->at(i) = emptymat6x1;
    strainratelastsurf_->at(i) = emptymat6x1;
    strainratecurr_->at(i) = emptymat6x1;
    strainratecurrsurf_->at(i) = emptymat6x1;

    etalast_->at(i) = emptymatempty;
    etacurr_->at(i) = emptymatempty;

    etalastphi_->at(i) = emptymat;
    etalastphisurf_->at(i) = emptymat;
    etacurrphi_->at(i) = emptymat;

    Nfil_ ->at(i) = 0.0;
    rate_ ->at(i) = 0.0;
    D_ ->at(i) = 0.0;
    ccurr_ ->at(i) = 0.0;
    clast_ ->at(i) = 0.0;
    nsfhor_->at(i) = 0.0;
    nsfver_->at(i) = 0.0;
    nsfdiagup_->at(i) = 0.0;
    nsfdiagdown_->at(i) = 0.0;
  }

  // Setup of passive material
  matpassive_ = Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Material::Factory(params_->idmatpassive_));
  matpassive_->Setup(numgp, linedef);

  isinit_ = true;
  return;

}   // Setup()


/*----------------------------------------------------------------------*
 |  ResetAll                                                rauch  01/16|
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::ResetAll(const int numgp)
{
  // initialize history variables
  histdefgrdlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  histdefgrdlastsurf_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  histdefgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  histdefgrdcurrsurf_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  strainratelast_= Teuchos::rcp( new std::vector<LINALG::Matrix<6,1> > );
  strainratelastsurf_= Teuchos::rcp( new std::vector<LINALG::Matrix<6,1> > );
  strainratecurr_= Teuchos::rcp( new std::vector<LINALG::Matrix<6,1> > );
  strainratecurrsurf_= Teuchos::rcp( new std::vector<LINALG::Matrix<6,1> > );
  etalast_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );
  etacurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );
  etalastphi_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );
  etacurrphi_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );
  etalastphisurf_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );
  etacurrphisurf_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );

  Nfil_ = Teuchos::rcp( new std::vector<double> );
  rate_ = Teuchos::rcp( new std::vector<double> );
  D_ = Teuchos::rcp( new std::vector<double> );
  ccurr_ = Teuchos::rcp( new std::vector<double> );
  clast_ = Teuchos::rcp( new std::vector<double> );
  nsfhor_=Teuchos::rcp( new std::vector<double> );
  nsfver_=Teuchos::rcp( new std::vector<double> );
  nsfdiagup_=Teuchos::rcp( new std::vector<double> );
  nsfdiagdown_=Teuchos::rcp( new std::vector<double> );

  // set all history variables to zero
  histdefgrdlast_->resize(numgp);
  histdefgrdcurr_->resize(numgp);
  histdefgrdlastsurf_->resize(numgp);
  histdefgrdcurrsurf_->resize(numgp);
  etalast_->resize(numgp);
  etacurr_->resize(numgp);
  etalastphi_->resize(numgp);
  etacurrphi_->resize(numgp);
  etalastphisurf_->resize(numgp);
  etacurrphisurf_->resize(numgp);
  Nfil_->resize(numgp);
  rate_->resize(numgp);
  D_->resize(numgp);
  ccurr_->resize(numgp);
  clast_->resize(numgp);
  nsfhor_->resize(numgp);
  nsfver_->resize(numgp);
  nsfdiagup_->resize(numgp);
  nsfdiagdown_->resize(numgp);

//  etahat_->resize(numgp);
//  etahor_->resize(numgp);
//  etaver_->resize(numgp);
//  etadiag_->resize(numgp);


  LINALG::Matrix<6,1> emptymat6x1(true);
  LINALG::Matrix<3,3> emptymat3x3(true);
  LINALG::Matrix<3,3> emptymatempty(true);
  LINALG::Matrix<nphi,1> emptymat(true);
  for (int i=0; i<3; i++)
    emptymat3x3(i,i) = 1.0;

  for (int i=0; i<numgp; i++)
  {
    histdefgrdlast_->at(i) = emptymat3x3;
    histdefgrdcurr_->at(i) = emptymat3x3;
    histdefgrdlastsurf_->at(i) = emptymat3x3;
    histdefgrdcurrsurf_->at(i) = emptymat3x3;
    strainratelast_->at(i) = emptymat6x1;
    strainratelastsurf_->at(i) = emptymat6x1;
    strainratecurr_->at(i) = emptymat6x1;
    strainratecurrsurf_->at(i) = emptymat6x1;
    etalast_->at(i) = emptymatempty;
    etacurr_->at(i) = emptymatempty;
    etalastphi_->at(i) =  emptymat;
    etacurrphi_->at(i) =  emptymat;
    etalastphisurf_->at(i) =  emptymat;
    etacurrphisurf_->at(i) =  emptymat;
    Nfil_ ->at(i) = 0.0;
    rate_ ->at(i) = 0.0;
    D_ ->at(i) = 0.0;
    ccurr_ ->at(i) = 0.0;
    clast_ ->at(i) = 0.0;
    nsfhor_->at(i) = 0.0;
    nsfver_->at(i) = 0.0;
    nsfdiagup_->at(i) = 0.0;
    nsfdiagdown_->at(i) = 0.0;
  }

  matpassive_->ResetAll(numgp);

  isinit_ = true;
  return;

}   // ResetAll()


/*----------------------------------------------------------------------*
 |  Update internal variables                               rauch  01/16|
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::Update()
{

  // make current values at time step t_last+1 to values of last step t_last
  histdefgrdlast_ = histdefgrdcurr_;

  histdefgrdlastsurf_ = histdefgrdcurrsurf_;

  strainratelast_ =     strainratecurr_;

  strainratelastsurf_ =     strainratecurrsurf_;

  etalast_ = etacurr_;
  etalastphi_ = etacurrphi_;
  etalastphisurf_ = etacurrphisurf_;
  clast_=ccurr_;


  //std::cout<<"test \n"<<std::endl;
  // empty vectors of current data
  histdefgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  histdefgrdcurrsurf_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  strainratecurr_ =     Teuchos::rcp(new std::vector<LINALG::Matrix<6,1> >);
  strainratecurrsurf_ =   Teuchos::rcp(new std::vector<LINALG::Matrix<6,1> >);

  etacurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );
  etacurrphi_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );
  etacurrphisurf_ = Teuchos::rcp( new std::vector<LINALG::Matrix<nphi,1> > );

  ccurr_=Teuchos::rcp( new std::vector<double> );

  // get the size of the vector
  // (use the last vector, because it includes latest results, current is empty)
  const int histsize = histdefgrdlast_->size();

  histdefgrdcurr_->resize(histsize);
  histdefgrdcurrsurf_->resize(histsize);
  strainratecurr_ ->resize(histsize);
  strainratecurrsurf_ ->resize(histsize);
  etacurr_->resize(histsize);
  etacurrphi_->resize(histsize);
  etacurrphisurf_->resize(histsize);
  ccurr_->resize(histsize);
  //sigmaomegaphicurr_->resize(histsize);

  LINALG::Matrix<6,1> emptymat6x1(true);
  LINALG::Matrix<3,3> emptymat3x3(true);
  LINALG::Matrix<3,3> emptymatempty(true);
  LINALG::Matrix<nphi,1> emptymat(true);
  double tmp_scalar = 0.0;
  for (int i=0; i<3; i++)
    emptymat3x3(i,i) = 1.0;
  for (int i=0; i<histsize; i++)
  {
    histdefgrdcurr_->at(i) = emptymat3x3;
    histdefgrdcurrsurf_->at(i) = emptymat3x3;
    strainratecurr_ ->at(i) = emptymat6x1;
    strainratecurrsurf_ ->at(i) = emptymat6x1;
    etacurr_->at(i) = emptymatempty;
    etacurrphi_->at(i) = emptymat;
    etacurrphisurf_->at(i) = emptymat;
    ccurr_->at(i) = tmp_scalar;
  }

  matpassive_->Update();

  return;

}   // Update()


/*----------------------------------------------------------------------*
 |  Reset internal variables                                rauch  01/16|
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::ResetStep()
{
  matpassive_->ResetStep();

}  //ResetStep()


/*----------------------------------------------------------------------*
 |  Evaluate Material                                       rauch  01/16|
 *----------------------------------------------------------------------*
 The stress response is decomposed into a passive and an active part:
     \sigma = \sigma_{passive} + \sigma_{active}
 */
void MAT::BioChemoMechanoCellActiveFiber::Evaluate(const LINALG::Matrix<3,3>* defgrd,
                           const LINALG::Matrix<6,1>* glstrain,
                           Teuchos::ParameterList& params,
                           LINALG::Matrix<6,1>* stress,
                           LINALG::Matrix<6,6>* cmat,
                           const int eleGID)
{
  const Teuchos::ParameterList& strucdynparams =  DRT::Problem::Instance()->CellMigrationParams().sublist("STRUCTURAL DYNAMIC");

  // Set source constant for Scatra Boundary Surface stress dependent calculation
  double sourceconst =   params_->sourceconst_;
  params.set<double>("source const",sourceconst);

  int Scatra = params.get<int>("FromSactraBoundary",0);

  double Theta = -1.0;
  const Teuchos::ParameterList& OneStepThetaParams = strucdynparams.get<Teuchos::ParameterList>("ONESTEPTHETA");

  Theta = OneStepThetaParams.get<double>("THETA");
  if (Theta == -1.0)   dserror("no Theta for time integration provided in material");

  double time = params.get<double>("total time",-1.0);
  if (time == -1.0)   dserror("no time provided in material");


  // Get gauss point number
  const int gp = params.get<int>("gp",-1);
  if (gp == -1)   dserror("no Gauss point number provided in material");
  // Get time algorithmic parameters
  double dt = params.get<double>("delta time",-1.0);
  if (dt == -1.0) dserror("no time step size provided in material");

  bool analyticalmaterialtangent = params_->analyticalmaterialtangent_;
  if (analyticalmaterialtangent)
    dserror("no analytical material tangent calculation possible for this material at the moment.");


//  Teuchos::RCP<const Epetra_MultiVector> testphi = DRT::Problem::Instance()->GetDis("cell")->GetState("phinp");


  DRT::ImmersedFieldExchangeManager* immersedmanager = DRT::ImmersedFieldExchangeManager::Instance();
  Teuchos::RCP<Epetra_MultiVector> testphi = immersedmanager->GetPointerToPhinps();
  if (time>0.0){
  if(testphi == Teuchos::null)
    dserror("testphi = Teuchos::null");
  }

  LINALG::Matrix<1,8> myphinp;
  DRT::Element* scatraele = DRT::Problem::Instance()->GetDis("cellscatra")->gElement(eleGID);
  if(scatraele == NULL)
    dserror("could not get 'cellscatra' element with gid %d",eleGID);

  const Teuchos::RCP<DRT::Discretization> discretization = DRT::Problem::Instance()->GetDis("cellscatra");
  DRT::Element::LocationArray la(1);
  scatraele->LocationVector(*discretization, la,false);
  const std::vector<int>&      lm = la[0].lm_;
  std::vector<double> myphinp2(lm.size());

  const int numofnodes= scatraele->NumNode();
  const int numofphis = scatraele->NumDofPerNode(*scatraele->Nodes()[0]);
  LINALG::Matrix<8,1> cROCKnodes(true);
  LINALG::Matrix<8,1> cActinnodes(true);
  if ( time >0.0){
    DRT::UTILS::ExtractMyValues(*testphi,myphinp2,lm);

    const int numofActin = 8;
    const int numofcrock = 6; // 7-1=6
    for (int node=0; node<numofnodes; ++node)
    {
      const int dof = node*numofphis+numofcrock;
      const int dofActin = node*numofphis+numofActin;
      cROCKnodes(node) = myphinp2[dof];
      cActinnodes(node) = myphinp2[dofActin];
    }
  }


   LINALG::Matrix<3,1> xsi(true);
   LINALG::Matrix<NUMNOD_SOH8,1> shapefcts;


   if (Scatra ==1)
     xsi = params.get<LINALG::Matrix<3,1> >("xsi");
   else
   {
     DRT::UTILS::IntPointsAndWeights<NUMDIM_SOH8> intpoints(DRT::UTILS::intrule_hex_8point);
     for (int i=0;i<3;i++){
       xsi(i)=(intpoints.IP().qxg)[gp][i];
     }
   }

   DRT::UTILS::shape_function<DRT::Element::hex8>(xsi,shapefcts);

   //cRock at GP
   double cROCK=0.0;
   double cActin =0.0;

   for (int i=0;i<numofnodes;i++){
     cROCK+=shapefcts(i)*cROCKnodes(i);
     cActin +=shapefcts(i)*cActinnodes(i);
   }


  //******************
  // PASSIVE PART
  //******************
  // Initialize passive stress and elasticity tensor
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatpassive(true);
  LINALG::Matrix<NUM_STRESS_3D,1> Spassive(true);

  //if (cmat != NULL)
  //{
    // Evaluate passive PK2 stress Spassive and passive elasticity tensor cmatpassive
    matpassive_->Evaluate(defgrd,glstrain,params,&Spassive,&cmatpassive,eleGID);
  //}

  //******************
  // ACTIVE PART
  //******************
  // Parameters for active constitutive model

  // Non-dimensional parameter governing the rate of formation of stress fibers
  double kforwards =   params_->kforwards_;
  // Non-dimensional parameter governing the rate of dissociation of stress fibers
  double kbackwards =  params_->kbackwards_;
  double ratemax =   params_->ratemax_;
  double Nmax = params_->nmax_;
  double kstress = params_->kstress_;

  // Setup inverse of deformation gradient
  LINALG::Matrix<3,3> invdefgrd(*defgrd);
  invdefgrd.Invert();

  // Setup deformation gradient rate, rotation tensor, strain rate and rotation rate
  // \dot{F} = \frac {F^n - F^{n-1}} {\Delta t}
  LINALG::Matrix<3,3> defgrdrate(true);
  // R = F * U^{-1}
  LINALG::Matrix<3,3> R(true);
  // \dot{\epsilon} = d = 0.5 * (\dot{F}F^{-1} + (\dot{F}F^{-1})^{T}
  LINALG::Matrix<6,1> strainrate(true);
  // \dot{R} = \frac {R^n - R^{n-1}} {\Delta t}
  LINALG::Matrix<3,3> rotationrate(true);

  // e = F^{-T} * E * F^{-1}
  LINALG::Matrix<3,3> eastrain(true);

  GLtoEA(*glstrain,invdefgrd,eastrain);


  SetupRates(*defgrd,invdefgrd,params,defgrdrate,R,strainrate,rotationrate,gp,dt);
  LINALG::Matrix<3,3> diagstrainrate(true);
  LINALG::Matrix<3,3> trafo(true);
  LINALG::Matrix<3,3> strainratefull;

  strainratefull(0,0)=strainrate(0);
  strainratefull(1,1)=strainrate(1);
  strainratefull(2,2)=strainrate(2);
  strainratefull(0,1)=strainrate(3);
  strainratefull(1,0)=strainrate(3);
  strainratefull(1,2)=strainrate(4);
  strainratefull(2,1)=strainrate(4);
  strainratefull(0,2)=strainrate(5);
  strainratefull(2,0)=strainrate(5);


  // Evaluate active stress
  //LINALG::Matrix<NUM_STRESS_3D,1> sigma(true);//6x1
  LINALG::Matrix<3,3> cauchystress(true);//3x3

  // Setup history variables

  LINALG::Matrix<3,3> etalast = etalast_->at(gp);



  LINALG::Matrix<nphi,1> etalastphi;
  if (Scatra==1){
    LINALG::Matrix<nphi,1> nodeetaphilast(true);
    EtaAtBrdyGP(etalastphi,params);
  }
  else{
    etalastphi = etalastphi_->at(gp);
  }


  LINALG::Matrix<3,3> etanew(true);
  LINALG::Matrix<3,1> etanewii(true);


  ////////////////////////////////////////////////////////////////////
  // Calculation of the active material    //
  ////////////////////////////////////////////////////////////////////

#ifdef option1
 double fac1 = kforwards*cROCK/Nmax;
 double Dxx=0.;
 double Dyy=0.;
 double Dzz=0.;

  // A x = b ----> x = A^(-1) b
  LINALG::Matrix<3,3> A(true);
  LINALG::Matrix<3,1> b(true);

  A(0,0)= 1.0/dt+ fac1;
  A(1,1)= 1.0/dt+ fac1;
  A(2,2)= 1.0/dt+ fac1;
  A(0,1)= fac1;
  A(0,2)= fac1;
  A(1,0)= fac1;
  A(2,0)= fac1;
  A(1,2)= fac1;
  A(2,1)= fac1;
  Dissociation(strainratefull(0,0),Dxx,kbackwards,ratemax);
  Dissociation(strainratefull(1,1),Dyy,kbackwards,ratemax);
  Dissociation(strainratefull(2,2),Dzz,kbackwards,ratemax);
//  Dxx=0.;
//  Dyy=0.;
//  Dzz=0.;
  b(0)=etalast(0,0)/dt+kforwards*cROCK-Dxx;
  b(1)=etalast(1,1)/dt+kforwards*cROCK-Dyy;
  b(2)=etalast(2,2)/dt+kforwards*cROCK-Dzz;
  // compute inverse of A
  LINALG::Matrix<3,3> Apinv(true);
  double det = A.Determinant();
  if (det < 1.0e-13){
    // Compute Pseudo inverse of A
    LINALG::Matrix<3,3> U(true);
    LINALG::Matrix<3,3> Diag(true);
    LINALG::Matrix<3,3> Vt(true);
    //A = U * Diag * transpose(V)
    LINALG::SVD(A,U,Diag,Vt);

    LINALG::Matrix<3,3> OneByDiag(true);
    for (int i=0; i<3; i++)
    {
    if (Diag(i,i)<1.0e-13)
      OneByDiag(i,i)=0.;
    else
      OneByDiag(i,i)=1.0/Diag(i,i);
    }
    LINALG::Matrix<3,3> temp(true);
    temp.MultiplyNT(OneByDiag,U);
    Apinv.MultiplyTN(Vt,temp);
  }
  else{
  Apinv.Invert(A);
  }
  etanewii.MultiplyNN(Apinv,b);
  //  Check result:
  LINALG::Matrix<3,1> tempvector(true);
  tempvector.Clear();
  LINALG::Matrix<3,1> res(true);
  res.MultiplyNN(A,etanewii);
  res.Update(1.0,b,-1.0);
  double absres = abs(res(0))+abs(res(1))+abs(res(2));
  if (absres>1.0e-12)
      {
    std::cout<<"res "<<res<<std::endl;
    dserror("calculation of eta_(ii)^{n+1} went wrong");
      }
//  std::cout<<"A "<<A<<std::endl;
//  std::cout<<"U "<<U<<std::endl;
//  std::cout<<"Vt "<<Vt<<std::endl;
//  std::cout<<"Sigma "<<Sigma<<std::endl;
//  std::cout<<"OneBySigma "<<OneBySigma<<std::endl;
//  std::cout<<"Apinv "<<Apinv<<std::endl;
//  std::cout<<"b "<<b<<std::endl;
//  std::cout<<"etanewii "<<etanewii<<std::endl;
//    dserror("blub");
  double time = params.get<double>("total time",-1.0);
//  if (time<=dt)
//  {
//    for (int i=0; i<3; i++)
//    {
//     etanewii(i)=0.0;
//    }
//  }
  double Nfil=0.;
  for (int i=0; i<3; i++)
  {
   etanew(i,i)=etanewii(i);
   Nfil += etanewii(i);
   if (etanewii(i)<-1.0e-13)
   {
     std::cout<<"Dxx "<<Dxx<<std::endl;
     std::cout<<"Dyy "<<Dyy<<std::endl;
     std::cout<<"Dzz "<<Dzz<<std::endl;
     std::cout<<"vector eta_{ii} "<<etanewii<<std::endl;
     //dserror("eta_{ii}<0");
   }
  }
//  LINALG::Matrix<3,3> temp(true);
//  temp.Clear();
//  temp.MultiplyNT(etanew,trafo);
//  etanew.MultiplyNN(trafo,temp);

  double Dp;
  double Dn;

  for (int i=0; i<3;i++)
  {
    for (int j=i+1; j<3; j++)
    {
      //if (i != j)
      {
        Dp=0.;
        Dissociation(strainratefull(i,j),Dp,kbackwards,ratemax);
        Dn=0.;
        Dissociation((-1.0)*strainratefull(i,j),Dn,kbackwards,ratemax);
        etanew(i,j) = dt*(kforwards*cROCK*(1-Nfil/Nmax)) + etalast(i,j) - dt*0.5*(Dp+Dp);
        //etanew(i,j) = (-1.0)*etanew(i,j) ;
        //etanew(i,j)=0.;
        etanew(j,i) = etanew(i,j);
        //printf("(%d,%d)", i, j);
      }
    }
  }
#endif

#ifdef option2
  double fac1 = kforwards*cROCK/Nmax;
  LINALG::Matrix<6,6> A(true);
  LINALG::Matrix<6,1> b(true);
  LINALG::Matrix<6,1> c(true);
  LINALG::Matrix<6,1> etanewvec(true);
  LINALG::Matrix<6,1> etalastvec(true);


  etalastvec(0) = etalast(0,0);
  etalastvec(1) = etalast(1,1);
  etalastvec(2) = etalast(2,2);
  etalastvec(3) = etalast(0,1);
  etalastvec(4) = etalast(1,2);
  etalastvec(5) = etalast(0,2);

  double D;
  for (int i=0;i<6;i++){
    A(i,i)= 1/dt+fac1;
    for (int j=i+1;j<6;j++){
        A(i,j)=fac1;
        A(j,i)=fac1;
    }
    D=0.;
    Dissociation(strainrate(i),D,kbackwards,ratemax);
    //if (i>2){
      b(i)=etalastvec(i)/dt+(kforwards*cROCK-D);
    //}
    //else{
    //  b(i)=etalastvec(i)+(kforwards*cROCK-D);
    //}
  }
  //LINALG::Matrix<6,6> Apinv(true);
  //Apinv.Invert(A);
  LINALG::Matrix<6,6> Ainv(true);
  Ainv.Update(A);
  LINALG::FixedSizeSerialDenseSolver<6,6> Ainversion;
  Ainversion.SetMatrix(Ainv);
  int err = Ainversion.Invert();
  if(err != 0){
    std::cout<<"A "<<A<<std::endl;
    dserror("Inversion of 6x6 matrix failed with errorcode %d",err);
  }



  etanewvec.MultiplyNN(Ainv,b);
  //  Check result:
  LINALG::Matrix<6,1> tempvector(true);
  //tempvector.Clear();
  LINALG::Matrix<6,1> res(true);
  res.MultiplyNN(A,etanewvec);
  res.Update(1.0,b,-1.0);
  double absres = abs(res(0))+abs(res(1))+abs(res(2))+abs(res(3))+abs(res(4))+abs(res(5));
//  if (absres>(1.0e-12)/3)
//      {
//    std::cout<<"A "<<A<<std::endl;
//    std::cout<<"Ainv "<<Ainv<<std::endl;
//    std::cout<<"b "<<b<<std::endl;
//    std::cout<<"res "<<res<<std::endl;
//    dserror("calculation of eta_(ij)^(n+1) went wrong");
//      }
  etanew(0,0)=etanewvec(0);
  etanew(1,1)=etanewvec(1);
  etanew(2,2)=etanewvec(2);
  etanew(0,1)=etanewvec(3);
  etanew(1,0)=etanewvec(3);
  etanew(1,2)=etanewvec(4);
  etanew(2,1)=etanewvec(4);
  etanew(0,2)=etanewvec(5);
  etanew(2,0)=etanewvec(5);

  double Nfil=0;
  for (int i=0;i<6;i++){
    if (i<3){
      Nfil+=etanewvec(i);
    }
    Nfil += 2*etanewvec(i);
  }
  //Testing
//  sigma(0)=strainrate(0);
//  sigma(1)=strainrate(1);
//  sigma(2)=strainrate(2);
//  sigma(3)=strainrate(3);
//  sigma(4)=strainrate(4);
//  sigma(5)=strainrate(5);

//  if (time<=dt){
//      sigma(0)=0.2;
//      sigma(1)=0.1;
//      sigma(2)=0.0;
//      sigma(3)=0.0;
//      sigma(4)=0.0;
//      sigma(5)=0.0;
//
//  }

  double etadiff= (abs(etanew(0,1)-etanew(1,0))+abs(etanew(0,2)-etanew(2,0))+abs(etanew(1,2)-etanew(2,1)));
      if (etadiff>1.0e-12)
        dserror("eta not diagonal");
  // Convert stress like 6x1-Voigt vector to 3x3 matrix

      LINALG::Matrix<6,1> singstrain(true);
      for (int i=0;i<6;i++){
        if (strainrate(i)<0.0){
          singstrain(i)=-1.0;
        }
        else
          singstrain(i)=1.0;
      }
#endif

// Newton Raphson loop (inkl. TraceBack) for finding cROCKused
#ifdef option5
      LINALG::Matrix<nphi,nphi> J(true);
      LINALG::Matrix<nphi,1> f(true);
      LINALG::Matrix<nphi,1> b(true);
      LINALG::Matrix<nphi,1> etanewphi(true);
      etanewphi.Update(1.0,etalastphi);
      LINALG::Matrix<nphi,1> etanewphioldguess(true);
      LINALG::Matrix<nphi,1> deltaetaphi(true);
      LINALG::Matrix<nphi,1> strainratephi(true);
      double rateaverage=0.0;
      double cROCKusednew=0.0;
      double cROCKusedold=0.0;
      if (Scatra==1){
          cROCKusedAtBrdyGP(cROCKusedold, params);
//            if (time>3.0){
//                std::cout<<"cROCKusedold "<<cROCKusedold<<std::endl;
//            for (int gp =0;gp<8;gp++){
//              cROCKusedold = clast_->at(gp);
//               std::cout<<"cROCKusedold "<<cROCKusedold<<std::endl;
//            }
//            dserror("stop");}
      }
      else{
        cROCKusedold = clast_->at(gp);
      }
      double deltaphi = 2.0*M_PI/((double)(nphi));

      double kROCKeta = params_->kRockEta_;
      double Nfilold=0.0;
      for (int i=0; i<nphi; i++){
        Nfilold+=etalastphi(i) * deltaphi;
      }


      double Nfil=0;
      //double kforwards =   params_->kforwards_;
      double  Dcouple = params_->DissCoupling_;
      double Daverage=0.0;

      double phi;
      LINALG::Matrix<nphi,1> Diss(true);
      LINALG::Matrix<2,1> m(true);
      for (int i=0; i<nphi; i++)
      {
          phi  = deltaphi*(double)(i);
          m(0) = cos(phi);
          m(1) = sin(phi);

      // Transform strain rate at each point to fiber strain rate in (phi) direction
      strainratephi(i) =              strainrate(0) * m(0)*m(0)
                                 +    strainrate(1) * m(1)*m(1)
                                 + 2.*strainrate(3) * m(0)*m(1);

      Dissociation(strainratephi(i),Diss(i),kbackwards,ratemax);
      Daverage+=Diss(i)*deltaphi;
      }
      Daverage=Daverage/(2.0*M_PI);
      int iter=0;

      Nfil=0.0;
      for (int i=0; i<nphi; i++){
        Nfil+=etanewphi(i) * deltaphi;
      }

      for (int i =0; i<nphi;i++){
      f(i)=(-1.0)*(etanewphi(i)-etalastphi(i))/dt
          + kforwards*(cROCK - (cROCKusedold + 1.0 / kROCKeta *(Nfil-Nfilold)+ Dcouple * Daverage))*(1-Nfil/Nmax)-Diss(i);
      }
//      if (time>0.0){
//        std::cout<<"f "<< f <<std::endl;
//      }
   double itertraceback=1.0;
   double absvalue=f.Norm2();
   while(absvalue>1.0e-9){
      iter+=1;
      //if (abs(etanewphi.MaxValue()-etanewphi.MinValue())>1.0e-10){
      double fac = kforwards*(-1.0 / kROCKeta *deltaphi )*(1-Nfil/Nmax)
                    + kforwards*(cROCK - (cROCKusedold + 1/kROCKeta*(Nfil-Nfilold) + Dcouple * Daverage) )*(-1.0)*deltaphi/Nmax;
      for (int i=0; i<nphi; i++){
        for (int j=i;j<nphi;j++){
            J(i,j)= fac;
            J(j,i)= fac;
        }
        J(i,i)=-1.0/dt + fac;
      }// end for loop over phi (i.e. phi_i)


      LINALG::Matrix<nphi,nphi> Jinv(true);
      Jinv.Update(J);
      LINALG::FixedSizeSerialDenseSolver<nphi,nphi> Jinversion;
      Jinversion.SetMatrix(Jinv);
      int err = Jinversion.Invert();
      if(err != 0){
//        LINALG::Matrix<nphi,nphi> U(true);
//        LINALG::Matrix<nphi,nphi> Diag(true);
//        LINALG::Matrix<nphi,nphi> Vt(true);
//        //A = U * Diag * transpose(V)
//        LINALG::SVD(J,U,Diag,Vt);
//
//        LINALG::Matrix<nphi,nphi> OneByDiag(true);
//        for (int i=0; i<nphi; i++)
//        {
//        if (Diag(i,i)<1.0e-13)
//          OneByDiag(i,i)=0.;
//        else
//          OneByDiag(i,i)=1.0/Diag(i,i);
//        }
//        LINALG::Matrix<nphi,nphi> temp(true);
//        temp.MultiplyNT(OneByDiag,U);
//        Jinv.MultiplyTN(Vt,temp);

        std::cout<<"J "<< J <<std::endl;
//        std::cout<<"absvalue "<< absvalue <<std::endl;
//        std::cout<<"f "<< f <<std::endl;
        dserror("Inversion of (nphi)x(nphi) matrix failed with errorcode %d",err);
      }
      deltaetaphi.MultiplyNN(Jinv,f);
      deltaetaphi.Scale(-1.0);
//      }
//      else{
//        std::cout<<"in scalar" <<std::endl;
//        double Jsc=0;
//        Jsc=-etanewphi(0)/dt + kforwards*(-1.0 / kROCKeta *2.0*M_PI )*(1-Nfil/Nmax)
//            + kforwards*(cROCK - (cROCKusedold + 1/kROCKeta*(Nfil-Nfilold)))*(-1.0)*2.0*M_PI/Nmax;
//        for (int i=0; i<nphi; i++){
//          deltaetaphi(i)=-1.0/Jsc * f(i);
//        }
//      }
      etanewphioldguess(true);
      etanewphioldguess.Update(1.0,etanewphi);
//      if (time>0.0){
//        std::cout<<"etanewphioldguess"<< etanewphioldguess <<std::endl;
//        std::cout<<"etanewphi"<< etanewphi <<std::endl;
//      }
      etanewphi.Update(1.0,deltaetaphi,1.0);
      Nfil=0.0;
      for (int i=0; i<nphi; i++){
        Nfil+=etanewphi(i) * deltaphi;
      }
      for (int i=0; i<nphi; i++){
        f(i)=(-1.0)*(etanewphi(i)-etalastphi(i))/dt
            + kforwards*(cROCK - (cROCKusedold + 1.0 / kROCKeta *(Nfil-Nfilold)+ Dcouple * Daverage))*(1-Nfil/Nmax)-Diss(i);
      }

      itertraceback=1.0;
      while(absvalue < f.Norm2()) {
        etanewphi.Update(1.0,etanewphioldguess,0.0);

        itertraceback=itertraceback/2.0;
        etanewphi.Update(itertraceback,deltaetaphi,1.0);
        if (time>0.0){
          //std::cout<<"f in traceback"<< f <<std::endl;
          //std::cout<<"etanewphi"<< etanewphi <<std::endl;
        }
        Nfil=0.0;
        for (int i=0; i<nphi; i++){
          Nfil+=etanewphi(i) * deltaphi;
        }
        for (int i=0; i<nphi; i++){
          f(i)=(-1.0)*(etanewphi(i)-etalastphi(i))/dt
              + kforwards*(cROCK - (cROCKusedold + 1.0 / kROCKeta *(Nfil-Nfilold) + Dcouple * Daverage ))*(1-Nfil/Nmax)-Diss(i);
        }
        if (itertraceback<1.0e-4){
          std::cout<<"deltaetaphi "<< deltaetaphi <<std::endl;
          std::cout<<"f "<< f <<std::endl;
          dserror("local loop trace back did not find smaller abs(f)");
        }
      }
      absvalue = f.Norm2();

      if (iter>12){
        std::cout<<"f "<< f <<std::endl;
        std::cout<<"etanewphi "<< etanewphi <<std::endl;
        dserror("local loop did not converge");
      }
//      if (time>0.0){
//        std::cout<<"f "<< f <<std::endl;
//        //std::cout<<"abs(etanewphi.MaxValue()-etanewphi.MinValue()) "<< abs(etanewphi.MaxValue()-etanewphi.MinValue()) <<std::endl;
//        //std::cout<<"etanewphi "<< etanewphi <<std::endl;
//        //std::cout<<"Diss "<< Diss <<std::endl;
//        //std::cout<<"strainratephi "<< strainratephi<<std::endl;
//      }
      }


   cROCKusednew = cROCKusedold + 1/kROCKeta * (Nfil -Nfilold) + Dcouple * Daverage;

   for (int i=0; i<nphi; i++){
     phi  = deltaphi*(double)(i);
     m(0) = cos(phi);
     m(1) = sin(phi);
     etanew(0,0) += etanewphi(i) * deltaphi * m(0)*m(0); // eta_11
     etanew(1,1) += etanewphi(i) * deltaphi * m(1)*m(1); // eta_22
     etanew(0,1) += etanewphi(i) * deltaphi * m(0)*m(1); // eta_12 = eta_21
   }
   etanew(1,0)= etanew(0,1);
   etanew.Scale(1.0/M_PI);
   if (time>0.0){
    // std::cout<<"iter "<< iter<<std::endl;
   //dserror("stop");
   }
      rateaverage=kforwards*(cROCK-cROCKusednew)*(1-Nfil/Nmax)-Daverage;

#endif

#ifdef option3
  const std::string MyIntMethod =   params_->myintmethod_;

  LINALG::Matrix<nphi,nphi> A(true);
  LINALG::Matrix<nphi,1> b(true);
  LINALG::Matrix<nphi,1> etanewphi(true);
  LINALG::Matrix<nphi,1> strainratephi(true);
  LINALG::Matrix<nphi,1> strainratephilast(true);
  LINALG::Matrix<6,1> strainratelast(true); // For Theta \neq 1


  double Nfil=0;
  double Nfilold=0;
  double Daverage=0.0;
  double phi;
  double newfac1;
  double D;
  double Dlast=0.0;
  double cROCKlast;
  double contractionrateaverage=0;
  //double Theta=1.0;
  double thetafac = 0.0;
  if (Scatra==1){
    cROCKlastAtBrdyGP(cROCKlast, params);
  }
  else{
    cROCKlast = clast_->at(gp);
  }


  if (MyIntMethod == "2DTrapez"){
    double deltaphi = 2.0*M_PI/((double)(nphi));
    Nfilold=0.0;
    for (int i=0; i<nphi; i++){
      Nfilold+=etalastphi(i) * deltaphi;
    }
    if (Theta!=1.0){
      strainratelast = strainratelast_->at(gp);
      thetafac = kforwards*cROCKlast*(1-Nfilold/(0.5*M_PI * Nmax));
    }
    newfac1= kforwards*cActin*cROCK*deltaphi/(0.5*M_PI*Nmax);

    LINALG::Matrix<2,1> m(true);
    for (int i=0; i<nphi; i++){
      phi  = deltaphi*(double)(i);
      m(0) = cos(phi);
      m(1) = sin(phi);

      // Transform strain rate at each point to fiber strain rate in (phi) direction
      strainratephi(i) =       strainrate(0) * m(0)*m(0)
                          +    strainrate(1) * m(1)*m(1)
                          + 2.*strainrate(3) * m(0)*m(1);
      if (Theta != 1.0){
        strainratephilast(i) =       strainratelast(0) * m(0)*m(0)
                                +    strainratelast(1) * m(1)*m(1)
                                + 2.*strainratelast(3) * m(0)*m(1);
      }


      for (int j=i+1;j<nphi;j++){
        A(i,j)=Theta*newfac1;
        A(j,i)=Theta*newfac1;
      }
      D=0.;
      Dissociation(strainratephi(i),D,kbackwards,ratemax);
      A(i,i)= 1/dt+Theta*(newfac1+D);
      if (Theta != 1.0){
        Dlast = 0.0;
        Dissociation(strainratephilast(i),Dlast,kbackwards,ratemax);
      }

      //b(i)=etalastphi(i)/dt+Theta*(kforwards*(cROCK)-D) + (1.0-Theta)*(thetafac-Dlast);
      b(i)=etalastphi(i)/dt+Theta*(kforwards*(cROCK)*cActin) + (1.0-Theta)*(thetafac-Dlast);
      contractionrateaverage +=strainratephi(i)*deltaphi;
      //Daverage+=D*deltaphi;
    }// end for loop over phi (i.e. phi_i)
    contractionrateaverage=contractionrateaverage/(2.0*M_PI);

    //Daverage=Daverage/(2.0*M_PI);

    LINALG::Matrix<nphi,1> etanewphisolver(true);
    // solve A.X=B
    LINALG::FixedSizeSerialDenseSolver<nphi,nphi,1> solver;
    solver.SetMatrix(A);              // set A
    solver.SetVectors(etanewphisolver, b);           // set X=etanewphi, B=b
    solver.FactorWithEquilibration(true); //
    //int err2 = solver.Factor();           //
    int err3 = solver.Solve();             // X = A^-1 B
    //if ((err3!=0) || (err2!=0))
    if (err3!=0)
      dserror("solving linear system in cell contraction for new eta failed");

    etanewphi.Update(etanewphisolver);


//    #ifdef DEBUG
//      //  Check result:
//      LINALG::Matrix<nphi,1> tempvector(true);
//      LINALG::Matrix<nphi,1> res(true);
//      res.MultiplyNN(A,etanewphi);
//      res.Update(1.0,b,-1.0);
//      double absres = abs(res(0))+abs(res(1))+abs(res(2))+abs(res(3))+abs(res(4))+abs(res(5));
//      if (absres>(1.0e-12)/3){
//        std::cout<<"A "<<A<<std::endl;
//        //std::cout<<"Ainv "<<Ainv<<std::endl;
//        std::cout<<"b "<<b<<std::endl;
//        std::cout<<"res "<<res<<std::endl;
//        dserror("calculation of eta_(ij)^(n+1) went wrong");
//      }
//    #endif

    Nfil=0.0;
    Daverage=0.0;
    for (int i=0; i<nphi; i++){
      phi  = deltaphi*(double)(i);
      m(0) = cos(phi);
      m(1) = sin(phi);
      etanew(0,0) += etanewphi(i) * deltaphi * m(0)*m(0);
      etanew(1,1) += etanewphi(i) * deltaphi * m(1)*m(1);
      etanew(0,1) += etanewphi(i) * deltaphi * m(0)*m(1);
      Nfil+=etanewphi(i) * deltaphi;
      D=0.;
      Dissociation(strainratephi(i),D,kbackwards,ratemax);
      Daverage+=D*deltaphi*etanewphi(i);
    }
    Daverage=Daverage/(2.0*M_PI);
    etanew.Scale(1.0/M_PI); // Vol Average Factor for 2D
    etanew(1,0)= etanew(0,1);


    }// End Trapez Sum
  //2D GAUSSIANINTEGRAL:
  else if(MyIntMethod=="2DGauss"){
    // Setting up gauss quadrature
    const DRT::UTILS::IntegrationPoints1D gausspoints(DRT::UTILS::intrule_line_20point);
    if (nphi==20)
      const DRT::UTILS::IntegrationPoints1D gausspoints(DRT::UTILS::intrule_line_20point);
    else if (nphi==16)
      const DRT::UTILS::IntegrationPoints1D gausspoints(DRT::UTILS::intrule_line_16point);
    else if (nphi==10)
      const DRT::UTILS::IntegrationPoints1D gausspoints(DRT::UTILS::intrule_line_10point);
    else if (nphi==9)
      const DRT::UTILS::IntegrationPoints1D gausspoints(DRT::UTILS::intrule_line_9point);
    else if (nphi==8)
      const DRT::UTILS::IntegrationPoints1D gausspoints(DRT::UTILS::intrule_line_8point);
    else if (nphi==7)
      const DRT::UTILS::IntegrationPoints1D gausspoints(DRT::UTILS::intrule_line_7point);
    else if (nphi==6)
      const DRT::UTILS::IntegrationPoints1D gausspoints(DRT::UTILS::intrule_line_6point);
    else if (nphi==5)
      const DRT::UTILS::IntegrationPoints1D gausspoints(DRT::UTILS::intrule_line_5point);
    else {
      std::cout<<"For 2D Gauss integration the following numbers of integration points are available: 20, 16, 10, 9, 8, 7 ,6, 5"<<std::endl;
      dserror("Fix number of integration points (nphi) in cellcontraction.h");
    }

    Nfilold=0.0;
    for (int i=0; i<nphi; i++){
      Nfilold+=M_PI*etalastphi(i) * gausspoints.qwgt[i];
    }

    if (Theta!=1.0){
      strainratelast = strainratelast_->at(gp);
      thetafac = kforwards*cROCKlast*(1-Nfilold/(M_PI*0.5 * Nmax));
    }
    newfac1= kforwards*cROCK/(0.5*Nmax);
    LINALG::Matrix<2,1> m(true);
    for (int i=0; i<nphi; i++){
      // Setting up gauss quadrature
      phi  = M_PI * gausspoints.qxg[i][0] + M_PI;
      m(0) = cos(phi);
      m(1) = sin(phi);

      // Transform strain rate at each point to fiber strain rate in (phi) direction
      strainratephi(i) =       strainrate(0) * m(0)*m(0)
                          +    strainrate(1) * m(1)*m(1)
                          + 2.*strainrate(3) * m(0)*m(1);
      if (Theta != 1.0){
        strainratephilast(i) =       strainratelast(0) * m(0)*m(0)
                                +    strainratelast(1) * m(1)*m(1)
                                + 2.*strainratelast(3) * m(0)*m(1);
      }


      for (int j=0;j<nphi;j++){
        A(i,j)=Theta*newfac1*gausspoints.qwgt[j];
        //A(j,i)=Theta*newfac1*gausspoints.qwgt[j];
      }
      D=0.;
      Dissociation(strainratephi(i),D,kbackwards,ratemax);
      //A(i,i)(i,i)= 1/dt+Theta*gausspoints.qwgt[i]*newfac1
      A(i,i)= 1/dt+Theta*gausspoints.qwgt[i]*newfac1 + D;

      if (Theta != 1.0){
        Dlast = 0.0;
        Dissociation(strainratephilast(i),Dlast,kbackwards,ratemax);
      }

      //b(i)=etalastphi(i)/dt+Theta*(kforwards*(cROCK)-D) + (1.0-Theta)*(thetafac-Dlast);
      b(i)=etalastphi(i)/dt+Theta*(kforwards*(cROCK)) + (1.0-Theta)*(thetafac-Dlast);
      Daverage+=D*gausspoints.qwgt[i]*M_PI;
    }// end for loop over phi (i.e. phi_i)

    Daverage=Daverage/(2.0*M_PI);

    LINALG::Matrix<nphi,1> etanewphisolver(true);
    // solve A.X=B
    LINALG::FixedSizeSerialDenseSolver<nphi,nphi,1> solver;
    solver.SetMatrix(A);              // set A
    solver.SetVectors(etanewphisolver, b);           // set X=etanewphi, B=b
    solver.FactorWithEquilibration(true); //
    //int err2 = solver.Factor();           //
    int err3 = solver.Solve();             // X = A^-1 B
    //if ((err3!=0) || (err2!=0))
    if (err3!=0)
      dserror("solving linear system in cell contraction for new eta failed");

    etanewphi.Update(etanewphisolver);


//    #ifdef DEBUG
//      //  Check result:
//      LINALG::Matrix<nphi,1> tempvector(true);
//      LINALG::Matrix<nphi,1> res(true);
//      res.MultiplyNN(A,etanewphi);
//      res.Update(1.0,b,-1.0);
//      double absres = abs(res(0))+abs(res(1))+abs(res(2))+abs(res(3))+abs(res(4))+abs(res(5));
//      if (absres>(1.0e-12)/3){
//        std::cout<<"A "<<A<<std::endl;
//        std::cout<<"Ainv "<<Ainv<<std::endl;
//        std::cout<<"b "<<b<<std::endl;
//        std::cout<<"res "<<res<<std::endl;
//        dserror("calculation of eta_(ij)^(n+1) went wrong");
//      }
//    #endif

    for (int i=0; i<nphi; i++){
      phi  = M_PI * gausspoints.qxg[i][0] + M_PI;
      m(0) = cos(phi);
      m(1) = sin(phi);
      etanew(0,0) += M_PI * etanewphi(i) * gausspoints.qwgt[i] * m(0)*m(0);
      etanew(1,1) += M_PI * etanewphi(i) * gausspoints.qwgt[i] * m(1)*m(1);
      etanew(0,1) += M_PI * etanewphi(i) * gausspoints.qwgt[i] * m(0)*m(1);
      Nfil+=M_PI*etanewphi(i) * gausspoints.qwgt[i];
    }
    etanew.Scale(1.0/M_PI); // Vol Average Factor for 2D
    etanew(1,0)= etanew(0,1);
    // Calculate amount of stress fibers in certain direction

  }// End Gauss Quadratur 2D
  //#endif
  else if (MyIntMethod=="3DGauss"){
    // Choose 1D integration rule for each dimension
    // N.B. nphi is the total number of integration points in 3DGauss routine.
    // nphi = n_{azimuthal} * n_{polar}
    int nphiphi = nphi/ntheta; // azimuthal angle phi

    // Integration points in phi direction
    const DRT::UTILS::IntegrationPoints1D gausspointsphi(DRT::UTILS::intrule_line_20point);
    if (nphiphi==20)
      const DRT::UTILS::IntegrationPoints1D gausspointsphi(DRT::UTILS::intrule_line_20point);
    else if (nphiphi==16)
      const DRT::UTILS::IntegrationPoints1D gausspointsphi(DRT::UTILS::intrule_line_16point);
    else if (nphiphi==10)
      const DRT::UTILS::IntegrationPoints1D gausspointsphi(DRT::UTILS::intrule_line_10point);
    else if (nphiphi==9)
      const DRT::UTILS::IntegrationPoints1D gausspointsphi(DRT::UTILS::intrule_line_9point);
    else if (nphiphi==8)
      const DRT::UTILS::IntegrationPoints1D gausspointsphi(DRT::UTILS::intrule_line_8point);
    else if (nphiphi==7)
      const DRT::UTILS::IntegrationPoints1D gausspointsphi(DRT::UTILS::intrule_line_7point);
    else if (nphiphi==6)
      const DRT::UTILS::IntegrationPoints1D gausspointsphi(DRT::UTILS::intrule_line_6point);
    else if (nphiphi==5)
      const DRT::UTILS::IntegrationPoints1D gausspointsphi(DRT::UTILS::intrule_line_5point);
    else if (nphiphi==4)
      const DRT::UTILS::IntegrationPoints1D gausspointsphi(DRT::UTILS::intrule_line_4point);
    else {
      std::cout<<"For 3D Gauss integration the following numbers of integration points are available: 20, 16, 10, 9, 8, 7 ,6, 5, 4"<<std::endl;
      dserror("Fix number of total integration points ('nphi=nphi*ntheata') in cellcontraction.h");
    }

    // Integration points in theta direction
    const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_20point);
    if (ntheta==20)
      const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_20point);
    else if (ntheta==16)
      const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_16point);
    else if (ntheta==10)
      const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_10point);
    else if (ntheta==9)
      const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_9point);
    else if (ntheta==8)
      const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_8point);
    else if (ntheta==7)
      const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_7point);
    else if (ntheta==6)
      const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_6point);
    else if (ntheta==5)
      const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_5point);
    else if (ntheta==4)
      const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_4point);
    else if (ntheta==3)
      const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_3point);
    else if (ntheta==2)
      const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_2point);
    else if (ntheta==1)
      const DRT::UTILS::IntegrationPoints1D gausspointstheta(DRT::UTILS::intrule_line_1point);
    else {
      std::cout<<"For 3D Gauss integration over theta the following numbers of integration points are available:"<<std::endl;
      std::cout<<"20, 16, 10, 9, 8, 7 ,6, 5, 4, 3, 2, 1"<<std::endl;
      dserror("Fix number of total integration points ('nphi=nphi*ntheata') in cellcontraction.h");
    }
    double theta;
    double thetatemp;
    Nfilold=0.0;
    for (int i=0; i<nphiphi; i++){
      for (int j = 0; j<ntheta; j++){
        theta  = M_PI /2.0 * gausspointstheta.qxg[j][0] + M_PI/2.0;
        Nfilold+=etalastphi(j+i*ntheta) *sin(theta)* gausspointsphi.qwgt[i]* gausspointstheta.qwgt[j];
      }
    }
    Nfilold = M_PI*M_PI/2.0*Nfilold;

    if (Theta!=1.0){
      strainratelast = strainratelast_->at(gp);
      thetafac = kforwards*cROCKlast*(1-Nfilold/(M_PI*2.0/3.0 * Nmax));
    }
    newfac1= 3.0/4.0*M_PI* kforwards*cROCK/(Nmax);
    LINALG::Matrix<3,1> m(true);
    for (int i=0; i<nphiphi; i++){
      for (int j = 0; j<ntheta; j++){
        LINALG::Matrix<3,1> m(true);
        theta  = M_PI /2.0 * gausspointstheta.qxg[j][0] + M_PI/2.0;
        phi  = M_PI * gausspointstheta.qxg[i][0] + M_PI;
        m(0) = sin(theta)*cos(phi);
        m(1) = sin(theta)*sin(phi);
        m(2) = cos(theta);
        // Transform strain rate at each point to fiber strain rate in (omega,phi) direction
        // \dot{\epsilon} = \dot{\epsilon}_{ij} m_{i} m_{j}
        strainrate(j+i*ntheta) =        strainrate(0) * m(0)*m(0)
                                   +    strainrate(1) * m(1)*m(1)
                                   +    strainrate(2) * m(2)*m(2)
                                   + 2.*strainrate(3) * m(0)*m(1)
                                   + 2.*strainrate(4) * m(1)*m(2)
                                   + 2.*strainrate(5) * m(0)*m(2);

        // Setting up gauss quadrature


      // Transform strain rate at each point to fiber strain rate in (phi) direction

      if (Theta != 1.0){
        strainratephilast(j+i*ntheta) =      strainratelast(0) * m(0)*m(0)
                                        +    strainratelast(1) * m(1)*m(1)
                                        +    strainratelast(2) * m(2)*m(2)
                                        + 2.*strainratelast(3) * m(0)*m(1)
                                        + 2.*strainratelast(4) * m(1)*m(2)
                                        + 2.*strainratelast(5) * m(0)*m(2);
      }

      //
      for (int k=0; k<nphiphi; k++){
        for (int l = 0; l<ntheta; l++){
          thetatemp  = M_PI /2.0 * gausspointstheta.qxg[l][0] + M_PI/2.0;
          A(j+i*ntheta,l+k*ntheta)=Theta*newfac1*gausspointsphi.qwgt[k]*gausspointstheta.qwgt[l]*sin(thetatemp);
        }
      }
      A(j+i*ntheta,j+i*ntheta)= 1/dt+Theta*gausspointsphi.qwgt[i]*gausspointstheta.qwgt[j]*newfac1;
      D=0.;
      Dissociation(strainratephi(j+i*ntheta),D,kbackwards,ratemax);
      if (Theta != 1.0){
        Dlast = 0.0;
        Dissociation(strainratephilast(j+i*ntheta),Dlast,kbackwards,ratemax);
      }

      b(j+i*ntheta)=etalastphi(j+i*ntheta)/dt+Theta*(kforwards*(cROCK)-D) + (1.0-Theta)*(thetafac-Dlast);

      Daverage+=D*gausspointsphi.qwgt[i]*gausspointstheta.qwgt[j]*sin(theta);
    }// end for loop over phi (i.e. phi_i)
    }
    Daverage=Daverage*M_PI/(8.0);


    LINALG::Matrix<nphi,1> etanewphisolver(true);
    // solve A.X=B
    LINALG::FixedSizeSerialDenseSolver<nphi,nphi,1> solver;
    solver.SetMatrix(A);              // set A
    solver.SetVectors(etanewphisolver, b);           // set X=etanewphi, B=b
    solver.FactorWithEquilibration(true); //
    //int err2 = solver.Factor();           //
    int err3 = solver.Solve();             // X = A^-1 B
    //if ((err3!=0) || (err2!=0))
    if (err3!=0)
      dserror("solving linear system in cell contraction for new eta failed");

    etanewphi.Update(etanewphisolver);

    for (int i=0; i<nphiphi; i++){
      for (int j = 0; j<ntheta; j++){
        theta  = M_PI /2.0 * gausspointstheta.qxg[j][0] + M_PI/2.0;
        phi  = M_PI * gausspointstheta.qxg[i][0] + M_PI;
        m(0) = sin(theta)*cos(phi);
        m(1) = sin(theta)*sin(phi);
        m(2) = cos(theta);
        etanew(0,0) += M_PI * M_PI / 2.0 * etanewphi(j+i*ntheta) * gausspointsphi.qwgt[i]* gausspointstheta.qwgt[j] * sin(theta) * m(0)*m(0);
        etanew(1,1) += M_PI * M_PI / 2.0 * etanewphi(j+i*ntheta) * gausspointsphi.qwgt[i]* gausspointstheta.qwgt[j] * sin(theta) * m(1)*m(1);
        etanew(2,2) += M_PI * M_PI / 2.0 * etanewphi(j+i*ntheta) * gausspointsphi.qwgt[i]* gausspointstheta.qwgt[j] * sin(theta) * m(2)*m(2);
        etanew(0,1) += M_PI * M_PI / 2.0 * etanewphi(j+i*ntheta) * gausspointsphi.qwgt[i]* gausspointstheta.qwgt[j] * sin(theta) * m(0)*m(1);
        etanew(1,2) += M_PI * M_PI / 2.0 * etanewphi(j+i*ntheta) * gausspointsphi.qwgt[i]* gausspointstheta.qwgt[j] * sin(theta) * m(1)*m(2);
        etanew(0,2) += M_PI * M_PI / 2.0 * etanewphi(j+i*ntheta) * gausspointsphi.qwgt[i]* gausspointstheta.qwgt[j] * sin(theta) * m(0)*m(2);
        Nfil += M_PI * M_PI / 2.0 * etanewphi(j+i*ntheta) * gausspointsphi.qwgt[i]* gausspointstheta.qwgt[j] * sin(theta);
      }
    }

    // Scale according to vol average
    etanew.Scale(1.0/(M_PI*4.0/3.0));
    etanew(1,0)=etanew(0,1);
    etanew(2,1)=etanew(1,2);
    etanew(2,0)=etanew(0,2);
  }
  else{
    dserror("Choose one of the following integration Methods for internal material law evaluation: 2DTrapez, 2DGauss, 3DGauss");
  }

#endif
double Nfilrate=0.0;
if (Scatra!=1){

Nfilrate=(Nfil-Nfilold)/dt;
Teuchos::RCP<Epetra_MultiVector> rates = immersedmanager->GetPointerToRates();
Teuchos::RCP<Epetra_MultiVector> ratesActin = immersedmanager->GetPointerToRatesActin();
if (time>0.0){
  if(rates == Teuchos::null)
    dserror("rates = Teuchos::null");
  if(ratesActin == Teuchos::null)
    dserror("ratesActin = Teuchos::null");
  double detF = defgrd->Determinant();
  double kROCKeta = params_->kRockEta_;
  double kActin = params_->kActin_;
  double value;
  double valueActin;
  value = (-1.0)*kROCKeta*Nfilrate*detF;
  valueActin = (-1.0)*kActin*Nfilrate*detF;
//  if (kROCKeta==1.0)
//    {value = -1.0e-10;}

  rates->ReplaceGlobalValue(eleGID,gp,value);
  ratesActin->ReplaceGlobalValue(eleGID,gp,valueActin);
  }
}// if Scatra!=1


  LINALG::Matrix<NUM_STRESS_3D,1> sigma(true);//6x1

  sigma(0)=etanew(0,0);
  sigma(1)=etanew(1,1);
  sigma(2)=etanew(2,2);
  sigma(3)=etanew(0,1);
  sigma(4)=etanew(1,2);
  sigma(5)=etanew(0,2);
  sigma.Scale(kstress);


  // Transform Cauchy stress to PK2 stress
  // S = J * F^{-1} sigma F^{-T}
  LINALG::Matrix<NUM_STRESS_3D,1> Sactive(true);//6x1
  CauchytoPK2(Sactive,cauchystress,*defgrd,invdefgrd,sigma);


  LINALG::Matrix<NUM_STRESS_3D,1> Stot(true);//6x1


  // Stress including active and passive part
  if(params.get<int>("iostress")==0)
  {
    Stot.Update(1.0,Sactive,0.0);
    Stot.Update(1.0,Spassive,1.0);
    stress ->Update(1.0,Stot,0.0);
  }
  // only active stress as output
  else
  {
    if (time >0){
      int error=0;
    for (int j=0;j<nphi;j++){
      if (etanewphi(j)<-1.0e-3){
        printf("Negative stress fiber nsf=%f concentration in direction phi(%d) \n",etanewphi(j),j);
        error+=1;
      }
    }
    if (error>0)
    std::cout<<"etanewphi"<<etanewphi<<std::endl;
    if (error>0)
      dserror("Negative stress fiber concentration in ele %d",eleGID);
    }
    //ToDo Adpated for other Integration methods
    double deltaphi = 2.0*M_PI/((double)(nphi));
    // Calculate percentage
    if(Nfil<0)
      dserror("Number of Filaments is smaller than zero NFil=%f",Nfil);

    if(Nfil>0)
    {
      nsfhor_->at(gp) = (etanewphi(0)*deltaphi+etanewphi(10)*deltaphi)/Nfil;
      nsfver_->at(gp) = (etanewphi(5)*deltaphi+etanewphi(15)*deltaphi)/Nfil;
      nsfdiagup_->at(gp) = deltaphi*(etanewphi(2)+etanewphi(3) +etanewphi(12)+ etanewphi(13))/Nfil;
      nsfdiagdown_->at(gp) = deltaphi*(etanewphi(7)+etanewphi(8) +etanewphi(17)+ etanewphi(18))/Nfil;
      // nsfdiag_->at(gp) = (etanewphi(5)*deltaphi+etanewphi(15)*deltaphi)/Nfil;
    }
    else
    {
      nsfhor_->at(gp) = 0.0;
      nsfver_->at(gp) = 0.0;
      nsfdiagup_->at(gp) = 0.0;
      nsfdiagdown_->at(gp) = 0.0;
      // nsfdiag_->at(gp) = (etanewphi(5)*deltaphi+etanewphi(15)*deltaphi)/Nfil;
    }


    Stot.Update(1.0,Sactive,0.0);
    stress ->Update(1.0,Stot,0.0);
  }


  // Update history only when accessed due to material evaluation NOT via Scatra Boundary Dependent Surface Integral
    if (Scatra !=1){
      //etacurr_->at(gp) = contractionrateaverage;
      etacurrphi_->at(gp) = etanewphi;
      //ccurr_->at(gp)=cROCKusednew;
      ccurr_->at(gp)=contractionrateaverage;
      Nfil_ ->at(gp) = Nfil;
      rate_->at(gp) = Nfilrate;
      D_ ->at(gp) = Daverage;
    }


#ifndef  MATERIALFDCHECK
  if (cmat != NULL and analyticalmaterialtangent)
  {
    dserror("no analytical material tangent calculation possible for this material at the moment");
//    // Setup active elasticity tensor cmatactive
//    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatactive(true);
//    SetupCmatActive(cmatactive,rotationrate,strainrate,*defgrd,defgrdrate,R,invdefgrd,etanew,sigmaomegaphinew,cauchystress,params,theta,Csignal);
//
//    // constitutive matrix including active and passive part
//    cmat->Update(1.0,cmatpassive,0.0);
//    cmat->Update(1.0,cmatactive,1.0);
  }
#else
  if (cmat != NULL)
  {
//    // Setup active elasticity tensor cmatactive
//    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatactive(true);
//    SetupCmatActive(cmatactive,rotationrate,strainrate,*defgrd,defgrdrate,R,invdefgrd,etanew,sigmaomegaphinew,cauchystress,params,theta,Csignal);
//
//    // constitutive matrix including active and passive part
//    cmat->Update(1.0,cmatpassive,0.0);
//    cmat->Update(1.0,cmatactive,1.0);
  }
#endif

} // end mat evaluate


///*----------------------------------------------------------------------*
// | Dummy Evaluate Material                                  rauch  01/16|
// *----------------------------------------------------------------------*
// The stress response is decomposed into a passive and an active part:
//     \sigma = \sigma_{passive} + \sigma_{active}
// */
//void MAT::BioChemoMechanoCellActiveFiber::Evaluate(const LINALG::Matrix<3,3>* defgrd,
//                           const LINALG::Matrix<6,1>* glstrain,
//                           Teuchos::ParameterList& params,
//                           LINALG::Matrix<6,1>* stress,
//                           LINALG::Matrix<6,6>* cmat,
//                           const int eleGID)
//{
//  //******************
//  // PASSIVE PART
//  //******************
//  // Initialize passive stress and elasticity tensor
//  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatpassive(true);
//  LINALG::Matrix<NUM_STRESS_3D,1> Spassive(true);
//
//  //if (cmat != NULL)
//  //{
//    // Evaluate passive PK2 stress Spassive and passive elasticity tensor cmatpassive
//    matpassive_->Evaluate(defgrd,glstrain,params,&Spassive,&cmatpassive,eleGID);
//  //}
//
//  //********************
//  // DUMMY ACTIVE PART
//  //********************
//
//  // Setup inverse of deformation gradient
//  LINALG::Matrix<3,3> invdefgrd(*defgrd);
//  invdefgrd.Invert();
//
//  LINALG::Matrix<NUM_STRESS_3D,1> sigma(true);//6x1
//  LINALG::Matrix<3,3> cauchystress(true);//3x3
//
//  sigma(0)=0.001;
//  sigma(1)=0.0;
//  sigma(2)=0.0;
//  sigma(3)=0.0;
//  sigma(4)=0.0;
//  sigma(5)=0.0;
//
//  // Transform Cauchy stress to PK2 stress
//  // S = J * F^{-1} sigma F^{-T}
//  LINALG::Matrix<NUM_STRESS_3D,1> Sactive(true);//6x1
//  CauchytoPK2(Sactive,cauchystress,*defgrd,invdefgrd,sigma);
//
//  LINALG::Matrix<NUM_STRESS_3D,1> Stot(true);//6x1
//
//  // Stress including active and passive part
//  if(params.get<int>("iostress")==0)
//  {
//    Stot.Update(1.0,Sactive,0.0);
//    Stot.Update(1.0,Spassive,1.0);
//    stress ->Update(1.0,Stot,0.0);
//  }
//  // only active stress as output
//  else
//  {
//    Stot.Update(1.0,Sactive,0.0);
//    stress ->Update(1.0,Stot,0.0);
//  }
//
//
//} // end mat evaluate


/*----------------------------------------------------------------------*
 | Calculate Dissociation                                   rauch 01/16 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::Dissociation(
  double rate,
  double& D,
  double k_B,
  double ratemax)
{
  /* D(rate)
     |
     |  \
     |   \
     |    \  line
     |     \
     |      \
     |       \
     |        \
     |          .  parabola
   0-|             .--------------------
     |_____|_______|__________|_________ rate
                ratemax     0
  */

  if (rate>ratemax/2.0)
      D=0.0;
  else if (ratemax<rate and rate<=ratemax/2.0)
      D=1.0/(-ratemax)*k_B*(rate-ratemax/2.0)*(rate-ratemax/2.0);
  else if (rate<=ratemax)
      D = -k_B*(rate -ratemax)-k_B*ratemax/4.0;
  else
    dserror("You should not end up here. Something went wrong in Dissociation().");

}  // Dissociation


/*----------------------------------------------------------------------*
 | Calculation of Nodal values from GP values               rauch 01/16 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::GPtoNodes(
    LINALG::Matrix<8,1> invalue,
    LINALG::Matrix<8,1>& outvalue)
{
    LINALG::Matrix<3,1> xsi(true);
    LINALG::Matrix<8,8> S(true);
    LINALG::Matrix<NUMNOD_SOH8,1> shapefcts;

    DRT::UTILS::IntPointsAndWeights<NUMDIM_SOH8> intpoints(DRT::UTILS::intrule_hex_8point);
    for (int gp=0;gp<8;gp++){
      for (int node=0;node<8;node++){
        for (int i=0;i<3;i++){
          xsi(i)=(intpoints.IP().qxg)[gp][i];
        }
        DRT::UTILS::shape_function<DRT::Element::hex8>(xsi,shapefcts);
        S(node,gp)=shapefcts(node);
      }
    }
    LINALG::Matrix<8,8> Sinv(true);
    Sinv.Update(S);
    LINALG::FixedSizeSerialDenseSolver<8,8> Sinversion;
    Sinversion.SetMatrix(Sinv);
    int err = Sinversion.Invert();
    if(err != 0){
      std::cout<<"A "<<S<<std::endl;
      dserror("Inversion of 8x8 matrix failed with errorcode %d",err);
    }

    outvalue.MultiplyNN(Sinv,invalue);

}


/*----------------------------------------------------------------------*
| Calculation of Nodal values from GP values for eta        rauch 01/16 |
*-----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::EtaGPtoNodes(
        LINALG::Matrix<nphi,1>& node1,
        LINALG::Matrix<nphi,1>& node2,
        LINALG::Matrix<nphi,1>& node3,
        LINALG::Matrix<nphi,1>& node4,
        LINALG::Matrix<nphi,1>& node5,
        LINALG::Matrix<nphi,1>& node6,
        LINALG::Matrix<nphi,1>& node7,
        LINALG::Matrix<nphi,1>& node8
        )
    {
      LINALG::Matrix<8,1> atgp(true);
      LINALG::Matrix<8,1> atnode(true);
   for (int entry=0;entry<nphi;entry++){
     for(int gp=0;gp<8;gp++){
       LINALG::Matrix<nphi,1> etalastphi = etalastphi_->at(gp);
       atgp(gp)= etalastphi(entry);

     }
     GPtoNodes(atgp,atnode);
     atgp.Clear();
     node1(entry)=atnode(0);
     node2(entry)=atnode(1);
     node3(entry)=atnode(2);
     node4(entry)=atnode(3);
     node5(entry)=atnode(4);
     node6(entry)=atnode(5);
     node7(entry)=atnode(6);
     node8(entry)=atnode(7);
     atnode.Clear();
   }
   }


/*----------------------------------------------------------------------*
| Calculation of eta at specific Boundary GP xsi from Scatra Brdy calc  |
*-----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::EtaAtBrdyGP(
        LINALG::Matrix<nphi,1>& outvalue,
        Teuchos::ParameterList& params)
    {
  LINALG::Matrix<nphi,1> node1(true);
  LINALG::Matrix<nphi,1> node2(true);
  LINALG::Matrix<nphi,1> node3(true);
  LINALG::Matrix<nphi,1> node4(true);
  LINALG::Matrix<nphi,1> node5(true);
  LINALG::Matrix<nphi,1> node6(true);
  LINALG::Matrix<nphi,1> node7(true);
  LINALG::Matrix<nphi,1> node8(true);
  EtaGPtoNodes(node1,node2,node3,node4,node5,node6,node7,node8);
  LINALG::Matrix<3,1> xsi(true);
  LINALG::Matrix<NUMNOD_SOH8,1> shapefcts;

    xsi = params.get<LINALG::Matrix<3,1> >("xsi");
    DRT::UTILS::shape_function<DRT::Element::hex8>(xsi,shapefcts);
    for (int entry=0;entry<nphi;entry++){
        outvalue(entry)=  shapefcts(0)*node1(entry) + shapefcts(1)*node2(entry)
                        + shapefcts(2)*node3(entry) + shapefcts(3)*node4(entry)
                        + shapefcts(4)*node5(entry) + shapefcts(5)*node6(entry)
                        + shapefcts(6)*node7(entry) + shapefcts(7)*node8(entry);
      }

}


/*----------------------------------------------------------------------*
| Calculation of Nodal values from GP values for Deformation Gradient   |
*-----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::DefGradGPtoNodes(
        LINALG::Matrix<6,1>& node1,
        LINALG::Matrix<6,1>& node2,
        LINALG::Matrix<6,1>& node3,
        LINALG::Matrix<6,1>& node4,
        LINALG::Matrix<6,1>& node5,
        LINALG::Matrix<6,1>& node6,
        LINALG::Matrix<6,1>& node7,
        LINALG::Matrix<6,1>& node8
        )
    {
      LINALG::Matrix<8,1> atgp(true);
      LINALG::Matrix<8,1> atnode(true);
      LINALG::Matrix<3,3> defgrdlast;
      LINALG::Matrix<6,1> defgrdlastvec;
   for (int entry=0;entry<6;entry++){
     for(int gp=0;gp<8;gp++){
        defgrdlast = histdefgrdlast_->at(gp);
        defgrdlastvec(0)= defgrdlast(0,0);
        defgrdlastvec(1)= defgrdlast(1,1);
        defgrdlastvec(2)= defgrdlast(2,2);
        defgrdlastvec(3)= defgrdlast(0,1);
        defgrdlastvec(4)= defgrdlast(1,2);
        defgrdlastvec(5)= defgrdlast(0,2);

       atgp(gp)= defgrdlastvec(entry);

     }
     GPtoNodes(atgp,atnode);
     atgp.Clear();
     node1(entry)=atnode(0);
     node2(entry)=atnode(1);
     node3(entry)=atnode(2);
     node4(entry)=atnode(3);
     node5(entry)=atnode(4);
     node6(entry)=atnode(5);
     node7(entry)=atnode(6);
     node8(entry)=atnode(7);
     atnode.Clear();
   }
   }


/*----------------------------------------------------------------------*
| Calculation of deformation gradient at specific Boundary GP xsi       |
*-----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::DefGradAtBrdyGP(
        LINALG::Matrix<6,1>& outvalue,
        Teuchos::ParameterList& params)
    {
  LINALG::Matrix<6,1> node1(true);
  LINALG::Matrix<6,1> node2(true);
  LINALG::Matrix<6,1> node3(true);
  LINALG::Matrix<6,1> node4(true);
  LINALG::Matrix<6,1> node5(true);
  LINALG::Matrix<6,1> node6(true);
  LINALG::Matrix<6,1> node7(true);
  LINALG::Matrix<6,1> node8(true);
  DefGradGPtoNodes(node1,node2,node3,node4,node5,node6,node7,node8);
  LINALG::Matrix<3,1> xsi(true);
  LINALG::Matrix<NUMNOD_SOH8,1> shapefcts;

    xsi = params.get<LINALG::Matrix<3,1> >("xsi");
    DRT::UTILS::shape_function<DRT::Element::hex8>(xsi,shapefcts);
    for (int entry=0;entry<6;entry++)
    {
        outvalue(entry)=  shapefcts(0)*node1(entry) + shapefcts(1)*node2(entry)
                        + shapefcts(2)*node3(entry) + shapefcts(3)*node4(entry)
                        + shapefcts(4)*node5(entry) + shapefcts(5)*node6(entry)
                        + shapefcts(6)*node7(entry) + shapefcts(7)*node8(entry);
    }

}


/*----------------------------------------------------------------------*
| Calculation of cROCKused at Boundary GP                               |
*----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::cROCKlastAtBrdyGP(
        double& cROCKlast,
         Teuchos::ParameterList& params)
    {
      LINALG::Matrix<8,1> cGP(true);
      LINALG::Matrix<8,1> cNode(true);

      for(int gp=0;gp<8;gp++){
        cGP(gp)= clast_->at(gp);
      }
      GPtoNodes(cGP,cNode);

  LINALG::Matrix<3,1> xsi(true);
  LINALG::Matrix<NUMNOD_SOH8,1> shapefcts;

    xsi = params.get<LINALG::Matrix<3,1> >("xsi");
    DRT::UTILS::shape_function<DRT::Element::hex8>(xsi,shapefcts);

    for (int i=0;i<8;i++){
      cROCKlast+=shapefcts(i)*cNode(i);
    }

}


/*----------------------------------------------------------------------*
| Calculation of cROCKlast at Boundary GP                               |
*-----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::cROCKusedAtBrdyGP(
        double& cROCKused,
         Teuchos::ParameterList& params)
    {
      LINALG::Matrix<8,1> cGP(true);
      LINALG::Matrix<8,1> cNode(true);

      for(int gp=0;gp<8;gp++){
        cGP(gp)= clast_->at(gp);
      }
      GPtoNodes(cGP,cNode);

  LINALG::Matrix<3,1> xsi(true);
  LINALG::Matrix<NUMNOD_SOH8,1> shapefcts;

    xsi = params.get<LINALG::Matrix<3,1> >("xsi");
    DRT::UTILS::shape_function<DRT::Element::hex8>(xsi,shapefcts);

    for (int i=0;i<8;i++){
      cROCKused+=shapefcts(i)*cNode(i);
    }

}


/*----------------------------------------------------------------------------------*
 | Calculation of deformation gradient rate, rotation tensor, strain rate and       |
 | rotation rate (finite difference scheme)                            rauch  01/16 |
 *----------------------------------------------------------------------------------*/
// see also: viscoelasthyper.cpp line 781 ff
void MAT::BioChemoMechanoCellActiveFiber::SetupRates(
    LINALG::Matrix<3,3> defgrd,
    LINALG::Matrix<3,3> invdefgrd,
    Teuchos::ParameterList& params,
    LINALG::Matrix<3,3>& defgrdrate,
    LINALG::Matrix<3,3>& R,
    LINALG::Matrix<6,1>& strainrate,
    LINALG::Matrix<3,3>& rotationrate,
    const int& gp,
    const double& dt)
{
  // Read history
  LINALG::Matrix<3,3> defgrdlast;
  int Scatra = params.get<int>("FromSactraBoundary",0);

  if (Scatra !=1)
  {
   defgrdlast = histdefgrdlast_->at(gp);
  }
  else
  {
    LINALG::Matrix<6,1> defgrdlastvec;
    DefGradAtBrdyGP(defgrdlastvec,params);
    defgrdlast(0,0)=defgrdlastvec(0);
    defgrdlast(1,1)=defgrdlastvec(1);
    defgrdlast(2,2)=defgrdlastvec(2);
    defgrdlast(0,1)=defgrdlastvec(3);
    defgrdlast(1,0)=defgrdlastvec(3);
    defgrdlast(2,1)=defgrdlastvec(4);
    defgrdlast(1,2)=defgrdlastvec(4);
    defgrdlast(0,2)=defgrdlastvec(5);
    defgrdlast(2,0)=defgrdlastvec(5);

  }

  // Rate of deformation gradient: \dot{F} = \frac {F^{n+1} - F^{n}} {\Delta t}
  defgrdrate.Update(1.0,defgrd,0.0);
  defgrdrate.Update(-1.0,defgrdlast,1.0);
  defgrdrate.Scale(1.0/dt);

  // Calculate velocity gradient l = \dot{F}.F^{-1}
  LINALG::Matrix<3,3> velgradient(true);
  velgradient.MultiplyNN(defgrdrate,invdefgrd);

  // Rate of strain/symmetric part of velocity gradient
  // d = 0.5 * (l + l^{T}) = 0.5 * (\dot{F}F^{-1} + (\dot{F}F^{-1})^{T})
  // Remark: strain-like 6-Voigt vector
  strainrate(0) = velgradient(0,0) + velgradient(0,0);
  strainrate(1) = velgradient(1,1) + velgradient(1,1);
  strainrate(2) = velgradient(2,2) + velgradient(2,2);
  strainrate(3) = velgradient(0,1) + velgradient(1,0);
  strainrate(4) = velgradient(1,2) + velgradient(2,1);
  strainrate(5) = velgradient(0,2) + velgradient(2,0);
  strainrate.Scale(0.5);


  // Rate of rotation tensor (!= skew part w of velocity gradient l, see Holzapfel S.99)
  // Determine rotation tensor R from F (F=R*U) -> polar decomposition of displacement based F
  LINALG::Matrix<3,3> Q(true);
  LINALG::Matrix<3,3> S(true);
  LINALG::Matrix<3,3> VT(true);

  // Calculate rotcurr from defgrd
  LINALG::SVD<3,3>(defgrd,Q,S,VT); // Singular Value Decomposition,analogously to micromaterial_evaluate.cpp lines 81ff
  R.MultiplyNN(Q,VT);

  // Update history of deformation gradient only via material evaluation
  // (NOT when evaluation is done via Scatra Surface Integral)
  if (Scatra !=1){
    histdefgrdcurr_->at(gp) = defgrd;
    strainratecurr_->at(gp) = strainrate;
   }
  else
  {
    strainratecurrsurf_->at(gp) = strainrate;
  }
}  // Evaluate()


/*----------------------------------------------------------------------*
 | pull back of spatial stresses                           rauch  01/16 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::CauchytoPK2(
  LINALG::Matrix<6,1>& Sactive,
  LINALG::Matrix<3,3>& cauchystress,
  LINALG::Matrix<3,3> defgrd,
  LINALG::Matrix<3,3> invdefgrd,
  LINALG::Matrix<6,1> sigma)
{
  // calculate the Jacobi-determinant
  double detF = defgrd.Determinant();   // const???

  // Convert stress like 6x1-Voigt vector to 3x3 matrix
  cauchystress(0,0) = sigma(0);
  cauchystress(0,1) = sigma(3);
  cauchystress(0,2) = sigma(5);
  cauchystress(1,0) = cauchystress(0,1);
  cauchystress(1,1) = sigma(1);
  cauchystress(1,2) = sigma(4);
  cauchystress(2,0) = cauchystress(0,2);
  cauchystress(2,1) = cauchystress(1,2);
  cauchystress(2,2) = sigma(2);

  // S = J * F^{-1} * sigma * F^{-T}
  LINALG::Matrix<3,3> temp(true);
  LINALG::Matrix<3,3> S(true);
  temp.MultiplyNN(invdefgrd,cauchystress);
  S.MultiplyNT(temp,invdefgrd);
  S.Scale(detF);

//#ifdef DEBUG
//  if(abs(S(1,2)-S(2,1))>1e-13 or abs(S(0,2)-S(2,0))>1e-13 or abs(S(0,1)-S(1,0))>1e-13)
//  {
//    std::cout<<S<<std::endl;
//    dserror("PK2 not symmetric!!");
//  }
//#endif


  // Sactive is stress like 6x1-Voigt vector
  Sactive(0) = S(0,0);
  Sactive(1) = S(1,1);
  Sactive(2) = S(2,2);
  Sactive(3) = S(0,1);
  Sactive(4) = S(1,2);
  Sactive(5) = S(0,2);

}  // CauchytoPK2()


/*----------------------------------------------------------------------*
 | push forward of material stresses                                    |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::PK2toCauchy(
  LINALG::Matrix<6,1> Sactive,
  LINALG::Matrix<3,3>& PK2stress,
  LINALG::Matrix<3,3> defgrd,
  LINALG::Matrix<6,1>& cauchystress)
{
  // calculate the Jacobi-determinant
  //const double detF = defgrd.Determinant();   // const???
  double detF = defgrd.Determinant();
  if (abs(detF)<1.0e-10)
  {
    std::cout<<"Error: detF in PK2toCauchy =  \n"<<detF<<std::endl;
  }
  // Convert stress like 6x1-Voigt vector to 3x3 matrix
  PK2stress(0,0) = Sactive(0);
  PK2stress(0,1) = Sactive(3);
  PK2stress(0,2) = Sactive(5);
  PK2stress(1,0) = PK2stress(0,1);
  PK2stress(1,1) = Sactive(1);
  PK2stress(1,2) = Sactive(4);
  PK2stress(2,0) = PK2stress(0,2);
  PK2stress(2,1) = PK2stress(1,2);
  PK2stress(2,2) = Sactive(2);

  // sigma = 1/J * F * sigma * F^{T}
  LINALG::Matrix<3,3> temp(true);
  LINALG::Matrix<3,3> sigma(true);
  temp.MultiplyNN(defgrd,PK2stress);
  sigma.MultiplyNT(temp,defgrd);
  sigma.Scale(1./detF);


  // Cauchy_passive is stress like 6x1-Voigt vector
  cauchystress(0) = sigma(0,0);
  cauchystress(1) = sigma(1,1);
  cauchystress(2) = sigma(2,2);
  cauchystress(3) = sigma(0,1);
  cauchystress(4) = sigma(1,2);
  cauchystress(5) = sigma(0,2);

}  // PK2toCauchy()


///*----------------------------------------------------------------------*
// | pull back of spatial stresses                           rauch  01/16 |
// *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::GLtoEA(
  LINALG::Matrix<6,1> glstrain,
  LINALG::Matrix<3,3> invdefgrd,
  LINALG::Matrix<3,3>& eastrain)    // muss noch initialisiert werden
{
  // Convert strain like 6x1-Voigt vector to 3x3 matrix
  LINALG::Matrix<3,3> greenlagrange(true);
  greenlagrange(0,0) = glstrain(0);
  greenlagrange(0,1) = 0.5*glstrain(3);            // 0.5???
  greenlagrange(0,2) = 0.5*glstrain(5);            // 0.5???
  greenlagrange(1,0) = greenlagrange(0,1);
  greenlagrange(1,1) = glstrain(1);
  greenlagrange(1,2) = 0.5*glstrain(4);            //0.5???
  greenlagrange(2,0) = greenlagrange(0,2);
  greenlagrange(2,1) = greenlagrange(1,2);
  greenlagrange(2,2) = glstrain(2);


  // e = F^{-T} * E * F^{-1}
  LINALG::Matrix<3,3> temp(true);
  temp.MultiplyTN(invdefgrd,greenlagrange);
  eastrain.MultiplyNT(temp,invdefgrd);

}  // GLtoEA


/*----------------------------------------------------------------------*
 |  Names of gp data to be visualized                      rauch  01/16 |
 *----------------------------------------------------------------------*/
void MAT::BioChemoMechanoCellActiveFiber::VisNames(std::map<std::string,int>& names)
{
  std::string cellcontr = "Nfil";
  names[cellcontr] = 1; // scalar
//  cellcontr  = "Contraction Rate";
//  names[cellcontr] = 1; // scalar
  cellcontr = "RateAverage";
  names[cellcontr] = 1; // scalar
  cellcontr = "Drate";
  names[cellcontr] = 1; // scalar
  cellcontr = "rate";
  names[cellcontr] = 1; // scalar
  cellcontr = "Nver";
  names[cellcontr] = 1; // scalar
  cellcontr = "Nhor";
  names[cellcontr] = 1; // scalar
  cellcontr = "Ndiagup";
  names[cellcontr] = 1; // scalar
  cellcontr = "Ndiagdown";
  names[cellcontr] = 1; // scalar
  matpassive_->VisNames(names);

} // VisNames()


/*----------------------------------------------------------------------*
 |  gp data to be visualized                               rauch  01/16 |
 *----------------------------------------------------------------------*/
bool MAT::BioChemoMechanoCellActiveFiber::VisData(const std::string& name, std::vector<double>& data, int numgp , int eleID)
{
    if (name == "Nfil")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<numgp; iter++)
        temp += Nfil_()->at(iter);
      data[0] = temp/numgp;
    }
  else if (name == "RateAverage")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<numgp; iter++)
        temp += clast_()->at(iter);
      data[0] = temp/numgp;
    }
  else if (name == "Drate")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<numgp; iter++)
        temp += D_()->at(iter);
      data[0] = temp/numgp;
    }
  else if (name == "rate")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<numgp; iter++)
        temp += rate_()->at(iter);
      data[0] = temp/numgp;
    }
  else if (name == "Nver")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<numgp; iter++)
        temp += nsfver_()->at(iter);
      data[0] = temp/numgp;
    }
  else if (name == "Nhor")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<numgp; iter++)
        temp += nsfhor_()->at(iter);
      data[0] = temp/numgp;
    }
  else if (name == "Ndiagup")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<numgp; iter++)
        temp += nsfdiagup_()->at(iter);
      data[0] = temp/numgp;
    }
  else if (name == "Ndiagdown")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<numgp; iter++)
        temp += nsfdiagdown_()->at(iter);
      data[0] = temp/numgp;
    }
  else
  {
    return matpassive_->VisData(name, data, numgp, eleID);
  }

  return true;
}
