/*----------------------------------------------------------------------*/
/*!
 \file matlist_bondreacs.cpp

 \brief material list for bond reactions.

 \level 2

 <pre>
   \maintainer Andreas Rauch
               rauch@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289 - 15240
 </pre>
 *----------------------------------------------------------------------*/


#include <vector>
#include "matlist_bondreacs.H"
#include "scatra_bondreac_mat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"



/*----------------------------------------------------------------------*
 | standard constructor                                     rauch 12/16 |
 *----------------------------------------------------------------------*/
MAT::PAR::MatListBondReacs::MatListBondReacs(
    Teuchos::RCP<MAT::PAR::Material> matdata
)
: MatList(matdata),
  MatListReactions(matdata)
{}


Teuchos::RCP<MAT::Material> MAT::PAR::MatListBondReacs::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MatListBondReacs(this));
}


MAT::MatListBondReacsType MAT::MatListBondReacsType::instance_;


DRT::ParObject* MAT::MatListBondReacsType::Create( const std::vector<char> & data )
{
  MAT::MatListBondReacs* MatListBondReacs = new MAT::MatListBondReacs();
  MatListBondReacs->Unpack(data);
  return MatListBondReacs;
}


/*----------------------------------------------------------------------*
 | construct empty material object                          rauch 12/16 |
 *----------------------------------------------------------------------*/
MAT::MatListBondReacs::MatListBondReacs()
: MatList(),
  MatListReactions(),
  paramsbondreac_(NULL)
{
}


/*----------------------------------------------------------------------*
 | construct the material object given material paramete    rauch 12/16 |
 *----------------------------------------------------------------------*/
MAT::MatListBondReacs::MatListBondReacs(MAT::PAR::MatListBondReacs* params)
: MatList(params),
  MatListReactions(params),
  paramsbondreac_(params)
{
  // setup of material map
  if (paramsbondreac_->local_)
  {
    SetupMatMap();
  }
}


/*----------------------------------------------------------------------*
 | setup of material map                                    rauch 12/16 |
 *----------------------------------------------------------------------*/
void MAT::MatListBondReacs::Initialize()
{

  if(paramsbondreac_!=NULL)
  {
    std::vector<int>::const_iterator m;
    for (m=paramsbondreac_->ReacIds()->begin(); m!=paramsbondreac_->ReacIds()->end(); ++m)
    {
      const int reacid = *m;
      Teuchos::RCP<MAT::Material> mat = MaterialById(reacid);
      if (mat == Teuchos::null) dserror("Failed to allocate this material");
      Teuchos::RCP<MAT::ScatraBondReacMat> reacmat =
          Teuchos::rcp_dynamic_cast<MAT::ScatraBondReacMat>(mat,true);
      reacmat->Initialize();
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | reset everything                                         rauch 12/16 |
 *----------------------------------------------------------------------*/
void MAT::MatListBondReacs::Clear()
{
  paramsbondreac_ = NULL;
  return;
}


/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class           rauch 12/16 |
 *----------------------------------------------------------------------*/
void MAT::MatListBondReacs::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (paramsbondreac_ != NULL) matid = paramsbondreac_->Id();  // in case we are in post-process mode

  AddtoPack(data,matid);

  // Pack base class material
  MAT::MatListReactions::Pack(data);
}


/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class           rauch 12/16 |
 *----------------------------------------------------------------------*/
void MAT::MatListBondReacs::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  Clear();

  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover paramsreac_
  int matid(-1);
  ExtractfromPack(position,data,matid);
  paramsbondreac_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        //Note: We need to do a dynamic_cast here since Bonds, Reaction, and Bond-reaction are in a diamond inheritance structure
        paramsbondreac_ = dynamic_cast<MAT::PAR::MatListBondReacs*>(mat);
      }
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  MAT::MatList::ExtractfromPack(position,data,basedata);
  MAT::MatListReactions::Unpack(basedata);

  // in the postprocessing mode, we do not unpack everything we have packed
  // -> position check cannot be done in this case
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}


/*----------------------------------------------------------------------*
 | calculate body force term                                rauch 12/16 |
 *----------------------------------------------------------------------*/
double MAT::MatListBondReacs::CalcReaBodyForceTerm(
    const int k,
    const std::vector<double>& phinp,
    const std::vector<double>& phin,
    const double violation,
    const double porosity,
    const double* gpcoord,
    const double scale
) const
{
  double bodyforcetermK=0.0;

  for (int condnum = 0;condnum < NumReac();condnum++)
  {
    const int reacid = ReacID(condnum);
    const Teuchos::RCP<const MAT::ScatraBondReacMat> reacmat = Teuchos::rcp_static_cast<const MAT::ScatraBondReacMat>(MaterialById(reacid));

    bodyforcetermK += reacmat->CalcReaBodyForceTerm(k, phinp, phin, violation, porosity, scale, gpcoord);
  }

  return bodyforcetermK;
}


/*----------------------------------------------------------------------*
 | calculate body force term derivative                     rauch 12/16 |
 *----------------------------------------------------------------------*/
void MAT::MatListBondReacs::CalcReaBodyForceDerivMatrix(
    const int k,
    std::vector<double>& derivs,
    const std::vector<double>& phinp,
    const std::vector<double>& phin,
    const double violation,
    const double porosity,
    const double* gpcoord,
    const double scale
) const
{

  for (int condnum = 0;condnum < NumReac();condnum++)
  {
    const int reacid = ReacID(condnum);
    const Teuchos::RCP<const MAT::ScatraBondReacMat> reacmat = Teuchos::rcp_static_cast<const MAT::ScatraBondReacMat>(MaterialById(reacid));

    reacmat->CalcReaBodyForceDerivMatrix(k, derivs, phinp, phin, violation, porosity, scale, gpcoord);
  }

  return;
}


/*----------------------------------------------------------------------*
 | calculate body force term                                rauch 12/16 |
 *----------------------------------------------------------------------*/
double MAT::MatListBondReacs::CalcReaBodyForceTerm(
    const int k,
    const std::vector<double>& phinp,
    const std::vector<double>& phin,
    const std::vector<std::pair<std::string,double> >& constants,
    const double violation,
    const double porosity,
    const double* gpcoord,
    const double scale
) const
{
  double bodyforcetermK=0.0;

  for (int condnum = 0;condnum < NumReac();condnum++)
  {
    const int reacid = ReacID(condnum);
    const Teuchos::RCP<const MAT::ScatraBondReacMat> reacmat = Teuchos::rcp_static_cast<const MAT::ScatraBondReacMat>(MaterialById(reacid));

    bodyforcetermK += reacmat->CalcReaBodyForceTerm(k, phinp, phin, violation, porosity, scale, gpcoord);
  }

  return bodyforcetermK;
}


/*----------------------------------------------------------------------*
 | calculate body force term derivative                     rauch 12/16 |
 *----------------------------------------------------------------------*/
void MAT::MatListBondReacs::CalcReaBodyForceDerivMatrix(
    const int k,
    std::vector<double>& derivs,
    const std::vector<double>& phinp,
    const std::vector<double>& phin,
    const std::vector<std::pair<std::string,double> >& constants,
    const double violation,
    const double porosity,
    const double* gpcoord,
    const double scale
) const
{

  for (int condnum = 0;condnum < NumReac();condnum++)
  {
    const int reacid = ReacID(condnum);
    const Teuchos::RCP<const MAT::ScatraBondReacMat> reacmat = Teuchos::rcp_static_cast<const MAT::ScatraBondReacMat>(MaterialById(reacid));

    reacmat->CalcReaBodyForceDerivMatrix(k, derivs, phinp, phin, violation, porosity, scale, gpcoord);
  }

  return;
}
