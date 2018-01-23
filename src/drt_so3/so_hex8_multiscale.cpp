/*!----------------------------------------------------------------------
\file so_hex8_multiscale.cpp
\brief
\level 2

<pre>
\maintainer Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/


#include "so_hex8.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "so_surface.H"



/*----------------------------------------------------------------------*
 |  homogenize material density (public)                        lw 07/07|
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void DRT::ELEMENTS::So_hex8::soh8_homog(Teuchos::ParameterList&  params)
{
  if(DRT::Problem::Instance(0)->GetNPGroup()->SubComm()->MyPID() == Owner())
  {
    double homogdens = 0.;
    const static std::vector<double> weights = soh8_weights();

    for (unsigned gp=0; gp<NUMGPT_SOH8; ++gp)
    {
      const double density = Material()->Density(gp);
      homogdens += detJ_[gp] * weights[gp] * density;
    }

    double homogdensity = params.get<double>("homogdens", 0.0);
    params.set("homogdens", homogdensity+homogdens);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Set EAS internal variables on the microscale (public)       lw 04/08|
 *----------------------------------------------------------------------*/
// the microscale internal EAS data have to be saved separately for every
// macroscopic Gauss point and set before the determination of microscale
// stiffness etc.

void DRT::ELEMENTS::So_hex8::soh8_set_eas_multi(Teuchos::ParameterList&  params)
{
  if (eastype_ != soh8_easnone)
  {
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > oldalpha =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("oldalpha",Teuchos::null);
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > oldfeas =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("oldfeas",Teuchos::null);
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > oldKaainv =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("oldKaainv",Teuchos::null);
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > oldKda =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("oldKda",Teuchos::null);

    if (oldalpha==Teuchos::null || oldfeas==Teuchos::null || oldKaainv==Teuchos::null || oldKda==Teuchos::null)
      dserror("Cannot get EAS internal data from parameter list for multi-scale problems");

    data_.Add("alpha", (*oldalpha)[Id()]);
    data_.Add("feas", (*oldfeas)[Id()]);
    data_.Add("invKaa", (*oldKaainv)[Id()]);
    data_.Add("Kda", (*oldKda)[Id()]);
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Initialize EAS internal variables on the microscale         lw 03/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_eas_init_multi(Teuchos::ParameterList&  params)
{
  if (eastype_ != soh8_easnone)
  {
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > lastalpha =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("lastalpha",Teuchos::null);
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > oldalpha =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("oldalpha",Teuchos::null);
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > oldfeas =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("oldfeas",Teuchos::null);
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > oldKaainv =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("oldKaainv",Teuchos::null);
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > oldKda =
      params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("oldKda",Teuchos::null);

    (*lastalpha)[Id()] = Teuchos::rcp(new Epetra_SerialDenseMatrix(neas_, 1));
    (*oldalpha)[Id()]  = Teuchos::rcp(new Epetra_SerialDenseMatrix(neas_, 1));
    (*oldfeas)[Id()]   = Teuchos::rcp(new Epetra_SerialDenseMatrix(neas_, 1));
    (*oldKaainv)[Id()] = Teuchos::rcp(new Epetra_SerialDenseMatrix(neas_, neas_));
    (*oldKda)[Id()]    = Teuchos::rcp(new Epetra_SerialDenseMatrix(neas_, NUMDOF_SOH8));
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Read restart on the microscale                              lw 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_read_restart_multi()
{
  Teuchos::RCP<MAT::Material> mat = Material();

  if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
  {
    MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
    int eleID = Id();
    bool eleowner = false;
    if (DRT::Problem::Instance()->GetDis("structure")->Comm().MyPID()==Owner()) eleowner = true;

    for (unsigned gp=0; gp<NUMGPT_SOH8; ++gp)
      micro->ReadRestart(gp, eleID, eleowner);
  }

  return;
}


