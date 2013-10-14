/*----------------------------------------------------------------------*/
/*!
\file str_statinvanalysis.cpp
\brief Statistical inverse analysis for structures

<pre>
Maintainer: Sebastian Kehl
            Kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de

</pre>
*/

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "str_statinvanalysis.H"
#include "../drt_inv_analysis/stat_inv_analysis.H"
#include "../drt_inpar/inpar_statinvanalysis.H"
#include "../drt_inv_analysis/stat_inv_ana_graddesc.H"
#include "../drt_inv_analysis/stat_inv_ana_mc.H"
#include "../drt_inv_analysis/stat_inv_ana_lbfgs.H"
#include "../drt_lib/drt_discret.H"

/*======================================================================*/
/* Statistical inverse analysis of structures */
void STR::statinvanalysis()
{
  // get input lists
  const Teuchos::ParameterList& statinvp = DRT::Problem::Instance()->StatInverseAnalysisParams();

  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis("structure");

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();
  if (!actdis->HaveDofs()) actdis->FillComplete();

  // context for output and restart
  //Teuchos::RCP<IO::DiscretizationWriter> output
  //  = Teuchos::rcp(new IO::DiscretizationWriter(actdis));

  switch(DRT::INPUT::IntegralValue<INPAR::STR::StatInvAnalysisType>(statinvp,"STAT_INV_ANALYSIS"))
  {
    case INPAR::STR::stat_inv_graddesc:
    {
      Teuchos::RCP<STR::INVANA::StatInvAnalysis>  ia = Teuchos::rcp(new STR::INVANA::StatInvAnaGradDesc(actdis));
      ia->Optimize();

      DRT::Problem::Instance()->AddFieldTest(ia->CreateFieldTest());
      DRT::Problem::Instance()->TestAll(actdis->Comm());
    }
    break;
    case INPAR::STR::stat_inv_lbfgs:
    {
      Teuchos::RCP<STR::INVANA::StatInvAnalysis>  ia = Teuchos::rcp(new STR::INVANA::StatInvAnaLBFGS(actdis));
      ia->Optimize();

      DRT::Problem::Instance()->AddFieldTest(ia->CreateFieldTest());
      DRT::Problem::Instance()->TestAll(actdis->Comm());
    }
    break;
    case INPAR::STR::stat_inv_mc:
    {
      Teuchos::RCP<STR::INVANA::StatInvAnalysis>  ia = Teuchos::rcp(new STR::INVANA::StatInvAnaMC(actdis));
      ia->Optimize();
    }
    break;
    default:
      dserror("Unknown type of statistical inverse analysis");
    break;
  }

  // done
  return;
} // end str_statinvanalysis()


/*----------------------------------------------------------------------*/
