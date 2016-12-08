/*!----------------------------------------------------------------------
\file cardiovascular0d_resulttest.cpp

\brief Monolithic coupling of 3D structural dynamics and 0D cardiovascular flow models

\level 2

<pre>
\maintainer Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>
*----------------------------------------------------------------------*/

#include <string>
#include "cardiovascular0d_resulttest.H"
#include "cardiovascular0d_manager.H"
#include "cardiovascular0d.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Cardiovascular0DResultTest::Cardiovascular0DResultTest(UTILS::Cardiovascular0DManager& cardvasc0dman, Teuchos::RCP<DRT::Discretization> discr)
  : DRT::ResultTest("CARDIOVASCULAR0D"),
    actdisc_(discr),
    cardvasc0d_dof_(cardvasc0dman.Get0D_dof_m()), // cardiovascular 0D dofs at generalized mid-point t_{n+\theta}
    havecardio_4elementwindkessel_(cardvasc0dman.GetCardvasc0D4ElementWindkessel()->HaveCardiovascular0D()),
    havecardio_arterialproxdist_(cardvasc0dman.GetCardvasc0DArterialProxDist()->HaveCardiovascular0D()),
    havecardio_syspulcirculation_(cardvasc0dman.GetCardvasc0DSysPulCirculation()->HaveCardiovascular0D())
{
  //empty
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Cardiovascular0DResultTest::TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  std::string quantity;
  res.ExtractString("QUANTITY",quantity);
  bool unknownquantity = true; // make sure the result value std::string can be handled
  double result = 0.0;    // will hold the actual result of run


  const Epetra_BlockMap& cardvasc0dmap = cardvasc0d_dof_->Map();
  const int offset = cardvasc0dmap.MinAllGID();

  bool havegid = false;

  if (havecardio_4elementwindkessel_)
    dserror("Testing not implemented for 4ElementWindkessel model!");

  if (havecardio_arterialproxdist_)
    dserror("Testing not implemented for ArterialProxDist model!");

  if (havecardio_syspulcirculation_)
  {
    // test for left atrial pressure
    if ( quantity == "p_at_l" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+0)];
      havegid = cardvasc0dmap.MyGID(offset+0);
    }
    // test for left ventricular in-flux
    if ( quantity == "q_vin_l" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+1)];
      havegid = cardvasc0dmap.MyGID(offset+1);
    }
    // test for left ventricular out-flux
    if ( quantity == "q_vout_l" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+2)];
      havegid = cardvasc0dmap.MyGID(offset+2);
    }
    // test for left ventricular pressure
    if ( quantity == "p_v_l" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+3)];
      havegid = cardvasc0dmap.MyGID(offset+3);
    }
    // test for systemic arterial pressure
    if ( quantity == "p_ar_sys" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+4)];
      havegid = cardvasc0dmap.MyGID(offset+4);
    }
    // test for systemic arterial flux
    if ( quantity == "q_ar_sys" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+5)];
      havegid = cardvasc0dmap.MyGID(offset+5);
    }
    // test for systemic venous pressure
    if ( quantity == "p_ven_sys" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+6)];
      havegid = cardvasc0dmap.MyGID(offset+6);
    }
    // test for systemic venous flux
    if ( quantity == "q_ven_sys" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+7)];
      havegid = cardvasc0dmap.MyGID(offset+7);
    }


    // test for right atrial pressure
    if ( quantity == "p_at_r" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+8)];
      havegid = cardvasc0dmap.MyGID(offset+8);
    }
    // test for right ventricular in-flux
    if ( quantity == "q_vin_r" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+9)];
      havegid = cardvasc0dmap.MyGID(offset+9);
    }
    // test for right ventricular out-flux
    if ( quantity == "q_vout_r" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+10)];
      havegid = cardvasc0dmap.MyGID(offset+10);
    }
    // test for right ventricular pressure
    if ( quantity == "p_v_r" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+11)];
      havegid = cardvasc0dmap.MyGID(offset+11);
    }
    // test for pulmonary arterial pressure
    if ( quantity == "p_ar_pul" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+12)];
      havegid = cardvasc0dmap.MyGID(offset+12);
    }
    // test for pulmonary arterial flux
    if ( quantity == "q_ar_pul" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+13)];
      havegid = cardvasc0dmap.MyGID(offset+13);
    }
    // test for pulmonary venous pressure
    if ( quantity == "p_ven_pul" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+14)];
      havegid = cardvasc0dmap.MyGID(offset+14);
    }
    // test for pulmonary venous flux
    if ( quantity == "q_ven_pul" )
    {
      unknownquantity = false;
      result = (*cardvasc0d_dof_)[cardvasc0dmap.LID(offset+15)];
      havegid = cardvasc0dmap.MyGID(offset+15);
    }

    // catch quantity strings, which are not handled by cardiovascular 0D result test
    if ( unknownquantity )
      dserror("Quantity '%s' not supported in cardiovascular 0D testing", quantity.c_str());

    if(havegid)
    {
      // compare values
      const int err = CompareValues(result, "SPECIAL", res);
      nerr += err;
      test_count++;
    }
  }

  return;
}
