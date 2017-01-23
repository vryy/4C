/*----------------------------------------------------------------------------*/
/*!
\file beam3contactinterface.cpp

\brief interface class for beam contact

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "beam3contactinterface.H"

#include "../drt_beaminteraction/beam3contact_defines.H"
#include "beam3contactnew.H"
#include "beam3contact.H"
#include "../drt_inpar/inpar_beamcontact.H"

Teuchos::RCP<CONTACT::Beam3contactinterface> CONTACT::Beam3contactinterface::Impl( const int numnodes,
                                                                      const int numnodalvalues,
                                                                      const DRT::Discretization& pdiscret,
                                                                      const DRT::Discretization& cdiscret,
                                                                      const std::map<int,int>& dofoffsetmap,
                                                                      DRT::Element* element1,
                                                                      DRT::Element* element2,
                                                                      Teuchos::ParameterList& beamcontactparams)
{
  //Decide, if beam contact with subsegment creation (beam3contact) or pure element based beam contact (beam3contactnew) should be applied
  bool beamssegcon = DRT::INPUT::IntegralValue<int>(beamcontactparams,"BEAMS_SEGCON");

  // note: numnodes is to be interpreted as number of nodes used for centerline interpolation.

  if(!beamssegcon)
  {
    switch (numnodalvalues)
    {
      case 1:
      {
        switch (numnodes)
        {
          case 2:
          {
            return Teuchos::rcp (new CONTACT::Beam3contactnew<2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
          }
          case 3:
          {
            return Teuchos::rcp (new CONTACT::Beam3contactnew<3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
          }
          case 4:
          {
            return Teuchos::rcp (new CONTACT::Beam3contactnew<4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
          }
          case 5:
          {
            return Teuchos::rcp (new CONTACT::Beam3contactnew<5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
          }
          default:
            dserror("No valid template parameter for the number of nodes (numnodes = 2,3,4,5 for Reissner beams) available!");
            break;
        }
        break;
      }
      case 2:
      {
        switch (numnodes)
        {
          case 2:
          {
            return Teuchos::rcp (new CONTACT::Beam3contactnew<2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
          }
          default:
            dserror("No valid template parameter combination for the number of nodes and number of types of nodal DoFs"
                    "(only numnodes = 2 in combination with numnodalvalues=2 possible so far, i.e. 3rd order Hermite interpolation)!");
            break;
        }
        break;
      }
      default:
        dserror("No valid template parameter for the number of types of nodal DoFs used for centerline interpolation!\n"
                "(numnodalvalues = 1, i.e. positions              for Lagrange interpolation,\n"
                " numnodalvalues = 2, i.e. positions AND tangents for Hermite interpolation)");
        break;
    }
  }
  else
  {
    switch (numnodalvalues)
    {
      case 1:
      {
        switch (numnodes)
        {
          case 2:
          {
            return Teuchos::rcp (new CONTACT::Beam3contact<2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
          }
          case 3:
          {
            return Teuchos::rcp (new CONTACT::Beam3contact<3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
          }
          case 4:
          {
            return Teuchos::rcp (new CONTACT::Beam3contact<4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
          }
          case 5:
          {
            return Teuchos::rcp (new CONTACT::Beam3contact<5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
          }
          default:
            dserror("No valid template parameter for the number of nodes (numnodes = 2,3,4,5 for Reissner beams) available!");
            break;
        }
        break;
      }
      case 2:
      {
        switch (numnodes)
        {
          case 2:
          {
            return Teuchos::rcp (new CONTACT::Beam3contact<2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
          }
          default:
            dserror("No valid template parameter combination for the number of nodes and number of types of nodal DoFs"
                    "(only numnodes = 2 in combination with numnodalvalues=2 possible so far, i.e. 3rd order Hermite interpolation)!");
            break;
        }
        break;
      }
      default:
        dserror("No valid template parameter for the number of types of nodal DoFs used for centerline interpolation!\n"
                "(numnodalvalues = 1, i.e. positions              for Lagrange interpolation,\n"
                " numnodalvalues = 2, i.e. positions AND tangents for Hermite interpolation)");
        break;
    }
  }
  return Teuchos::null;

}
