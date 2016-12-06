/*----------------------------------------------------------------------*/
/*!
\file inpar_beamcontact.cpp

\brief input parameter definitions for beam contact

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_beamcontact.H"



void INPAR::BEAMCONTACT::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& beamcontact = list->sublist("BEAM CONTACT",false,"");

  setStringToIntegralParameter<int>("BEAMS_STRATEGY","None","Type of employed solving strategy",
        tuple<std::string>("None","none",
                           "Penalty", "penalty",
                           "Gmshonly","gmshonly"),
        tuple<int>(
                bstr_none, bstr_none,
                bstr_penalty, bstr_penalty,
                bstr_gmshonly,bstr_gmshonly),
        &beamcontact);

  setStringToIntegralParameter<int>("BEAMS_NEWGAP","No","choose between original or enhanced gapfunction",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_SEGCON","No","choose between beam contact with and without subsegment generation",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_DEBUG","No","This flag can be used for testing purposes. When it is switched on, some sanity checks are not performed!",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_INACTIVESTIFF","No","Always apply contact stiffness in first Newton step for pairs which have active in last time step",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_BTSOLMT","No","decide, if also meshtying between beams and solids is possible",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_BTSOL","No","decide, if also the contact between beams and solids is possible",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_BTSPH","No","decide, if also the contact between beams and spheres is possible",
                               yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_ENDPOINTPENALTY","No","Additional consideration of endpoint-line and endpoint-endpoint contacts",
                                 yesnotuple,yesnovalue,&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_SMOOTHING","None","Application of smoothed tangent field",
       tuple<std::string>("None","none",
                          "Cpp", "cpp"),
       tuple<int>(
                  bsm_none,bsm_none,
                  bsm_cpp,bsm_cpp),
       &beamcontact);

  setStringToIntegralParameter<int>("BEAMS_DAMPING","No","Application of a contact damping force",
       tuple<std::string>("No","no",
                          "Yes", "yes"),
       tuple<int>(
                  bd_no,bd_no,
                  bd_yes,bd_yes),
       &beamcontact);

  DoubleParameter("BEAMS_BTBPENALTYPARAM",0.0,"Penalty parameter for beam-to-beam point contact",&beamcontact);
  DoubleParameter("BEAMS_BTBLINEPENALTYPARAM",-1.0,"Penalty parameter per unit length for beam-to-beam line contact",&beamcontact);
  DoubleParameter("BEAMS_BTSMTPENALTYPARAM",0.0,"Penalty parameter for beam-to-solid meshtying",&beamcontact);
  DoubleParameter("BEAMS_BTSPENALTYPARAM",0.0,"Penalty parameter for beam-to-solid contact",&beamcontact);
  DoubleParameter("BEAMS_BTSPH_PENALTYPARAM",0.0,"Penalty parameter for beam-to-rigidsphere penalty",&beamcontact);
  DoubleParameter("BEAMS_DAMPINGPARAM",-1000.0,"Damping parameter for contact damping force",&beamcontact);
  DoubleParameter("BEAMS_DAMPREGPARAM1",-1000.0,"First (at gap1, with gap1>gap2) regularization parameter for contact damping force",&beamcontact);
  DoubleParameter("BEAMS_DAMPREGPARAM2",-1000.0,"Second (at gap2, with gap1>gap2) regularization parameter for contact damping force",&beamcontact);
  DoubleParameter("BEAMS_MAXDISISCALEFAC",-1.0,"Scale factor in order to limit maximal iterative displacement increment (resiudal displacement)",&beamcontact);
  DoubleParameter("BEAMS_MAXDELTADISSCALEFAC",1.0,"Scale factor in order to limit maximal displacement per time step",&beamcontact);

  DoubleParameter("BEAMS_PERPSHIFTANGLE1",-1.0,"Lower shift angle (in degrees) for penalty scaling of large-angle-contact",&beamcontact);
  DoubleParameter("BEAMS_PERPSHIFTANGLE2",-1.0,"Upper shift angle (in degrees) for penalty scaling of large-angle-contact",&beamcontact);
  DoubleParameter("BEAMS_PARSHIFTANGLE1",-1.0,"Lower shift angle (in degrees) for penalty scaling of small-angle-contact",&beamcontact);
  DoubleParameter("BEAMS_PARSHIFTANGLE2",-1.0,"Upper shift angle (in degrees) for penalty scaling of small-angle-contact",&beamcontact);
  DoubleParameter("BEAMS_SEGANGLE",-1.0,"Maximal angle deviation allowed for contact search segmentation",&beamcontact);
  IntParameter("BEAMS_NUMINTEGRATIONINTERVAL",1,"Number of integration intervals per element",&beamcontact);

  setStringToIntegralParameter<int>("BEAMS_PENALTYLAW","LinPen","Applied Penalty Law",
       tuple<std::string>("LinPen",
                          "QuadPen",
                          "LinNegQuadPen",
                          "LinPosQuadPen",
                          "LinPosCubPen",
                          "LinPosDoubleQuadPen",
                          "LinPosExpPen"),
       tuple<int>(
                  pl_lp,
                  pl_qp,
                  pl_lnqp,
                  pl_lpqp,
                  pl_lpcp,
                  pl_lpdqp,
                  pl_lpep),
       &beamcontact);

  DoubleParameter("BEAMS_PENREGPARAM_G0",-1.0,"First penalty regularization parameter G0 >=0: For gap<G0 contact is active!",&beamcontact);
  DoubleParameter("BEAMS_PENREGPARAM_F0",-1.0,"Second penalty regularization parameter F0 >=0: F0 represents the force at the transition point between regularized and linear force law!",&beamcontact);
  DoubleParameter("BEAMS_PENREGPARAM_C0",-1.0,"Third penalty regularization parameter C0 >=0: C0 has different physical meanings for the different penalty laws!",&beamcontact);
  DoubleParameter("BEAMS_GAPSHIFTPARAM",0.0,"Parameter to shift penalty law!",&beamcontact);
  DoubleParameter("BEAMS_BASICSTIFFGAP",-1.0,"For gaps > -BEAMS_BASICSTIFFGAP, only the basic part of the contact linearization is applied!",&beamcontact);

  // enable octree search and determine type of bounding box (aabb = axis aligned, cobb = cylindrical oriented)
  setStringToIntegralParameter<int>("BEAMS_OCTREE","None","octree and bounding box type for octree search routine",
       tuple<std::string>("None","none","octree_axisaligned","octree_cylorient","octree_spherical"),
       tuple<int>(boct_none,boct_none,
                  boct_aabb,boct_cobb,boct_spbb),
       &beamcontact);

  setStringToIntegralParameter<int>("BEAMS_ADDITEXT","Yes","Switch between No==multiplicative extrusion factor and Yes==additive extrusion factor",
                               yesnotuple,yesnovalue,&beamcontact);
  setNumericStringParameter("BEAMS_EXTVAL","-1.0", "extrusion value(s) of the bounding box, Depending on BEAMS_ADDITIVEEXTFAC is either additive or multiplicative. Give one or two values.",&beamcontact);
  IntParameter("BEAMS_TREEDEPTH",6,"max, tree depth of the octree",&beamcontact);
  IntParameter("BEAMS_BOXESINOCT",8,"max number of bounding boxes in any leaf octant",&beamcontact);


}
