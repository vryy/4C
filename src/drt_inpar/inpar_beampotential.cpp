/*----------------------------------------------------------------------*/
/*!
\file inpar_beampotential.cpp
\maintainer Christoph Meier

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/


#include "drt_validparameters.H"
#include "inpar_beampotential.H"
#include "inpar_beamcontact.H"
#include "inpar_structure.H"
#include "inpar_tsi.H"
#include "inpar_parameterlist_utils.H"



void INPAR::BEAMPOTENTIAL::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  /* parameters for potential-based beam interaction */
  Teuchos::ParameterList& beampotential = list->sublist("BEAM POTENTIAL",false,"");

  setNumericStringParameter("POT_LAW_EXPONENT","1.0", "negative(!) exponent(s) m_i of potential law Phi(r) = sum_i (k_i * r^(-m_i)).",&beampotential);
  setNumericStringParameter("POT_LAW_PREFACTOR","0.0", "prefactor(s) k_i of potential law Phi(r) = sum_i (k_i * r^(-m_i)).",&beampotential);
  DoubleParameter("CUTOFFRADIUS",-1.0,"cutoff radius for search of potential-based interaction pairs",&beampotential);

  setStringToIntegralParameter<int>("BEAMPOTENTIAL_TYPE","Surface","Type of potential interaction: surface (default) or volume potential",
       tuple<std::string>("Surface","surface",
                          "Volume", "volume"),
       tuple<int>(
                  beampot_surf,beampot_surf,
                  beampot_vol,beampot_vol),
       &beampotential);

  setStringToIntegralParameter<int>("BEAMPOT_BTSOL","No","decide, whether potential-based interaction between beams and solids is considered",
                               yesnotuple,yesnovalue,&beampotential);

  setStringToIntegralParameter<int>("BEAMPOT_BTSPH","No","decide, whether potential-based interaction between beams and spheres is considered",
                               yesnotuple,yesnovalue,&beampotential);

  // enable octree search and determine type of bounding box (aabb = axis aligned, spbb = spherical)
  setStringToIntegralParameter<int>("BEAMPOT_OCTREE","None","octree and bounding box type for octree search routine",
       tuple<std::string>("None","none","octree_axisaligned","octree_cylorient","octree_spherical"),
       tuple<int>(INPAR::BEAMCONTACT::boct_none,INPAR::BEAMCONTACT::boct_none,
                  INPAR::BEAMCONTACT::boct_aabb,INPAR::BEAMCONTACT::boct_cobb,INPAR::BEAMCONTACT::boct_spbb),
       &beampotential);

  IntParameter("BEAMPOT_TREEDEPTH",6,"max, tree depth of the octree",&beampotential);
  IntParameter("BEAMPOT_BOXESINOCT",8,"max number of bounding boxes in any leaf octant",&beampotential);
}
