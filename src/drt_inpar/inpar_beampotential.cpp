/*----------------------------------------------------------------------*/
/*!
\file inpar_beampotential.cpp

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

  /*----------------------------------------------------------------------*/
  /* parameters for semi-smooth Newton plasticity algorithm */
  Teuchos::ParameterList& iplast = list->sublist("SEMI-SMOOTH PLASTICITY",false,"");

  DoubleParameter("SEMI_SMOOTH_CPL",1.0,"Weighting factor cpl for semi-smooth PDASS",&iplast);
  DoubleParameter("STABILIZATION_S",1.0,"Stabilization factor s for semi-smooth PDASS",&iplast);

  // solver convergence test parameters for semi-smooth plasticity formulation
  setStringToIntegralParameter<int>("NORMCOMBI_RESFPLASTCONSTR","And",
    "binary operator to combine plasticity constraints and residual force values",
    tuple<std::string>(
      "And",
      "Or"),
    tuple<int>(
      INPAR::STR::bop_and,
      INPAR::STR::bop_or),
    &iplast
    );

  setStringToIntegralParameter<int>("NORMCOMBI_DISPPLASTINCR","And",
      "binary operator to combine displacement increments and plastic flow (Delta Lp) increment values",
      tuple<std::string>(
        "And",
        "Or"),
      tuple<int>(
        INPAR::STR::bop_and,
        INPAR::STR::bop_or),
      &iplast
      );

  DoubleParameter("TOLPLASTCONSTR",1.0E-8,
                  "tolerance in the plastic constraint norm for the newton iteration",
                  &iplast);
  DoubleParameter("TOLDELTALP",1.0E-8,
                  "tolerance in the plastic flow (Delta Lp) norm for the Newton iteration",
                  &iplast);

  setStringToIntegralParameter<int>("NORMCOMBI_EASRES","And",
    "binary operator to combine EAS-residual and residual force values",
    tuple<std::string>(
      "And",
      "Or"),
    tuple<int>(
      INPAR::STR::bop_and,
      INPAR::STR::bop_or),
    &iplast
    );

  setStringToIntegralParameter<int>("NORMCOMBI_EASINCR","And",
      "binary operator to combine displacement increments and EAS increment values",
      tuple<std::string>(
        "And",
        "Or"),
      tuple<int>(
        INPAR::STR::bop_and,
        INPAR::STR::bop_or),
      &iplast
      );

  DoubleParameter("TOLEASRES",1.0E-8,
                  "tolerance in the EAS residual norm for the newton iteration",
                  &iplast);
  DoubleParameter("TOLEASINCR",1.0E-8,
                  "tolerance in the EAS increment norm for the Newton iteration",
                  &iplast);

  setStringToIntegralParameter<int>("DISSIPATION_MODE","pl_multiplier",
      "method to calculate the plastic dissipation",
      tuple<std::string>(
        "pl_multiplier",
        "pl_flow"),
      tuple<int>(
        INPAR::TSI::pl_multiplier,
        INPAR::TSI::pl_flow),
      &iplast
      );

}
