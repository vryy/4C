/*----------------------------------------------------------------------*/
/*! \file
\brief Unit tests for input parameters of structure field.

\level 0
*/
/*----------------------------------------------------------------------*/
#ifndef UNIT_INPAR_STRUCTURE_H
#define UNIT_INPAR_STRUCTURE_H

#include "src/common/unit_cxx_test_wrapper.H"

#include "src/drt_inpar/inpar_structure.H"

// Declaration
namespace INPAR
{
  namespace STR
  {
    class InparStructure_TestSuite;
  }
}  // namespace INPAR

// Implementation
class INPAR::STR::InparStructure_TestSuite : public BACICxxTestWrapper
{
 public:
  void TestString2ModelType()
  {
    TS_ASSERT_EQUALS(INPAR::STR::String2ModelType("Structure"), INPAR::STR::model_structure);
    TS_ASSERT_EQUALS(INPAR::STR::String2ModelType("Contact"), INPAR::STR::model_contact);
    TS_ASSERT_EQUALS(INPAR::STR::String2ModelType("Meshtying"), INPAR::STR::model_meshtying);
    TS_ASSERT_EQUALS(
        INPAR::STR::String2ModelType("SpringDashpot"), INPAR::STR::model_springdashpot);
    TS_ASSERT_EQUALS(
        INPAR::STR::String2ModelType("Lag-Pen-Constraint"), INPAR::STR::model_lag_pen_constraint);
    TS_ASSERT_EQUALS(INPAR::STR::String2ModelType("Partitioned-Coupling"),
        INPAR::STR::model_partitioned_coupling);
    TS_ASSERT_EQUALS(
        INPAR::STR::String2ModelType("BeamInteractionOld"), INPAR::STR::model_beam_interaction_old);
    TS_ASSERT_EQUALS(INPAR::STR::String2ModelType("BrownianDyn"), INPAR::STR::model_browniandyn);
    TS_ASSERT_EQUALS(
        INPAR::STR::String2ModelType("BeamInteraction"), INPAR::STR::model_beaminteraction);
  }

  void TestModelTypeString()
  {
    TS_ASSERT_EQUALS(INPAR::STR::ModelTypeString(INPAR::STR::model_structure), "Structure");
    TS_ASSERT_EQUALS(INPAR::STR::ModelTypeString(INPAR::STR::model_contact), "Contact");
    TS_ASSERT_EQUALS(INPAR::STR::ModelTypeString(INPAR::STR::model_meshtying), "Meshtying");
    TS_ASSERT_EQUALS(INPAR::STR::ModelTypeString(INPAR::STR::model_springdashpot), "SpringDashpot");
    TS_ASSERT_EQUALS(
        INPAR::STR::ModelTypeString(INPAR::STR::model_lag_pen_constraint), "Lag-Pen-Constraint");
    TS_ASSERT_EQUALS(INPAR::STR::ModelTypeString(INPAR::STR::model_partitioned_coupling),
        "Partitioned-Coupling");
    TS_ASSERT_EQUALS(
        INPAR::STR::ModelTypeString(INPAR::STR::model_beam_interaction_old), "BeamInteractionOld");
    TS_ASSERT_EQUALS(INPAR::STR::ModelTypeString(INPAR::STR::model_browniandyn), "BrownianDyn");
    TS_ASSERT_EQUALS(
        INPAR::STR::ModelTypeString(INPAR::STR::model_beaminteraction), "BeamInteraction");
  }
};

#endif /* UNIT_INPAR_STRUCTURE_H */
