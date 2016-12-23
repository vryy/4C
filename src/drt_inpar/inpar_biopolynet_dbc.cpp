/*-----------------------------------------------------------*/
/*!
\file inpar_biopolynet_dbc.cpp

\brief input parameter for statistical mechanic problem

\maintainer Jonas Eichinger, Maximilian Grill

\level 2

*/
/*-----------------------------------------------------------*/


#include "drt_validparameters.H"
#include "inpar_biopolynet_dbc.H"


void INPAR::BIOPOLYNETDBC::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& biopolynetdbc = list->sublist("BIOPOLYNET DBC",false,"");



  //Reading which kind of Dirichlet boundary condition should be applied
  setStringToIntegralParameter<int>("DBCTYPE","std","Dirichlet BC type applied",
                                 //listing possible std::strings in input file in category DBCTYPE
                                 tuple<std::string>("none",
                                                    "std",
                                                    "shearfixed",
                                                    "shearfixeddel",
                                                    "sheartrans",
                                                    "pinnodes" ,
                                                    "affineshear",
                                                    "affinesheardel",
                                                    "movablesupport1d"),
                                 //translating input std::strings into BACI input parameters
                                 tuple<int>(dbctype_none,
                                            dbctype_std,
                                            dbctype_shearfixed,
                                            dbctype_shearfixeddel,
                                            dbctype_sheartrans,
                                            dbctype_pinnodes,
                                            dbctype_affineshear,
                                            dbctype_affinesheardel,
                                            dbctype_movablesupport1d),
                                            &biopolynetdbc);

   //Reading double parameter for shear flow field
  DoubleParameter("SHEARAMPLITUDE",0.0,"Shear amplitude of flow in z-direction; note: not amplitude of displacement, but of shear strain!",&biopolynetdbc);
  //Reading direction of oscillatory motion that DBC nodes are subjected to (we need this when using periodic BCs)
  IntParameter("DBCDISPDIR",0,"Global spatial direction of oscillatory motion by Dirichlet BCs",&biopolynetdbc);
  //Reading time curve number for Dirichlet boundary conditions
  IntParameter("CURVENUMBER",0,"Specifies Time Curve number of imposed Dirichlet BCs",&biopolynetdbc);


}
