/*----------------------------------------------------------------------------*/
/*! \file
\brief shell8

\level 1

\maintainer Christoph Meier

*/
/*---------------------------------------------------------------------------*/
#ifdef D_SHELL8

#include "shell8.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Shell8::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  linedef->ExtractDouble("THICK", thickness_);

  std::vector<int> gp;
  linedef->ExtractIntVector("GP", gp);
  std::copy(gp.begin(), gp.end(), ngp_);

  linedef->ExtractInt("GP_TRI", ngptri_);

  std::string buffer;
  linedef->ExtractString("FORCES", buffer);

  if (buffer == "XYZ")
    forcetype_ = s8_xyz;
  else if (buffer == "RST")
    forcetype_ = s8_rst;
  else if (buffer == "RST_ortho")
    forcetype_ = s8_rst_ortho;
  else
    dserror("Reading of SHELL8 element failed");

  linedef->ExtractString("EAS", buffer);
  if (buffer == "none")
    eas_[0] = 0;
  else if (buffer == "N4_1")
    eas_[0] = 1;
  else if (buffer == "N4_2")
    eas_[0] = 2;
  else if (buffer == "N4_3")
    eas_[0] = 3;
  else if (buffer == "N4_4")
    eas_[0] = 4;
  else if (buffer == "N4_5")
    eas_[0] = 5;
  else if (buffer == "N4_7")
    eas_[0] = 7;
  else if (buffer == "N9_7")
    eas_[0] = 7;
  else if (buffer == "N9_9")
    eas_[0] = 9;
  else if (buffer == "N9_11")
    eas_[0] = 11;
  else
    dserror("Illegal eas parameter '%s'", buffer.c_str());

  linedef->ExtractString("EAS2", buffer);
  if (buffer == "none")
    eas_[1] = 0;
  else if (buffer == "N4_4")
    eas_[1] = 4;
  else if (buffer == "N4_5")
    eas_[1] = 5;
  else if (buffer == "N4_6")
    eas_[1] = 6;
  else if (buffer == "N4_7")
    eas_[1] = 7;
  else if (buffer == "N9_9")
    eas_[1] = 9;
  else if (buffer == "N9_11")
    eas_[1] = 11;
  else
    dserror("Illegal eas parameter '%s'", buffer.c_str());

  linedef->ExtractString("EAS3", buffer);
  if (buffer == "none")
    eas_[2] = 0;
  else if (buffer == "N_1")
    eas_[2] = 1;
  else if (buffer == "N_3")
    eas_[2] = 3;
  else if (buffer == "N_4")
    eas_[2] = 4;
  else if (buffer == "N_6")
    eas_[2] = 6;
  else if (buffer == "N_8")
    eas_[2] = 8;
  else if (buffer == "N_9")
    eas_[2] = 9;
  else
    dserror("Illegal eas parameter '%s'", buffer.c_str());

  linedef->ExtractString("EAS4", buffer);
  if (buffer == "none")
    eas_[3] = 0;
  else if (buffer == "N4_2")
    eas_[3] = 2;
  else if (buffer == "N4_4")
    eas_[3] = 4;
  else if (buffer == "N9_2")
    eas_[3] = 2;
  else if (buffer == "N9_4")
    eas_[3] = 4;
  else if (buffer == "N9_6")
    eas_[3] = 6;
  else
    dserror("Illegal eas parameter '%s'", buffer.c_str());

  linedef->ExtractString("EAS5", buffer);
  if (buffer == "none")
    eas_[4] = 0;
  else if (buffer == "N4_2")
    eas_[4] = 2;
  else if (buffer == "N4_4")
    eas_[4] = 4;
  else if (buffer == "N9_2")
    eas_[4] = 2;
  else if (buffer == "N9_4")
    eas_[4] = 4;
  else if (buffer == "N9_6")
    eas_[4] = 6;
  else
    dserror("Illegal eas parameter '%s'", buffer.c_str());

  // count no. eas parameters
  nhyb_ = 0;
  for (int i = 0; i < 5; ++i) nhyb_ += eas_[i];

  // create arrays alfa, alfa_inc, Dtildinv, Lt, Rtild in data_
  std::vector<double> alfa(nhyb_);
  std::vector<double> alfao(nhyb_);
  std::vector<double> alfa_inc(nhyb_);
  std::vector<double> Rtild(nhyb_);
  std::fill(alfa.begin(), alfa.end(), 0);
  std::fill(alfao.begin(), alfao.end(), 0);
  std::fill(Rtild.begin(), Rtild.end(), 0);

  Epetra_SerialDenseMatrix Dtildinv;
  Epetra_SerialDenseMatrix Lt;
  Dtildinv.Shape(nhyb_, nhyb_);
  Lt.Shape(nhyb_, NumNode() * 6);

  data_.Add("alfa", alfa);
  data_.Add("alfao", alfao);
  data_.Add("alfa_inc", alfa_inc);
  data_.Add("Rtild", Rtild);
  data_.Add("Dtildinv", Dtildinv);
  data_.Add("Lt", Lt);

  // read ANS
  linedef->ExtractString("ANS", buffer);
  if (buffer == "none")
    ans_ = 0;
  else if (buffer == "Q")
    ans_ = 1;
  else if (buffer == "T")
    ans_ = 2;
  else if (buffer == "QT")
    ans_ = 3;
  else if (buffer == "TQ")
    ans_ = 3;
  else
    dserror("Illegal ans parameter '%s'", buffer.c_str());

  // read SDC
  linedef->ExtractDouble("SDC", sdc_);

  return true;
}

#endif  // #ifdef D_SHELL8
