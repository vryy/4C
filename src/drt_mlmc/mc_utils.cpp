/*----------------------------------------------------------------------*/
/*! \file
\brief Some utility functions for Monte Carlo analysis

\maintainer Jonas Nitzler

\level 3
*/
/*----------------------------------------------------------------------*/

#include "mlmc.H"
#include "mc_utils.H"
#include "randomfield.H"
#include "randomfield_fourier.H"
#include "randomfield_spectral.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_control.H"

#include "../drt_inpar/inpar_material.H"
#include "../drt_inpar/inpar_mlmc.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../linalg/linalg_utils.H"

#ifdef HAVE_FFTW

void UQ::MLMC::WriteStdVectorToFile(std::vector<double> myvector, std::string FileNameWithPath)
{
  const char* c = FileNameWithPath.c_str();

  if (discret_->Comm().MyPID() == 0)
  {
    std::ofstream File;
    File.open(c, std::ios::out);
    int size = myvector.size();
    for (int i = 0; i < size; i++)
    {
      File << myvector[i] << std::endl;
    }
    File.close();
  }
}


void UQ::ComputePeakAndQuantile(Teuchos::RCP<std::vector<double>> values1,
    Teuchos::RCP<std::vector<double>> values2, Teuchos::RCP<std::vector<int>> ele_ids,
    Teuchos::RCP<std::pair<double, std::pair<int, double>>> peak,
    Teuchos::RCP<std::pair<double, std::pair<int, double>>> quantile99, int tot_num_elements,
    std::string mode, const Epetra_Comm& mycomm)
{
  // Build a list< pair< double sort_value, pair< int EleId, second_value > >
  // and sort it
  std::list<std::pair<double, std::pair<int, double>>> my_quantity;

  // loop over all element ids
  std::pair<double, std::pair<int, double>> thispair;
  for (unsigned int i = 0; i < ele_ids->size(); i++)
  {
    thispair.first = values1->at(i);
    thispair.second.first = ele_ids->at(i);
    thispair.second.second = values2->at(i);
    my_quantity.push_back(thispair);
  }

  // sort the the list in ascending order by the estimated distance
  // this is the STL sorting, which is pretty fast

  if (mode == "min")  // minimum search
  {
    my_quantity.sort(MyComparePairsMin);
  }
  else  // maximum search
  {
    my_quantity.sort(MyComparePairsMax);
  }

  int num_top_ten = (int)(tot_num_elements * 0.1 + 0.5);
  if (num_top_ten < 1)  // to avoid errors when calculating very small models
  {
    num_top_ten = 1;
  }

  if (mode == "min")  // minimum search
  {
    thispair.first = 10000000;
    thispair.second.first = -1;
    thispair.second.second = -10000000;
  }
  else  // maximum search
  {
    thispair.first = -10000000;
    thispair.second.first = -1;
    thispair.second.second = -10000000;
  }
  // resize to top ten fill up the remainder if needed with small number and negative ele id
  my_quantity.resize(num_top_ten, thispair);

  // gather does not like lists, thus
  std::vector<int> my_quantity_id_vec(0);
  std::vector<double> my_quantity_value1_vec(0);
  std::vector<double> my_quantity_value2_vec(0);

  // put data in separate vectors because Gather doesnt like vectors of pairs
  for (std::list<std::pair<double, std::pair<int, double>>>::iterator it = my_quantity.begin();
       it != my_quantity.end(); it++)
  {
    my_quantity_value1_vec.push_back(it->first);
    my_quantity_id_vec.push_back(it->second.first);
    my_quantity_value2_vec.push_back(it->second.second);
  }

  // vector to gather all the data in
  std::vector<double> my_quantity_values1_gathered_vec(0);
  std::vector<int> my_quantity_id_gathered_vec(0);
  std::vector<double> my_quantity_values2_gathered_vec(0);

  // gather information across all processors
  const int my_numprocs = mycomm.NumProc();

  // information how many processors participate in total
  std::vector<int> allproc(mycomm.NumProc());
  for (int i = 0; i < mycomm.NumProc(); ++i) allproc[i] = i;

  LINALG::Gather<int>(
      my_quantity_id_vec, my_quantity_id_gathered_vec, my_numprocs, &allproc[0], mycomm);
  LINALG::Gather<double>(
      my_quantity_value1_vec, my_quantity_values1_gathered_vec, my_numprocs, &allproc[0], mycomm);
  LINALG::Gather<double>(
      my_quantity_value2_vec, my_quantity_values2_gathered_vec, my_numprocs, &allproc[0], mycomm);

  // put data back into list for sorting
  std::list<std::pair<double, std::pair<int, double>>> my_quantity_gathered;
  for (unsigned int i = 0; i < my_quantity_id_gathered_vec.size(); i++)
  {
    std::pair<double, std::pair<int, double>> mytemp;
    mytemp.first = my_quantity_values1_gathered_vec[i];
    mytemp.second.first = my_quantity_id_gathered_vec[i];
    mytemp.second.second = my_quantity_values2_gathered_vec[i];
    my_quantity_gathered.push_back(mytemp);
  }
  // sort the list
  if (mode == "min")  // minimum search
  {
    my_quantity_gathered.sort(MyComparePairsMin);
  }
  else  // maximum search
  {
    my_quantity_gathered.sort(MyComparePairsMax);
  }

  peak->first = my_quantity_gathered.begin()->first;
  peak->second.first = my_quantity_gathered.begin()->second.first;
  peak->second.second = my_quantity_gathered.begin()->second.second;

  std::list<std::pair<double, std::pair<int, double>>>::iterator it = my_quantity_gathered.begin();
  std::advance(it, (int)(tot_num_elements * 0.01));
  quantile99->first = it->first;
  quantile99->second.first = it->second.first;
  quantile99->second.second = it->second.second;
}


#endif
