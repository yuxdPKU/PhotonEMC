#include <G4_ActsGeom.C>
#include <G4_Magnet.C>
#include <GlobalVariables.C>
#include <G4_Global.C>

#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <cdbobjects/CDBTTree.h>

#include <calobase/RawTowerDefs.h>

#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawTowerCalibration.h>

#include <photonemc/PhotonEMC.h>

#include <stdio.h>
#include <iostream>
#include <filesystem>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libphotonemc.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libg4calo.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4eval.so)

using namespace std;

namespace fs = std::filesystem;

std::string GetFirstLine(std::string listname);
bool is_directory_empty(const fs::path& dir_path);

void Fun4All_CaloANA(
    const int nEvents = 10,
    vector<string> myInputLists = {
        "dst_calo_cluster.list",
        "g4hits.list"}, 
    bool doEMcalRadiusCorr = false)
{
  int verbosity = 0;

  std::cout << "Including " << myInputLists.size() << " files." << std::endl;
  std::string firstfile = GetFirstLine(myInputLists[0]);
  if (*firstfile.c_str() == '\0') return;
  std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(firstfile);
  int runnumber = runseg.first;
  int segment = runseg.second;

  //The next set of lines figures out folder revisions, file numbers etc
  string outDir = "/sphenix/u/xyu3/workarea/PhotonEMC/macro";

  string outputAnaFileName = "TrackCalo_" + to_string(segment) + "_ana.root";

  string outputRecoDir = outDir + "/inReconstruction/" + to_string(runnumber) + "/";
  string makeDirectory = "mkdir -p " + outputRecoDir;
  system(makeDirectory.c_str());
  string outputAnaFile = outputRecoDir + outputAnaFileName;
  std::cout << "Reco ANA file: " << outputAnaFile << std::endl;

  //Create the server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(verbosity);

  //Add all required input files
  for (unsigned int i = 0; i < myInputLists.size(); ++i)
  {
    Fun4AllInputManager *infile = new Fun4AllDstInputManager("DSTin_" + to_string(i));
    std::cout << "Including file " << myInputLists[i] << std::endl;
    //infile->AddFile(myInputLists[i]);
    infile->AddListFile(myInputLists[i]);
    se->registerInputManager(infile);
  }

  auto rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  G4MAGNET::magfield_tracking = G4MAGNET::magfield;
  G4MAGNET::magfield_rescale = 1;

  Global_Reco();

  std::cout << "Building clusters" << std::endl;
  RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  ClusterBuilder->Detector("CEMC");
  ClusterBuilder->set_threshold_energy(0.030);  // for when using basic calibration
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  ClusterBuilder->LoadProfile(emc_prof);
  ClusterBuilder->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
  se->registerSubsystem(ClusterBuilder);

  float new_cemc_rad = 100.70;//(1-(-0.077))*93.5 recommended cemc radius
  PhotonEMC *photonemc = new PhotonEMC("Single_photon_emc_reco", outputAnaFile);
  photonemc->EMcalRadiusUser(doEMcalRadiusCorr);
  photonemc->setEMcalRadius(new_cemc_rad);
  se->registerSubsystem(photonemc);

  se->run(nEvents);
  se->End();
  se->PrintTimer();

  ifstream file_ana(outputAnaFile.c_str(), ios::binary | ios::ate);
  if (file_ana.good() && (file_ana.tellg() > 100))
  {
    string outputRecoDirMove = outDir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputRecoDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + outputAnaFile + " " + outDir + "/Reconstructed/" + to_string(runnumber);
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  delete se;
  std::cout << "All done" << std::endl;
  gSystem->Exit(0);

  return;
}

std::string GetFirstLine(std::string listname)
{
  std::ifstream file(listname);

  std::string firstLine = "";
  if (file.is_open()) {
      if (std::getline(file, firstLine)) {
          std::cout << "First Line: " << firstLine << std::endl;
      } else {
          std::cerr << "Unable to read first line of file" << std::endl;
      }
      file.close();
  } else {
      std::cerr << "Unable to open file" << std::endl;
  }
  return firstLine;
}

bool is_directory_empty(const fs::path& dir_path) {
    if (fs::exists(dir_path) && fs::is_directory(dir_path)) {
        return fs::directory_iterator(dir_path) == fs::directory_iterator();
    }
    return false;
}
