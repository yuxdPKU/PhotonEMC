// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 *  \file               PhotonEMC.h
 *  \brief              Record Photon response in EMCal, and Truth info
 *  \author Xudong Yu <xyu3@bnl.gov>
 */

#ifndef PHOTONEMC_H
#define PHOTONEMC_H

#include <fun4all/SubsysReco.h>

#include <globalvertex/GlobalVertexMap.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawTowerGeomContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4VtxPoint.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TH1;
class TH2;
class TFile;
class TTree;

class PhotonEMC : public SubsysReco
{
 public:

  PhotonEMC(const std::string &name = "PhotonEMC", const std::string &file = "output.root");

  ~PhotonEMC() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void ResetTreeVectors();

  void fillTree();

  void createBranches();

  void EMcalRadiusUser(bool use) {m_use_emcal_radius = use;}
  void IHcalRadiusUser(bool use) {m_use_ihcal_radius = use;}
  void OHcalRadiusUser(bool use) {m_use_ohcal_radius = use;}
  void setEMcalRadius(float r) {m_emcal_radius_user = r;}
  void setIHcalRadius(float r) {m_ihcal_radius_user = r;}
  void setOHcalRadius(float r) {m_ohcal_radius_user = r;}

  void setRawClusContEMName(std::string name) {m_RawClusCont_EM_name = name;}

  void setEmcalELowCut(float e) {m_emcal_e_low_cut = e;}

 private:
   int cnt = 0;
   bool m_use_emcal_radius = false;
   bool m_use_ihcal_radius = false;
   bool m_use_ohcal_radius = false;
   float m_emcal_radius_user = 93.5;
   float m_ihcal_radius_user = 117;
   float m_ohcal_radius_user = 177.423;
   std::string _outfilename;
   TFile *_outfile = nullptr;
   TTree *_tree = nullptr;
   TTree *_tree_KFP = nullptr;

   std::string m_RawClusCont_EM_name = "CLUSTERINFO_CEMC";

   int _runNumber;
   int _eventNumber;

   std::vector<int> _emcal_id;
   std::vector<float> _emcal_phi;
   std::vector<float> _emcal_eta;
   std::vector<float> _emcal_x;
   std::vector<float> _emcal_y;
   std::vector<float> _emcal_z;
   std::vector<float> _emcal_e;
   std::vector<float> _emcal_ecore;
   std::vector<float> _emcal_chi2;
   std::vector<float> _emcal_prob;

   std::vector<float> _emcal_corr_phi;
   std::vector<float> _emcal_corr_eta;
   std::vector<float> _emcal_corr_x;
   std::vector<float> _emcal_corr_y;
   std::vector<float> _emcal_corr_z;
   std::vector<float> _emcal_corr_e;
   std::vector<float> _emcal_corr_pt;

   std::vector<float> _emcal_corr_corephi;
   std::vector<float> _emcal_corr_coreeta;
   std::vector<float> _emcal_corr_corex;
   std::vector<float> _emcal_corr_corey;
   std::vector<float> _emcal_corr_corez;
   std::vector<float> _emcal_corr_ecore;
   std::vector<float> _emcal_corr_corept;

   std::vector<int> _emcal_tower_cluster_id;
   std::vector<float> _emcal_tower_e;
   std::vector<float> _emcal_tower_phi;
   std::vector<float> _emcal_tower_eta;
   std::vector<int> _emcal_tower_status;

   float _truth_vtx_x;
   float _truth_vtx_y;
   float _truth_vtx_z;

   std::vector<float> _truth_px;
   std::vector<float> _truth_py;
   std::vector<float> _truth_pz;
   std::vector<float> _truth_e;
   std::vector<float> _truth_pt;
   std::vector<float> _truth_eta;
   std::vector<float> _truth_phi;
   std::vector<int> _truth_pid;

   std::vector<float> _CEMC_Hit_Edep;
   std::vector<float> _CEMC_Hit_Evis;
   std::vector<float> _CEMC_Hit_ch;
   std::vector<float> _CEMC_Hit_x;
   std::vector<float> _CEMC_Hit_y;
   std::vector<float> _CEMC_Hit_z;
   std::vector<float> _CEMC_Hit_t;
   std::vector<float> _CEMC_Hit_particle_x;
   std::vector<float> _CEMC_Hit_particle_y;
   std::vector<float> _CEMC_Hit_particle_z;

   GlobalVertexMap *vertexmap = nullptr;
   RawClusterContainer *clustersEM = nullptr;
   RawClusterContainer *EMCAL_RawClusters = nullptr;
   TowerInfoContainer *EMCAL_Container = nullptr;
   RawTowerGeomContainer *EMCalGeo = nullptr;
   RawTowerGeomContainer *IHCalGeo = nullptr;
   RawTowerGeomContainer *OHCalGeo = nullptr;
   PHG4HitContainer *hits_CEMC = nullptr;
   PHG4TruthInfoContainer *truthinfo = nullptr;

   float m_emcal_e_low_cut = 0.2;
};

#endif // PHOTONEMC_H
