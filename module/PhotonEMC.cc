/*!
 *  \file               PhotonEMC.cc
 *  \brief              Record Photon response in EMCal, and Truth info
 *  \author Xudong Yu <xyu3@bnl.gov>
 */
#include "PhotonEMC.h"

#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>

#include <ffaobjects/EventHeaderv1.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <globalvertex/SvtxVertex.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>

#include <CLHEP/Vector/ThreeVector.h>
#include <math.h>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>

//____________________________________________________________________________..
PhotonEMC::PhotonEMC(const std::string &name, const std::string &file):
 SubsysReco(name),
 _outfilename(file),
 _outfile(nullptr),
 _tree(nullptr)
{
  std::cout << "PhotonEMC::PhotonEMC(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
PhotonEMC::~PhotonEMC()
{
  std::cout << "PhotonEMC::~PhotonEMC() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int PhotonEMC::Init(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  std::cout << "PhotonEMC::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  delete _outfile;
  _outfile = new TFile(_outfilename.c_str(), "RECREATE");
  createBranches();

  cnt=0;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void PhotonEMC::createBranches()
{
  delete _tree;
  _tree = new TTree("tree", "A tree with calo reco/truth info");
  _tree->Branch("_runNumber", &_runNumber);
  _tree->Branch("_eventNumber", &_eventNumber);
  _tree->Branch("_emcal_id", &_emcal_id);
  _tree->Branch("_emcal_phi", &_emcal_phi);
  _tree->Branch("_emcal_eta", &_emcal_eta);
  _tree->Branch("_emcal_x", &_emcal_x);
  _tree->Branch("_emcal_y", &_emcal_y);
  _tree->Branch("_emcal_z", &_emcal_z);
  _tree->Branch("_emcal_e", &_emcal_e);
  _tree->Branch("_emcal_ecore", &_emcal_ecore);
  _tree->Branch("_emcal_chi2", &_emcal_chi2);
  _tree->Branch("_emcal_prob", &_emcal_prob);
  _tree->Branch("_emcal_corr_phi", &_emcal_corr_phi);
  _tree->Branch("_emcal_corr_eta", &_emcal_corr_eta);
  _tree->Branch("_emcal_corr_x", &_emcal_corr_x);
  _tree->Branch("_emcal_corr_y", &_emcal_corr_y);
  _tree->Branch("_emcal_corr_z", &_emcal_corr_z);
  _tree->Branch("_emcal_corr_e", &_emcal_corr_e);
  _tree->Branch("_emcal_corr_pt", &_emcal_corr_pt);
  _tree->Branch("_emcal_corr_corephi", &_emcal_corr_corephi);
  _tree->Branch("_emcal_corr_coreeta", &_emcal_corr_coreeta);
  _tree->Branch("_emcal_corr_corex", &_emcal_corr_corex);
  _tree->Branch("_emcal_corr_corey", &_emcal_corr_corey);
  _tree->Branch("_emcal_corr_corez", &_emcal_corr_corez);
  _tree->Branch("_emcal_corr_ecore", &_emcal_corr_ecore);
  _tree->Branch("_emcal_corr_corept", &_emcal_corr_corept);
  _tree->Branch("_emcal_tower_cluster_id", &_emcal_tower_cluster_id);
  _tree->Branch("_emcal_tower_e", &_emcal_tower_e);
  _tree->Branch("_emcal_tower_phi", &_emcal_tower_phi);
  _tree->Branch("_emcal_tower_eta", &_emcal_tower_eta);
  _tree->Branch("_emcal_tower_status", &_emcal_tower_status);
  _tree->Branch("_truth_vtx_x", &_truth_vtx_x);
  _tree->Branch("_truth_vtx_y", &_truth_vtx_y);
  _tree->Branch("_truth_vtx_z", &_truth_vtx_z);

  //truth info
  _tree->Branch("_truth_px", &_truth_px);
  _tree->Branch("_truth_py", &_truth_py);
  _tree->Branch("_truth_pz", &_truth_pz);
  _tree->Branch("_truth_e", &_truth_e);
  _tree->Branch("_truth_pt", &_truth_pt);
  _tree->Branch("_truth_eta", &_truth_eta);
  _tree->Branch("_truth_phi", &_truth_phi);
  _tree->Branch("_truth_pid", &_truth_pid);

  _tree->Branch("_CEMC_Hit_Evis", &_CEMC_Hit_Evis);
  _tree->Branch("_CEMC_Hit_Edep", &_CEMC_Hit_Edep);
  _tree->Branch("_CEMC_Hit_ch", &_CEMC_Hit_ch);
  _tree->Branch("_CEMC_Hit_x", &_CEMC_Hit_x);
  _tree->Branch("_CEMC_Hit_y", &_CEMC_Hit_y);
  _tree->Branch("_CEMC_Hit_z", &_CEMC_Hit_z);
  _tree->Branch("_CEMC_Hit_t", &_CEMC_Hit_t);
  _tree->Branch("_CEMC_Hit_particle_x", &_CEMC_Hit_particle_x);
  _tree->Branch("_CEMC_Hit_particle_y", &_CEMC_Hit_particle_y);
  _tree->Branch("_CEMC_Hit_particle_z", &_CEMC_Hit_particle_z);

}

//____________________________________________________________________________..
int PhotonEMC::process_event(PHCompositeNode *topNode)
{
  std::cout<<"PhotonEMC::process_event event "<<cnt<<std::endl;
  cnt++;
  PHNodeIterator nodeIter(topNode);
  PHNode* evtNode = dynamic_cast<PHNode*>(nodeIter.findFirst("EventHeader"));
  if (evtNode)
  {
    EventHeaderv1* evtHeader = findNode::getClass<EventHeaderv1>(topNode, "EventHeader");
    std::cout<<"runNumber = "<<evtHeader->get_RunNumber()<<" , m_evtNumber = "<<evtHeader->get_EvtSequence()<<std::endl;
    _runNumber = evtHeader->get_RunNumber();
    _eventNumber = evtHeader->get_EvtSequence();
  }
  else
  {
    _runNumber = 0;
    _eventNumber = -1;
  }

  if (!clustersEM)
  {
    clustersEM = findNode::getClass<RawClusterContainer>(topNode, m_RawClusCont_EM_name);
    if (!clustersEM)
    {
      std::cout << "PhotonEMC::process_event: cannot find cluster container " << m_RawClusCont_EM_name << std::endl;
    }
  }

  if(!EMCAL_Container)
  {
    EMCAL_Container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
    if(!EMCAL_Container)
    {
      std::cout << "PhotonEMC::process_event: TOWERINFO_CALIB_CEMC not found!!!" << std::endl;
    }
  }

  if(!EMCalGeo)
  {
    EMCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if(!EMCalGeo)
    {
      std::cout << "PhotonEMC::process_event: TOWERGEOM_CEMC not found!!!" << std::endl;
    }
  }

  if(!IHCalGeo)
  {
    IHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    if(!IHCalGeo)
    {
      std::cout << "PhotonEMC::process_event: TOWERGEOM_HCALIN not found!!!" << std::endl;
    }
  }

  if(!OHCalGeo)
  {
    OHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    if(!OHCalGeo)
    {
      std::cout << "PhotonEMC::process_event: TOWERGEOM_HCALOUT not found!!!" << std::endl;
    }
  }

  if(!vertexmap)
  {
    vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if(!vertexmap)
    {
      std::cout << "PhotonEMC::process_event: GlobalVertexMap not found!!! (but not necessary)" << std::endl;
    }
  }

  if (!truthinfo)
  {
    truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if(!truthinfo)
    {
      std::cout << "PhotonEMC::process_event: G4TruthInfo not found!!!" << std::endl;
    }
  }

  if (!hits_CEMC)
  {
    hits_CEMC = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_CEMC");
    if(!hits_CEMC)
    {
      std::cout << "PhotonEMC::process_event: G4HIT_CEMC not found!!!" << std::endl;
    }
  }


  ResetTreeVectors();
  fillTree();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void PhotonEMC::fillTree()
{
  if (!clustersEM || !EMCAL_Container || !EMCalGeo || !IHCalGeo || !OHCalGeo || !hits_CEMC || !truthinfo)
  {
    std::cout << PHWHERE << "missing node trees, can't continue with track calo matching"
              << std::endl;
    return;
  }

  CLHEP::Hep3Vector vertex(0., 0., 0.);

  PHG4TruthInfoContainer::VtxRange vtxrange = truthinfo->GetPrimaryVtxRange();
  for (PHG4TruthInfoContainer::ConstVtxIterator iter = vtxrange.first; iter != vtxrange.second; ++iter)
  {
    PHG4VtxPoint *vtx = iter->second;
    _truth_vtx_x = vtx->get_x();
    _truth_vtx_y = vtx->get_y();
    _truth_vtx_z = vtx->get_z();
    vertex.setX(_truth_vtx_x);
    vertex.setY(_truth_vtx_y);
    vertex.setZ(_truth_vtx_z);
//std::cout << " _truth_vtx_x = " << _truth_vtx_x << " _truth_vtx_y = " << _truth_vtx_y << " _truth_vtx_z = " << _truth_vtx_z << std::endl;
  }

  double caloRadiusEMCal;
  double caloRadiusIHCal;
  double caloRadiusOHCal;
  if (m_use_emcal_radius)
  {
    caloRadiusEMCal = m_emcal_radius_user;
  }
  else
  {
    caloRadiusEMCal = EMCalGeo->get_radius();
  }
  if (m_use_ihcal_radius)
  {
    caloRadiusIHCal = m_ihcal_radius_user;
  }
  else
  {
    caloRadiusIHCal = IHCalGeo->get_radius();
  }
  if (m_use_ohcal_radius)
  {
    caloRadiusOHCal = m_ohcal_radius_user;
  }
  else
  {
    caloRadiusOHCal = OHCalGeo->get_radius();
  }

  RawCluster *cluster = nullptr;

  RawClusterContainer::Range begin_end_EMC = clustersEM->getClusters();
  RawClusterContainer::Iterator clusIter_EMC;

  /// Loop over the EMCal clusters
  for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
  {
    cluster = clusIter_EMC->second;
    if(cluster->get_energy() < m_emcal_e_low_cut) continue;

    _emcal_id.push_back(clusIter_EMC->first);
    _emcal_e.push_back(cluster->get_energy());
    _emcal_phi.push_back(RawClusterUtility::GetAzimuthAngle(*cluster, vertex));
    _emcal_eta.push_back(RawClusterUtility::GetPseudorapidity(*cluster, vertex));
    _emcal_x.push_back(cluster->get_x());
    _emcal_y.push_back(cluster->get_y());
    _emcal_z.push_back(cluster->get_z());
    _emcal_ecore.push_back(cluster->get_ecore());
    _emcal_chi2.push_back(cluster->get_chi2());
    _emcal_prob.push_back(cluster->get_prob());

    //std::cout<<"cluster->get_energy() = "<<cluster->get_energy()<<" , phi = "<<RawClusterUtility::GetAzimuthAngle(*cluster, vertex)<<" , eta = "<<RawClusterUtility::GetPseudorapidity(*cluster, vertex)<<" , x = "<<cluster->get_x()<<" , y = "<<cluster->get_y()<<" , z = "<<cluster->get_z()<<" , ecore = "<<cluster->get_ecore()<<" , chi2 = "<<cluster->get_chi2()<<" , prob = "<<cluster->get_prob()<<std::endl;

    CLHEP::Hep3Vector E_vec_cluster_Full = RawClusterUtility::GetEVec(*cluster, vertex);
    float clusE = E_vec_cluster_Full.mag();
    float clus_eta = E_vec_cluster_Full.pseudoRapidity();
    float clus_phi = E_vec_cluster_Full.phi();
    float clus_x = E_vec_cluster_Full.x();
    float clus_y = E_vec_cluster_Full.y();
    float clus_z = E_vec_cluster_Full.z();
    float clus_pt = E_vec_cluster_Full.perp();
    _emcal_corr_e.push_back(clusE);
    _emcal_corr_eta.push_back(clus_eta);
    _emcal_corr_phi.push_back(clus_phi);
    _emcal_corr_x.push_back(clus_x);
    _emcal_corr_y.push_back(clus_y);
    _emcal_corr_z.push_back(clus_z);
    _emcal_corr_pt.push_back(clus_pt);

    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);
    float clusEcore = E_vec_cluster.mag();
    float clus_coreeta = E_vec_cluster.pseudoRapidity();
    float clus_corephi = E_vec_cluster.phi();
    float clus_corex = E_vec_cluster.x();
    float clus_corey = E_vec_cluster.y();
    float clus_corez = E_vec_cluster.z();
    float clus_corept = E_vec_cluster.perp();
    _emcal_corr_ecore.push_back(clusEcore);
    _emcal_corr_coreeta.push_back(clus_coreeta);
    _emcal_corr_corephi.push_back(clus_corephi);
    _emcal_corr_corex.push_back(clus_corex);
    _emcal_corr_corey.push_back(clus_corey);
    _emcal_corr_corez.push_back(clus_corez);
    _emcal_corr_corept.push_back(clus_corept);

//std::cout<<"clusE = "<<clusE<<" , clus_eta = "<<clus_eta<<" , clus_phi = "<<clus_phi<<" , clus_pt = "<<clus_pt<<" , x = "<<E_vec_cluster_Full.x()<<" , y = "<<E_vec_cluster_Full.y()<<" , z = "<<E_vec_cluster_Full.z()<<std::endl;
//std::cout<<"clusEcore = "<<clusEcore<<" , clus_coreeta = "<<clus_coreeta<<" , clus_corephi = "<<clus_corephi<<" , clus_corept = "<<clus_corept<<" , x = "<<E_vec_cluster.x()<<" , y = "<<E_vec_cluster.y()<<" , z = "<<E_vec_cluster.z()<<std::endl;

    RawCluster::TowerConstRange towers = cluster->get_towers();
    RawCluster::TowerConstIterator toweriter;

    TowerInfo *towerInfo = nullptr;

    for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
    {
      _emcal_tower_cluster_id.push_back(clusIter_EMC->first);
      RawTowerGeom *tower_geom = EMCalGeo->get_tower_geometry(toweriter->first);
      _emcal_tower_phi.push_back(tower_geom->get_phi());
      _emcal_tower_eta.push_back(tower_geom->get_eta());
      _emcal_tower_e.push_back(toweriter->second);
      unsigned int key = TowerInfoDefs::encode_emcal(tower_geom->get_bineta(), tower_geom->get_binphi());
      towerInfo = EMCAL_Container->get_tower_at_key(key);
      _emcal_tower_status.push_back(towerInfo->get_status());

    }
  }

  // loop over truth primary particles
  PHG4TruthInfoContainer::ConstRange range = truthinfo->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first; truth_itr != range.second; ++truth_itr)
  {
    PHG4Particle *truth = truth_itr->second;
    _truth_px.push_back(truth->get_px());
    _truth_py.push_back(truth->get_py());
    _truth_pz.push_back(truth->get_pz());
    _truth_e.push_back(truth->get_e());
    _truth_pt.push_back(sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py()));
    _truth_eta.push_back(atanh(truth->get_pz() / sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz())));
    _truth_phi.push_back(atan2(truth->get_py(),truth->get_px()));
    _truth_pid.push_back(truth->get_pid());
//std::cout<<"pid = "<<truth->get_pid()<<std::endl;
  }

  // map between track id and charge deposit
  std::map<int, float> trkmap;

  TLorentzVector v4;
  PHG4HitContainer::ConstRange hit_range = hits_CEMC->getHits();
  for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
  {
    PHG4Hit *this_hit = hit_iter->second;
    float light_yield = hit_iter->second->get_light_yield();
    float edep = hit_iter->second->get_edep();
    int ch = hit_iter->second->get_layer();
    float x = hit_iter->second->get_x(0);
    float y = hit_iter->second->get_y(0);
    float z = hit_iter->second->get_z(0);
    float t = hit_iter->second->get_t(0);
    int trkid = hit_iter->second->get_trkid();
    int truthtrkid = trkid;
    PHG4Particle *part = truthinfo->GetParticle(trkid);

    v4.SetPxPyPzE(part->get_px(), part->get_py(), part->get_pz(), part->get_e());
    int pid = part->get_pid();
//std::cout<<"g4hit with mother pid = "<<pid<<std::endl;

    int vtxid = part->get_vtx_id();
    PHG4VtxPoint *vtx = truthinfo->GetVtx(vtxid);

    // add trkid to a set
    _CEMC_Hit_Evis.push_back(light_yield);
    _CEMC_Hit_Edep.push_back(edep);
    _CEMC_Hit_ch.push_back(ch);
    _CEMC_Hit_x.push_back(x);
    _CEMC_Hit_y.push_back(y);
    _CEMC_Hit_z.push_back(z);
    _CEMC_Hit_t.push_back(t);
    _CEMC_Hit_particle_x.push_back(vtx->get_x());
    _CEMC_Hit_particle_y.push_back(vtx->get_y());
    _CEMC_Hit_particle_z.push_back(vtx->get_z());

    // print all info
    // std::cout << "CEMC hit: " << ch << " " << x << " " << y << " " << z << " " << light_yield << " " << edep << " " << trkid << std::endl;
  }

  _tree->Fill();

}

//____________________________________________________________________________..
int PhotonEMC::End(PHCompositeNode *topNode)
{
  std::cout << topNode << std::endl;
  _outfile->cd();
  _outfile->Write();
  _outfile->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

void PhotonEMC::ResetTreeVectors()
{
  _emcal_id.clear();
  _emcal_phi.clear();
  _emcal_eta.clear();
  _emcal_x.clear();
  _emcal_y.clear();
  _emcal_z.clear();
  _emcal_e.clear();
  _emcal_ecore.clear();
  _emcal_prob.clear();
  _emcal_chi2.clear();
  _emcal_corr_phi.clear();
  _emcal_corr_eta.clear();
  _emcal_corr_x.clear();
  _emcal_corr_y.clear();
  _emcal_corr_z.clear();
  _emcal_corr_e.clear();
  _emcal_corr_corephi.clear();
  _emcal_corr_coreeta.clear();
  _emcal_corr_corex.clear();
  _emcal_corr_corey.clear();
  _emcal_corr_corez.clear();
  _emcal_corr_ecore.clear();
  _emcal_tower_cluster_id.clear();
  _emcal_tower_e.clear();
  _emcal_tower_phi.clear();
  _emcal_tower_eta.clear();
  _emcal_tower_status.clear();

  _truth_vtx_x = -999;
  _truth_vtx_y = -999;
  _truth_vtx_z = -999;

  _truth_px.clear();
  _truth_py.clear();
  _truth_pz.clear();
  _truth_e.clear();
  _truth_pt.clear();
  _truth_eta.clear();
  _truth_phi.clear();
  _truth_pid.clear();

  _CEMC_Hit_Evis.clear();
  _CEMC_Hit_Edep.clear();
  _CEMC_Hit_ch.clear();
  _CEMC_Hit_x.clear();
  _CEMC_Hit_y.clear();
  _CEMC_Hit_z.clear();
  _CEMC_Hit_t.clear();
  _CEMC_Hit_particle_x.clear();
  _CEMC_Hit_particle_y.clear();
  _CEMC_Hit_particle_z.clear();

}
