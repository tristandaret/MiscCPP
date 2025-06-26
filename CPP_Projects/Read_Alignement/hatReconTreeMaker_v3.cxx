#include <iostream>
#include <memory>
#include <numeric>
#include <vector>
#include <map>

#include <TND280Event.hxx>
#include <TND280Log.hxx>
#include <nd280EventLoop.hxx>
#include <TOARuntimeParams.hxx>
#include <TParametersOptionManager.hxx>
#include <TReconTrack.hxx>
#include <TReconPID.hxx>
#include <TAlgorithmResult.hxx>
#include <TIntegerDatum.hxx>
#include <THandle.hxx>
#include <TComboHit.hxx>
#include <THitSelection.hxx>
#include <TRealDatum.hxx>
#include <THATUtils.hxx>
#include "TFieldManager.hxx"
#include "HEPUnits.hxx"

#include "TGeomInfo.hxx"
#include "hatReconTreeMaker.hxx"

#include "TTree.h"
#include "TFile.h"

namespace ND {
/// An exception when the TAlgorithm is missing from the loaded library
OA_EXCEPTION(EMissingAlgorithm, EoaCore);
/// An exception when the TData is missing from the event
OA_EXCEPTION(EMissingDatum, EoaCore);
} // namespace ND

#include <TG4Trajectory.hxx>
#include <TrackTruthInfo.hxx>
#include <TrackingUtils.hxx>
#include <fstream>

ND::TG4TrajectoryPoint GetCloserG4Point(ND::TG4Trajectory &traj,
                                        ND::THandle<ND::TTrackState> front, double &dist)
{

   std::vector<ND::TG4TrajectoryPoint> points = traj.GetTrajectoryPoints();

   double x0 = front->GetPosition()[0];
   double y0 = front->GetPosition()[1];
   double z0 = front->GetPosition()[2];

   double dmin = DBL_MAX;
   int icloser = -1;

   for (unsigned int i = 0; i < points.size(); i++) {
      double deltaz = (points[i].GetPosition().Z() - z0);
      double deltay = (points[i].GetPosition().Y() - y0);
      double deltax = (points[i].GetPosition().X() - x0);
      double distance = TMath::Sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz);
      if (distance < dmin) {
         dmin = distance;
         icloser = i;
      }
   }

   if (icloser >= 0) {
      // std::cout << " X " << points[icloser].GetPosition().X() << "  " <<
      // front->GetPosition()[0] << std::endl; std::cout << " Y " <<
      // points[icloser].GetPosition().Y() << "  " << front->GetPosition()[1] <<
      // std::endl; std::cout << " Z " << points[icloser].GetPosition().Z() << "  " <<
      // front->GetPosition()[2] << std::endl;

      dist = dmin;
      return points[icloser];
   } else {
      std::cerr << " No Point found " << std::endl;
      return {};
   }
}

class THATReconTreeMakerLoop : public ND::TND280EventLoopFunction {

public:
   THATReconTreeMakerLoop()
      : fEvent(0), fMomentum(0), fNclusters(0), fCurvature(0), fDirNorm(0)
   {
      fPos[0] = 0;
      fPos[1] = 0;
      fPos[2] = 0;
      fDir[0] = 0;
      fDir[1] = 0;
      fDir[2] = 0;

      fOutputTree = new TTree("outTree", "HAT tracks parameters");
      fOutputTree->Branch("event", &fEvent, "event/I");
      fOutputTree->Branch("mom", &fMomentum, "mom/D");
      fOutputTree->Branch("ncl", &fNclusters, "ncl/I");
      fOutputTree->Branch("resolution", &fResTrack, "resolution/D");
      fOutputTree->Branch("starttime", &fStartTime, "starttime/D");
      fOutputTree->Branch("endtime", &fEndTime, "endtime/D");
      // fOutputTree->Branch("mean_mult", &fMeanMultiplicity,"mean_mult/D");
      fOutputTree->Branch("chi2", &fChi2, "chi2/D");
      fOutputTree->Branch("curv", &fCurvature, "curv/D");
      fOutputTree->Branch("pos", &fPos, "pos[3]/D");
      fOutputTree->Branch("dir", &fDir, "dir[3]/D");
      fOutputTree->Branch("dir_norm", &fDirNorm, "dir_norm/D");
      fOutputTree->Branch("true_momentum", &fTrueMomentum, "true_momentum/D");
      fOutputTree->Branch("pid", &fPID, "pid/I");
      fOutputTree->Branch("pull_prot", &fPullProt, "pull_prot/D");
      fOutputTree->Branch("pull_ele", &fPullElec, "pull_ele/D");
      fOutputTree->Branch("pull_muon", &fPullMuon, "pull_muon/D");
      fOutputTree->Branch("dEdx", &fdEdx, "dEdx/D");
      fOutputTree->Branch("EndPlate_ID", &fEndPlate_ID, "EndPlate_ID/I");
      fOutputTree->Branch("mean_time", &fMeanTime, "mean_time/D");

      treeT = new TTree("treeT", "treeT");
      treeT->Branch("events", "events", &runT, 32000, 6);
      m_counter_Evt = 0;

      fres_his_ep0 = new TH1F("res_ep0", "res_ep0", 200, -10, 10);
      fres_his_ep1 = new TH1F("res_ep1", "res_ep1", 200, -10, 10);
      residuals_col =
         new TH2D("residuals_col", "residuals_col", 1152, 0, 1152, 200, -10, 10);
      residuals_row =
         new TH2D("residuals_row", "residuals_row", 1024, 0, 1024, 200, -10, 10);

      fres_his_ep0->SetTitle("cluster_residuals_EP0;counts;residuals [mm]");
      fres_his_ep1->SetTitle("cluster_residuals_EP1;counts;residuals [mm]");
      residuals_row->SetTitle(
         "residuals vs row; col + 36 * mm + 36 * 8 * EP; residuals [mm]");
      residuals_col->SetTitle(
         "residuals vs col; row + 32 * mm + 32 * 8 * EP; residuals [mm]");
   };

   ~THATReconTreeMakerLoop() override = default;

   TTree *FindOrClone(const char *name, TTree *src)
   {
      auto *tree = dynamic_cast<TTree *>(fRootFile->Get(name));
      if (tree != nullptr)
         return tree;
      tree = src->CloneTree(0);
      tree->SetName(name);
      return tree;
   }

   void Initialize() override
   {
      if (fOutputFileName.empty()) {
         ND280Error("No output filename provided via '-O outfile=XXX.root' option");
         exit(1);
      }
      fRootFile = new TFile(fOutputFileName.c_str(), "RECREATE");
      if (fRootFile == nullptr or fRootFile->IsZombie()) {
         ND280Error("Cannot open output file: " << fOutputFileName);
         exit(1);
      }
   }

   bool operator()(ND::TND280Event &event) override
   {
      ND280Log("Open: Event " << m_counter_Evt);
      fRootFile->cd();
      auto output_tree = FindOrClone("outTree", fOutputTree);

      ND280Log("Event " << event.GetContext());
      auto fit = event.GetFit("THATRecon");

      EventT eventT;
      eventT.eventid = m_counter_Evt;
      eventT.tracks.clear();

      if (!fit) {
         ND280Warn("No HATRECON fit found; next!");
         return false;
      }

      if (fit->GetResultsContainer("THATPID")) {
         ND280Log("We have PID");
      } else if (fit->GetResultsContainer("THATTrackFitter")) {
         ND280Error("We have THATTrackFitter; need THATPID!");
         return false;
      }

      auto fusedHitSelection = fit->GetHitSelection("used");
      if (!fusedHitSelection) {
         ND280Log("Used hits selection not found");
         return false;
      }
      ND280Log("Used hits selection size: " << fusedHitSelection->size());

      uint iPattern = 0;
      int i_trk = 0;
      for (auto const &pattern : *fit->GetResultsContainer("THATPID")) {
         ND::THandle<ND::TReconPID> pid_result = pattern;
         if (!pid_result) {
            ND280Log("Not a pid");
            continue;
         }
         iPattern++;

         for (auto const &constituent : *pid_result->GetConstituents()) {
            ND::THandle<ND::TReconTrack> trk = constituent;
            if (!trk) {
               ND280Log("Not a track");
               continue;
            }
            if (!trk->GetHits()) {
               ND280Log("No hits");
               continue;
            }
            if (trk->Get<ND::TIntegerDatum>("PathId"))
               ND280Verbose("Pattern id: "
                            << pattern->Get<ND::TIntegerDatum>("PathId")->GetValue());
            ND280Verbose("\tNumber of hits " << trk->GetHits()->size());

            fCurvature = trk->GetCurvature();

            TVector3 Pos(trk->GetPosition()[0], trk->GetPosition()[1],
                         trk->GetPosition()[2]);
            TVector3 Dir(trk->GetDirection()[0], trk->GetDirection()[1],
                         trk->GetDirection()[2]);

            double B = 0.2;
            // project into the bending plane
            // double B = ND::TFieldManager::GetFieldValue(Pos).X() / unit::tesla;
            double factor = -(0.3 * B) / sqrt(1. - Dir.X() * Dir.X());
            double p, q;

            p = fabs(factor / fCurvature);
            q = -fCurvature / fabs(fCurvature);
            fMomentum = p * q;

            fNclusters = 0;
            for (const auto &hit : *trk->GetHits()) {
               ND::THandle<ND::TComboHit> comboHit = hit;
               if (!comboHit) {
                  ND280Warn("Not a comboHit");
                  continue;
               }
               fNclusters++;
            }
            fPos[0] = trk->GetPosition()[0];
            fPos[1] = trk->GetPosition()[1];
            fPos[2] = trk->GetPosition()[2];
            fDir[0] = trk->GetDirection()[0];
            fDir[1] = trk->GetDirection()[1];
            fDir[2] = trk->GetDirection()[2];

            fStartTime = trk->Get<ND::TRealDatum>("minTime")->GetValue();
            fEndTime = trk->Get<ND::TRealDatum>("maxTime")->GetValue();
            fcenter_circle_X = trk->Get<ND::TRealDatum>("circleX")->GetValue();
            fcenter_circle_Y = trk->Get<ND::TRealDatum>("circleY")->GetValue();

            fEndPlate_ID = GetEndPlate_ID(trk);
            ND280Log("Endplate ID: " << fEndPlate_ID);

            // ---------------------------------------------------------------//
            // example on how to read clusters information
            ComputeResolution(trk); // Need fEndPlate_ID being defined
            // ---------------------------------------------------------------//

            // ---------------------------------------------------------------//
            // example on how to read pad hits Waveforms
            fMeanTime = GetMeanTimeOfMax_WF(trk, fusedHitSelection);
            // ---------------------------------------------------------------//

            TrackingUtils::Curvature_to_MomentumAndCharge(Pos, Dir, fCurvature, p, q);
            // TODO we could also extract the particle ID and dEdx for each PID
            fPID = pid_result->GetParticleId();
            // fPullMuon = pid_result->GetPIDWeight();
            fdEdx = pid_result->Get<ND::TRealDatum>("dEdx_PID")->GetValue();
            // std::cout<<fPID<<" "<<fdEdx<<" pull mu
            // "<<pid_result->GetPIDWeight()<<std::endl;

            ND::TReconObjectContainer::const_iterator it;
            for (it = pid_result->GetAlternates().begin();
                 it != pid_result->GetAlternates().end(); ++it) {
               ND::THandle<ND::TReconPID> alter = *it;
               if (ND::TReconPID::ConvertParticleId(alter->GetParticleId()) ==
                   "Electron") {
                  // std::cout << " name " <<
                  // ND::TReconPID::ConvertParticleId(alter->GetParticleId()) << " pull "
                  //           << alter->GetPIDWeight() << std::endl;
                  fPullElec = alter->GetPIDWeight();
               }
               if (ND::TReconPID::ConvertParticleId(alter->GetParticleId()) == "Proton") {
                  // std::cout << " name " <<
                  // ND::TReconPID::ConvertParticleId(alter->GetParticleId()) << " pull "
                  //           << alter->GetPIDWeight() << std::endl;
                  fPullProt = alter->GetPIDWeight();
               }
               if (ND::TReconPID::ConvertParticleId(alter->GetParticleId()) == "Muon") {
                  // std::cout << " name " <<
                  // ND::TReconPID::ConvertParticleId(alter->GetParticleId()) << " pull "
                  //           << alter->GetPIDWeight() << std::endl;
                  fPullMuon = alter->GetPIDWeight();
               }
            }

            output_tree->Fill();

            ///////////////////////////////////////////////////////////////////////////////////////////
            ///Start Tristan

            TrackT track;
            track.trackid = i_trk;
            std::map<int, SubtrackT> moduleSubtracks;

            track.params.push_back(fPos[0]);
            track.params.push_back(fPos[1]);
            track.params.push_back(fPos[2]);
            track.params.push_back(fDir[0]);
            track.params.push_back(fDir[1]);
            track.params.push_back(fDir[2]);
            track.params.push_back(trk->GetCurvature());

            // Clusters
            for (auto clus : *trk->GetConstituents()) {
               ClusterT cluster;
               ND::THandle<ND::TReconCluster> clt = clus;
               TVector3 global;
               if (clt->GetPosition().X() > 0)
                  global.SetXYZ(980.59, clt->GetPosition().Y(),
                                clt->GetPosition().Z()); // X is the drift
               else
                  global.SetXYZ(-980.59, clt->GetPosition().Y(), clt->GetPosition().Z());
               if (clt->GetPosition().Y() > 1200 or clt->GetPosition().Y() < 400 or
                   clt->GetPosition().Z() > -1000 or clt->GetPosition().Z() < -2900)
                  continue;
               ND::TGeometryId id;
               bool is_inside = ND::TGeomInfo::HAT().GlobalXYZToGeomId(global, id);
               if (!is_inside) {
                  ND280Log("ID not found");
                  continue;
               }

               auto mms = ND::TGeomInfo::HAT().GeomIdToMM(id);
               auto half = ND::TGeomInfo::HAT().GeomIdToHalf(id);
               int modID = mms + 8 * half;
               cluster.x = clt->GetPosition().X();
               cluster.y = clt->GetPosition().Y();
               cluster.z = clt->GetPosition().Z();
               cluster.ex = clt->GetPositionVariance().X();
               cluster.ey = clt->GetPositionVariance().Y();
               cluster.ez = clt->GetPositionVariance().Z();

               if (moduleSubtracks.find(modID) ==
                   moduleSubtracks.end()) { // Create a new SubTrack for this module
                  SubtrackT subtrack;
                  subtrack.moduleid = modID;
                  subtrack.clusters.push_back(cluster);
                  moduleSubtracks[modID] = subtrack;
               } else { // Add the cluster to the existing SubTrack with this module
                  moduleSubtracks[modID].clusters.push_back(cluster);
               }
            }
            // Add all SubTracks to the Track
            for (const auto &pair : moduleSubtracks) {
               if (!pair.second.clusters.empty())
                  track.subtracks.push_back(pair.second);
               track.n_subtracks = track.subtracks.size();
            }

            eventT.tracks.push_back(track);
            i_trk++;
         }
      }
      runT.events.push_back(eventT);

      // Dump ////////////////////
      int nsubtrk = 0;
      for (const TrackT &track : eventT.tracks)
         nsubtrk += track.n_subtracks;
      ND280Log("  Subtracks: " << nsubtrk);
      int i_subtrk = 0;
      for (const TrackT &track : eventT.tracks) {
         ND280Log("  Track parameters: X = "
                  << track.params[0] << " | Y = " << track.params[1]
                  << " | Z = " << track.params[2] << " | Vx = " << track.params[3]
                  << " | Vy = " << track.params[4] << " | Vz = " << track.params[5]
                  << " | curve = " << track.params[6]);
         for (const SubtrackT &subtrack : track.subtracks) {
            ND280Log("    Open:  Subtracks " << i_subtrk + 1 << "/" << nsubtrk << " (MM #"
                                             << subtrack.moduleid << "):");
            int iC = 0;
            for (const ClusterT &cluster : subtrack.clusters) {
               ND280Log("      Cluster "
                        << std::setw(3) << iC << ": "
                        << "X = " << cluster.x << " ± " << cluster.ex << " | "
                        << "Y = " << cluster.y << " ± " << cluster.ey << " | "
                        << "Z = " << cluster.z << " ± " << 0.3);
               iC++;
            }
            ND280Log("    Close: Subtracks " << i_subtrk + 1 << "/" << nsubtrk << " (MM #"
                                             << subtrack.moduleid << "):");
            i_subtrk++;
         }
      }
      ND280Log("Close: Event " << m_counter_Evt);
      m_counter_Evt++;
      std::cout << std::endl;
      ///////////////////////////////////////////////////////////////////////////////////////////
      ///End Tristan

      ND280Log("Finished event");
      fEvent++;
      fRootFile->cd();
      output_tree->Write("", TObject::kOverwrite);
      fres_his_ep0->Write("", TObject::kOverwrite);
      fres_his_ep1->Write("", TObject::kOverwrite);
      residuals_col->Write("", TObject::kOverwrite);
      residuals_row->Write("", TObject::kOverwrite);
      return true;
   }

   void Finalize(ND::TND280Output *const file) override
   {
      fRootFile->cd();
      std::cout << runT.events.size() << " entries in runT" << std::endl;
      treeT->Fill();
      treeT->Write("", TObject::kOverwrite);
      treeT->Print();
      std::cout << treeT->GetEntries() << " entries in treeT" << std::endl;
      fRootFile->Close();
      // treeT->Fill();
      // fileT->cd();
      // fileT->Write();
      // fileT->Close();
   }

   bool SetOption(std::string option, std::string value) override
   {
      ND280Log("Option " << option << " = " << value);
      if (option == "outfile") {
         fOutputFileName = value;
      } else if (fParOptMan.IsRelevantOption(option)) {
         fParOptMan.UseRelevantOption(value);
      } else
         return false;
      return true;
   }

   void Usage() override
   {
      std::cout << "Generate TTree from reconstructed track parameters." << std::endl
                << std::endl;
      fParOptMan.Usage();
   }

   double GetMeanTimeOfMax_WF(ND::THandle<ND::TReconTrack> trk,
                              ND::THandle<ND::THitSelection> fusedHitSelection)
   {

      int i = 0;
      double mean = 0.;

      for (const auto &hit : *(trk->GetHits())) { // loop over hits

         ND::THandle<ND::TComboHit> comboHit = hit; // choosing only combo hits
         if (!comboHit) {
            ND280Warn("Not a comboHit");
            continue;
         }

         for (const auto &subhit : comboHit->GetHits()) { // loop over pad hits in combo

            auto mhit = THATUtils::FindWaveform(fusedHitSelection, subhit);
            if (not mhit) {
               ND280Log("Cannot find TMultiHit");
               continue;
            }

            TVector3 global_pad;
            ND::TGeomInfo::HAT().GeomIdToGlobalXYZ(subhit->GetGeomId(), global_pad);

            double q_max = -1e10;
            double t_max = -1e10;
            for (const auto &tmphit : *mhit) { // loop over WF time stamps
               if (tmphit->GetCharge() > q_max) {
                  q_max = tmphit->GetCharge();
                  t_max = tmphit->GetTime();
               }
            }
            mean += t_max;
            i++;
         }
      }
      return mean / double(i);
   }

   void ComputeResolution(ND::THandle<ND::TReconTrack> trk)
   {

      ND280Log("computing res");
      bool EP0 = false, EP1 = false;
      auto R = 1. / abs(fCurvature);
      int n_clus = 0;
      fResTrack = 0.;

      for (auto clus : *trk->GetConstituents()) {
         ND::THandle<ND::TReconCluster> clt = clus;
         if (clt) {

            TVector3 global;
            if (clt->GetPosition().X() > 0)
               global.SetXYZ(980.59, clt->GetPosition().Y(), clt->GetPosition().Z());
            else
               global.SetXYZ(-980.59, clt->GetPosition().Y(), clt->GetPosition().Z());

            if (clt->GetPosition().Y() > 1200 or clt->GetPosition().Y() < 400 or
                clt->GetPosition().Z() > -1000 or clt->GetPosition().Z() < -2900)
               continue;

            ND::TGeometryId id;
            bool is_inside = ND::TGeomInfo::HAT().GlobalXYZToGeomId(global, id);
            if (!is_inside) {
               ND280Log("ID not found");
               continue;
            }
            auto EP_num = ND::TGeomInfo::HAT().GeomIdToHalf(id);
            auto row = ND::TGeomInfo::HAT().GeomIdToRow(id);
            auto col = ND::TGeomInfo::HAT().GeomIdToColumn(id);
            auto mms = ND::TGeomInfo::HAT().GeomIdToMM(id);

            auto gen_col = col + 36 * mms + 36 * 8 * EP_num;
            auto gen_row = row + 32 * mms + 32 * 8 * EP_num;

            if (EP_num == 0)
               EP0 = true;
            else if (EP_num == 1)
               EP1 = true;
            auto r = sqrt((clt->GetPosition().Z() - fcenter_circle_X) *
                             (clt->GetPosition().Z() - fcenter_circle_X) +
                          (clt->GetPosition().Y() - fcenter_circle_Y) *
                             (clt->GetPosition().Y() - fcenter_circle_Y)) -
                     R;

            n_clus++;
            fResTrack += abs(r);
            residuals_col->Fill(gen_col, r);
            residuals_row->Fill(gen_row, r);

            if (EP0)
               fres_his_ep0->Fill(r);
            else if (EP1)
               fres_his_ep1->Fill(r);
         }
         EP0 = false;
         EP1 = false;
      }
      fResTrack /= double(n_clus);
   }

   int GetEndPlate_ID(ND::THandle<ND::TReconTrack> trk)
   {

      bool EP0 = false, EP1 = false;

      for (const auto &hit : *trk->GetHits()) {
         ND::THandle<ND::TComboHit> comboHit = hit;
         if (!comboHit) {
            ND280Warn("Not a comboHit");
            continue;
         }
         auto EP_num = ND::TGeomInfo::HAT().GeomIdToHalf(comboHit->GetGeomId());
         if (EP_num == 0)
            EP0 = true;
         else if (EP_num == 1)
            EP1 = true;
      }

      if (EP0 and EP1)
         return 12;
      else if (EP0 and !EP1)
         return 8;
      else if (EP1 and !EP0)
         return 4;
      else
         return 0;
   }

private:
   ND::TParametersOptionManager fParOptMan;
   std::string fOutputFileName;
   TFile *fRootFile;
   TTree *fOutputTree;
   Int_t fEvent;
   Double_t fMomentum;
   Int_t fNclusters;
   Double_t fResolution{};
   Double_t fChi2{};
   // Double_t fMeanMultiplicity{};
   Double_t fTrueMomentum{};

   Double_t fResTrack;

   Double_t fStartTime;
   Double_t fEndTime;
   Double_t fMeanTime;

   Double_t fPullElec;
   Double_t fPullProt;
   Double_t fPullMuon;
   Double_t fdEdx;
   Int_t fEndPlate_ID;

   Double_t fcenter_circle_X;
   Double_t fcenter_circle_Y;
   Double_t fCurvature;
   Double_t fPos[3]{};
   Double_t fDir[3]{};
   Double_t fDirNorm;
   Int_t fPID{};

   TH1F *fres_his_ep0;
   TH1F *fres_his_ep1;
   TH2D *residuals_col;
   TH2D *residuals_row;

   int m_counter_Evt;
   // TFile *fileT = new TFile("outputT.root", "RECREATE");
   TTree *treeT;
   RunT runT;
};

int main(int argc, char **argv)
{
   THATReconTreeMakerLoop userCode;
   return nd280EventLoop(argc, argv, userCode);
}
