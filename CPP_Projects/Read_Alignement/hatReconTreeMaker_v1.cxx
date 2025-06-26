#include <iostream>
#include <memory>
#include <numeric>

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
      ND280Log("Beacon: THATReconTreeMakerLoop()");
      m_counter_Evt = 0;
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
      ND280Log("Beacon: FindOrClone()");
      auto *tree = dynamic_cast<TTree *>(fRootFile->Get(name));
      if (tree != nullptr)
         return tree;
      tree = src->CloneTree(0);
      tree->SetName(name);
      return tree;
   }

   void Initialize() override
   {
      ND280Log("Beacon: Initialize()");
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
      ND280Log("Beacon: operator()");
      ND280Log("Open: Event " << m_counter_Evt);
      fRootFile->cd();
      auto output_tree = FindOrClone("outTree", fOutputTree);

      ND280Log("Event " << event.GetContext());
      auto fit = event.GetFit("THATRecon");

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
      int counter_Trk = 0;
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

            /*
            double complete = -999.;
            double clean = -999.;
            ND::THandle<ND::TG4Trajectory> g4traj = TrackTruthInfo::GetG4Trajectory(trk,
            complete, clean); if (!g4traj) { ND280Warn("Didn't find G4 traj!"); continue;
            }
            double dist;
            ND::THandle<ND::TTrackState> trackState = trk->GetState();
            ND::TG4TrajectoryPoint G4Point2 = GetCloserG4Point(*g4traj, trackState, dist);
            fTrueMomentum = G4Point2.GetMomentum().Mag();
            */
            fCurvature = trk->GetCurvature();

            TVector3 Pos(trk->GetPosition()[0], trk->GetPosition()[1],
                         trk->GetPosition()[2]);
            TVector3 Dir(trk->GetDirection()[0], trk->GetDirection()[1],
                         trk->GetDirection()[2]);

            ND280Log("  Open:  Trk " << counter_Trk << " | " << Pos[0] << " | " << Pos[1]
                                     << " | " << Pos[2]);

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

            // int n_clus = 0;

            // Clusters (working)
            std::vector<std::vector<float>> arr_subtracks_x;
            std::vector<std::vector<float>> arr_subtracks_y;
            std::vector<std::vector<float>> arr_subtracks_z;
            std::vector<std::vector<float>> arr_subtracks_ex;
            std::vector<std::vector<float>> arr_subtracks_ey;
            std::vector<std::vector<float>> arr_subtracks_ez;
            for (int k = 0; k < 32; k++) {
               std::vector<float> v_subtracks_x;
               std::vector<float> v_subtracks_y;
               std::vector<float> v_subtracks_z;
               std::vector<float> v_subtracks_ex;
               std::vector<float> v_subtracks_ey;
               std::vector<float> v_subtracks_ez;
               arr_subtracks_x.push_back(v_subtracks_x);
               arr_subtracks_y.push_back(v_subtracks_y);
               arr_subtracks_z.push_back(v_subtracks_z);
               arr_subtracks_ex.push_back(v_subtracks_ex);
               arr_subtracks_ey.push_back(v_subtracks_ey);
               arr_subtracks_ez.push_back(v_subtracks_ez);
            }
            int list_subtrack[32] = {};
            int n_clus[32] = {};
            for (auto clus : *trk->GetConstituents()) {
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
               auto row = ND::TGeomInfo::HAT().GeomIdToRow(id);
               auto col = ND::TGeomInfo::HAT().GeomIdToColumn(id);
               auto mms = ND::TGeomInfo::HAT().GeomIdToMM(id);
               auto half = ND::TGeomInfo::HAT().GeomIdToHalf(id);
               // arr_subtracks[mms+8*half].push_back(clt) ;
               arr_subtracks_x[mms + 8 * half].push_back(clt->GetPosition().X());
               arr_subtracks_y[mms + 8 * half].push_back(clt->GetPosition().Y());
               arr_subtracks_z[mms + 8 * half].push_back(clt->GetPosition().Z());
               arr_subtracks_ex[mms + 8 * half].push_back(clt->GetPositionVariance().X());
               arr_subtracks_ey[mms + 8 * half].push_back(clt->GetPositionVariance().Y());
               arr_subtracks_ez[mms + 8 * half].push_back(clt->GetPositionVariance().Z());
               list_subtrack[mms + 8 * half] = 1;
               n_clus[mms + 8 * half]++;
               // ND280Log("    Cluster " << std::setw(3) << n_clus <<  ": " << "MM " <<
               // std::setw(2) << mms+8*half << " | col: " << std::setw(2) << col << " |
               // row: " << std::setw(2) << row << " | " <<
               //             "X = " << clt->GetPosition().X() << " ± " <<
               //             clt->GetPositionVariance().X() << " | " << "Y = " <<
               //             clt->GetPosition().Y() << " ± " <<
               //             clt->GetPositionVariance().Y() << " | " << "Z = " <<
               //             clt->GetPosition().Z() << " ± " <<
               //             clt->GetPositionVariance().Z() );
               // n_clus++;
            }
            int n_subtracks =
               std::accumulate(std::begin(list_subtrack), std::end(list_subtrack), 0);
            ND280Log("    Subtracks: " << n_subtracks);
            int i_track = 0;
            for (int iMM = 0; iMM < 32; iMM++) {
               if (list_subtrack[iMM] == 0)
                  continue;
               ND280Log("    Open:  Subtracks " << i_track + 1 << "/" << n_subtracks
                                                << " (MM #" << iMM << "):");
               for (int iC = 0; iC < (int)n_clus[iMM]; iC++) {
                  ND280Log("      Cluster " << std::setw(3) << iC << ": "
                                            << "X = " << arr_subtracks_x[iMM][iC] << " ± "
                                            << arr_subtracks_ex[iMM][iC] << " | "
                                            << "Y = " << arr_subtracks_y[iMM][iC] << " ± "
                                            << arr_subtracks_ey[iMM][iC] << " | "
                                            << "Z = " << arr_subtracks_z[iMM][iC] << " ± "
                                            << arr_subtracks_ez[iMM][iC]);
               }
               ND280Log("    Close: Subtracks " << i_track + 1 << "/" << n_subtracks
                                                << " (MM #" << iMM << "):");
               i_track++;
            }

            // for (const auto& hit: *(trk->GetHits())) { //loop over hits
            //     ND::THandle<ND::TComboHit> chit = hit; // choosing only combo hits
            //     if (!chit) { ND280Warn("Not a comboHit"); continue; }
            //     TVector3 global_clus;
            //     if( chit->GetPosition().X() > 0 )
            //     global_clus.SetXYZ(chit->GetPosition().X(), chit->GetPosition().Y(),
            //     chit->GetPosition().Z() ); else
            //     global_clus.SetXYZ(chit->GetPosition().X(), chit->GetPosition().Y(),
            //     chit->GetPosition().Z() ); if (chit->GetPosition().Y() > 1200 or
            //     chit->GetPosition().Y() < 400 or chit->GetPosition().Z() > -1000 or
            //     chit->GetPosition().Z() < -2900) continue; ND::TGeometryId id; bool
            //     is_inside = ND::TGeomInfo::HAT().GlobalXYZToGeomId(global_clus, id); if
            //     (!is_inside) { ND280Log("ID not found"); continue; } auto row =
            //     ND::TGeomInfo::HAT().GeomIdToRow(id); auto col =
            //     ND::TGeomInfo::HAT().GeomIdToColumn(id); auto mms =
            //     ND::TGeomInfo::HAT().GeomIdToMM(id); ND280Log(       "Cluster " <<
            //     std::setw(3) << n_clus <<  ": " <<
            //                     "X = " << chit->GetPosition().X()   << " ± " <<
            //                     std::setw(5) << std::setprecision(4)
            //                     << chit->GetUncertainty().X() << " | " << "Y = " <<
            //                     chit->GetPosition().Y()   << " ± " << std::setw(5) <<
            //                     std::setprecision(4) << chit->GetUncertainty().Y() << "
            //                     | " << "Z = " << chit->GetPosition().Z()   << " ± " <<
            //                     std::setw(5) << std::setprecision(4)
            //                     << chit->GetUncertainty().Z() << " | " << "Endplate: "
            //                     << fEndPlate_ID << " | MM " << mms << " | col: " << col
            //                     << " | row: " << row);

            //     int n_pad = 0;
            //     for (const auto &subhit: chit->GetHits()) {
            //         ND280Log(   "    Pad " << std::setw(3) << n_pad << ": " <<
            //                     "X = " << subhit->GetPosition().X() << std::setw(11) <<
            //                     " | " << "Y = " << subhit->GetPosition().Y() <<
            //                     std::setw(10) << " | " << "Z = " <<
            //                     subhit->GetPosition().Z()  );
            //         n_pad++;
            //     }
            //     n_clus++;
            // }

            // // Cluster and pads (not working)
            // auto hit_selection = trk->Get<ND::THitSelection>("clusters");
            // for (auto clus: *hit_selection){
            //     ND::THandle<ND::TComboHit> chit = clus;
            //     if (not chit) continue; // this is not a Combo hit
            //     // Clusters (ComboHit)
            //     TVector3 global_clus;
            //     if( chit->GetPosition().X() > 0 )   global_clus.SetXYZ( 980.59,
            //     chit->GetPosition().Y(), chit->GetPosition().Z() ); else
            //     global_clus.SetXYZ(-980.59, chit->GetPosition().Y(),
            //     chit->GetPosition().Z() ); if (chit->GetPosition().Y() > 1200 or
            //     chit->GetPosition().Y() < 400 or chit->GetPosition().Z() > -1000 or
            //     chit->GetPosition().Z() < -2900) continue; ND::TGeometryId id; bool
            //     is_inside = ND::TGeomInfo::HAT().GlobalXYZToGeomId(global_clus, id); if
            //     (!is_inside) { ND280Log("ID not found"); continue; } auto row =
            //     ND::TGeomInfo::HAT().GeomIdToRow(id); auto col =
            //     ND::TGeomInfo::HAT().GeomIdToColumn(id); auto mms =
            //     ND::TGeomInfo::HAT().GeomIdToMM(id); ND280Log("Cluster " << n_clus << "
            //     | MM " << mms << " | Endplate: " << fEndPlate_ID     << " | col: " <<
            //     col << " | row: " << row << " | " <<
            //                 "X = " << chit->GetPosition().X() << " ± " <<
            //                 chit->GetUncertainty().X() << " | " << "Y = " <<
            //                 chit->GetPosition().Y() << " ± " <<
            //                 chit->GetUncertainty().Y() << " | " << "Z = " <<
            //                 chit->GetPosition().Z() << " ± " <<
            //                 chit->GetUncertainty().Z() );
            //     n_clus++;
            //     // Pads (Hit)
            //     int n_pad = 0;
            //     for (auto hit: chit->GetHits()){
            //         TVector3 global_pad;
            //         if (hit->GetPosition().X() > 0) global_pad.SetXYZ( 980.59,
            //         hit->GetPosition().Y(), hit->GetPosition().Z()); else
            //         global_pad.SetXYZ(-980.59, hit->GetPosition().Y(),
            //         hit->GetPosition().Z()); if (hit->GetPosition().Y() > 1200 or
            //         hit->GetPosition().Y() < 400 or
            //             hit->GetPosition().Z() > -1000 or hit->GetPosition().Z() <
            //             -2900) continue;
            //         ND280Log(   "Pad " << n_pad << ": " <<
            //                     "X = " << hit->GetPosition().X() << " ± " <<
            //                     hit->GetUncertainty().X() << " | " << "Y = " <<
            //                     hit->GetPosition().Y() << " ± " <<
            //                     hit->GetUncertainty().Y() << " | " << "Z = " <<
            //                     hit->GetPosition().Z() << " ± " <<
            //                     hit->GetUncertainty().Z() );
            //         n_pad++;
            //     }
            // }

            // // Mathieu's proposal (not working)
            // auto hit_selection = trk->Get<ND::THitSelection>("clusters");
            // hit_selection.ls() ;
            // for (auto clus: *hit_selection) {
            //     ND::THandle<ND::TComboHit> chit = clus;
            //     if (not chit) continue; // this is not a Combo hit
            //     for (auto hit: chit->GetHits()) {
            //         TVector3 global;
            //         if (hit->GetPosition().X() > 0)
            //             global.SetXYZ(980.59, hit->GetPosition().Y(),
            //             hit->GetPosition().Z());
            //         else global.SetXYZ(-980.59, hit->GetPosition().Y(),
            //         hit->GetPosition().Z()); if (hit->GetPosition().Y() > 1200 or
            //         hit->GetPosition().Y() < 400 or
            //             hit->GetPosition().Z() > -1000 or hit->GetPosition().Z() <
            //             -2900) continue;
            //         ND::TGeometryId id;
            //         bool is_inside = ND::TGeomInfo::HAT().GlobalXYZToGeomId(global,
            //         id); if (!is_inside) {
            //             ND280Log("ID not found");
            //             continue;
            //         }
            //         auto mms = ND::TGeomInfo::HAT().GeomIdToMM(id);
            //     }
            // }

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
            ND280Log("  Close: Trk " << counter_Trk);
            counter_Trk++;
         }
      }

      ND280Log("Close: Event " << m_counter_Evt);
      m_counter_Evt++;
      std::cout << std::endl;
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
      ND280Log("Beacon: Finalize()");
      fRootFile->Close();
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
            double t_max;
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
};

int main(int argc, char **argv)
{
   THATReconTreeMakerLoop userCode;
   return nd280EventLoop(argc, argv, userCode);
}
