#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconCluster.hxx>
#include <CubeReconNode.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>

#include <ToolPrimaryId.hxx>
#include <ToolG4Hits.hxx>
#include <ToolMainTrajectory.hxx>
#include <ToolTrueDirection.hxx>
#include <ToolContained.hxx>
#include <ToolRecon.hxx>

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>

#include <iostream>

Cube::Event* event = nullptr;
double EvtVtx[4]; 
double StdHepP4[1000][4]; 
int StdHepStatus[1000]; 
int StdHepPdg[1000]; 
int StdHepN; 
const double THRESHOLD = 20;

TH1D tof("","(true-reco)/true tof",100,-2,2);
TH1D leverarm("","(true-reco)/true lever arm",100,-2,2);

int main2()
{
    int eventNum = 0;
    TH2D true_reco_KE("","true vs reco neutron KE;reco;true",50,0,500,50,0,500);
    for (int i = 100; i < 200; i++) {
        TFile inputFile(Form("/Users/gwon/Analysis/datafiles/full3DST.antineutrino.%d.cuberecon_latest.root",i));
        TFile inputGenieFile(Form("/Users/gwon/CubeAnalysis/datafiles/latest/full3DST.antineutrino.%d.rootracker.root",i+1));
        TTree* inputChain = (TTree*)inputFile.Get("CubeEvents");
        TTree* inputGenieTree = (TTree*)inputGenieFile.Get("gRooTracker");
        inputChain->SetBranchAddress("Event", &event);
        inputGenieTree->SetBranchAddress("EvtVtx", &EvtVtx);
        inputGenieTree->SetBranchAddress("StdHepP4", &StdHepP4);
        inputGenieTree->SetBranchAddress("StdHepStatus", &StdHepStatus);
        inputGenieTree->SetBranchAddress("StdHepPdg", &StdHepPdg);
        inputGenieTree->SetBranchAddress("StdHepN", &StdHepN);
        for (int i = 0; i < inputChain->GetEntries(); i++) {
            try {
                std::cout << "event: " << eventNum << std::endl;
                eventNum++;
                inputChain->GetEntry(i);
                inputGenieTree->GetEntry(i);

                Cube::Event::G4TrajectoryContainer& trajectories = event->G4Trajectories;
                Cube::Handle<Cube::ReconObjectContainer> objects = event->GetObjectContainer();

                int numPrimaryAntiMuon = 0;
                Cube::Handle<Cube::ReconTrack> muonObject;

                for (auto& o : *objects) {
                    Cube::Handle<Cube::ReconTrack> track = o;
                    if (!track) 
                        continue;
                    int mainTraj = Cube::Tool::MainTrajectory(*event, *track);
                    Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
                    if (!traj || traj->GetPDGCode() != -13 || traj->GetParentId() != -1)
                        continue;
                    numPrimaryAntiMuon++;
                    muonObject = track;
                }
                if (numPrimaryAntiMuon != 1 || !muonObject)
                    continue;

                int numberOfMuonAssociated = 0;
                double earliestTime = 1E+8;
                Cube::Handle<Cube::ReconObject> earliestObject;
                for (auto& o : *objects) {
                    if (Cube::Tool::AreNeighboringObjects(*muonObject, *o)) {
                        numberOfMuonAssociated++;
                        continue;
                    }
                    double objTime = -1.0;
                    double objEDep = -1.0;
                    Cube::Handle<Cube::ReconTrack> track = o;
                    Cube::Handle<Cube::ReconCluster> cluster = o;
                    if (track) {
                        objTime = track->GetMedian().T();
                        objEDep = track->GetEDeposit();
                    } else if (cluster) {
                        objTime = cluster->GetMedian().T();
                        objEDep = cluster->GetEDeposit();
                    } else 
                        continue;
                    if (objTime < earliestTime && objEDep > THRESHOLD) {
                        earliestTime = objTime;
                        earliestObject = o;
                    }
                }
                if (numberOfMuonAssociated != 1) continue;
                if (!earliestObject) continue; 
                Cube::Handle<Cube::ReconTrack> earliestTrack = earliestObject;
                Cube::Handle<Cube::ReconCluster> earliestCluster = earliestObject;

                if (!earliestCluster)
                    continue;

                int mainTraj = Cube::Tool::MainTrajectory(*event,*earliestObject);
                Cube::Handle<Cube::G4Trajectory> earliestTraj = trajectories[mainTraj];
                if (!earliestTraj) continue;
                int mainTraj2 = Cube::Tool::MainTrajectory(*event,*muonObject);
                Cube::Handle<Cube::G4Trajectory> muonTraj = trajectories[mainTraj2];
                if (!muonTraj) continue;

                std::vector<Cube::Handle<Cube::G4Hit>> earliestSegs
                    = Cube::Tool::ObjectG4Hits(*event,*earliestObject);
                double earliestTruth = 1E+8;
                TVector3 earliestVector;
                for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator
                        t = earliestSegs.begin();
                        t != earliestSegs.end(); ++t) {
                    if ((*t)->GetStart().T() < earliestTruth) {
                        earliestTruth = (*t)->GetStart().T();
                        earliestVector.SetX((*t)->GetStart().X());
                        earliestVector.SetY((*t)->GetStart().Y());
                        earliestVector.SetZ((*t)->GetStart().Z());
                    }
                }

                std::vector<Cube::Handle<Cube::G4Hit>> muonSegs
                    = Cube::Tool::ObjectG4Hits(*event,*muonObject);
                double muonTruth = 1E+8;
                TVector3 muonVector;
                for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator
                        t = muonSegs.begin();
                        t != muonSegs.end(); ++t) {
                    if ((*t)->GetStart().T() < muonTruth) {
                        muonTruth = (*t)->GetStart().T();
                        muonVector.SetX((*t)->GetStart().X());
                        muonVector.SetY((*t)->GetStart().Y());
                        muonVector.SetZ((*t)->GetStart().Z());
                    }
                }

                std::cout << mainTraj << ", " << mainTraj2 << std::endl;

                double trueLeverArm = (earliestVector - muonVector).Mag();
                double recoLeverArm = (earliestCluster->GetPosition().Vect() - muonObject->GetPosition().Vect()).Mag();
                double trueTof = (earliestTruth - muonTruth);
                double recoTof = (earliestCluster->GetMedian().T() - muonObject->GetPosition().T());
                double trueBeta = (trueLeverArm/trueTof)/300;
                std::cout << "true earliestTraj T: " << earliestTraj->GetInitialPosition().T() << std::endl;
                std::cout << "true muon T: " << muonTraj->GetInitialPosition().T() << std::endl;
                std::cout << "true tof: " << trueTof << std::endl;
                std::cout << "true lever arm: " << trueLeverArm << std::endl;
                std::cout << "trueBeta: " << trueBeta << std::endl;
                double recoBeta = (recoLeverArm/recoTof)/300;
                double trueKE = 939.565*(1./std::pow(1.-std::pow(trueBeta,2),0.5)-1.);
                double recoKE = 939.565*(1./std::pow(1.-std::pow(recoBeta,2),0.5)-1.);
                true_reco_KE.Fill(recoKE, trueKE);
                std::cout << "trueKE: " << trueKE << std::endl;
                std::cout << "recoKE: " << recoKE << std::endl;

                tof.Fill((trueTof - recoTof)/trueTof);
                leverarm.Fill((trueLeverArm - recoLeverArm)/trueLeverArm);
            } catch(...) {
                continue;
            }
        }
    }

    tof.SaveAs("tof.C");
    leverarm.SaveAs("leverarm.C");
    true_reco_KE.SaveAs("true_reco_KE.C");

    return 0;
}

int main()
{
    int eventNum = 0;
    TH2D true_reco_KE("","true vs reco neutron KE;reco;true",50,0,500,50,0,500);
    for (int i = 100; i < 200; i++) {
        TFile inputFile(Form("/Users/gwon/Analysis/datafiles/full3DST.antineutrino.%d.cuberecon_latest.root",i));
        TFile inputGenieFile(Form("/Users/gwon/CubeAnalysis/datafiles/latest/full3DST.antineutrino.%d.rootracker.root",i+1));
        TTree* inputChain = (TTree*)inputFile.Get("CubeEvents");
        TTree* inputGenieTree = (TTree*)inputGenieFile.Get("gRooTracker");
        inputChain->SetBranchAddress("Event", &event);
        inputGenieTree->SetBranchAddress("EvtVtx", &EvtVtx);
        inputGenieTree->SetBranchAddress("StdHepP4", &StdHepP4);
        inputGenieTree->SetBranchAddress("StdHepStatus", &StdHepStatus);
        inputGenieTree->SetBranchAddress("StdHepPdg", &StdHepPdg);
        inputGenieTree->SetBranchAddress("StdHepN", &StdHepN);
        for (int i = 0; i < inputChain->GetEntries(); i++) {
            try {
                std::cout << "event: " << eventNum << std::endl;
                eventNum++;
                inputChain->GetEntry(i);
                inputGenieTree->GetEntry(i);

                Cube::Event::G4TrajectoryContainer& trajectories = event->G4Trajectories;
                Cube::Handle<Cube::ReconObjectContainer> objects = event->GetObjectContainer();

                std::vector<int> muons;
                std::vector<int> antiMuons;
                std::vector<int> electrons;
                std::vector<int> gammas;
                std::vector<int> pions;
                std::vector<int> neutrons;
                std::vector<int> protons;
                std::vector<int> others;
                // Get the primaries for this event.
                std::vector<int> primId = Cube::Tool::AllPrimaries(*event);
                for (std::vector<int>::iterator tr = primId.begin();
                        tr != primId.end(); ++tr) {
                    Cube::Handle<Cube::G4Trajectory> traj = trajectories[*tr];
                    if (!traj) {
                        std::cout << "This cannot happen!";
                        throw std::runtime_error("Invalid primary id");
                    }
                    switch (traj->GetPDGCode()) {
                        case 11: case -11: electrons.push_back(*tr); break;
                        case 13:  muons.push_back(*tr); break;
                        case -13: antiMuons.push_back(*tr); break;
                        case 22: gammas.push_back(*tr); break;
                        case 111: case 211: case -211: pions.push_back(*tr); break;
                        case 2112: neutrons.push_back(*tr); break;
                        case 2212: protons.push_back(*tr); break;
                        default:
                                   others.push_back(*tr);
                                   break;
                    }
                }

                if (!electrons.empty()) continue;
                if (!muons.empty()) continue;
                if (antiMuons.empty()) continue;
                if (!pions.empty()) continue;
                if (!protons.empty()) continue;
                if (neutrons.size() != 1) continue;
                if (!gammas.empty()) continue;
                if (!others.empty()) continue;

                int numPrimaryAntiMuon = 0;
                Cube::Handle<Cube::ReconTrack> muonObject;

                for (auto& o : *objects) {
                    Cube::Handle<Cube::ReconTrack> track = o;
                    if (!track) 
                        continue;
                    int mainTraj = Cube::Tool::MainTrajectory(*event, *track);
                    Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
                    if (!traj || traj->GetPDGCode() != -13 || traj->GetParentId() != -1)
                        continue;
                    numPrimaryAntiMuon++;
                    muonObject = track;
                }
                if (numPrimaryAntiMuon != 1 || !muonObject)
                    continue;

                double earliestTime = 1E+8;
                Cube::Handle<Cube::ReconObject> earliestObject;
                for (auto& o : *objects) {
                    if (Cube::Tool::AreNeighboringObjects(*muonObject, *o)) {
                        continue;
                    }
                    double objTime = -1.0;
                    double objEDep = -1.0;
                    Cube::Handle<Cube::ReconTrack> track = o;
                    Cube::Handle<Cube::ReconCluster> cluster = o;
                    if (track) {
                        objTime = track->GetMedian().T();
                        objEDep = track->GetEDeposit();
                    } else if (cluster) {
                        objTime = cluster->GetMedian().T();
                        objEDep = cluster->GetEDeposit();
                    } else 
                        continue;
                    if (objTime < earliestTime && objEDep > THRESHOLD) {
                        earliestTime = objTime;
                        earliestObject = o;
                    }
                }
                if (!earliestObject) continue; 
                Cube::Handle<Cube::ReconTrack> earliestTrack = earliestObject;
                Cube::Handle<Cube::ReconCluster> earliestCluster = earliestObject;

                if (!earliestCluster)
                    continue;

                int mainTraj = Cube::Tool::MainTrajectory(*event,*earliestObject);
                Cube::Handle<Cube::G4Trajectory> earliestTraj = trajectories[mainTraj];
                if (!earliestTraj) continue;
                int mainTraj2 = Cube::Tool::MainTrajectory(*event,*muonObject);
                Cube::Handle<Cube::G4Trajectory> muonTraj = trajectories[mainTraj2];
                if (!muonTraj) continue;

                std::vector<Cube::Handle<Cube::G4Hit>> earliestSegs
                    = Cube::Tool::ObjectG4Hits(*event,*earliestObject);
                double earliestTruth = 1E+8;
                TVector3 earliestVector;
                for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator
                        t = earliestSegs.begin();
                        t != earliestSegs.end(); ++t) {
                    if ((*t)->GetStart().T() < earliestTruth) {
                        earliestTruth = (*t)->GetStart().T();
                        earliestVector.SetX((*t)->GetStart().X());
                        earliestVector.SetY((*t)->GetStart().Y());
                        earliestVector.SetZ((*t)->GetStart().Z());
                    }
                }

                std::vector<Cube::Handle<Cube::G4Hit>> muonSegs
                    = Cube::Tool::ObjectG4Hits(*event,*muonObject);
                double muonTruth = 1E+8;
                TVector3 muonVector;
                for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator
                        t = muonSegs.begin();
                        t != muonSegs.end(); ++t) {
                    if ((*t)->GetStart().T() < muonTruth) {
                        muonTruth = (*t)->GetStart().T();
                        muonVector.SetX((*t)->GetStart().X());
                        muonVector.SetY((*t)->GetStart().Y());
                        muonVector.SetZ((*t)->GetStart().Z());
                    }
                }

                std::cout << mainTraj << ", " << mainTraj2 << std::endl;

                double trueLeverArm = (earliestVector - muonVector).Mag();
                double recoLeverArm = (earliestCluster->GetPosition().Vect() - muonObject->GetPosition().Vect()).Mag();
                double trueTof = (earliestTruth - muonTruth);
                double recoTof = (earliestCluster->GetMedian().T() - muonObject->GetPosition().T());
                double trueBeta = (trueLeverArm/trueTof)/300;
                std::cout << "true earliestTraj T: " << earliestTraj->GetInitialPosition().T() << std::endl;
                std::cout << "true muon T: " << muonTraj->GetInitialPosition().T() << std::endl;
                std::cout << "true tof: " << trueTof << std::endl;
                std::cout << "true lever arm: " << trueLeverArm << std::endl;
                std::cout << "trueBeta: " << trueBeta << std::endl;
                double recoBeta = (recoLeverArm/recoTof)/300;
                double trueKE = 939.565*(1./std::pow(1.-std::pow(trueBeta,2),0.5)-1.);
                double recoKE = 939.565*(1./std::pow(1.-std::pow(recoBeta,2),0.5)-1.);
                true_reco_KE.Fill(recoKE, trueKE);
                std::cout << "trueKE: " << trueKE << std::endl;
                std::cout << "recoKE: " << recoKE << std::endl;

                tof.Fill((trueTof - recoTof)/trueTof);
                leverarm.Fill((trueLeverArm - recoLeverArm)/trueLeverArm);
            } catch(...) {
                continue;
            }
        }
    }

    tof.SaveAs("tof_True_Channel.C");
    leverarm.SaveAs("leverarm_True_Channel.C");
    true_reco_KE.SaveAs("true_reco_KE_trueChannel.C");

    return 0;
}
