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
#include <sstream>
#include <map>
#include <set>
#include <getopt.h>

Cube::Event* event = NULL;

TH1D lever_arm_sig("","lever arm",20,0,1000);
TH1D lever_arm_sig_cluster("","lever arm",20,0,1000);
TH1D lever_arm_bkg("","lever arm",20,0,1000);
TH1D lever_arm_bkg_cluster("","lever arm",20,0,1000);

TH1D eDep_sig("","edep",20,0,20000);
TH1D eDep_bkg("","edep",20,0,20000);
TH1D eDep_sig_cluster("","edep",20,0,2000);
TH1D eDep_bkg_cluster("","edep",20,0,2000);

TH1D trackLength_sig("","track length",20,0,500);
TH1D trackLength_bkg("","track length",20,0,500);

TH1D angle_sig("","angle",20,-1,1);
TH1D angle_bkg("","angle",20,-1,1);
TH1D angle_sig_cluster("","angle",20,-1,1);
TH1D angle_bkg_cluster("","angle",20,-1,1);

TH1D nbhdist_sig("","neighbor distance",20,0,2000);
TH1D nbhdist_bkg("","neighbor distance",20,0,2000);
TH1D nbhdist_sig_cluster("","neighbor distance",20,0,2000);
TH1D nbhdist_bkg_cluster("","neighbor distance",20,0,2000);

double leverArm;
double eDep;
double trackLength;
double angle;
double nbhdist;
auto outputFile  = std::make_shared<TFile> ("variableOutput.root","RECREATE");
auto outputTree  = std::make_shared<TTree> ("tree", "tree");

int category; //0: sig track, 1: sig cluster, 2: bkg track, 3: bkg cluster

const double THRESHOLD = 20;

void Analysis(Cube::Event* event);

int main(int argc, char** argv) {
    std::string inputFileName;

    while (optind < argc) {
        inputFileName = argv[optind];
        optind++;
    }

    std::unique_ptr<TChain> inputChain = std::make_unique<TChain> ("CubeEvents");
    inputChain->Add(inputFileName.c_str());
    inputChain->SetBranchAddress("Event", &event);

    outputTree->Branch("leverArm", &leverArm, "lever arm/D");
    outputTree->Branch("eDep", &eDep, "edeposit/D");
    outputTree->Branch("trackLength", &trackLength, "track length/D");
    outputTree->Branch("angle", &angle, "cos(angle)/D");
    outputTree->Branch("nbhdist", &nbhdist, "neighbor distance/D");
    outputTree->Branch("category", &category, "0: sig track, 1: sig cluster, 2: bkg track, 3: bkg cluster/I");

    for (int i = 0; i < inputChain->GetEntries(); i++) {
        std::cout << "event: " << i << std::endl;
        inputChain->GetEntry(i);
        Analysis(event);
    }
    outputFile->Write();
    outputFile->Close();

    TCanvas canLeverArm;
    lever_arm_sig.Scale(1/lever_arm_sig.Integral(),"nosw2");
    lever_arm_sig.GetYaxis()->SetRangeUser(0,1);
    lever_arm_bkg.Scale(1/lever_arm_bkg.Integral(),"nosw2");
    lever_arm_bkg.GetYaxis()->SetRangeUser(0,1);
    lever_arm_sig.SetLineColor(2);
    lever_arm_bkg.SetLineColor(4);
    lever_arm_sig.Draw();
    lever_arm_bkg.Draw("same");
    canLeverArm.SaveAs("leverArm.pdf");

    TCanvas caneDep;
    eDep_sig.Scale(1/eDep_sig.Integral(),"nosw2");
    eDep_sig.GetYaxis()->SetRangeUser(0,1);
    eDep_bkg.Scale(1/eDep_bkg.Integral(),"nosw2");
    eDep_bkg.GetYaxis()->SetRangeUser(0,1);
    eDep_sig.SetLineColor(2);
    eDep_bkg.SetLineColor(4);
    eDep_sig.Draw();
    eDep_bkg.Draw("same");
    caneDep.SaveAs("eDep.pdf");

    TCanvas cannbhdist;
    nbhdist_sig.Scale(1/nbhdist_sig.Integral(),"nosw2");
    nbhdist_sig.GetYaxis()->SetRangeUser(0,1);
    nbhdist_bkg.Scale(1/nbhdist_bkg.Integral(),"nosw2");
    nbhdist_bkg.GetYaxis()->SetRangeUser(0,1);
    nbhdist_sig.SetLineColor(2);
    nbhdist_bkg.SetLineColor(4);
    nbhdist_sig.Draw();
    nbhdist_bkg.Draw("same");
    cannbhdist.SaveAs("nbhdist.pdf");

    TCanvas cantrackLength;
    trackLength_sig.Scale(1/trackLength_sig.Integral(),"nosw2");
    trackLength_sig.GetYaxis()->SetRangeUser(0,1);
    trackLength_bkg.Scale(1/trackLength_bkg.Integral(),"nosw2");
    trackLength_bkg.GetYaxis()->SetRangeUser(0,1);
    trackLength_sig.SetLineColor(2);
    trackLength_bkg.SetLineColor(4);
    trackLength_sig.Draw();
    trackLength_bkg.Draw("same");
    cantrackLength.SaveAs("trackLength.pdf");

    TCanvas canangle;
    angle_sig.Scale(1/angle_sig.Integral(),"nosw2");
    angle_sig.GetYaxis()->SetRangeUser(0,1);
    angle_bkg.Scale(1/angle_bkg.Integral(),"nosw2");
    angle_bkg.GetYaxis()->SetRangeUser(0,1);
    angle_sig.SetLineColor(2);
    angle_bkg.SetLineColor(4);
    angle_sig.Draw();
    angle_bkg.Draw("same");
    canangle.SaveAs("angle.pdf");

    TCanvas canLeverArm_cluster;
    lever_arm_sig_cluster.Scale(1/lever_arm_sig_cluster.Integral(),"nosw2");
    lever_arm_sig_cluster.GetYaxis()->SetRangeUser(0,1);
    lever_arm_bkg_cluster.Scale(1/lever_arm_bkg_cluster.Integral(),"nosw2");
    lever_arm_bkg_cluster.GetYaxis()->SetRangeUser(0,1);
    lever_arm_sig_cluster.SetLineColor(2);
    lever_arm_bkg_cluster.SetLineColor(4);
    lever_arm_sig_cluster.Draw();
    lever_arm_bkg_cluster.Draw("same");
    canLeverArm_cluster.SaveAs("leverArm_cluster.pdf");

    TCanvas caneDep_cluster;
    eDep_sig_cluster.Scale(1/eDep_sig_cluster.Integral(),"nosw2");
    eDep_sig_cluster.GetYaxis()->SetRangeUser(0,1);
    eDep_bkg_cluster.Scale(1/eDep_bkg_cluster.Integral(),"nosw2");
    eDep_bkg_cluster.GetYaxis()->SetRangeUser(0,1);
    eDep_sig_cluster.SetLineColor(2);
    eDep_bkg_cluster.SetLineColor(4);
    eDep_sig_cluster.Draw();
    eDep_bkg_cluster.Draw("same");
    caneDep_cluster.SaveAs("eDep_cluster.pdf");

    TCanvas cannbhdist_cluster;
    nbhdist_sig_cluster.Scale(1/nbhdist_sig_cluster.Integral(),"nosw2");
    nbhdist_sig_cluster.GetYaxis()->SetRangeUser(0,1);
    nbhdist_bkg_cluster.Scale(1/nbhdist_bkg_cluster.Integral(),"nosw2");
    nbhdist_bkg_cluster.GetYaxis()->SetRangeUser(0,1);
    nbhdist_sig_cluster.SetLineColor(2);
    nbhdist_bkg_cluster.SetLineColor(4);
    nbhdist_sig_cluster.Draw();
    nbhdist_bkg_cluster.Draw("same");
    cannbhdist_cluster.SaveAs("nbhdist_cluster.pdf");

    TCanvas canangle_cluster;
    angle_sig_cluster.Scale(1/angle_sig_cluster.Integral(),"nosw2");
    angle_sig_cluster.GetYaxis()->SetRangeUser(0,1);
    angle_bkg_cluster.Scale(1/angle_bkg_cluster.Integral(),"nosw2");
    angle_bkg_cluster.GetYaxis()->SetRangeUser(0,1);
    angle_sig_cluster.SetLineColor(2);
    angle_bkg_cluster.SetLineColor(4);
    angle_sig_cluster.Draw();
    angle_bkg_cluster.Draw("same");
    canangle_cluster.SaveAs("angle_cluster.pdf");

    return 0;
}

void Analysis(Cube::Event* event) {

    leverArm = -1;
    eDep = -1;
    trackLength = -1;
    angle = -1;
    nbhdist = -1;

    Cube::Event::G4TrajectoryContainer trajectories = event->G4Trajectories;
    Cube::Handle<Cube::ReconObjectContainer> objects = event->GetObjectContainer();
    if (!objects) 
        return;

    int muonLike = 0;
    Cube::Handle<Cube::ReconTrack> muonObject;
    for (Cube::ReconObjectContainer::iterator o = objects->begin(); o != objects->end(); ++o) {
        Cube::Handle<Cube::ReconTrack> track = *o;
        if (!track) 
            continue;
        int mainTraj = Cube::Tool::MainTrajectory(*event, *track);
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
        if (!traj || traj->GetPDGCode() != -13 || traj->GetParentId() != -1) 
            continue;
        muonLike++;
        muonObject = track;
    }
    if (muonLike != 1 || !muonObject) return;

    int numberOfMuonAssociated = 0;
    double earliestTime = 1E+8;
    Cube::Handle<Cube::ReconObject> earliestObject;
    for (Cube::ReconObjectContainer::iterator o = objects->begin(); o != objects->end(); ++o) {
        if (Cube::Tool::AreNeighboringObjects(*muonObject, *(*o))) {
            numberOfMuonAssociated++;
            continue;
        }
        double objTime = -1.0;
        double objEDep = -1.0;
        Cube::Handle<Cube::ReconTrack> track = *o;
        Cube::Handle<Cube::ReconCluster> cluster = *o;
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
            earliestObject = *o;
        }
    }
    if (numberOfMuonAssociated != 1) return;
    if (!earliestObject) return; 
    int mainTraj = Cube::Tool::MainTrajectory(*event,*earliestObject);
    Cube::Handle<Cube::G4Trajectory> earliestTraj = trajectories[mainTraj];
    if (!earliestTraj) return;

    Cube::Handle<Cube::ReconTrack> earliestTrack = earliestObject;
    Cube::Handle<Cube::ReconCluster> earliestCluster = earliestObject;
    TVector3 beamDirection(0,0,1);
    double neighborDistance = 1E+8;
    for (Cube::ReconObjectContainer::iterator o = objects->begin(); o != objects->end(); ++o) {
        if (*o == muonObject) continue;
        if (*o == earliestObject) continue;
        double distance = Cube::Tool::DistanceBetweenObjects(*earliestObject, *(*o));
        if (distance < neighborDistance) {
            neighborDistance = distance;
        }
    }
    if (neighborDistance == 1E+8) {
        neighborDistance = -1;
    }
    nbhdist = neighborDistance;

    if (earliestTrack) {
        leverArm = (muonObject->GetPosition().Vect() - earliestTrack->GetPosition().Vect()).Mag();
        std::cout << "leverArm: " << leverArm << std::endl;
        eDep = earliestTrack->GetEDeposit();
        std::cout << "eDep: " << eDep << std::endl;
        Cube::ReconNodeContainer::iterator n = earliestTrack->GetNodes().begin();
        Cube::Handle<Cube::TrackState> lastState = (*(n++))->GetState();
        while (n != earliestTrack->GetNodes().end()) {
            Cube::Handle<Cube::TrackState> nodeState = (*(n++))->GetState();
            trackLength += (nodeState->GetPosition().Vect()
                    - lastState->GetPosition().Vect()).Mag();
            lastState = nodeState;
        }
        std::cout << "trackLength: " << trackLength << std::endl;
        angle = (earliestTrack->GetPosition() - muonObject->GetPosition()).Angle(beamDirection);
        std::cout << "angle: " << TMath::Cos(angle) << std::endl;
        std::cout << "neighborDistance: " << neighborDistance << std::endl;
        if (trajectories[mainTraj]->GetPDGCode() == 2112) {
            lever_arm_sig.Fill(leverArm);
            eDep_sig.Fill(eDep);
            trackLength_sig.Fill(trackLength);
            angle_sig.Fill(angle);
            nbhdist_sig.Fill(neighborDistance);
            category = 0;
        }
        else {
            lever_arm_bkg.Fill(leverArm);
            eDep_bkg.Fill(eDep);
            trackLength_bkg.Fill(trackLength);
            angle_bkg.Fill(angle);
            nbhdist_bkg.Fill(neighborDistance);
            category = 2;
        }
    }
    if (earliestCluster) {
        leverArm = (muonObject->GetPosition().Vect() - earliestCluster->GetPosition().Vect()).Mag();
        std::cout << "leverArm: " << leverArm << std::endl;
        eDep = earliestCluster->GetEDeposit();
        std::cout << "eDep: " << eDep << std::endl;
        std::cout << "trackLength: " << trackLength << std::endl;
        angle = (earliestCluster->GetPosition() - muonObject->GetPosition()).Angle(beamDirection);
        std::cout << "angle: " << TMath::Cos(angle) << std::endl;
        std::cout << "neighborDistance: " << neighborDistance << std::endl;
        if (trajectories[mainTraj]->GetPDGCode() == 2112) {
            lever_arm_sig_cluster.Fill(leverArm);
            eDep_sig_cluster.Fill(eDep);
            angle_sig_cluster.Fill(angle);
            nbhdist_sig_cluster.Fill(neighborDistance);
            category = 1;
        }
        else {
            lever_arm_bkg_cluster.Fill(leverArm);
            eDep_bkg_cluster.Fill(eDep);
            angle_bkg_cluster.Fill(angle);
            nbhdist_bkg_cluster.Fill(neighborDistance);
            category = 3;
        }
    }
    outputTree->Fill();
}
