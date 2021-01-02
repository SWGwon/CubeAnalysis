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
double leverArm;
double eDep;
double trackLength;
double angle;
double nbhdist;
double neutrinoEBeforeSelection;
double neutrinoEAfterSelection;
int category; //0: sig track, 1: sig cluster, 2: bkg track, 3: bkg cluster
auto outputFile  = std::make_shared<TFile> ("variableOutput.root","RECREATE");
auto outputTree  = std::make_shared<TTree> ("tree", "tree");

double EvtVtx[4]; 
double StdHepP4[1000][4]; 
int StdHepStatus[1000]; 
int StdHepPdg[1000]; 
int StdHepN; 

TH1D nuEBeforeSelection("","neutrinoE before selection",100,0,10);
TH1D nuEAfterSelection("","neutrinoE after selection",100,0,10);
TH1D nuEAfterSelectionForDivide("","neutrinoE after selection",100,0,10);

const double THRESHOLD = 20;

void Analysis(Cube::Event* event);

int main(int argc, char** argv) {
    std::string inputFileName;

    while (optind < argc) {
        inputFileName = argv[optind];
        optind++;
    }
    outputTree->Branch("leverArm", &leverArm, "lever arm/D");
    outputTree->Branch("eDep", &eDep, "edeposit/D");
    outputTree->Branch("trackLength", &trackLength, "track length/D");
    outputTree->Branch("angle", &angle, "cos(angle)/D");
    outputTree->Branch("nbhdist", &nbhdist, "neighbor distance/D");
    outputTree->Branch("neutrinoEBeforeSelection", &neutrinoEBeforeSelection, "neutirnoEe/D");
    outputTree->Branch("neutrinoEAfterSelection", &neutrinoEAfterSelection, "neutrinoE/D");
    outputTree->Branch("category", &category, "category/I");

    int eventNum = 0;

    //std::unique_ptr<TChain> inputChain = std::make_unique<TChain> ("CubeEvents");

    for (int j = 100; j < 150; j++) {
        TFile inputFile(Form("/Users/gwon/Analysis/datafiles/full3DST.antineutrino.%d.cuberecon_latest.root",j));
        TFile inputGenieFile(Form("/Users/gwon/CubeAnalysis/datafiles/latest/full3DST.antineutrino.%d.rootracker.root",j+1));
        TTree* inputChain = (TTree*)inputFile.Get("CubeEvents");
        TTree* inputGenieTree = (TTree*)inputGenieFile.Get("gRooTracker");
        inputGenieTree->SetBranchAddress("EvtVtx", &EvtVtx);
        inputGenieTree->SetBranchAddress("StdHepP4", &StdHepP4);
        inputGenieTree->SetBranchAddress("StdHepStatus", &StdHepStatus);
        inputGenieTree->SetBranchAddress("StdHepPdg", &StdHepPdg);
        inputGenieTree->SetBranchAddress("StdHepN", &StdHepN);
        inputChain->SetBranchAddress("Event", &event);
        for (int i = 0; i < inputChain->GetEntries(); i++) {
            std::cout << "event: " << eventNum << std::endl;
            eventNum++;
            inputChain->GetEntry(i);
            inputGenieTree->GetEntry(i);
            nuEBeforeSelection.Fill(StdHepP4[0][3]);
            Analysis(event);
        }
    }
    outputFile->Write();
    outputFile->Close();

    TCanvas nuE;
    nuEBeforeSelection.Draw();
    nuE.SaveAs("nuEBeforeSelection.pdf");
    nuE.SaveAs("nuEBeforeSelection.C");
    nuE.Clear();
    nuEAfterSelection.Draw();
    nuE.SaveAs("nuEAfterSelection.pdf");
    nuE.SaveAs("nuEAfterSelection.C");
    nuE.Clear();
    nuEAfterSelectionForDivide.Divide(&nuEBeforeSelection);
    nuEAfterSelectionForDivide.Draw();
    nuE.SaveAs("nuEAfterSelectionForDivide.pdf");
    nuE.SaveAs("nuEAfterSelectionForDivide.C");


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
    std::cout << "reco vertex: " << muonObject->GetPosition().X() << ", " << muonObject->GetPosition().Y() << ", " << muonObject->GetPosition().Z() << std::endl;
    std::cout << "true vertex: " << EvtVtx[0]*1000 << ", " << EvtVtx[1]*1000 << ", " << EvtVtx[2]*1000 << std::endl;

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
            category = 0;
        }
        else {
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
            category = 1;
        }
        else {
            category = 3;
        }
        if (leverArm < 100)
        {
            nuEAfterSelection.Fill(StdHepP4[0][3]);
            nuEAfterSelectionForDivide.Fill(StdHepP4[0][3]);
        }
    }
    outputTree->Fill();
}
