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

#include <TRandom.h>
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

TVector3 beamDirection(0,0,1);
Cube::Event* event = NULL;
float leverArm;
float eDep;
float trackLength;
float angle;
float nbhdist;
float trueNeutrinoE;
float recoNeutrinoE;
float recoNeutronKE;
float trueNeutronKE;
float genieNu;
float recoNuTrueTof;
float recoNuRecoTof;
int category; //0: sig track, 1: sig cluster, 2: bkg track, 3: bkg cluster
auto outputFile = std::make_shared<TFile> ("variableOutput.root","RECREATE");
auto outputTree = std::make_shared<TTree> ("tree", "tree");

float muonE;

double EvtVtx[4]; 
double StdHepP4[1000][4]; 
int StdHepStatus[1000]; 
int StdHepPdg[1000]; 
int StdHepN; 

TH1D nuEBeforeSelection("","neutrinoE before selection",100,0,10);
TH1D nuEAfterSelection("","neutrinoE after selection",100,0,10);
TH1D nuEAfterSelectionForDivide("","neutrinoE after selection",100,0,10);

TH2D true_reco_nu("","true vs reco nu;reco nu;true nu",50,0,1000,50,0,1000);

const double THRESHOLD = 20;

void Analysis(Cube::Event* event);

int main(int argc, char** argv) {
    std::string inputFileName;

    while (optind < argc) {
        inputFileName = argv[optind];
        optind++;
    }
    outputTree->Branch("leverArm", &leverArm, "lever arm/F");
    outputTree->Branch("eDep", &eDep, "edeposit/F");
    outputTree->Branch("trackLength", &trackLength, "track length/F");
    outputTree->Branch("angle", &angle, "cos(angle)/F");
    outputTree->Branch("nbhdist", &nbhdist, "neighbor distance/F");
    outputTree->Branch("trueNeutrinoE", &trueNeutrinoE, "neutirnoEe/F");
    outputTree->Branch("recoNeutrinoE", &recoNeutrinoE, "neutrinoE/F");
    outputTree->Branch("recoNeutronKE", &recoNeutronKE, "recoNeutronKE/F");
    outputTree->Branch("trueNeutronKE", &trueNeutronKE, "trueNeutronKE/F");
    outputTree->Branch("genieNu", &genieNu, "genieNu/F");
    outputTree->Branch("recoNuRecoTof", &recoNuRecoTof, "recoNuRecoTof/F");
    outputTree->Branch("recoNuTrueTof", &recoNuTrueTof, "recoNuTrueTof/F");
    outputTree->Branch("category", &category, "category/I");

    int eventNum = 0;
    int fileNum = 0;
    std::cout << "file num?" << std::endl;
    std::cin >> fileNum;

    //std::unique_ptr<TChain> inputChain = std::make_unique<TChain> ("CubeEvents");

    for (int j = 100; j < fileNum+100; j++) {
        TFile inputFile(Form("/Users/gwon/Analysis/datafiles/full3DST.antineutrino.%d.cuberecon_latest.root",j));
        TFile inputGenieFile(Form("/Users/gwon/CubeAnalysis/datafiles/latest/full3DST.antineutrino.%d.rootracker.root",j+1));
        if (!inputFile.IsOpen() || !inputGenieFile.IsOpen())
            continue;
        if (inputFile.TestBit(TFile::kRecovered) || inputGenieFile.TestBit(TFile::kRecovered))
            continue;
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
    nuE.Clear();
    nuE.SetLogz();
    true_reco_nu.Draw("colz");
    nuE.SaveAs("true_reco_nu.pdf");
    nuE.SaveAs("true_reco_nu.C");

    return 0;
}

Cube::Handle<Cube::ReconTrack> GetMuonObject(Cube::Event* event, 
                                             Cube::Event::G4TrajectoryContainer& trajectories, 
                                             Cube::Handle<Cube::ReconObjectContainer>& objects) {
    int muonLike = 0;
    Cube::Handle<Cube::ReconTrack> muonObject;
    for (auto& o : *objects) {
        Cube::Handle<Cube::ReconTrack> track = o;
        if (!track) 
            continue;
        int mainTraj = Cube::Tool::MainTrajectory(*event, *track);
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
        if (!traj || traj->GetPDGCode() != -13 || traj->GetParentId() != -1) 
            continue;
        muonLike++;
        muonE = traj->GetInitialMomentum().E();
        muonObject = track;
    }

    if (muonLike > 1) {
        throw std::runtime_error("more than 1 anti muon track");
    }
    if (!muonObject) {
        throw std::runtime_error("no anti muon track");
    }

    return muonObject;
}

int NumberOfAssociated(const Cube::Handle<Cube::ReconTrack>& muonObject, 
                       Cube::Handle<Cube::ReconObjectContainer>& objects) {
    int numberOfMuonAssociated = 0;
    for (auto& o : *objects) {
        if (Cube::Tool::AreNeighboringObjects(*muonObject, *o)) {
            numberOfMuonAssociated++;
            continue;
        }
        //TVector3 muonVertex;
        //muonVertex.SetX(muonObject->GetPosition().X());
        //muonVertex.SetY(muonObject->GetPosition().Y());
        //muonVertex.SetZ(muonObject->GetPosition().Z());
        //TVector3 objectPosition;
        //Cube::Handle<Cube::ReconTrack> tempTrack = o;
        //Cube::Handle<Cube::ReconCluster> tempCluster = o;
        //if (!tempCluster && !tempTrack)
        //    continue;
        //objectPosition.SetX(tempTrack? tempTrack->GetPosition().X() : tempCluster->GetPosition().X());
        //objectPosition.SetY(tempTrack? tempTrack->GetPosition().Y() : tempCluster->GetPosition().Y());
        //objectPosition.SetZ(tempTrack? tempTrack->GetPosition().Z() : tempCluster->GetPosition().Z());
        //if ((objectPosition - muonVertex).Mag() < 40) {
        //    numberOfMuonAssociated++;
        //    continue;
        //}
    }
    return numberOfMuonAssociated;
}

Cube::Handle<Cube::ReconObject> GetEarliestObject(const Cube::Handle<Cube::ReconTrack>& muonObject,
                                                  Cube::Handle<Cube::ReconObjectContainer>& objects,
                                                  double threshold) {
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
        if (objTime < earliestTime && objEDep > threshold) {
            earliestTime = objTime;
            earliestObject = o;
        }
    }
    if (!earliestObject) {
        throw std::runtime_error("no earliest object candidate");
    }
    return earliestObject;
}

Cube::Handle<Cube::G4Trajectory> GetObjectTrajectory(Cube::Event* event,  
                                                     Cube::Handle<Cube::ReconObject>& object,
                                                     Cube::Event::G4TrajectoryContainer& trajectories) {
    int mainTraj = Cube::Tool::MainTrajectory(*event, *object);
    Cube::Handle<Cube::G4Trajectory> Traj = trajectories[mainTraj];
    if (!Traj) {
        throw std::runtime_error("no Traj corresponding object");
    } else {
        return Traj;
    }
}

double GetNeighborDistance(Cube::Handle<Cube::ReconTrack>& muonObject, 
                           Cube::Handle<Cube::ReconObject>& earliestObject,
                           Cube::Handle<Cube::ReconObjectContainer>& objects) {
    double neighborDistance = 1E+8;
    for (auto& o : *objects) {
        if (o == muonObject) continue;
        if (o == earliestObject) continue;
        double distance = Cube::Tool::DistanceBetweenObjects(*earliestObject, *o);
        if (distance < neighborDistance) {
            neighborDistance = distance;
        }
    }
    if (neighborDistance == 1E+8) {
        neighborDistance = -1;
    }
    return neighborDistance;
}

double GetTrackLength(Cube::Handle<Cube::ReconTrack>& earliestTrack) {
    double tempTrackLength = -1;
    Cube::ReconNodeContainer::iterator n = earliestTrack->GetNodes().begin();
    Cube::Handle<Cube::TrackState> lastState = (*(n++))->GetState();
    while (n != earliestTrack->GetNodes().end()) {
        Cube::Handle<Cube::TrackState> nodeState = (*(n++))->GetState();
        tempTrackLength += (nodeState->GetPosition().Vect()
                - lastState->GetPosition().Vect()).Mag();
        lastState = nodeState;
    }
    return tempTrackLength;
}

void Analysis(Cube::Event* event) {

    leverArm = -1;
    eDep = -1;
    trackLength = -1;
    angle = -1;
    nbhdist = -1;
    trueNeutrinoE = -1;
    recoNeutrinoE = -1;
    recoNeutronKE = -1;
    trueNeutronKE = -1;

    Cube::Event::G4TrajectoryContainer trajectories = event->G4Trajectories;
    Cube::Handle<Cube::ReconObjectContainer> objects = event->GetObjectContainer();
    if (!objects) 
        return;

    //single muon track selection
    Cube::Handle<Cube::ReconTrack> muonObject;
    try {
        muonObject = GetMuonObject(event, trajectories, objects);
    } catch (std::runtime_error e) {
        std::cout << e.what() << std::endl;
        return;
    }

    //single track selectrion
    if (NumberOfAssociated(muonObject, objects) != 1)
        return;

    //earliestObject selection
    Cube::Handle<Cube::ReconObject> earliestObject;
    try {
        earliestObject = GetEarliestObject(muonObject, objects, 20);
    } catch (std::runtime_error e) {
        std::cout << e.what() << std::endl;
        return;
    }

    //earliestObject trajectory
    Cube::Handle<Cube::G4Trajectory> earliestTraj;
    try {
        earliestTraj = GetObjectTrajectory(event, earliestObject, trajectories);
    } catch (std::runtime_error e) {
        std::cout << e.what() << std::endl;
        return;
    }

    Cube::Handle<Cube::ReconTrack> earliestTrack = earliestObject;
    Cube::Handle<Cube::ReconCluster> earliestCluster = earliestObject;
    if (!earliestCluster && !earliestTrack)
        return;

    double muonTime = muonObject->GetPosition().T();
    double recoTof = (earliestTrack ? 
                      earliestTrack->GetMedian().T() - muonTime : 
                      earliestCluster->GetMedian().T() - muonTime);

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
    double trueNu = 0;
    for (int k = 0; k < StdHepN; k++)
    {
        if (StdHepPdg[k] == -13)
        {
            trueNu = StdHepP4[0][3] - StdHepP4[k][3];
            break;
        }
    }

    double trueTof = earliestTruth - muonTruth; 
    double trueLeverArm = (earliestVector - muonVector).Mag();
    leverArm = (earliestTrack? 
                (earliestTrack->GetPosition().Vect() - muonObject->GetPosition().Vect()).Mag() :
                (earliestCluster->GetPosition().Vect() - muonObject->GetPosition().Vect()).Mag());
    eDep = (earliestTrack?
            earliestTrack->GetEDeposit() :
            earliestCluster->GetEDeposit());
    angle = (earliestTrack?
             TMath::Cos((earliestTrack->GetPosition() - muonObject->GetPosition()).Angle(beamDirection)) :
             TMath::Cos((earliestCluster->GetPosition() - muonObject->GetPosition()).Angle(beamDirection)));
    nbhdist = GetNeighborDistance(muonObject, earliestObject, objects);

    double recoBeta = (leverArm/recoTof)/300;
    double trueBeta = (trueLeverArm/trueTof)/300;
    
    recoNeutronKE = 939.565*(1./std::pow(1.-std::pow(recoBeta,2),0.5)-1.); 
    trueNeutronKE = 939.565*(1./std::pow(1.-std::pow(trueBeta,2),0.5)-1.); 
    trueNeutrinoE = StdHepP4[0][3]*1000.;
    recoNeutrinoE = recoNeutronKE + muonE * gRandom->Gaus(1,0.04) + 40;
    nuEAfterSelection.Fill(StdHepP4[0][3]);
    nuEAfterSelectionForDivide.Fill(StdHepP4[0][3]);
    genieNu = trueNu*1000. - 40;
    recoNuRecoTof = recoNeutronKE; 
    recoNuTrueTof = trueNeutronKE; 

    int parentId = earliestTraj->GetParentId();
    if (parentId > trajectories.size())
        return;
    int parentPdg = trajectories[parentId]->GetPDGCode();

    if (earliestTrack) {
        trackLength = GetTrackLength(earliestTrack);
        if (earliestTraj->GetPDGCode() == 2112 || parentPdg == 2112) {
            category = 0;
        }
        else {
            category = 2;
        }
    }

    if (earliestCluster) {
        if (earliestTraj->GetPDGCode() == 2112 || parentPdg == 2112) {
            category = 1;
        } else {
            category = 3;
        }
        nuEBeforeSelection.Fill(StdHepP4[0][3]);
        
        //selection
        //if (eDep > 600 && nbhdist < 130 && angle < 0.1) {
            std::cout << "leverArm: " << leverArm << std::endl;
            std::cout << "eDep: " << eDep << std::endl;
            std::cout << "trackLength: " << trackLength << std::endl;
            std::cout << "angle: " << TMath::Cos(angle) << std::endl;
            std::cout << "neighborDistance: " << nbhdist << std::endl;
            std::cout << " asdasdasdasdasda " << std::endl;
            true_reco_nu.Fill(recoNeutronKE, trueNu*1000);
        //}
    }
            outputTree->Fill();
}
