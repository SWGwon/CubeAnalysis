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
Cube::Event* event = nullptr;
float trueNeutrinoE;
float recoNeutrinoE;
float trueNu;
float recoNu;
float trueMuonE;
float recoMuonE;
float trueNeutronKE;
float recoNeutronKE;
float leverArm;
float tof;
float trackLength;
float angle;
float eDep;
float neighborDistance;
int numberOfBranches;
int primaryPDG;
int parentPDG;
int trackNum;
int objectPDG;
//0: sig track, 1: sig cluster, 2: bkg track, 3: bkg cluster
int category;
float leptonAngle;
float leptonMomentum;
float Q2;
float Q32;

auto outputFile = std::make_shared<TFile> ("variableOutput.root","RECREATE");
auto outputTree = std::make_shared<TTree> ("tree", "tree");

float muonE;

double EvtVtx[4]; 
double StdHepP4[1000][4]; 
int StdHepStatus[1000]; 
int StdHepPdg[1000]; 
int StdHepN; 

const double THRESHOLD = 20;

void Analysis(Cube::Event* event);

int main(int argc, char** argv) {
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    outputTree->Branch("trueNeutrinoE", &trueNeutrinoE, "trueNeutrinoE/F");
    outputTree->Branch("recoNeutrinoE", &recoNeutrinoE, "recoNeutrinoE/F");
    outputTree->Branch("trueNu", &trueNu, "trueNu/F");
    outputTree->Branch("recoNu", &recoNu, "recoNu/F");
    outputTree->Branch("trueMuonE", &trueMuonE, "trueMuonE/F");
    outputTree->Branch("recoMuonE", &recoMuonE, "recoMuonE/F");
    outputTree->Branch("trueNeutronKE", &trueNeutronKE, "trueNeutronKE/F");
    outputTree->Branch("recoNeutronKE", &recoNeutronKE, "recoNeutronKE/F");
    outputTree->Branch("leverArm", &leverArm, "leverArm/F");
    outputTree->Branch("tof", &tof, "tof/F");
    outputTree->Branch("trackLength", &trackLength, "trackLength/F");
    outputTree->Branch("angle", &angle, "angle/F");
    outputTree->Branch("eDep", &eDep, "eDep/F");
    outputTree->Branch("neighborDistance", &neighborDistance, "neighborDistance/F");
    outputTree->Branch("trackNum", &trackNum, "trackNum/I");
    outputTree->Branch("category", &category, "category/I");
    outputTree->Branch("primaryPDG", &primaryPDG, "primaryPDG/I");
    outputTree->Branch("parentPDG", &parentPDG, "parentPDG/I");
    outputTree->Branch("objectPDG", &objectPDG, "objectPDG/I");
    outputTree->Branch("numberOfBranches", &numberOfBranches, "numberOfBranches/I");
    outputTree->Branch("leptonAngle", &leptonAngle, "leptonAngle/F");
    outputTree->Branch("leptonMomentum", &leptonMomentum, "leptonMomentum/F");
    outputTree->Branch("Q2", &Q2, "Q2/F");
    outputTree->Branch("Q32", &Q32, "Q32/F");

    int eventNum = 0;
    int fileNum = 0;
    std::cout << "file num?" << std::endl;
    std::cin >> fileNum;
    double avgTime = 0;

    //std::unique_ptr<TChain> inputChain = std::make_unique<TChain> ("CubeEvents");

    for (int j = 0; j < fileNum; j++) {
        //TFile inputFile(Form("/Users/gwon/Analysis/datafiles/newfile/full3DST.antineutrino.%d.cuberecon_Jan1_2021.root",j));
        TFile inputFile(Form("/Users/gwon/Analysis/datafiles/CC/full3DST.antineutrino.%d.cuberecon_Jan1_2021_CCSelected.root",j));
        TFile inputGenieFile(Form("/Users/gwon/Analysis/datafiles/CC/full3DST.antineutrino.%d.rootracker_CCSelected.root",j+1));
        //TFile inputGenieFile(Form("/Users/gwon/CubeAnalysis/datafiles/latest/full3DST.antineutrino.%d.rootracker.root",j+1));
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

    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::cout << "took " << sec.count() << "sec" << std::endl;

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
        if (traj && traj->GetPDGCode() == -13 && traj->GetParentId() == -1) {
            muonLike++;
            muonE = traj->GetInitialMomentum().E();
            muonObject = track;
        }
        if (muonLike > 1) {
            throw std::runtime_error("more than 1 anti muon track");
        }
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
        //if (Cube::Tool::AreNeighboringObjects(*muonObject, *o)) {
        //    numberOfMuonAssociated++;
        //    continue;
        //}
        TVector3 muonVertex;
        muonVertex.SetX(muonObject->GetPosition().X());
        muonVertex.SetY(muonObject->GetPosition().Y());
        muonVertex.SetZ(muonObject->GetPosition().Z());
        TVector3 objectPosition;
        Cube::Handle<Cube::ReconTrack> tempTrack = o;
    //    Cube::Handle<Cube::ReconCluster> tempCluster = o;
        if (!tempTrack)
            continue;
        //objectPosition.SetX(tempTrack? tempTrack->GetPosition().X() : tempCluster->GetPosition().X());
        //objectPosition.SetY(tempTrack? tempTrack->GetPosition().Y() : tempCluster->GetPosition().Y());
        //objectPosition.SetZ(tempTrack? tempTrack->GetPosition().Z() : tempCluster->GetPosition().Z());
        objectPosition.SetX(tempTrack->GetPosition().X());
        objectPosition.SetY(tempTrack->GetPosition().Y());
        objectPosition.SetZ(tempTrack->GetPosition().Z());
        //if ((objectPosition - muonVertex).Mag() < 40) {
        if ((objectPosition - muonVertex).Mag() < 40) {
        //if ((objectPosition - muonVertex).Mag() < 20) {
            numberOfMuonAssociated++;
            continue;
        }
    }
    return numberOfMuonAssociated;
}

Cube::Handle<Cube::ReconObject> GetEarliestObject(const Cube::Handle<Cube::ReconTrack>& muonObject,
                                                  Cube::Handle<Cube::ReconObjectContainer>& objects,
                                                  double threshold) {
    double earliestTime = 1E+8;
    Cube::Handle<Cube::ReconObject> earliestObject;
    for (auto& o : *objects) {
        if (Cube::Tool::AreNeighboringObjects(*muonObject, *o, 30)) {
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

    trueNeutrinoE = -10;
    recoNeutrinoE = -10;
    trueNu = -10;
    recoNu = -10;
    trueMuonE = -10;
    recoMuonE = -10;
    trueNeutronKE = -10;
    recoNeutronKE = -10;
    leverArm = -10;
    tof = -10;
    trackLength = -10;
    angle = -10;
    eDep = -10;
    neighborDistance = -10;
    trackNum = -10;
    category = -10;
    primaryPDG = -10;
    parentPDG = -10;
    objectPDG = -10;
    numberOfBranches = -10;
    leptonAngle = -10;
    leptonMomentum = -10;
    Q2 = -10;
    Q32 = -10;

    for (int i = 0; i < StdHepN; ++i) {
        if (StdHepStatus[i] == 1 && StdHepPdg[i] == -13) {
            TVector3 muonMomentum;
            muonMomentum.SetX(StdHepP4[i][0]);
            muonMomentum.SetY(StdHepP4[i][1]);
            muonMomentum.SetZ(StdHepP4[i][2]);
            leptonAngle = muonMomentum.Angle(beamDirection);
            leptonMomentum = muonMomentum.Mag();
            Q2 = 4*StdHepP4[i][3]*StdHepP4[0][3]*std::pow(TMath::Sin(leptonAngle/2),2);
            Q32 = std::pow(std::pow(StdHepP4[0][0],2) + std::pow(StdHepP4[0][1],2) + std::pow(StdHepP4[0][2],2),0.5) - leptonMomentum;
        }
    }

    Cube::Event::G4TrajectoryContainer trajectories = event->G4Trajectories;
    Cube::Handle<Cube::ReconObjectContainer> objects = event->GetObjectContainer();
    if (!objects) 
        return;

    //single muon track selection
    //There should be only one muon track object
    Cube::Handle<Cube::ReconTrack> muonObject;
    try {
        muonObject = GetMuonObject(event, trajectories, objects);
    } catch (std::runtime_error e) {
        std::cout << e.what() << std::endl;
        return;
    }

    //single track selection
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

    //earliestObject check
    Cube::Handle<Cube::ReconTrack> earliestTrack = earliestObject;
    Cube::Handle<Cube::ReconCluster> earliestCluster = earliestObject;
    if (!earliestCluster && !earliestTrack)
        return;

    double tempAvgNeighborDistance = 1E+8;
    int tempNumberOfBranches = 0;
    for (auto& o : *objects) {
        if (o == muonObject) continue;
        if (o == earliestObject) continue;
        double distance = 10000;
        distance = Cube::Tool::DistanceBetweenObjects(*earliestObject, *o);
        if (distance < 18)
            tempNumberOfBranches++;
    }
    numberOfBranches = tempNumberOfBranches;

    double muonTime = muonObject->GetPosition().T();

    //tof;
    tof = (earliestTrack ? 
            earliestTrack->GetMedian().T() - muonTime : 
            earliestCluster->GetMedian().T() - muonTime);
    //if (tof < 0)
        //return;
    //leverArm;
    leverArm = (earliestTrack? 
            (earliestTrack->GetPosition().Vect() - muonObject->GetPosition().Vect()).Mag() :
            (earliestCluster->GetPosition().Vect() - muonObject->GetPosition().Vect()).Mag());
    //angle;
    angle = (earliestTrack?
            TMath::Cos((earliestTrack->GetPosition() - muonObject->GetPosition()).Angle(beamDirection)) :
            TMath::Cos((earliestCluster->GetPosition() - muonObject->GetPosition()).Angle(beamDirection)));
    //eDep;
    eDep = (earliestTrack?
            earliestTrack->GetEDeposit() :
            earliestCluster->GetEDeposit());
    //neighborDistance;
    neighborDistance = GetNeighborDistance(muonObject, earliestObject, objects);
    //trueMuonE 
    trueMuonE = muonE;
    //recoMuonE
    recoMuonE = muonE * gRandom->Gaus(1, 0.04);
    //recoNeutronKE
    double recoBeta = (leverArm/tof)/300;
    recoNeutronKE = 939.565*(1./std::pow(1.-std::pow(recoBeta,2),0.5)-1.); 
    //trueNeutrinoE
    trueNeutrinoE = StdHepP4[0][3]*1000.;
    //recoNeutrinoE
    recoNeutrinoE = recoNeutronKE + recoMuonE + 40;
    //trueNu
    for (int k = 0; k < StdHepN; k++)
    {
        if (StdHepPdg[k] == -13)
        {
            trueNu = StdHepP4[0][3]*1000. - StdHepP4[k][3]*1000.;
            break;
        }
    }
    //recoNu
    recoNu = recoNeutronKE;
    //trueNeutronKE
    trueNeutronKE = -10;
    trackNum = NumberOfAssociated(muonObject, objects);

    int earliestTrajID = Cube::Tool::MainTrajectory(*event,*earliestObject);
    int earliestPrim = Cube::Tool::PrimaryId(*event,earliestTrajID);

    //parentPDG
    int parentId = earliestTraj->GetParentId();
    int parentPdg = 0;
    if (parentId == -1) {
        parentPdg = earliestTraj->GetPDGCode();
        //std::cout << "parentPdg: " << parentPdg << std::endl;
        //std::cout << "trajectories[earliestPrim]->GetPDGCode()" << trajectories[earliestPrim]->GetPDGCode() << std::endl;
    } else {
        parentPdg = trajectories[parentId]->GetPDGCode();
    }
    objectPDG = earliestTraj->GetPDGCode();
    parentPDG = parentPdg;
    primaryPDG = trajectories[earliestPrim]->GetPDGCode();

    //0: sig track, 1: sig cluster, 2: bkg track, 3: bkg cluster
    if (earliestTrack) {
        trackLength = GetTrackLength(earliestTrack);
        if (trajectories[earliestPrim]->GetPDGCode() == 2112 && parentPdg == 2112) {
            category = 0;
            trueNeutronKE = trajectories[earliestPrim]->GetInitialMomentum().E() - 939.565;
            std::cout << "trueNeutronKE: " << trueNeutronKE << std::endl;
        }
        else {
            category = 2;
        }
    }

    if (earliestCluster) {
        if (trajectories[earliestPrim]->GetPDGCode() == 2112 && parentPdg == 2112) {
            category = 1;
            trueNeutronKE = trajectories[earliestPrim]->GetInitialMomentum().E() - 939.565;
            std::cout << "trueNeutronKE: " << trueNeutronKE << std::endl;
        } else {
            category = 3;
        }
    }
    outputTree->Fill();
}
