#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconCluster.hxx>
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

TH2D histDeltaT("","true #DeltaT vs reco #DeltaT;reco #DeltaT;true #DeltaT",100,-10,10,50,0,10);
TH2D histMedianDeltaT("","true #DeltaMedianT vs reco #DeltaMedianT;reco #DeltaMedianT;true #DeltaT",100,-10,10,50,0,10);

void DeltaT(Cube::Event* event);
void DeltaMedianT(Cube::Event* event);

int main(int argc, char* argv[])
{
    std::string inputNamesCubeReconFile;

    auto neutronTrueKE = std::make_unique<TH1D> ("","neutron True KE",100,0,2000);

    if (argc <= optind) throw std::runtime_error("Missing input file");

    while (optind < argc) {
        inputNamesCubeReconFile = argv[1];
        optind++;
    }

    Cube::Event* event = NULL;
    std::unique_ptr<TChain> inputChainCubeRecon = std::make_unique<TChain> ("CubeEvents");
    inputChainCubeRecon->Add(inputNamesCubeReconFile.c_str());
    inputChainCubeRecon->SetBranchAddress("Event",&event);

    for (int i = 0; i < inputChainCubeRecon->GetEntries(); i++)
    {
        std::cout << "event: " << i << std::endl;
        inputChainCubeRecon->GetEntry(i);
        try
        {
            DeltaT(event);
        }
        catch(...)
        {
            continue;
        }
        try
        {
            DeltaMedianT(event);
        }
        catch(...)
        {
            continue;
        }
    }

    TCanvas can;
    can.SetLogz();
    histDeltaT.SetStats(false);
    histDeltaT.Draw("colz");
    can.SaveAs("histDeltaT.C");

    TCanvas can1;
    can1.SetLogz();
    histMedianDeltaT.SetStats(false);
    histMedianDeltaT.Draw("colz");
    can1.SaveAs("histMedianDeltaT.C");

    return 0;
}

void DeltaT(Cube::Event* event)
{
    Cube::Event::G4TrajectoryContainer& trajectories = event->G4Trajectories;

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

    if (!electrons.empty()) return;
    if (!muons.empty()) return;
    if (antiMuons.empty()) return;
    if (!pions.empty()) return;
    if (!protons.empty()) return;
    if (neutrons.size() != 1) return;
    if (!gammas.empty()) return;
    if (!others.empty()) return;

    Cube::Handle<Cube::ReconObjectContainer> objects
        = event->GetObjectContainer();
    if (!objects) return;

    // Make sure there's a single muon track.
    int muonLike = 0;
    Cube::Handle<Cube::ReconTrack> muonObject;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
            o != objects->end(); ++o) {
        Cube::Handle<Cube::ReconTrack> track = *o;
        if (!track) {
            continue;
        }
        int mainTraj = Cube::Tool::MainTrajectory(*event,*track);
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
        if (!traj) {
            continue;
        }
        if (std::abs(traj->GetPDGCode()) != 13) {
            continue;
        }
        ++muonLike;
        muonObject = track;
    }

    if (muonLike > 1) return;
    if (!muonObject) {
        return;
    }

    // Collect the earliest objects and the muon.
    double muonTime = muonObject->GetPosition().T();
    double muonMedianTime = muonObject->GetMedian().T();
    double earliestTrajTime = 1E+8;
    double earliestTime = 1E+8;
    Cube::Handle<Cube::ReconObject> earliestObject;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
            o != objects->end(); ++o) {
        if (Cube::Tool::AreNeighboringObjects(*muonObject,*(*o))) continue;
        Cube::Handle<Cube::ReconTrack> track = *o;
        double objTime = -1.0;
        if (track) {
            if (track == muonObject) continue;
            objTime = track->GetPosition().T();
        }
        Cube::Handle<Cube::ReconCluster> cluster = *o;
        if (cluster) {
            objTime = cluster->GetPosition().T();
        }
        if (objTime < 0) continue;
        if (objTime < earliestTime) {
            earliestTime = objTime;
            int mainTraj = Cube::Tool::MainTrajectory(*event,*(*o));
            Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
            earliestTrajTime = traj->GetInitialPosition().T();
            earliestObject = *o;
        }

#ifdef DUMP_STACK
        int mainTraj = Cube::Tool::MainTrajectory(event,*(*o));
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
        int parentTraj = mainTraj;
        std::cout << "Object stack: "
            while (parentTraj >= 0) {
                traj = trajectories[parentTraj];
                std::cout << parentTraj << "(" << traj->GetPDGCode()
                    << "," << traj->GetInitialPosition().T()
                    << "," << (traj->GetInitialMomentum().E()
                            - traj->GetInitialMomentum().M()) << ") ";
                parentTraj = traj->GetParentId();
            }
        std::cout << std::endl;
#endif
    }
    if (!earliestObject) return;

    double earliestTrueTime = 1E+8;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
            o != objects->end(); ++o) 
    {
        if (*o == muonObject)
        {
            continue;
        }
        std::vector<Cube::Handle<Cube::G4Hit>> objectSegs
            = Cube::Tool::ObjectG4Hits(*event,*(*o));
        for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator
                t = objectSegs.begin();
                t != objectSegs.end(); ++t) {
            if ((*t)->GetStart().T() < earliestTrueTime) {
                earliestTrueTime = (*t)->GetStart().T();
            }
        }
    }
    int earliestTraj = Cube::Tool::MainTrajectory(*event,*earliestObject);
    int earliestPrim = Cube::Tool::PrimaryId(*event,earliestTraj);
    std::cout << "earliestPrim: " << earliestPrim << std::endl;
    if (earliestPrim != 1)
    {
        return;
    }

    double recoDeltaT = earliestTime - muonTime;
    double trueDeltaT = earliestTrueTime - 1;

    std::cout << "recoDeltaT: " << recoDeltaT << std::endl;
    std::cout << "trueDeltaT: " << trueDeltaT << std::endl;

    histDeltaT.Fill(recoDeltaT, trueDeltaT);
}

void DeltaMedianT(Cube::Event* event)
{
    Cube::Event::G4TrajectoryContainer& trajectories = event->G4Trajectories;

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

    if (!electrons.empty()) return;
    if (!muons.empty()) return;
    if (antiMuons.empty()) return;
    if (!pions.empty()) return;
    if (!protons.empty()) return;
    if (neutrons.size() != 1) return;
    if (!gammas.empty()) return;
    if (!others.empty()) return;

    Cube::Handle<Cube::ReconObjectContainer> objects
        = event->GetObjectContainer();
    if (!objects) return;

    // Make sure there's a single muon track.
    int muonLike = 0;
    Cube::Handle<Cube::ReconTrack> muonObject;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
            o != objects->end(); ++o) {
        Cube::Handle<Cube::ReconTrack> track = *o;
        if (!track) {
            continue;
        }
        int mainTraj = Cube::Tool::MainTrajectory(*event,*track);
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
        if (!traj) {
            continue;
        }
        if (std::abs(traj->GetPDGCode()) != 13) {
            continue;
        }
        ++muonLike;
        muonObject = track;
    }

    if (muonLike > 1) return;
    if (!muonObject) {
        return;
    }

    // Collect the earliest objects and the muon.
    double muonTime = muonObject->GetPosition().T();
    double earliestTrajTime = 1E+8;
    double earliestTime = 1E+8;
    Cube::Handle<Cube::ReconObject> earliestObject;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
            o != objects->end(); ++o) {
        if (Cube::Tool::AreNeighboringObjects(*muonObject,*(*o))) continue;
        Cube::Handle<Cube::ReconTrack> track = *o;
        double objTime = -1.0;
        if (track) {
            if (track == muonObject) continue;
            objTime = track->GetMedian().T();
        }
        Cube::Handle<Cube::ReconCluster> cluster = *o;
        if (cluster) {
            objTime = cluster->GetMedian().T();
        }
        if (objTime < 0) continue;
        if (objTime < earliestTime) {
            earliestTime = objTime;
            int mainTraj = Cube::Tool::MainTrajectory(*event,*(*o));
            Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
            earliestTrajTime = traj->GetInitialPosition().T();
            earliestObject = *o;
        }
    }
    if (!earliestObject) return;

    double earliestTrueTime = 1E+8;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
            o != objects->end(); ++o) 
    {
        if (*o == muonObject)
        {
            continue;
        }
        std::vector<Cube::Handle<Cube::G4Hit>> objectSegs
            = Cube::Tool::ObjectG4Hits(*event,*(*o));
        for (std::vector<Cube::Handle<Cube::G4Hit>>::iterator
                t = objectSegs.begin();
                t != objectSegs.end(); ++t) {
            if ((*t)->GetStart().T() < earliestTrueTime) {
                earliestTrueTime = (*t)->GetStart().T();
            }
        }
    }
    int earliestTraj = Cube::Tool::MainTrajectory(*event,*earliestObject);
    int earliestPrim = Cube::Tool::PrimaryId(*event,earliestTraj);
    std::cout << "earliestPrim: " << earliestPrim << std::endl;
    if (earliestPrim != 1)
    {
        return;
    }

    double recoDeltaT = earliestTime - muonTime;
    double trueDeltaT = earliestTrueTime - 1;

    std::cout << "recoDeltaMedianT: " << recoDeltaT << std::endl;
    std::cout << "trueDeltaMedianT: " << trueDeltaT << std::endl;

    histMedianDeltaT.Fill(recoDeltaT, trueDeltaT);
}
