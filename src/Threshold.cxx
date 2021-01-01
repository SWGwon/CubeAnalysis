//threshold : 20pe , 100 pe
//ratio signal:background
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

TH2D histMedianDeltaT_20("","true #DeltaT vs reco #DeltaMedianT, threshold: 20;reco #DeltaMedianT;true #DeltaT",100,-10,10,50,0,10);
TH2D histMedianDeltaT_50("","true #DeltaT vs reco #DeltaMedianT, threshold: 50;reco #DeltaMedianT;true #DeltaT",100,-10,10,50,0,10);
TH2D histMedianDeltaT_200("","true #DeltaT vs reco #DeltaMedianT, threshold: 200;reco #DeltaMedianT;true #DeltaT",100,-10,10,50,0,10);

Cube::Event* event = NULL;
double eff_20 = 0;
double eff_50 = 0;
double eff_200 = 0;
double numSignal_20 = 0;
double numSignal_50 = 0;
double numSignal_200 = 0;
double numTotalSignal_20 = 0;
double numTotalSignal_50 = 0;
double numTotalSignal_200 = 0;
double numBackground = 0;
const int THRESHOLD = 100;

void Analysis(Cube::Event* event, int threshold);

int NumberOfTrack();

int main(int argc, char** argv)
{
    std::string inputFileName;

    while (optind < argc)
    {
        inputFileName = argv[optind];
        optind++;
    }
    std::unique_ptr<TChain> inputChain = std::make_unique<TChain> ("CubeEvents");
    inputChain->Add(inputFileName.c_str());

    inputChain->SetBranchAddress("Event", &event);

    for (int i = 0; i < inputChain->GetEntries(); i++)
    {
        std::cout << "event: " << i << std::endl;
        inputChain->GetEntry(i);
        try
        {
            Analysis(event, 20);
        }
        catch(...)
        {
            continue;
        }
    }
    for (int i = 0; i < inputChain->GetEntries(); i++)
    {
        std::cout << "event: " << i << std::endl;
        inputChain->GetEntry(i);
        try
        {
            Analysis(event, 50);
        }
        catch(...)
        {
            continue;
        }
    }
    for (int i = 0; i < inputChain->GetEntries(); i++)
    {
        std::cout << "event: " << i << std::endl;
        inputChain->GetEntry(i);
        try
        {
            Analysis(event, 200);
        }
        catch(...)
        {
            continue;
        }
    }

    TCanvas can1;
    can1.SetLogz();
    histMedianDeltaT_20.SetStats(false);
    histMedianDeltaT_20.Draw("colz");
    can1.SaveAs(Form("histMedianDeltaT_with_threshold_%d.C",20));
    can1.SaveAs(Form("histMedianDeltaT_with_threshold_%d.pdf",20));
    TCanvas can2;
    can2.SetLogz();
    histMedianDeltaT_50.SetStats(false);
    histMedianDeltaT_50.Draw("colz");
    can2.SaveAs(Form("histMedianDeltaT_with_threshold_%d.C",50));
    can2.SaveAs(Form("histMedianDeltaT_with_threshold_%d.pdf",50));
    TCanvas can3;
    can3.SetLogz();
    histMedianDeltaT_200.SetStats(false);
    histMedianDeltaT_200.Draw("colz");
    can3.SaveAs(Form("histMedianDeltaT_with_threshold_%d.C",200));
    can3.SaveAs(Form("histMedianDeltaT_with_threshold_%d.pdf",200));
    eff_20 = numSignal_20/numTotalSignal_20;
    eff_50 = numSignal_50/numTotalSignal_50;
    eff_200 = numSignal_200/numTotalSignal_200;
    std::cout << "eff threshold 20 : " << eff_20 << std::endl;
    std::cout << "eff threshold 50 : " << eff_50 << std::endl;
    std::cout << "eff threshold 200 : " << eff_200 << std::endl;


    return 0;
}

void Analysis(Cube::Event* event, int threshold)
{
    Cube::Event::G4TrajectoryContainer& trajectories = event->G4Trajectories;
    Cube::Handle<Cube::ReconObjectContainer> objects;
    try{
        objects = event->GetObjectContainer();
    } catch (...) {
        return;
    }
    if (!objects) return;

    // Make sure there's a single muon track.
    int muonLike = 0;
    Cube::Handle<Cube::ReconTrack> muonObject;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
            o != objects->end(); ++o) {
        Cube::Handle<Cube::ReconTrack> track;
        try{
            track= *o;
        } catch(...) {
            continue;
        }
        if (!track) {
            continue;
        }
        int mainTraj = Cube::Tool::MainTrajectory(*event,*track);
        Cube::Handle<Cube::G4Trajectory> traj;
        try{
            traj = trajectories[mainTraj];
        } catch (...) {
            continue;
        }
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

    int numberOfMuonAssociated = 0;
    // Collect the earliest objects and the muon.
    double muonTime = muonObject->GetPosition().T();
    double earliestTrajTime = 1E+8;
    double earliestTime = 1E+8;

    Cube::Handle<Cube::ReconObject> earliestObject;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
            o != objects->end(); ++o) 
    {
        if (Cube::Tool::AreNeighboringObjects(*muonObject,*(*o))) 
        {
            numberOfMuonAssociated++;
            continue;
        }
        Cube::Handle<Cube::ReconTrack> track = *o;
        double objTime = -1.0;
        double objEDep = -1.0;
        if (track) 
        {
            if (track == muonObject) 
                continue;
            objTime = track->GetMedian().T();
            objEDep = track->GetEDeposit();
        }
        Cube::Handle<Cube::ReconCluster> cluster = *o;
        if (cluster) 
        {
            objTime = cluster->GetMedian().T();
            objEDep = cluster->GetEDeposit();
        }
        if (objTime < 0) continue;
        if (objTime < earliestTime && objEDep > threshold) 
        {
            earliestTime = objTime;
            int mainTraj = Cube::Tool::MainTrajectory(*event,*(*o));
            Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
            earliestTrajTime = traj->GetInitialPosition().T();
            earliestObject = *o;
        }
    }
    if (!earliestObject) return;
    if (numberOfMuonAssociated != 1) return;
    {
        double earliestTime = 1E+8;
        Cube::Handle<Cube::ReconObject> earliestObject;
        for (Cube::ReconObjectContainer::iterator o = objects->begin();
                o != objects->end(); ++o) 
        {
            if (Cube::Tool::AreNeighboringObjects(*muonObject,*(*o))) 
            {
                continue;
            }
            Cube::Handle<Cube::ReconTrack> track = *o;
            double objTime = -1.0;
            double objEDep = -1.0;
            if (track) 
            {
                if (track == muonObject) 
                    continue;
                objTime = track->GetMedian().T();
                objEDep = track->GetEDeposit();
            }
            Cube::Handle<Cube::ReconCluster> cluster = *o;
            if (cluster) 
            {
                objTime = cluster->GetMedian().T();
                objEDep = cluster->GetEDeposit();
            }
            if (objTime < 0) continue;
            if (objTime < earliestTime) 
            {
                earliestTime = objTime;
                int mainTraj = Cube::Tool::MainTrajectory(*event,*(*o));
                Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
                earliestTrajTime = traj->GetInitialPosition().T();
                earliestObject = *o;
            }
        }
        if (earliestObject)
        {
            int earliestTraj = Cube::Tool::MainTrajectory(*event,*earliestObject);
            if (trajectories[earliestTraj]->GetPDGCode() == 2112)
            {
                if (threshold == 20)
                numTotalSignal_20++;
                if (threshold == 50)
                numTotalSignal_50++;
                if (threshold == 200)
                numTotalSignal_200++;
            }
        }
    }

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

    if (trajectories[earliestTraj]->GetPDGCode() == 2112)
    {
        double recoDeltaT = earliestTime - muonTime;
        double trueDeltaT = earliestTrueTime - 1;

        std::cout << "recoDeltaMedianT: " << recoDeltaT << std::endl;
        std::cout << "trueDeltaMedianT: " << trueDeltaT << std::endl;

        if (threshold == 20)
        {
            numSignal_20++;
            histMedianDeltaT_20.Fill(recoDeltaT, trueDeltaT);
        }
        if (threshold == 50)
        {
            numSignal_50++;
            histMedianDeltaT_50.Fill(recoDeltaT, trueDeltaT);
        }
        if (threshold == 200)
        {
            numSignal_200++;
            histMedianDeltaT_200.Fill(recoDeltaT, trueDeltaT);
        }
    }
    else
    {
        numBackground++;
    }
}
