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
        Analysis(event, 20);
        Analysis(event, 50);
        Analysis(event, 200);
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
    Cube::Handle<Cube::ReconObjectContainer> objects = event->GetObjectContainer();
    if (!objects) 
        return;

    // Make sure there's a single muon track.
    int muonLike = 0;
    Cube::Handle<Cube::ReconTrack> muonObject;
    for (Cube::ReconObjectContainer::iterator o = objects->begin();
            o != objects->end(); ++o) 
    {
        Cube::Handle<Cube::ReconTrack> track = *o;
        if (!track) 
            continue;
        int mainTraj = Cube::Tool::MainTrajectory(*event, *track);
        if (mainTraj > trajectories.size()) 
            continue;
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[mainTraj];
        if (traj->GetPDGCode() != -13 || traj->GetParentId() != -1) 
            continue;
        muonLike++;
        muonObject = track;
    }
    if (muonLike != 1 || !muonObject) return;

    int numberOfMuonAssociated = 0;
    // Collect the earliest objects and the muon.
    double muonTime = muonObject->GetPosition().T();
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
        double objTime = -1.0;
        double objEDep = -1.0;
        Cube::Handle<Cube::ReconTrack> track = *o;
        Cube::Handle<Cube::ReconCluster> cluster = *o;
        if (track) 
        {
            objTime = track->GetMedian().T();
            objEDep = track->GetEDeposit();
        }
        else if (cluster) 
        {
            objTime = cluster->GetMedian().T();
            objEDep = cluster->GetEDeposit();
        }
        else 
            continue;
        if (objTime < earliestTime && objEDep > threshold) 
        {
            earliestTime = objTime;
            earliestObject = *o;
        }
    }
    if (!earliestObject) return;
    if (numberOfMuonAssociated != 1) return;
    int earliestTraj = Cube::Tool::MainTrajectory(*event,*earliestObject);
    if (earliestTraj > trajectories.size())
        return;

    //no threshold
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
            Cube::Handle<Cube::ReconCluster> cluster = *o;
            double objTime = -1.0;
            double objEDep = -1.0;
            if (track) 
            {
                objTime = track->GetMedian().T();
                objEDep = track->GetEDeposit();
            }
            else if (cluster) 
            {
                objTime = cluster->GetMedian().T();
                objEDep = cluster->GetEDeposit();
            }
            else
                continue;
            if (objTime < earliestTime) 
            {
                earliestTime = objTime;
                earliestObject = *o;
            }
        }
        if (earliestObject)
        {
            int earliestTraj = Cube::Tool::MainTrajectory(*event,*earliestObject);
            if (earliestTraj > trajectories.size())
                return;
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
