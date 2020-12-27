#include <CubeEvent.hxx>
#include <CubeG4Hit.hxx>
#include <CubeG4Trajectory.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>

#include <ToolG4Hits.hxx>
#include <ToolTrueDirection.hxx>
#include <ToolContained.hxx>
#include <ToolRecon.hxx>
#include <ToolPrimaryId.hxx>

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>

#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <getopt.h>
#include <memory>

#include "Event.hxx"

int main(int argc, char* argv[])
{
    auto deltaTNeutron = std::make_unique<TH1D> ("","#Delta T, neutron", 100, -10 ,10);
    auto deltaTOther = std::make_unique<TH1D> ("","#Delta T, other", 100, -10 ,10);
    auto recoNuVsRealNu = std::make_unique<TH2D> ("","reco nu vs. real nu;real nu;reco nu",100,0,1,100,0,1);

    std::string inputNamesCubeReconFile;
    std::string inputNamesGenieFile;

    if (argc <= optind) throw std::runtime_error("Missing input file");

    while (optind < argc) {
        inputNamesCubeReconFile = argv[1];
        inputNamesGenieFile = argv[2];
        optind++;
    }

    int numSingleTrackEventSignal = 0;
    int nuBelow500MeV = 0;
    int notLownu = 0;

    int totalCCEvent = 0;
    int numSingleTrackEvent = 0;

    for (int fileNum = 0; fileNum < 5; fileNum++)
    {
    std::unique_ptr<TChain> inputChainCubeRecon = std::make_unique<TChain> ("CubeEvents");
    //inputChainCubeRecon->Add(inputNamesCubeReconFile.c_str());
    inputChainCubeRecon->Add(Form("~/CubeAnalysis/datafiles/latest/full3DST.antineutrino.%d.cuberecon_latest.root",fileNum));

    std::unique_ptr<TChain> inputChainGenie = std::make_unique<TChain> ("gRooTracker");
    inputChainGenie->Add(Form("~/CubeAnalysis/datafiles/latest/full3DST.antineutrino.%d.rootracker.root",fileNum+1));
    //inputChainGenie->Add(inputNamesGenieFile.c_str());

    std::cout << inputChainCubeRecon->GetEntries() << std::endl;
    std::cout << "genie: " << inputChainGenie->GetEntries() << std::endl;
    double t_StdHepP4[1000][4];
    int t_StdHepPdg[1000];
    int t_StdHepStatus[1000];
    int t_StdHepN;
    double t_EvtVtx[4];
    inputChainGenie->SetBranchAddress("StdHepP4", &t_StdHepP4);
    inputChainGenie->SetBranchAddress("StdHepPdg", &t_StdHepPdg);
    inputChainGenie->SetBranchAddress("StdHepStatus", &t_StdHepStatus);
    inputChainGenie->SetBranchAddress("StdHepN", &t_StdHepN);
    inputChainGenie->SetBranchAddress("EvtVtx", &t_EvtVtx);

    inputChainCubeRecon->SetBranchAddress("Event",&event);

    for (int i = 0; i < inputChainCubeRecon->GetEntries(); i++)
    {
        inputChainCubeRecon->GetEntry(i);
        inputChainGenie->GetEntry(i);
        std::unique_ptr<Event> testEvent = NULL;
        try
        {
            testEvent = std::make_unique<Event> ();
        }
        catch(std::runtime_error& e)
        {
            std::cout << "event: " << i;
            std::cout << ", " << e.what() << std::endl;
            continue;
        }
        catch(...)
        {
            std::cout << "event: " << i;
            std::cout << ", making Event failed" << std::endl;
            continue;
        }
        /*
        std::cout << "event: " << i << std::endl;
        for (int i : testEvent->GetPrimaryParticles().mTrackIdNeutrons)
        {
            std::cout << "testEvent->GetPrimaryParticles().mTrackIdNeutrons): " << i << std::endl;
        }
        for (int i : testEvent->GetPrimaryParticles().mTrackIdPions)
        {
            std::cout <<  "testEvent->GetPrimaryParticles().mTrackIdPions): " << i << std::endl;
        }
        for (int i : testEvent->GetPrimaryParticles().mTrackIdGammas)
        {
            std::cout <<  "testEvent->GetPrimaryParticles().mTrackIdGammas): " << i << std::endl;
        }
        for (int i : testEvent->GetPrimaryParticles().mTrackIdMuons)
        {
            std::cout <<  "testEvent->GetPrimaryParticles().mTrackIdMuons): " << i << std::endl;
        }
        for (int i : testEvent->GetPrimaryParticles().mTrackIdElectrons)
        {
            std::cout <<  "testEvent->GetPrimaryParticles().mTrackIdElectrons): " << i << std::endl;
        }
        for (int i : testEvent->GetPrimaryParticles().mTrackIdOthers)
        {
            std::cout <<  "testEvent->GetPrimaryParticles().mTrackIdOthers): " << i << std::endl;
        }
        for (int i : testEvent->GetPrimaryParticles().mTrackIdAntiMuons)
        {
            std::cout <<  "testEvent->GetPrimaryParticles().mTrackIdAntiMuons): " << i << std::endl;
        }
        for (int i : testEvent->GetPrimaryParticles().mTrackIdProtons)
        {
            std::cout <<  "testEvent->GetPrimaryParticles().mTrackIdProtons): " << i << std::endl;
        }
        */
        //CC0pi
        //if (testEvent->GetPrimaryParticles().GetNumberOfMuon() != 0) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfAntiMuon() != 1) continue;
        //if (testEvent->GetPrimaryParticles().GetNumberOfElectron() != 0) continue;
        //if (testEvent->GetPrimaryParticles().GetNumberOfGamma() != 0) continue;
        //if (testEvent->GetPrimaryParticles().GetNumberOfPion() != 0) continue;
        //if (testEvent->GetPrimaryParticles().GetNumberOfNeutron() != 1) continue;
        //if (testEvent->GetPrimaryParticles().GetNumberOfProton() != 0) continue;
        //if (testEvent->GetPrimaryParticles().GetNumberOfOther() != 0) continue;
        //
        totalCCEvent++;

        int numberOfPrimaryAntiMuon = 0;
        //select single primary anti muon object event
        for (auto o : testEvent->GetObjects())
        {
            if (o.GetPdg() == -13)
            {
                numberOfPrimaryAntiMuon++;
            }
        }
        if (numberOfPrimaryAntiMuon != 1)
        {
            continue;
        }

        //single track event
        if (testEvent->GetNumberOfVertexAssociated() != 1)
        {
            continue;
        }

        Object firstObject = testEvent->GetFirstObject();

        double earliestTime = testEvent->GetFirstObject().GetPosition().T();
        double muonTime = testEvent->GetVertex().GetPosition().T();

        Cube::Handle<Cube::ReconTrack> tempTrack;
        Cube::Handle<Cube::ReconCluster> tempCluster;
        int earliestTraj;
        if (firstObject.IsTrack)
        {
            tempTrack = firstObject.GetTrack();
            earliestTraj = Cube::Tool::MainTrajectory(*event,*tempTrack);
        }
        if (firstObject.IsCluster)
        {
            tempCluster = firstObject.GetCluster();
            earliestTraj= Cube::Tool::MainTrajectory(*event,*tempCluster);
        }
        int earliestPrim = Cube::Tool::PrimaryId(*event,earliestTraj);
        //if (testEvent->GetFirstObject().GetPdg() == 2112)
        if (earliestPrim == 1)
        {
            deltaTNeutron->Fill(earliestTime - muonTime);
        }
        else
        {
            deltaTOther->Fill(earliestTime - muonTime);
        }
        if (earliestTime - muonTime < 0)
        {
            continue;
        }
        numSingleTrackEvent++;
        if (testEvent->GetFirstObject().GetPdg() != 2112)
        {
            continue;
        }
        double realMuonE = 0;
        for (int j = 0; j < t_StdHepN; j++)
        {
            if (t_StdHepPdg[j] == -13 && t_StdHepStatus[j] == 1)
            {
                realMuonE = t_StdHepP4[j][3];
            }
        }
        numSingleTrackEventSignal++;
        double tof = earliestTime - muonTime;
        double leverArm = std::pow(std::pow(testEvent->GetFirstObject().GetPosition().X() - testEvent->GetVertex().GetPosition().X(),2)
                + std::pow(testEvent->GetFirstObject().GetPosition().Y() - testEvent->GetVertex().GetPosition().Y(),2)
                + std::pow(testEvent->GetFirstObject().GetPosition().Z() - testEvent->GetVertex().GetPosition().Z(),2),0.5);
        std::cout << "event: " << i << std::endl;
        std::cout << "genie vertex: " << t_EvtVtx[0]*1000. << ", " << t_EvtVtx[1]*1000. << ", " << t_EvtVtx[2]*1000. << std::endl;
        std::cout << "reco vetex: " << testEvent->GetVertex().GetPosition().X() << ", " << testEvent->GetVertex().GetPosition().Y() << ", " << testEvent->GetVertex().GetPosition().Z() << std::endl;
        std::cout << "anti enutrino energy: " << t_StdHepP4[0][3] << std::endl;
        std::cout << "lever arm: " << leverArm << std::endl;
        std::cout << "tof: " << tof << std::endl;
        double beta = (leverArm/tof)/300.;
        std::cout << "beta: " << beta << std::endl;
        double recoNu = 939.565*(1/std::pow(1-std::pow(beta,2),0.5)-1.);
        std::cout << "recoNu: " << recoNu/1000. << " GeV" << std::endl;
        std::cout << "realNu: " << t_StdHepP4[0][3] - realMuonE << std::endl;
        std::cout << "-------------" << std::endl;
        recoNuVsRealNu->Fill(t_StdHepP4[0][3] - realMuonE, recoNu/1000.);
        if (recoNu < 500)
        {
            nuBelow500MeV++;
        }
        else
        {
            notLownu++;
        }
    }
        //testEvent->Show();
    }
    TCanvas can1;
    deltaTNeutron->Draw();
    can1.SaveAs("deltaTNeutron.pdf");

    TCanvas can2;
    deltaTOther->Draw();
    can2.SaveAs("deltaTOther.pdf");

    TCanvas can3;
    recoNuVsRealNu->Draw("colz");
    can3.SaveAs("recoNuVsRealNu.pdf");

    std::cout << "totalCCEvent: " << totalCCEvent << std::endl;
    std::cout << "numSingleTrackEvent: " << numSingleTrackEvent << std::endl;
    std::cout << "numSingleTrackEventSignal: " << numSingleTrackEventSignal << std::endl;
    std::cout << "nuBelow500MeV: " << nuBelow500MeV << std::endl;
    std::cout << "not low nu: " << notLownu << std::endl;

    return 0;
}
