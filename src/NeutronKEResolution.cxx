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
#include <TH1.h>
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
    std::string inputNamesCubeReconFile;

    auto neutronTrueKE = std::make_unique<TH1D> ("","neutron True KE",100,0,2000);
    TH2D histDeltaT("","true #DeltaT vs reco #DeltaT;reco #DeltaT;true #DeltaT",100,-10,10,100,-10,10);
    TH1D histTrue_Reco("","true-reco/true",100,-10,10);
    TH2D histTrueRecoNeutronKE("","true vs reco neutron KE;reco;true",100,0,500,100,0,500);

    if (argc <= optind) throw std::runtime_error("Missing input file");

    while (optind < argc) {
        inputNamesCubeReconFile = argv[1];
        optind++;
    }

    std::unique_ptr<TChain> inputChainCubeRecon = std::make_unique<TChain> ("CubeEvents");
    inputChainCubeRecon->Add(inputNamesCubeReconFile.c_str());
    inputChainCubeRecon->SetBranchAddress("Event",&event);
    for (int i = 0; i < inputChainCubeRecon->GetEntries(); i++)
    {
        inputChainCubeRecon->GetEntry(i);
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
        if (testEvent->GetPrimaryParticles().GetNumberOfMuon() != 0) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfAntiMuon() != 1) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfElectron() != 0) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfGamma() != 0) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfPion() != 0) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfNeutron() != 1) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfProton() != 0) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfOther() != 0) continue;
        //

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

        Object firstObject = testEvent->GetFirstObject();
        double earliestTime = testEvent->GetFirstObject().GetPosition().T();
        double muonTime = testEvent->GetVertex().GetPosition().T();
        double recoDeltaT = earliestTime - muonTime;

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
        double trueMuonT = 0;
        for (auto t : testEvent->mTtrajectories)
        {
            if (t.second->GetPDGCode() == -13)
            {
                trueMuonT = t.second->GetInitialPosition().T();
                break;
            }
        }
        double trueDeltaT = 0;

        std::cout << "recoDeltaT: " << recoDeltaT << std::endl;
        try
        {
            trueDeltaT = testEvent->mTtrajectories[earliestTraj]->GetInitialPosition().T() - trueMuonT;
        }
        catch(...)
        {
        }
        std::cout << "trueMuonT: " << trueMuonT << std::endl;
        std::cout << "trueDeltaT: " << trueDeltaT << std::endl;
        histDeltaT.Fill(recoDeltaT, trueDeltaT);
        
        //negative tof: skip
        if (earliestTime - muonTime < 0)
        {
            continue;
        }

        //signal
        if (testEvent->GetFirstObject().GetPdg() != 2112)
        {
            continue;
        }

        double tof = earliestTime - muonTime;
        double leverArm = std::pow(std::pow(testEvent->GetFirstObject().GetPosition().X() - testEvent->GetVertex().GetPosition().X(),2)
                + std::pow(testEvent->GetFirstObject().GetPosition().Y() - testEvent->GetVertex().GetPosition().Y(),2)
                + std::pow(testEvent->GetFirstObject().GetPosition().Z() - testEvent->GetVertex().GetPosition().Z(),2),0.5);
        int histLeverArm = leverArm/10.; //cm
        int histTof = tof;
        std::cout << "event: " << i << std::endl;
        std::cout << "reco vetex: " << testEvent->GetVertex().GetPosition().X() << ", " << testEvent->GetVertex().GetPosition().Y() << ", " << testEvent->GetVertex().GetPosition().Z() << std::endl;
        std::cout << "lever arm: " << leverArm << std::endl;
        std::cout << "hist lever arm: " << histLeverArm << std::endl;
        std::cout << "tof: " << tof << std::endl;
        double beta = (leverArm/tof)/300.;
        std::cout << "beta: " << beta << std::endl;
        double recoNu = 939.565*(1./std::pow(1.-std::pow(beta,2),0.5)-1.);
        std::cout << "reco neutron KE: " << recoNu << " MeV" << std::endl;
        std::cout << "testEvent->mTtrajectories[earliestTraj]->GetInitialMomentum().E() - 939.565: " << testEvent->mTtrajectories[earliestTraj]->GetInitialMomentum().E() - 939.565 << std::endl;
        std::cout << "testEvent->mTtrajectories[earliestTraj]->GetInitialMomentum().Beta(): " << testEvent->mTtrajectories[earliestTraj]->GetInitialMomentum().Beta() << std::endl;
        neutronTrueKE->Fill(testEvent->mTtrajectories[earliestTraj]->GetInitialMomentum().E() - 939.565);
        double trueNeutronKE = testEvent->mTtrajectories[earliestTraj]->GetInitialMomentum().E() - 939.565;
        std::cout << "testEvent->mTtrajectories[earliestTraj]->GetPDGCode(): " << testEvent->mTtrajectories[earliestTraj]->GetPDGCode() << std::endl;
        histTrue_Reco.Fill((trueNeutronKE - recoNu)/trueNeutronKE);
        histTrueRecoNeutronKE.Fill(recoNu, trueNeutronKE);
    }

    neutronTrueKE->SaveAs("neutronTrueKE.C");
    histTrue_Reco.SaveAs("histTrue_Reco.C");
    TCanvas can;
    can.SetLogz();
    histTrueRecoNeutronKE.SetStats(false);
    histTrueRecoNeutronKE.Draw("colz");
    can.SaveAs("histTrueRecoNeutronKE.C");

    return 0;
}
