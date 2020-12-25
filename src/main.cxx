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
#include <memory>

#include "Event.hxx"

int main(int argc, char* argv[])
{
    auto deltaTNeutron = std::make_unique<TH1D> ("","#Delta T, neutron", 100, -10 ,10);
    auto deltaTOther = std::make_unique<TH1D> ("","#Delta T, other", 100, -10 ,10);
    std::string inputNames;
    if (argc <= optind) throw std::runtime_error("Missing input file");
    while (optind < argc) {
        inputNames = argv[optind++];
    }

    std::unique_ptr<TChain> inputChain = std::make_unique<TChain> ("CubeEvents");
    inputChain->Add(inputNames.c_str());

    std::cout << inputChain->GetEntries() << std::endl;
    inputChain->SetBranchAddress("Event",&event);

    for (int i = 0; i < inputChain->GetEntries(); i++)
    {
        inputChain->GetEntry(i);
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
        //CC0pi
        if (testEvent->GetPrimaryParticles().GetNumberOfMuon() != 0) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfAntiMuon() != 1) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfElectron() != 0) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfGamma() != 0) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfPion() != 0) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfNeutron() != 1) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfProton() != 0) continue;
        if (testEvent->GetPrimaryParticles().GetNumberOfOther() != 0) continue;

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

        double earliestTime = testEvent->GetFirstObject().GetPosition().T();
        double muonTime = testEvent->GetVertex().GetPosition().T();

        if (testEvent->GetFirstObject().GetPdg() == 2112)
        {
            deltaTNeutron->Fill(earliestTime - muonTime);
        }
        else
        {
            deltaTOther->Fill(earliestTime - muonTime);
        }
        std::cout << "event: " << i << std::endl;
        testEvent->Show();
    }
    TCanvas can1;
    deltaTNeutron->Draw();
    can1.SaveAs("deltaTNeutron.pdf");

    TCanvas can2;
    deltaTOther->Draw();
    can2.SaveAs("deltaTOther.pdf");

    return 0;
}
