void cut()
{
    TH1D* afterCutNeutronKECluster = new TH1D("","efficiency of selection, cluster;reco neutron KE, MeV;efficiency",20,0,100);
    TH1D* beforeCutNeutronKECluster = new TH1D("","efficiency, cluster;reco neutron KE",20,0,100);
    TH1D* afterCutNeutronKETrack = new TH1D("","efficiency of selection, track;reco neutron KE, MeV;efficiency",20,0,100);
    TH1D* beforeCutNeutronKETrack = new TH1D("","efficiency, track;reco neutron KE",20,0,100);

    TH1D* afterCutLeptonAngleCluster = new TH1D("","efficiency of selection, cluster;lepton angle;efficiency",20,0,1.5);
    TH1D* beforeCutLeptonAngleCluster = new TH1D("","efficiency, cluster;lepton angle",20,0,1.5);
    TH1D* afterCutLeptonAngleTrack = new TH1D("","efficiency of selection, track;lepton angle;efficiency",20,0,1.5);
    TH1D* beforeCutLeptonAngleTrack = new TH1D("","efficiency, track;lepton angle",20,0,1.5);

    TH1D* afterCutLeptonMomentumCluster = new TH1D("","efficiency of selection, cluster;lepton momentum;efficiency",20,0,10);
    TH1D* beforeCutLeptonMomentumCluster = new TH1D("","efficiency, cluster;lepton momentum",20,0,10);
    TH1D* afterCutLeptonMomentumTrack = new TH1D("","efficiency of selection, track;lepton momentum;efficiency",20,0,10);
    TH1D* beforeCutLeptonMomentumTrack = new TH1D("","efficiency, track;lepton momentum",20,0,10);

    TH1D* afterCutQ2Cluster = new TH1D("","efficiency of selection, cluster;Q2;efficiency",20,0,5);
    TH1D* beforeCutQ2Cluster = new TH1D("","efficiency, cluster;Q2",20,0,5);
    TH1D* afterCutQ2Track = new TH1D("","efficiency of selection, track;Q2;efficiency",20,0,5);
    TH1D* beforeCutQ2Track = new TH1D("","efficiency, track;Q2",20,0,5);

    TH1D* afterCutQ32Cluster = new TH1D("","efficiency of selection, cluster;Q32;efficiency",20,0,5);
    TH1D* beforeCutQ32Cluster = new TH1D("","efficiency, cluster;Q32",20,0,5);
    TH1D* afterCutQ32Track = new TH1D("","efficiency of selection, track;Q32;efficiency",20,0,5);
    TH1D* beforeCutQ32Track = new TH1D("","efficiency, track;Q32",20,0,5);

    //TFile* inputFile = new TFile("variableOutput_noinf.root");
    TFile* inputFile = new TFile("variableOutput.root");
    TTree* inputTree = (TTree*)inputFile->Get("tree");
    int numberOfBranches; inputTree->SetBranchAddress("numberOfBranches", &numberOfBranches);
    //float angle; inputTree->SetBranchAddress("angle", &angle);
    float tof; inputTree->SetBranchAddress("tof", &tof);
    float eDep; inputTree->SetBranchAddress("eDep", &eDep);
    float trackLength; inputTree->SetBranchAddress("trackLength", &trackLength);
    float neighborDistance; inputTree->SetBranchAddress("neighborDistance", &neighborDistance);
    float recoNeutronKE; inputTree->SetBranchAddress("recoNeutronKE", &recoNeutronKE);
    float trueNeutronKE; inputTree->SetBranchAddress("trueNeutronKE", &trueNeutronKE);
    int category; inputTree->SetBranchAddress("category", &category);
    float leptonAngle; inputTree->SetBranchAddress("leptonAngle", &leptonAngle);
    float leptonMomentum; inputTree->SetBranchAddress("leptonMomentum", &leptonMomentum);
    float Q2; inputTree->SetBranchAddress("Q2", &Q2);
    float Q32; inputTree->SetBranchAddress("Q32", &Q32);

    auto outputFile = std::make_unique<TFile> ("variableOutputAfterCut.root","RECREATE");
    TTree* cloneTree = inputTree->CloneTree(0);

    double trackSigBeforeCuts = 0;
    double trackSigAfterCuts = 0;
    double clusterSigBeforeCuts = 0;
    double clusterSigAfterCuts = 0;
    double trackBkgBeforeCuts = 0;
    double trackBkgAfterCuts = 0;
    double clusterBkgBeforeCuts = 0;
    double clusterBkgAfterCuts = 0;

    float eDepMultiplier = 1;
    float eDepMultiplierCluster = 0.5;
    int eDepTime = 2000;
    float ArrclusterSigBeforeCuts[2000][10] = {};
    float ArrclusterSigAfterCuts[2000][10] = {};
    float ArrclusterBkgBeforeCuts[2000][10] = {};
    float ArrclusterBkgAfterCuts[2000][10] = {};
    float ArrtrackSigBeforeCuts[2000][10][10] = {};
    float ArrtrackSigAfterCuts[2000][10][10] = {};
    float ArrtrackBkgBeforeCuts[2000][10][10] = {};
    float ArrtrackBkgAfterCuts[2000][10][10] = {};

    for (int i = 0; i < inputTree->GetEntries(); i++) {
        inputTree->GetEntry(i);
        for (int j = 0; j < eDepTime; j++) {
            for (int nBran = 0; nBran < 10; nBran++) {
                for (int trackL = 0; trackL < 10; trackL++) {
                    if (category == 0) {//sig track
                        ArrtrackSigBeforeCuts[j][nBran][trackL]++;
                        if (eDep > eDepMultiplier*j && numberOfBranches < nBran && trackLength < 10*trackL && tof > 0) {
                            ArrtrackSigAfterCuts[j][nBran][trackL]++;
                        }
                    }
                    if (category == 2) {//bkg track
                        ArrtrackBkgBeforeCuts[j][nBran][trackL]++;
                        if (eDep > eDepMultiplier*j && numberOfBranches < nBran && trackLength < 10*trackL && tof > 0) {
                            ArrtrackBkgAfterCuts[j][nBran][trackL]++;
                        }
                    }
                }
            }

            for (int ndist = 0; ndist < 10; ++ndist) {
                if (category == 1) {//sig cluster
                    ArrclusterSigBeforeCuts[j][ndist]++;
                    if (eDep > eDepMultiplierCluster*j && neighborDistance > ndist*10 && tof > 0) {
                        ArrclusterSigAfterCuts[j][ndist]++;
                    }
                }
                if (category == 3) {//bkg cluster
                    ArrclusterBkgBeforeCuts[j][ndist]++;
                    if (eDep > eDepMultiplierCluster*j && neighborDistance > ndist*10 && tof > 0) {
                        ArrclusterBkgAfterCuts[j][ndist]++;
                    }
                }
            }
        }
    }
    int bestEDep = 0;
    int bestndist = 0;
    int bestEDepTrack = 0;
    int bestnBran = 0;
    int besttrackL = 0;
    float bestPurity = 0;
    float bestEfficiency = 0;
    float bestPurityTrack = 0;
    float bestEfficiencyTrack = 0;
    float score = 0;
    float scoreTrack = 0;
    for (int j = 0; j < eDepTime; j++) {
        for (int nBran = 0; nBran < 10; nBran++) {
            for (int trackL = 0; trackL < 10; trackL++) {
                double trackPurity = ArrtrackSigAfterCuts[j][nBran][trackL]/(ArrtrackSigAfterCuts[j][nBran][trackL] + ArrtrackBkgAfterCuts[j][nBran][trackL]);
                double trackEfficiency = ArrtrackSigAfterCuts[j][nBran][trackL]/ArrtrackSigBeforeCuts[j][nBran][trackL];
                if (trackPurity > 0.95 && trackPurity*trackEfficiency > scoreTrack) {
                //if (trackPurity > scoreTrack) {
                    //scoreTrack = trackPurity;
                    scoreTrack = trackPurity * trackEfficiency;
                    bestEDepTrack = j;
                    bestnBran = nBran;
                    besttrackL = trackL;
                    bestPurityTrack = trackPurity;
                    bestEfficiencyTrack = trackEfficiency;
                }
            }
        }
        for (int ndist = 0; ndist < 10; ndist++) {
            double clusterPurity = ArrclusterSigAfterCuts[j][ndist]/(ArrclusterSigAfterCuts[j][ndist] + ArrclusterBkgAfterCuts[j][ndist]);
            double clusterEfficiency = ArrclusterSigAfterCuts[j][ndist]/ArrclusterSigBeforeCuts[j][ndist];
            if (clusterPurity > 0.95 && clusterPurity*clusterEfficiency > score) {
            //if (clusterPurity > score) {
                //score = clusterPurity;
                score = clusterPurity * clusterEfficiency;
                bestEDep = j;
                bestndist = ndist;
                bestPurity = clusterPurity;
                bestEfficiency = clusterEfficiency;
            }
        }
    }
    std::cout << "------------------" << std::endl;
    std::cout << "cluster best cut" << std::endl;
    std::cout << "cuts: eDep > " << eDepMultiplierCluster*bestEDep << ", neighborDistance > " << 10* bestndist << std::endl;
    std::cout << "clusterPurity: " << bestPurity << std::endl;
    std::cout << "clusterEfficiency: " << bestEfficiency << std::endl;
    std::cout << "eff*purity: " << score << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "track best cut" << std::endl;
    std::cout << "cuts: eDep > " << eDepMultiplier*bestEDepTrack << ", numberOfBranches < " << bestnBran << ", trackLength < " << 10*besttrackL << std::endl;
    std::cout << "trackPurity: " << bestPurityTrack << std::endl;
    std::cout << "trackEfficiency: " << bestEfficiencyTrack << std::endl;
    std::cout << "eff*purity: " << scoreTrack << std::endl;
    std::cout << "------------------" << std::endl;
    
    for (int i = 0; i < inputTree->GetEntries(); i++) {
        inputTree->GetEntry(i);
        if (category == 0) {//sig track
            trackSigBeforeCuts++;
            beforeCutNeutronKETrack->Fill(trueNeutronKE);
            beforeCutLeptonAngleTrack->Fill(leptonAngle);
            beforeCutLeptonMomentumTrack->Fill(leptonMomentum);
            beforeCutQ2Track->Fill(Q2);
            beforeCutQ32Track->Fill(Q32);
            if (eDep > eDepMultiplier*bestEDepTrack && numberOfBranches < bestnBran && trackLength < 10*besttrackL && tof > 0) {
                afterCutNeutronKETrack->Fill(trueNeutronKE);
                afterCutLeptonAngleTrack->Fill(leptonAngle);
                afterCutLeptonMomentumTrack->Fill(leptonMomentum);
                afterCutQ2Track->Fill(Q2);
                afterCutQ32Track->Fill(Q32);
                trackSigAfterCuts++;
                cloneTree->Fill();
            }
        }
        if (category == 2) {//bkg track
            trackBkgBeforeCuts++;
            if (eDep > eDepMultiplier*bestEDepTrack && numberOfBranches < bestnBran && trackLength < 10*besttrackL && tof > 0) {
                trackBkgAfterCuts++;
                cloneTree->Fill();
            }
        }
        if (category == 1) {//sig cluster
            clusterSigBeforeCuts++;
            beforeCutNeutronKECluster->Fill(trueNeutronKE);
            beforeCutLeptonAngleCluster->Fill(leptonAngle);
            beforeCutLeptonMomentumCluster->Fill(leptonMomentum);
            beforeCutQ2Cluster->Fill(Q2);
            beforeCutQ32Cluster->Fill(Q32);
            if (eDep > eDepMultiplierCluster*bestEDep && neighborDistance > 10*bestndist && tof > 0) {
                clusterSigAfterCuts++;
                afterCutNeutronKECluster->Fill(trueNeutronKE);
                afterCutLeptonAngleCluster->Fill(leptonAngle);
                afterCutLeptonMomentumCluster->Fill(leptonMomentum);
                afterCutQ2Cluster->Fill(Q2);
                afterCutQ32Cluster->Fill(Q32);
                cloneTree->Fill();
            }
        }
        if (category == 3) {//bkg cluster
            clusterBkgBeforeCuts++;
            if (eDep > eDepMultiplierCluster*bestEDep && neighborDistance > 10*bestndist && tof > 0) {
                clusterBkgAfterCuts++;
                cloneTree->Fill();
            }
        }
    }
    TCanvas can;
    afterCutNeutronKECluster->SetStats(false);
    afterCutNeutronKECluster->GetYaxis()->SetRangeUser(0,1);
    afterCutNeutronKECluster->Divide(beforeCutNeutronKECluster);
    afterCutNeutronKECluster->Draw();
    can.SaveAs("efficiencyNeutronKECluster.pdf");
    TCanvas can1;
    afterCutNeutronKETrack->SetStats(false);
    afterCutNeutronKETrack->GetYaxis()->SetRangeUser(0,1);
    afterCutNeutronKETrack->Divide(beforeCutNeutronKETrack);
    afterCutNeutronKETrack->Draw();
    can1.SaveAs("efficiencyNeutronKETrack.pdf");

    TCanvas can2;
    afterCutLeptonAngleCluster->SetStats(false);
    afterCutLeptonAngleCluster->GetYaxis()->SetRangeUser(0,1);
    afterCutLeptonAngleCluster->Divide(beforeCutLeptonAngleCluster);
    afterCutLeptonAngleCluster->Draw();
    can2.SaveAs("efficiencyLeptonAngleCluster.pdf");
    TCanvas can3;
    afterCutLeptonAngleTrack->SetStats(false);
    afterCutLeptonAngleTrack->GetYaxis()->SetRangeUser(0,1);
    afterCutLeptonAngleTrack->Divide(beforeCutLeptonAngleTrack);
    afterCutLeptonAngleTrack->Draw();
    can3.SaveAs("efficiencyLeptonAngleTrack.pdf");
    
    TCanvas can4;
    afterCutLeptonMomentumCluster->SetStats(false);
    afterCutLeptonMomentumCluster->GetYaxis()->SetRangeUser(0,1);
    afterCutLeptonMomentumCluster->Divide(beforeCutLeptonMomentumCluster);
    afterCutLeptonMomentumCluster->Draw();
    can4.SaveAs("efficiencyLeptonMomentumCluster.pdf");
    TCanvas can5;
    afterCutLeptonMomentumTrack->SetStats(false);
    afterCutLeptonMomentumTrack->GetYaxis()->SetRangeUser(0,1);
    afterCutLeptonMomentumTrack->Divide(beforeCutLeptonMomentumTrack);
    afterCutLeptonMomentumTrack->Draw();
    can5.SaveAs("efficiencyLeptonMomentumTrack.pdf");

    TCanvas can6;
    afterCutQ2Cluster->SetStats(false);
    afterCutQ2Cluster->GetYaxis()->SetRangeUser(0,1);
    afterCutQ2Cluster->Divide(beforeCutQ2Cluster);
    afterCutQ2Cluster->Draw();
    can6.SaveAs("efficiencyQ2Cluster.pdf");
    TCanvas can7;
    afterCutQ2Track->SetStats(false);
    afterCutQ2Track->GetYaxis()->SetRangeUser(0,1);
    afterCutQ2Track->Divide(beforeCutQ2Track);
    afterCutQ2Track->Draw();
    can7.SaveAs("efficiencyQ2Track.pdf");

    TCanvas can8;
    afterCutQ32Cluster->SetStats(false);
    afterCutQ32Cluster->GetYaxis()->SetRangeUser(0,1);
    afterCutQ32Cluster->Divide(beforeCutQ32Cluster);
    afterCutQ32Cluster->Draw();
    can8.SaveAs("efficiencyQ32Cluster.pdf");
    TCanvas can9;
    afterCutQ32Track->SetStats(false);
    afterCutQ32Track->GetYaxis()->SetRangeUser(0,1);
    afterCutQ32Track->Divide(beforeCutQ32Track);
    afterCutQ32Track->Draw();
    can9.SaveAs("efficiencyQ32Track.pdf");

    double trackPurity = trackSigAfterCuts/(trackSigAfterCuts + trackBkgAfterCuts);
    double trackEfficiency = trackSigAfterCuts/trackSigBeforeCuts;
    double clusterPurity = clusterSigAfterCuts/(clusterSigAfterCuts + clusterBkgAfterCuts);
    double clusterEfficiency = clusterSigAfterCuts/clusterSigBeforeCuts;

    std::cout << "trackSigBeforeCuts: " << trackSigBeforeCuts << std::endl;
    std::cout << "trackBkgBeforeCuts: " << trackBkgBeforeCuts << std::endl;
    std::cout << "clusterSigBeforeCuts: " << clusterSigBeforeCuts << std::endl;
    std::cout << "clusterBkgBeforeCuts: " << clusterBkgBeforeCuts << std::endl;
    std::cout << "trackSigAfterCuts: " << trackSigAfterCuts << std::endl;
    std::cout << "trackBkgAfterCuts: " << trackBkgAfterCuts << std::endl;
    std::cout << "clusterSigAfterCuts: " << clusterSigAfterCuts << std::endl;
    std::cout << "clusterBkgAfterCuts: " << clusterBkgAfterCuts << std::endl;
    std::cout << "------------------" << std::endl;
    std::cout << "trackPurity: " << trackPurity << std::endl;
    std::cout << "trackEfficiency: " << trackEfficiency << std::endl;
    std::cout << "clusterPurity: " << clusterPurity << std::endl;
    std::cout << "clusterEfficiency: " << clusterEfficiency << std::endl;

    outputFile->Write();
    outputFile->Close();

    exit(0);
}
