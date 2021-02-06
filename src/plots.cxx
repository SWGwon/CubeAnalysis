void Normalize(TH1D& inHisto)
{
    inHisto.SetStats(false);
    inHisto.Scale(1/inHisto.Integral(),"nosw2");
    inHisto.GetYaxis()->SetRangeUser(0,1);
    inHisto.GetYaxis()->SetTitle("Normalized fraction");
}

void plots()
{
    
    //histograms
    //00 signal cluster
    //01 signal track
    //10 bkg cluster
    //11 bkg track
    //
    //{
    TH1D histLeverArm[2][2];
    for (int i = 0; i < 2; i++)
    {
        histLeverArm[i][0] = TH1D("","lever arm, cluster;mm;Normalized fraction", 30, 0, 2000);
        histLeverArm[i][1] = TH1D("","lever arm, track;mm;Normalized fraction", 30, 0, 2000);
        for (int j = 0; j < 2; j++)
        {
            histLeverArm[i][j].SetLineColor(2*i+2);
            histLeverArm[i][j].SetStats(false);
        }
    }

    TH1D histTrackLength[2];
    for (int i = 0; i < 2; i++)
    {
        histTrackLength[i] = TH1D("","track length, track;mm;Normalized fraction",30,0,200);
        histTrackLength[i].SetLineColor(2*i+2);
        histTrackLength[i].SetStats(false);
    }

    TH1D histObjectE[2][2];
    for (int i = 0; i < 2; i++)
    {
        histObjectE[i][0] = TH1D("","clusterE;pe",60,0,3000);
        histObjectE[i][1] = TH1D("","trackE;pe",20,0,20000);
        for (int j = 0; j < 2; j++)
        {
            histObjectE[i][j].SetLineColor(2*i+2);
            histObjectE[i][j].SetStats(false);
        }
    }

    TH1D histAngle[2][2];
    for (int i = 0; i < 2; i++)
    {
        histAngle[i][0] = TH1D("","angle, cluster;cos(#theta)",20,-1,1);
        histAngle[i][1] = TH1D("","angle, track;cos(#theta)",20,-1,1);
        for (int j = 0; j < 2; j++)
        {
            histAngle[i][j].SetLineColor(2*i+2);
            histAngle[i][j].SetStats(false);
        }
    }

    TH1D histNeighborDistance[2][2];
    for (int i = 0; i < 2; i++)
    {
        histNeighborDistance[i][0] = TH1D("","neighborDistance, cluster;mm",10,0,400);
        histNeighborDistance[i][1] = TH1D("","number of branches, track;mm",10,0,10);
        for (int j = 0; j < 2; j++)
        {
            histNeighborDistance[i][j].SetLineColor(2*i+2);
            histNeighborDistance[i][j].SetStats(false);
        }
    }

    std::shared_ptr<TFile> inputFile = std::make_shared<TFile> ("variableOutput.root");
    TTree* inputTree = (TTree*)inputFile->Get("tree");
    float angle; inputTree->SetBranchAddress("angle", &angle);
    float leverArm; inputTree->SetBranchAddress("leverArm", &leverArm);
    float eDep; inputTree->SetBranchAddress("eDep", &eDep);
    float trackLength; inputTree->SetBranchAddress("trackLength", &trackLength);
    float neighborDistance; inputTree->SetBranchAddress("neighborDistance", &neighborDistance);
    int numberOfBranches; inputTree->SetBranchAddress("numberOfBranches", &numberOfBranches);
    float recoNeutronKE; inputTree->SetBranchAddress("recoNeutronKE", &recoNeutronKE);
    int category; inputTree->SetBranchAddress("category", &category);
    for (int i = 0; i < inputTree->GetEntries(); i++) {
        inputTree->GetEntry(i);
        if (category == 0) {//sig track
            histLeverArm[0][1].Fill(leverArm);
            histTrackLength[0].Fill(trackLength);
            histObjectE[0][1].Fill(eDep);
            histAngle[0][1].Fill(angle);
            histNeighborDistance[0][1].Fill(numberOfBranches);
        }
        if (category == 2) {//bkg track
            histLeverArm[1][1].Fill(leverArm);
            histTrackLength[1].Fill(trackLength);
            histObjectE[1][1].Fill(eDep);
            histAngle[1][1].Fill(angle);
            histNeighborDistance[1][1].Fill(numberOfBranches);
        }
        if (category == 1) {//sig cluster
            histLeverArm[0][0].Fill(leverArm);
            histObjectE[0][0].Fill(eDep);
            histAngle[0][0].Fill(angle);
            histNeighborDistance[0][0].Fill(neighborDistance);
        }
        if (category == 3) {//bkg cluster
            histLeverArm[1][0].Fill(leverArm);
            histObjectE[1][0].Fill(eDep);
            histAngle[1][0].Fill(angle);
            histNeighborDistance[1][0].Fill(neighborDistance);
        }
    }

    TLegend l(0.25,0.60,0.75,0.9);
    l.SetHeader("DUNE: Simulation", "C");
    l.AddEntry(&histLeverArm[0][0],"signal","l");
    l.AddEntry(&histLeverArm[1][0],"background","l");

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            Normalize(histLeverArm[i][j]);
            Normalize(histObjectE[i][j]);
            Normalize(histAngle[i][j]);
            Normalize(histNeighborDistance[i][j]);
        }
        Normalize(histTrackLength[i]);
    }

    TCanvas can1;
    histTrackLength[0].Draw();
    histTrackLength[1].Draw("same");
    l.Draw();
    can1.SaveAs("track_length.C");
    can1.SaveAs("track_length.pdf");

    TCanvas can2;
    can2.Divide(2,2);
    can2.cd(1);
    histLeverArm[0][0].Draw();
    histLeverArm[1][0].Draw("same");
    l.Draw();
    can2.cd(2);
    histObjectE[0][0].Draw();
    histObjectE[1][0].Draw("same");
    l.Draw();
    can2.cd(3);
    histAngle[0][0].Draw();
    histAngle[1][0].Draw("same");
    l.Draw();
    can2.cd(4);
    histNeighborDistance[0][0].Draw();
    histNeighborDistance[1][0].Draw("same");
    l.Draw();
    can2.SaveAs("cluster.C");
    can2.SaveAs("cluster.pdf");

    TCanvas can3;
    can3.Divide(2,2);
    can3.cd(1);
    histLeverArm[0][1].Draw();
    histLeverArm[1][1].Draw("same");
    l.Draw();
    can3.cd(2);
    histObjectE[0][1].Draw();
    histObjectE[1][1].Draw("same");
    l.Draw();
    can3.cd(3);
    histAngle[0][1].Draw();
    histAngle[1][1].Draw("same");
    l.Draw();
    can3.cd(4);
    histNeighborDistance[0][1].Draw();
    histNeighborDistance[1][1].Draw("same");
    l.Draw();
    can3.SaveAs("track.pdf");
    can3.SaveAs("track.C");

    //exit(0);
}
