void geineNu()
{
    TH1D* hist = new TH1D("","nuE - muonE",100,0,200);
    double EvtVtx[4]; 
    double StdHepP4[1000][4]; 
    int StdHepStatus[1000]; 
    int StdHepPdg[1000]; 
    int StdHepN; 
    for (int j = 100; j < 200; j++) {
        TFile inputGenieFile(Form("/Users/gwon/CubeAnalysis/datafiles/latest/full3DST.antineutrino.%d.rootracker.root",j+1));
        TTree* inputGenieTree = (TTree*)inputGenieFile.Get("gRooTracker");
        inputGenieTree->SetBranchAddress("EvtVtx", &EvtVtx);
        inputGenieTree->SetBranchAddress("StdHepP4", &StdHepP4);
        inputGenieTree->SetBranchAddress("StdHepStatus", &StdHepStatus);
        inputGenieTree->SetBranchAddress("StdHepPdg", &StdHepPdg);
        inputGenieTree->SetBranchAddress("StdHepN", &StdHepN);
        for (int i = 0; i < inputGenieTree->GetEntries(); i++) {
            inputGenieTree->GetEntry(i);
            if (StdHepPdg[0] != -14)
                continue;
            for (int k = 0; k < StdHepN; k++)
            {
                if (StdHepPdg[k] == -13)
                {
                    hist->Fill((StdHepP4[0][3] - StdHepP4[k][3])*1000.);
                    break;
                }
            }
        }
    }
    hist->Draw();
}
