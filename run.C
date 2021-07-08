 double deltaR(double eta1, double eta2, double phi1, double phi2)
{
  double deta = eta1 - eta2;
  double dphi = phi1 - phi2;
  if (dphi<-M_PI) dphi += 2*M_PI;
  else if (dphi>=M_PI) dphi -= 2*M_PI;
  return sqrt(deta*deta + dphi*dphi);
}

void bteff(const char *file_in,const char *file_out){
  const char* chtree = "bTag_AntiKt4EMPFlowJets";
  TChain* tt = new TChain(chtree);
  tt->Add(file_in);
  TTree* tree = (TTree*)gROOT->FindObject(chtree);

  // setup useful branches
  vector<double>* jet_mv2c10 = 0; 
  vector<float>* jet_ip3d_llr = 0;
  vector<int>* jet_isPU = 0;
  vector<float>* jet_pt = 0;
  vector<float>* jet_eta = 0;
  vector<float>* jet_phi = 0;
  vector<int>* jet_LabDr_HadF = 0;
  vector<int>* jet_truthMatch =0;
  vector<int>*  jet_aliveAfterOR=0;
  vector<int>*  jet_aliveAfterORmu=0;
  vector<float>* jet_JVT=0;
  float truth_PVx,truth_PVy,truth_PVz, PVz;
  vector<double>* jet_dl1_pb = 0;
  vector<double>* jet_dl1_pu = 0;
  vector<double>* jet_dl1_pc = 0;
  
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  vector<vector<float>>*  jet_bH_Lxy = 0;
  
  tree->SetBranchAddress("jet_mv2c10",&jet_mv2c10);
  tree->SetBranchAddress("jet_ip3d_llr",&jet_ip3d_llr);
  tree->SetBranchAddress("jet_isPU",&jet_isPU);
  tree->SetBranchAddress("jet_bH_Lxy",&jet_bH_Lxy);
  tree->SetBranchAddress("jet_pt",&jet_pt);
  tree->SetBranchAddress("jet_eta",&jet_eta);
  tree->SetBranchAddress("jet_phi",&jet_phi);
  tree->SetBranchAddress("jet_LabDr_HadF",&jet_LabDr_HadF);
  tree->SetBranchAddress("jet_truthMatch", &jet_truthMatch);
  tree->SetBranchAddress("jet_aliveAfterOR", &jet_aliveAfterOR);
  tree->SetBranchAddress("jet_aliveAfterORmu", &jet_aliveAfterORmu);
  tree->SetBranchAddress("jet_JVT", &jet_JVT);
  tree->SetBranchAddress("truth_PVx",&truth_PVx);
  tree->SetBranchAddress("truth_PVy",&truth_PVy);
  tree->SetBranchAddress("jet_dl1_pb", &jet_dl1_pb);
  tree->SetBranchAddress("jet_dl1_pc", &jet_dl1_pc);
  tree->SetBranchAddress("jet_dl1_pu", &jet_dl1_pu);
  
  int nmmc;
  tree->SetBranchAddress("nmmc",&nmmc);
  vector<int>* mymc_pdgId = 0;
  tree->SetBranchAddress("mymc_pdgId",&mymc_pdgId);
  vector<float>* mymc_decayVtx_x = 0;
  tree->SetBranchAddress("mymc_decayVtx_x",&mymc_decayVtx_x);
  vector<float>* mymc_decayVtx_y = 0;
  tree->SetBranchAddress("mymc_decayVtx_y",&mymc_decayVtx_y);

  vector<int>* mymc_ix1 = 0;
  tree->SetBranchAddress("mymc_ix1",&mymc_ix1);
  vector<int>* mymc_ix2 = 0;
  tree->SetBranchAddress("mymc_ix2",&mymc_ix2);
  vector<float>* mymc1_eta = 0;
  tree->SetBranchAddress("mymc1_eta",&mymc1_eta);
  vector<float>* mymc1_phi = 0;
  tree->SetBranchAddress("mymc1_phi",&mymc1_phi);
  
  TFile* fout = new TFile(file_out,"RECREATE");//file for writing all histograms
  fout->cd();

  // book histograms

  const int ntag = 3;
  const char* chtag[ntag] = {"mv2c10", "ip3d", "dl1"};
  double xmin[ntag] = {-1.,-20.,-25.};
  double xmax[ntag] = { 1., 50., 25.};

  const int nf = 3; // l,c,b
  TH1F* h_w[ntag*nf];
  for (int jf = 0; jf<nf; ++jf) {
    for (int itag = 0; itag<ntag; ++itag) {
      TString sweight = "w_"; sweight+=chtag[itag];
      sweight+="_"; sweight+=jf;
      h_w[itag*nf+jf] = new TH1F(sweight,"weight",2000,xmin[itag],xmax[itag]);
    }
  }
  
  //double xdiv[] = {0.01,0.02,0.05,0.1,0.2,0.5,1.,2.,5.,10.,20.,50.,100.,200.,500.,1000.};
  //int ndiv = sizeof(xdiv)/sizeof(double)-1;

  TH1* h_dl1_all = new TH1F("dl1_all","",nf,0.,nf);
  TH1* h_dl1_tag = (TH1*)h_dl1_all->Clone("dl1_tag");
  TH1* h_mv2c10_all = (TH1*)h_dl1_all->Clone("mv2c10_all");
  TH1* h_mv2c10_tag = (TH1*)h_dl1_all->Clone("mv2c10_tag");

  TH1* h_mv2c10_all_r2v[nf];
  TH1* h_mv2c10_tag_r2v[nf];
  TH1* h_dl1_all_r2v[nf];
  TH1* h_dl1_tag_r2v[nf];
  
  TH1* h_mv2c10_all_lxy[nf];
  TH1* h_mv2c10_tag_lxy[nf];
  TH1* h_dl1_all_lxy[nf];
  TH1* h_dl1_tag_lxy[nf];
  
  for (int jf = 0; jf<nf; ++jf) {
    h_mv2c10_all_r2v[jf] = new TH1F(TString("mv2c10_all_r2v_"+std::to_string(jf)),"",50,0.,50.);
    h_mv2c10_tag_r2v[jf] = new TH1F(TString("mv2c10_tag_r2v_"+std::to_string(jf)),"",50,0.,50.);
    h_dl1_all_r2v[jf] = new TH1F(TString("dl1_all_r2v_"+std::to_string(jf)),"",50,0.,50.);
    h_dl1_tag_r2v[jf] = new TH1F(TString("dl1_tag_r2v_"+std::to_string(jf)),"",50,0.,50.);
  }
  for (int jf = 0; jf<nf; ++jf) {
    h_mv2c10_all_lxy[jf] = new TH1F(TString("mv2c10_all_lxy_"+std::to_string(jf)),"",50,0.,50.);
    h_mv2c10_tag_lxy[jf] = new TH1F(TString("mv2c10_tag_lxy_"+std::to_string(jf)),"",50,0.,50.);
    h_dl1_all_lxy[jf] = new TH1F(TString("dl1_all_lxy_"+std::to_string(jf)),"",50,0.,50.);
    h_dl1_tag_lxy[jf] = new TH1F(TString("dl1_tag_lxy_"+std::to_string(jf)),"",50,0.,50.);
  }
  
  int dumpsz;
  vector<float> lxyval;
  vector<float> r2vval;
  TH2F *lxy_r2v = new TH2F("Lxy_r2v","Jet Wise lxy & r2v",100,0.,100.,100,0.,100.);
  
   vector<double> wtag;
   Long64_t nentries = tree->GetEntriesFast();
   for (Long64_t jentry=0; jentry<nentries; jentry++) {
        if (tree->LoadTree(jentry)<0) break;
        tree->GetEntry(jentry);

    // check if there are jets
        int njet = 0; 
        if (jet_pt) njet = jet_pt->size();
        dumpsz=njet;
        if (njet==0) continue;
        
    
        wtag.clear();
    for (int ijet = 0; ijet<njet; ++ijet) {
        wtag.push_back((*jet_mv2c10)[ijet]);
        wtag.push_back((*jet_ip3d_llr)[ijet]);
        double frac_c = 0.080;
        double a = (*jet_dl1_pb)[ijet];
        double b = frac_c*(*jet_dl1_pc)[ijet] + (1-frac_c)*(*jet_dl1_pu)[ijet];
        double w_dl1 = -999; if (a>0 && b>0) w_dl1 = log(a/b);
        wtag.push_back(w_dl1);
    }

    // find LLP
    for (int i = 0; i<nmmc; ++i) {
      if ((*mymc_pdgId)[i]==36) {
          //cout<<"PDG ID 36 Confirmed:"<<endl;
          // calculate the LLP decay vertex position
	      double r2v = sqrt(pow((*mymc_decayVtx_x)[i]-truth_PVx,2) + pow((*mymc_decayVtx_y)[i]-truth_PVy,2));
           //cout<<"r2v is:"<<r2v<<endl;
	       r2vval.push_back(r2v);
            // loop over LLP decay products
        for (int ix = (*mymc_ix1)[i]; ix<(*mymc_ix2)[i]; ++ix) {
            // get the decay product direction
            double etaq = (*mymc1_eta)[ix];
            double phiq = (*mymc1_phi)[ix];
                    // loop over jets
                for (int j = 0; j<njet; ++j) {
                    // only consider b-jets
                    int jf = (*jet_LabDr_HadF)[j];
                    if (jf==4) jf = 1; else if (jf==5) jf = 2; else if (jf!=0) continue;
                 
                    //if ((*jet_isPU)[j]==1){cout<<"Jet PU is 1:"<<endl; continue;}
                    //if ((*jet_truthMatch)[j]==0){cout<<"Truth Match==0:"<<endl; continue;}
                 if ((*jet_pt)[j]<120e3 && (*jet_JVT)[j]<0.59){/**cout<<"JVT<0.59"<<endl;*/  continue;}
                    if ((*jet_aliveAfterOR)[j]==0){/**cout<<"Alive==0"<<endl;*/  continue;}
                    if ((*jet_aliveAfterORmu)[j]==0){/**cout<<"Alivemu==0"<<endl;*/  continue;}
                    if ((*jet_pt)[j]<20e3) {/**cout<<"JetPT<20K"<<endl;*/  continue;}
                    if (fabs((*jet_eta)[j])>2.5){/**cout<<"JetEta>2.5"<<endl;*/  continue;}
                                    
                    for (int itag = 0; itag<ntag; ++itag){h_w[itag*nf+jf]->Fill(wtag[j*ntag+itag]);}
            
                    h_mv2c10_all->Fill(jf);
                    //h_mv2c10_all->Write();
                    if (wtag[j*ntag]>0.83) h_mv2c10_tag->Fill(jf);
                    //h_mv2c10_tag->Write();
            
                    h_dl1_all->Fill(jf);
                    //h_dl1_all->Write();
                    if (wtag[j*ntag+ntag-1]>2.01) h_dl1_tag->Fill(jf);
                    //get Lxy
                    double lxy = 0;
                    double lxy1 = 0; 
          lxy1 = (*jet_bH_Lxy)[j][0];
          lxyval.push_back(lxy1);
          lxy_r2v->Fill(r2v,lxy1,1.);
          //lxy_r2v->Fill(r2v,lxy1);
                   
		    //!(*jet_bH_Lxy)[j].empty() &&
		    if ((*jet_bH_Lxy)[j][0]>0){
			    	lxy = (*jet_bH_Lxy)[j][0];
                    
                    		h_mv2c10_all_lxy[jf]->Fill(lxy);
                    		if (wtag[j*ntag]>0.831) {h_mv2c10_tag_lxy[jf]->Fill(lxy);}
		
                    		h_dl1_all_lxy[jf]->Fill(lxy);
                    		if (wtag[j*ntag+ntag-1]>2.01){h_dl1_tag_lxy[jf]->Fill(lxy);}
		    }
                    //h_dl1_tag->Write();
                    // get the jet direction
                    double etaj = (*jet_eta)[j];
                    double phij = (*jet_phi)[j];
                    
                    // calculate dR between the jet and the LLP decay product
                    double dR = deltaR(etaq, etaj, phiq, phij);
                    //cout<<"dR is:"<<dR<<endl;
                    if (dR<0.3) {//cout<<"Hiyaa, Entered dR cut loop:"<<endl;
                        // all jets
                        // & tagged jets
                        h_mv2c10_all_r2v[jf]->Fill(r2v);
                        //h_mv2c10_all_r2v[jf]->Write();
                        if (wtag[j*ntag]>0.831) {h_mv2c10_tag_r2v[jf]->Fill(r2v);
                            //h_mv2c10_tag_r2v[jf]->Write();
                            }
		
                        h_dl1_all_r2v[jf]->Fill(r2v);
                        //h_dl1_all_r2v[jf]->Write();
                        if (wtag[j*ntag+ntag-1]>2.01){h_dl1_tag_r2v[jf]->Fill(r2v);
                            //h_dl1_tag_r2v[jf]->Write();
                            
                        }
            
                    }    
                            
         }
       }
     }
   }
   }
   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(0);
   std::string str0(file_out);
   std::string str1;
   //str1=str0 +"LEGO.png";
   //str1=str0 +"COLZ.png";
   str1=str0 +"LEGO2Z.png";
   const char *h1 = str1.c_str();
  TCanvas *c12 = new TCanvas("c12","c12;r2v;lxy",1200,800) ;
  //lxy_r2v->Draw("COLZ");
  lxy_r2v->Draw("LEGO2Z");
  //lxy_r2v->SetLineColor(3);
  //lxy_r2v->SetLineWidth(3);
  //lxy_r2v->SetMarkerStyle(20);
  lxy_r2v->SetTitle("Jet count(r2v vs lxy);r2v(mm);lxy(mm)");
  c12->SaveAs(h1);
  fout->Write();
  delete c12;
  delete fout;
  delete tree;
}

void mv2c10_plot_r2v(const char *dat_file1, const char *dat_file2 ,const char *dat_file3){
    TFile *f_in1 = new TFile(dat_file1,"READ");
    TFile *f_in2 = new TFile(dat_file2,"READ");
    TFile *f_in3 = new TFile(dat_file3,"READ");
    TFile* fout = new TFile("mv2c10_r2v.root","RECREATE");//file for writing all histograms
    fout->cd();

    //TCanvas *c1 = new TCanvas("c1","c1",1200,800) ;
  
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitleX(0.1f);
    gStyle->SetTitleW(0.8f);
    TLegend* l = new TLegend(0.60,.90,0.75,0.95);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    //l->SetMargin(0.4);
    
    std::string str1;
    std::string str2;
    str1="mv2c10_tag_r2v_2";
    str2="mv2c10_all_r2v_2";
    
    const char *h1 = str1.c_str();
    const char *h2 = str2.c_str();
    
    TH1F* gram1[2];
    TH1F* gram2[2];
    TH1F* gram3[2];     
    
    gram1[0] = (TH1F*)f_in1->Get(h1);
    gram1[1] = (TH1F*)f_in1->Get(h2);
    gram2[0] = (TH1F*)f_in2->Get(h1);
    gram2[1] = (TH1F*)f_in2->Get(h2);
    gram3[0] = (TH1F*)f_in3->Get(h1);
    gram3[1] = (TH1F*)f_in3->Get(h2);
    
    TH1F* h_ratio1 = (TH1F*)gram1[0]->Clone("ratio_25");
    h_ratio1->Divide(gram1[0],gram1[1],1,1,"B");
    h_ratio1->SetLineColor(1);
    h_ratio1->SetLineWidth(3);
    h_ratio1->SetMarkerStyle(20);
    h_ratio1->SetTitle("B-Tag Eff vs r2v(MV2c10:0.83);r2v(mm);Efficiency");
    //l->AddEntry(h_ratio1,"d0max=5.0mm, z0max=5.0mm","l");
    h_ratio1->Write();

    TH1F* h_ratio2 = (TH1F*)gram2[0]->Clone("ratio_50");
    h_ratio2->Divide(gram2[0],gram2[1],1,1,"B");
    h_ratio2->SetLineColor(2);
    h_ratio2->SetLineWidth(3);
    h_ratio2->SetMarkerStyle(20);
    h_ratio2->SetTitle("B-Tag Eff vs r2v(MV2c10:0.83);r2v(mm);Efficiency");
    //l->AddEntry(h_ratio2,"d0max=7.5mm, z0max=7.5mm","l");
    h_ratio2->Write();
    //l->Write("SAME");
 
    
    TH1F* h_ratio3 = (TH1F*)gram3[0]->Clone("ratio_75");
    h_ratio3->Divide(gram3[0],gram3[1],1,1,"B");
    h_ratio3->SetLineColor(3);
    h_ratio3->SetLineWidth(3);
    h_ratio3->SetMarkerStyle(20);
    h_ratio3->SetTitle("B-Tag Eff vs r2v(MV2c10:0.83);r2v(mm);Efficiency");
    //l->AddEntry(h_ratio3,"d0max=7.5mm, z0max=7.5mm","l");
    h_ratio3->Write();
    //l->Write("SAME");



    //c1->Update();
    //c1->SaveAs("ZH_20kEff_MV2c10_r2v.png");
   // delete c1;
    delete f_in1;
    delete f_in2;
    delete f_in3;
    delete fout;


}
void dl1_plot_r2v(const char *dat_file1, const char *dat_file2 ,const char *dat_file3){
    TFile *f_in1 = new TFile(dat_file1,"READ");
    TFile *f_in2 = new TFile(dat_file2,"READ");
    TFile *f_in3 = new TFile(dat_file3,"READ");
    TFile* fout = new TFile("dl1_r2v.root","RECREATE");//file for writing all histograms
     fout->cd();
    // TCanvas *c1 = new TCanvas("c1","c1",1200,800) ;
  
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitleX(0.1f);
    gStyle->SetTitleW(0.8f);
    TLegend* l = new TLegend(0.60,.90,0.75,0.95);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    //l->SetMargin(0.4);
    
    std::string str1;
    std::string str2;
    str1="dl1_tag_r2v_2";
    str2="dl1_all_r2v_2";
    const char *h1 = str1.c_str();
    const char *h2 = str2.c_str();
    
    TH1F* gram1[2];
    TH1F* gram2[2];
    TH1F* gram3[2];     
    gram1[0] = (TH1F*)f_in1->Get(h1);
    gram1[1] = (TH1F*)f_in1->Get(h2);
    gram2[0] = (TH1F*)f_in2->Get(h1);
    gram2[1] = (TH1F*)f_in2->Get(h2);
    gram3[0] = (TH1F*)f_in3->Get(h1);
    gram3[1] = (TH1F*)f_in3->Get(h2);
    
    TH1F* h_ratio1 = (TH1F*)gram1[0]->Clone("ratio_25");
    h_ratio1->Divide(gram1[0],gram1[1],1,1,"B");
    h_ratio1->SetLineColor(1);
    h_ratio1->SetLineWidth(3);
    h_ratio1->SetMarkerStyle(20);
    h_ratio1->SetTitle("B-Tag Eff vs r2v(DL1:2.02);r2v(mm);Efficiency");
    //l->AddEntry(h_ratio1,"d0max=5.0mm, z0max=5.0mm","l");
    h_ratio1->Write();
    
    TH1F* h_ratio2 = (TH1F*)gram2[0]->Clone("ratio_50");
    h_ratio2->Divide(gram2[0],gram2[1],1,1,"B");
    h_ratio2->SetLineColor(2);
    h_ratio2->SetLineWidth(3);
    h_ratio2->SetMarkerStyle(20);
    h_ratio2->SetTitle("B-Tag Eff vs r2v(DL1:2.02);r2v(mm);Efficiency");
    //l->AddEntry(h_ratio2,"d0max=7.5mm, z0max=7.5mm","l");
    h_ratio2->Write();
    //l->Write("SAME");
   
    TH1F* h_ratio3 = (TH1F*)gram3[0]->Clone("ratio_75");
    h_ratio3->Divide(gram3[0],gram3[1],1,1,"B");
    h_ratio3->SetLineColor(3);
    h_ratio3->SetLineWidth(3);
    h_ratio3->SetMarkerStyle(20);
    h_ratio3->SetTitle("B-Tag Eff vs r2v(DL1:2.02);r2v(mm);Efficiency");
    //l->AddEntry(h_ratio3,"d0max=7.5mm, z0max=7.5mm","l");
    h_ratio3->Write();
    //l->Write("SAME");
  
    //c1->Update();
    //c1->SaveAs("ZH_20kEff_DL1_r2v.png");
    //delete c1;
    delete f_in1;
    delete f_in2;
    delete f_in3;
    delete fout;
    
}

void mv2c10_plot_lxy(const char *dat_file1, const char *dat_file2 ,const char *dat_file3){
    TFile *f_in1 = new TFile(dat_file1,"READ");
    TFile *f_in2 = new TFile(dat_file2,"READ");
    TFile *f_in3 = new TFile(dat_file3,"READ");
     TFile* fout = new TFile("mv2c10_lxy.root","RECREATE");//file for writing all histograms
     fout->cd();
    //TCanvas *c1 = new TCanvas("c1","c1",1200,800) ;
  
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitleX(0.1f);
    gStyle->SetTitleW(0.8f);
    TLegend* l = new TLegend(0.60,.90,0.75,0.95);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    //l->SetMargin(0.4);
    
    std::string str1;
    std::string str2;
    str1="mv2c10_tag_lxy_2";
    str2="mv2c10_all_lxy_2";
    const char *h1 = str1.c_str();
    const char *h2 = str2.c_str();
    
    TH1F* gram1[2];
    TH1F* gram2[2];
    TH1F* gram3[2];     
    gram1[0] = (TH1F*)f_in1->Get(h1);
    gram1[1] = (TH1F*)f_in1->Get(h2);
    gram2[0] = (TH1F*)f_in2->Get(h1);
    gram2[1] = (TH1F*)f_in2->Get(h2);
    gram3[0] = (TH1F*)f_in3->Get(h1);
    gram3[1] = (TH1F*)f_in3->Get(h2);
    
    TH1F* h_ratio1 = (TH1F*)gram1[0]->Clone("ratio_25");
    h_ratio1->Divide(gram1[0],gram1[1],1,1,"B");
    h_ratio1->SetLineColor(1);
    h_ratio1->SetLineWidth(3);
    h_ratio1->SetMarkerStyle(20);
    h_ratio1->SetTitle("B-Tag Eff vs lxy(MV2c10:0.83);lxy(mm);Efficiency");
    //l->AddEntry(h_ratio1,"d0max=5.0mm, z0max=5.0mm","l");
    h_ratio1->Write();
    
    TH1F* h_ratio2 = (TH1F*)gram2[0]->Clone("ratio_50");
    h_ratio2->Divide(gram2[0],gram2[1],1,1,"B");
    h_ratio2->SetLineColor(2);
    h_ratio2->SetLineWidth(3);
    h_ratio2->SetMarkerStyle(20);
    h_ratio2->SetTitle("B-Tag Eff vs lxy(MV2c10:0.83);lxy(mm);Efficiency");
    //l->AddEntry(h_ratio2,"d0max=7.5mm, z0max=7.5mm","l");
    h_ratio2->Write();
    //l->Write("SAME");
 
    
    TH1F* h_ratio3 = (TH1F*)gram3[0]->Clone("ratio_75");
    h_ratio3->Divide(gram3[0],gram3[1],1,1,"B");
    h_ratio3->SetLineColor(3);
    h_ratio3->SetLineWidth(3);
    h_ratio3->SetMarkerStyle(20);
    h_ratio3->SetTitle("B-Tag Eff vs lxy(MV2c10:0.83);lxy(mm);Efficiency");
    //l->AddEntry(h_ratio3,"d0max=7.5mm, z0max=7.5mm","l");
    h_ratio3->Write();
    //l->Write("SAME");
    
    //c1->Update();
    //c1->SaveAs("ZH_20kEff_MV2c10_lxy.png");
    //delete c1;
    delete f_in1;
    delete f_in2;
    delete f_in3;
    delete fout;
    
}
void dl1_plot_lxy(const char *dat_file1, const char *dat_file2,const char *dat_file3){
    TFile *f_in1 = new TFile(dat_file1,"READ");
    TFile *f_in2 = new TFile(dat_file2,"READ");
    TFile *f_in3 = new TFile(dat_file3,"READ");
    TFile* fout = new TFile("dl1_lxy.root","RECREATE");//file for writing all histograms
     fout->cd();
    //TCanvas *c1 = new TCanvas("c1","c1",1200,800) ;
  
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitleX(0.1f);
    gStyle->SetTitleW(0.8f);
    TLegend* l = new TLegend(0.60,.90,0.75,0.95);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    //l->SetMargin(0.4);
    
    std::string str1;
    std::string str2;
    
    str1="dl1_tag_lxy_2";
    str2="dl1_all_lxy_2";
    
    const char *h1 = str1.c_str();
    const char *h2 = str2.c_str();
    
    TH1F* gram1[2];
    TH1F* gram2[2];
    TH1F* gram3[2];     
    
    gram1[0] = (TH1F*)f_in1->Get(h1);
    gram1[1] = (TH1F*)f_in1->Get(h2);
    gram2[0] = (TH1F*)f_in2->Get(h1);
    gram2[1] = (TH1F*)f_in2->Get(h2);
    gram3[0] = (TH1F*)f_in3->Get(h1);
    gram3[1] = (TH1F*)f_in3->Get(h2);
    
    TH1F* h_ratio1 = (TH1F*)gram1[0]->Clone("ratio_25");
    h_ratio1->Divide(gram1[0],gram1[1],1,1,"B");
    h_ratio1->SetLineColor(1);
    h_ratio1->SetLineWidth(3);
    h_ratio1->SetMarkerStyle(20);
    h_ratio1->SetTitle("B-Tag Eff vs lxy(DL1:2.02);lxy(mm);Efficiency");
    //l->AddEntry(h_ratio1,"d0max=5.0mm, z0max=5.0mm","l");
    h_ratio1->Write();
    
    TH1F* h_ratio2 = (TH1F*)gram2[0]->Clone("ratio_50");
    h_ratio2->Divide(gram2[0],gram2[1],1,1,"B");
    h_ratio2->SetLineColor(2);
    h_ratio2->SetLineWidth(3);
    h_ratio2->SetMarkerStyle(20);
    h_ratio2->SetTitle("B-Tag Eff vs lxy(DL1:2.02);lxy(mm);Efficiency");
    //l->AddEntry(h_ratio2,"d0max=7.5mm, z0max=7.5mm","l");
    h_ratio2->Write();
    //l->Write("SAME");
    
    TH1F* h_ratio3 = (TH1F*)gram3[0]->Clone("ratio_75");
    h_ratio3->Divide(gram3[0],gram3[1],1,1,"B");
    h_ratio3->SetLineColor(3);
    h_ratio3->SetLineWidth(3);
    h_ratio3->SetMarkerStyle(20);
    h_ratio3->SetTitle("B-Tag Eff vs lxy(DL1:2.02);lxy(mm);Efficiency");
    //l->AddEntry(h_ratio3,"d0max=7.5mm, z0max=7.5mm","l");
    h_ratio3->Write();
    //l->Write("SAME");
    
    //c1->Update();
    //c1->SaveAs("ZH_20kEff_DL1_lxy.png");
    //delete c1;
    delete f_in1;
    delete f_in2;
    delete f_in3;
    delete fout;
    
}

void dl1_plot_lxy_r2v(){
    std::string str1;
    std::string str2;
    std::string str3;
    
    str1="ratio_25";
    str2="ratio_50";
    str3="ratio_75";
    
    const char *h1 = str1.c_str();
    const char *h2 = str2.c_str();
    const char *h3 = str3.c_str();
    
    TFile *f_in1 = new TFile("dl1_lxy.root","READ");
    TFile *f_in2 = new TFile("dl1_r2v.root","READ");
        
    TCanvas *c1 = new TCanvas("c1","c1",1200,800) ;
  
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitleX(0.1f);
    gStyle->SetTitleW(0.8f);
    TLegend* l = new TLegend(0.60,.90,0.75,0.95);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
        
    TH1F* gram1[3];
    TH1F* gram2[3];
        
    gram1[0] = (TH1F*)f_in1->Get(h1);
    gram1[1] = (TH1F*)f_in1->Get(h2);
    gram1[2] = (TH1F*)f_in1->Get(h3);
    
    gram2[0] = (TH1F*)f_in2->Get(h1);
    gram2[1] = (TH1F*)f_in2->Get(h2);
    gram2[2] = (TH1F*)f_in2->Get(h3);
    
    const Int_t size = gram1[0]->GetXaxis()->GetNbins();
    float gram1A[size]; float gram1B[size]; float gram1C[size];
    float gram2A[size]; float gram2B[size]; float gram2C[size];
    
    for(Int_t i=0; i<gram1[0]->GetXaxis()->GetNbins();i++) {
        gram1A[i]=gram1[0]->GetBinContent(i);
        gram2A[i]=gram2[0]->GetBinContent(i);
        gram1B[i]=gram1[1]->GetBinContent(i);
        gram2B[i]=gram2[1]->GetBinContent(i);
        gram1C[i]=gram1[2]->GetBinContent(i);
        gram2C[i]=gram2[2]->GetBinContent(i);
        
    }

    TGraph *mg1 = new TGraph(size,gram1A,gram2A);
    TGraph *mg2 = new TGraph(size,gram1B,gram2B);
    TGraph *mg3 = new TGraph(size,gram1C,gram2C);
    
    //mg1->SetLineColor(2);mg2->SetLineColor(3);mg3->SetLineColor(4);
    mg1->SetMarkerColor(2);mg2->SetMarkerColor(3);mg3->SetMarkerColor(4);
    mg1->SetMarkerStyle(20); mg2->SetMarkerStyle(20); mg3->SetMarkerStyle(20);
    
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("r2v vs lxy(DL1:2.02)");
    mg->GetXaxis()->SetTitle("lxy Efficiency");
    mg->GetYaxis()->SetTitle("r2v Efficiency");
    
    mg->Add(mg1);
    mg->Add(mg2);
    mg->Add(mg3);
    
    l->AddEntry(mg1,"d0max=2.5mm, z0max=2.5mm","l");
    l->AddEntry(mg2,"d0max=5.0mm, z0max=5.0mm","l");
    l->AddEntry(mg3,"d0max=7.5mm, z0max=7.5mm","l");
    
    mg->Draw("AP");
    l->Draw("SAME");
   
    c1->Update();
    c1->SaveAs("DL1_r2v_vs_lxy.png");
    
    delete c1;
    delete f_in1;
    delete f_in2;
   
}

void mv2c10_plot_lxy_r2v(){
    std::string str1;
    std::string str2;
    std::string str3;
    
    str1="ratio_25";
    str2="ratio_50";
    str3="ratio_75";
    
    const char *h1 = str1.c_str();
    const char *h2 = str2.c_str();
    const char *h3 = str3.c_str();
    
    TFile *f_in1 = new TFile("mv2c10_lxy.root","READ");
    TFile *f_in2 = new TFile("mv2c10_r2v.root","READ");
        
    TCanvas *c1 = new TCanvas("c1","c1",1200,800) ;
  
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetTitleX(0.1f);
    gStyle->SetTitleW(0.8f);
    TLegend* l = new TLegend(0.60,.90,0.75,0.95);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
        
    TH1F* gram1[3];
    TH1F* gram2[3];
        
    gram1[0] = (TH1F*)f_in1->Get(h1);
    gram1[1] = (TH1F*)f_in1->Get(h2);
    gram1[2] = (TH1F*)f_in1->Get(h3);
    
    gram2[0] = (TH1F*)f_in2->Get(h1);
    gram2[1] = (TH1F*)f_in2->Get(h2);
    gram2[2] = (TH1F*)f_in2->Get(h3);
    
    const Int_t size = gram1[0]->GetXaxis()->GetNbins();
    float gram1A[size]; float gram1B[size]; float gram1C[size];
    float gram2A[size]; float gram2B[size]; float gram2C[size];
    
    for(Int_t i=0; i<gram1[0]->GetXaxis()->GetNbins();i++) {
        gram1A[i]=gram1[0]->GetBinContent(i);
        gram2A[i]=gram2[0]->GetBinContent(i);
        gram1B[i]=gram1[1]->GetBinContent(i);
        gram2B[i]=gram2[1]->GetBinContent(i);
        gram1C[i]=gram1[2]->GetBinContent(i);
        gram2C[i]=gram2[2]->GetBinContent(i);
        
    }

    TGraph *mg1 = new TGraph(size,gram1A,gram2A);
    TGraph *mg2 = new TGraph(size,gram1B,gram2B);
    TGraph *mg3 = new TGraph(size,gram1C,gram2C);
    
    //mg1->SetLineColor(2);mg2->SetLineColor(3);mg3->SetLineColor(4);
    mg1->SetMarkerColor(2);mg2->SetMarkerColor(3);mg3->SetMarkerColor(4);
    mg1->SetMarkerStyle(20); mg2->SetMarkerStyle(20); mg3->SetMarkerStyle(20);
    
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("r2v vs lxy(MV2c10:0.83)");
    mg->GetXaxis()->SetTitle("lxy Efficiency");
    mg->GetYaxis()->SetTitle("r2v Efficiency");
    
    mg->Add(mg1);
    mg->Add(mg2);
    mg->Add(mg3);
    
    l->AddEntry(mg1,"d0max=2.5mm, z0max=2.5mm","l");
    l->AddEntry(mg2,"d0max=5.0mm, z0max=5.0mm","l");
    l->AddEntry(mg3,"d0max=7.5mm, z0max=7.5mm","l");
    
    mg->Draw("AP");
    l->Draw("SAME");
   
    c1->Update();
    c1->SaveAs("MV2c10_r2v_vs_lxy.png");
    
    delete c1;
    delete f_in1;
    delete f_in2;
   
   
}

void run(){
    
    std::string str1;
    std::string str2;
    std::string str3;
    std::string str4;
    std::string str5;
    std::string str6;
    str1="ZH25.root";
    str2="ZH50.root";
    str3="ZH75.root";
    str4="ZH25_ana.root";
    str5="ZH50_ana.root";
    str6="ZH75_ana.root";
    const char *a1 = str1.c_str();
    const char *a2 = str2.c_str();
    const char *a3 = str3.c_str();
    const char *a4 = str4.c_str();
    const char *a5 = str5.c_str();
    const char *a6 = str6.c_str();
    
    bteff(a1,a4);
    bteff(a2,a5);
    bteff(a3,a6);
    mv2c10_plot_r2v(a4,a5,a6);
    dl1_plot_r2v(a4,a5,a6);
    mv2c10_plot_lxy(a4,a5,a6);
    dl1_plot_lxy(a4,a5,a6);
    mv2c10_plot_lxy_r2v();
    dl1_plot_lxy_r2v();
    
}
    
