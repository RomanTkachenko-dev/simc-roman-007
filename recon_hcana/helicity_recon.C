/*
 * Description:
 * ================================================================
 * Time-stamp: "2024-06-06 16:29:53 trottar"
 * ================================================================
 *
 * Author:  Richard L. Trotta III <trotta@cua.edu>, Carlos Yero <cyero002@fiu.edu, cyero@jlab.org>
 *
 * Copyright (c) trottar
 */

#include<iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TCutG.h"
#include "recon_hcana.h"

using namespace std;

recon_hcana::recon_hcana(TString filename, TString reaction_str) {
  
  ReadReaction(reaction_str);
  
  buildFileName(filename);

  InSIMCHistname = InSIMCFilename + ".hist";
  InSIMCRootname = InSIMCFilename + ".root";
  
  cout << "InSIMCFilename: " << InSIMCFilename << endl;
  cout << "InSIMCHistname: " << InSIMCHistname << endl;
  cout << "InSIMCRootname: " << InSIMCRootname << endl;

  simc_nevents = stod(split(FindString("Ngen",InSIMCHistname)[0], '=')[1]);
  simc_normfactor = stod(split(FindString("normfac",InSIMCHistname)[0], '=')[1]);
  Ein = stod(split(FindString("Ebeam",InSIMCHistname)[0], '=')[1]);
  kf0 = num_split(split(FindString("momentum",InSIMCHistname)[0], '=')[1])[0];
  e_th = -num_split(split(FindString("angle",InSIMCHistname)[0], '=')[1])[0]; // Negative sign for HMS spec
  Pf0 = num_split(split(FindString("momentum",InSIMCHistname)[0], '=')[1])[1];
  h_th = num_split(split(FindString("angle",InSIMCHistname)[0], '=')[1])[1];

  cout << "Ngen: " << simc_nevents << endl;
  cout << "normfac: " << simc_normfactor << endl;
  cout << "Ein: " << Ein << endl;
  cout << "kf: " << kf0 << endl;
  cout << "e_th: " << e_th << endl;
  cout << "Pf: " << Pf0 << endl;
  cout << "h_th: " << h_th << endl;  

  //-----If H(e,e'p)
  if(reaction=="heep"){
    
    HeepReadTree();    
  
  }
  
  //-----Base production from SIMC
  else if(reaction=="production"){
    
    ProductionReadTree();
  
  }
  
  else{

    cerr << "ERROR: Invalid reaction of type '" << reaction << "'!\n\tMust be 'heep' or 'production'" << endl;
    exit(EXIT_FAILURE);
  }

  EventLoop();
  WriteHist();
  
}

void recon_hcana::ProductionReadTree(){
  
  cout << "Calling ProductionReadTree() . . . " << endl;

  f = new TFile(InSIMCRootname,"UPDATE");
  if (!f || f->IsZombie()) {
      cerr << "Error: failed to open file " << InSIMCRootname << endl;
      return;
   }
  
  tree = (TTree*)f->Get("h10");
  if (!tree) {
    cerr << "Error: failed to get tree h10 from file " << InSIMCRootname << endl;
    f->Close();
    return;
   }

  newTree = new TTree("h10", "A modified version of the original tree");
  
  //tree->GetListOfBranches()->Print();

  nentries = tree->GetEntries();
  
  tree->SetBranchAddress("hsdelta", &hsdelta);
  tree->SetBranchAddress("hsyptar", &hsyptar);
  tree->SetBranchAddress("hsxptar", &hsxptar);
  tree->SetBranchAddress("hsytar", &hsytar);
  tree->SetBranchAddress("hsxfp", &hsxfp);
  tree->SetBranchAddress("hsxpfp", &hsxpfp);
  tree->SetBranchAddress("hsyfp", &hsyfp);
  tree->SetBranchAddress("hsypfp", &hsypfp);
  tree->SetBranchAddress("hsdeltai", &hsdeltai);
  tree->SetBranchAddress("hsyptari", &hsyptari);
  tree->SetBranchAddress("hsxptari", &hsxptari);
  tree->SetBranchAddress("hsytari", &hsytari);
  tree->SetBranchAddress("ssdelta", &ssdelta);
  tree->SetBranchAddress("ssyptar", &ssyptar);
  tree->SetBranchAddress("ssxptar", &ssxptar);
  tree->SetBranchAddress("ssytar", &ssytar);
  tree->SetBranchAddress("ssxfp", &ssxfp);
  tree->SetBranchAddress("ssxpfp", &ssxpfp);
  tree->SetBranchAddress("ssyfp", &ssyfp);
  tree->SetBranchAddress("ssypfp", &ssypfp);
  tree->SetBranchAddress("ssdeltai", &ssdeltai);
  tree->SetBranchAddress("ssyptari", &ssyptari);
  tree->SetBranchAddress("ssxptari", &ssxptari);
  tree->SetBranchAddress("ssytari", &ssytari);
  tree->SetBranchAddress("q", &q);
  tree->SetBranchAddress("nu", &nu);
  tree->SetBranchAddress("Q2", &Q2);
  tree->SetBranchAddress("W", &W);
  tree->SetBranchAddress("epsilon", &epsilon);
  tree->SetBranchAddress("epscm", &epscm);  
  tree->SetBranchAddress("Em", &Em);
  tree->SetBranchAddress("Pm", &Pm);
  tree->SetBranchAddress("thetapq", &thetapq);
  tree->SetBranchAddress("thetacm", &thetacm);
  tree->SetBranchAddress("phipq", &phipq);
  tree->SetBranchAddress("missmass", &missmass);
  tree->SetBranchAddress("mmnuc", &mmnuc);
  tree->SetBranchAddress("phad", &phad);
  tree->SetBranchAddress("t", &t);
  tree->SetBranchAddress("pmpar", &pmpar);
  tree->SetBranchAddress("pmper", &pmper);
  tree->SetBranchAddress("pmoop", &pmoop);
  tree->SetBranchAddress("fry", &fry);
  tree->SetBranchAddress("radphot", &radphot);
  tree->SetBranchAddress("pfermi", &pfermi);
  tree->SetBranchAddress("siglab", &siglab);
  tree->SetBranchAddress("sigcm", &sigcm);
  tree->SetBranchAddress("Weight", &Weight);
  tree->SetBranchAddress("decdist", &decdist);
  tree->SetBranchAddress("Mhadron", &Mhadron);
  tree->SetBranchAddress("pdotqhat", &pdotqhat);
  tree->SetBranchAddress("Q2i", &Q2i);
  tree->SetBranchAddress("Wi", &Wi);
  tree->SetBranchAddress("ti", &ti);
  tree->SetBranchAddress("phipqi", &phipqi);
  tree->SetBranchAddress("saghai", &saghai);
  tree->SetBranchAddress("factor", &factor);
  tree->SetBranchAddress("helicity", &helicity);
  
  newTree->Branch("hsdelta", &hsdelta, "hsdelta/F");
  newTree->Branch("hsyptar", &hsyptar, "hsyptar/F");
  newTree->Branch("hsxptar", &hsxptar, "hsxptar/F");
  newTree->Branch("hsytar", &hsytar, "hsytar/F");
  newTree->Branch("hsxfp", &hsxfp, "hsxfp/F");
  newTree->Branch("hsxpfp", &hsxpfp, "hsxpfp/F");
  newTree->Branch("hsyfp", &hsyfp, "hsyfp/F");
  newTree->Branch("hsypfp", &hsypfp, "hsypfp/F");
  newTree->Branch("hsdeltai", &hsdeltai, "hsdeltai/F");
  newTree->Branch("hsyptari", &hsyptari, "hsyptari/F");
  newTree->Branch("hsxptari", &hsxptari, "hsxptari/F");
  newTree->Branch("hsytari", &hsytari, "hsytari/F");
  newTree->Branch("ssdelta", &ssdelta, "ssdelta/F");
  newTree->Branch("ssyptar", &ssyptar, "ssyptar/F");
  newTree->Branch("ssxptar", &ssxptar, "ssxptar/F");
  newTree->Branch("ssytar", &ssytar, "ssytar/F");
  newTree->Branch("ssxfp", &ssxfp, "ssxfp/F");
  newTree->Branch("ssxpfp", &ssxpfp, "ssxpfp/F");
  newTree->Branch("ssyfp", &ssyfp, "ssyfp/F");
  newTree->Branch("ssypfp", &ssypfp, "ssypfp/F");
  newTree->Branch("ssdeltai", &ssdeltai, "ssdeltai/F");
  newTree->Branch("ssyptari", &ssyptari, "ssyptari/F");
  newTree->Branch("ssxptari", &ssxptari, "ssxptari/F");
  newTree->Branch("ssytari", &ssytari, "ssytari/F");
  newTree->Branch("q", &q, "q/F");
  newTree->Branch("nu", &nu, "nu/F");
  newTree->Branch("Q2", &Q2, "Q2/F");
  newTree->Branch("W", &W, "W/F");
  newTree->Branch("epsilon", &epsilon, "epsilon/F");
  newTree->Branch("epscm", &epscm, "epscm/F");  
  newTree->Branch("Em", &Em, "Em/F");
  newTree->Branch("Pm", &Pm, "Pm/F");
  newTree->Branch("thetapq", &thetapq, "thetapq/F");
  newTree->Branch("thetacm", &thetacm, "thetacm/F");
  newTree->Branch("phipq", &phipq, "phipq/F");
  newTree->Branch("missmass", &missmass, "missmass/F");
  newTree->Branch("mmnuc", &mmnuc, "mmnuc/F");
  newTree->Branch("phad", &phad, "phad/F");
  newTree->Branch("t", &t, "t/F");
  newTree->Branch("pmpar", &pmpar, "pmpar/F");
  newTree->Branch("pmper", &pmper, "pmper/F");
  newTree->Branch("pmoop", &pmoop, "pmoop/F");
  newTree->Branch("fry", &fry, "fry/F");
  newTree->Branch("radphot", &radphot, "radphot/F");
  newTree->Branch("pfermi", &pfermi, "pfermi/F");
  newTree->Branch("siglab", &siglab, "siglab/F");
  newTree->Branch("sigcm", &sigcm, "sigcm/F");
  newTree->Branch("Weight", &Weight, "Weight/F");
  newTree->Branch("decdist", &decdist, "decdist/F");
  newTree->Branch("Mhadron", &Mhadron, "Mhadron/F");
  newTree->Branch("pdotqhat", &pdotqhat, "pdotqhat/F");
  newTree->Branch("Q2i", &Q2i, "Q2i/F");
  newTree->Branch("Wi", &Wi, "Wi/F");
  newTree->Branch("ti", &ti, "ti/F");
  newTree->Branch("phipqi", &phipqi, "phipqi/F");
  newTree->Branch("saghai", &saghai, "saghai/F");
  newTree->Branch("factor", &factor, "factor/F");
  newTree->Branch("paero_z_det", &paero_z_det, "paero_z_det/F");
  newTree->Branch("paero_x_det", &paero_x_det, "paero_x_det/F");
  newTree->Branch("paero_y_det", &paero_y_det, "paero_y_det/F");  
  newTree->Branch("phgcer_z_det", &phgcer_z_det, "phgcer_z_det/F");
  newTree->Branch("phgcer_x_det", &phgcer_x_det, "phgcer_x_det/F");
  newTree->Branch("phgcer_y_det", &phgcer_y_det, "phgcer_y_det/F");
  newTree->Branch("pend_z_det", &pend_z_det, "pend_z_det/F");
  newTree->Branch("pend_x_det", &pend_x_det, "pend_x_det/F");
  newTree->Branch("pend_y_det", &pend_y_det, "pend_y_det/F");
  newTree->Branch("helicity", &helicity, "helicity/F");
  
  //newTree = tree->CloneTree();
  
  cout << "Ending ProductionReadTree() . . . " << endl;

}

void recon_hcana::HeepReadTree(){
  
  cout << "Calling HeepReadTree() . . . " << endl;

  f = new TFile(InSIMCRootname,"UPDATE");
  if (!f || f->IsZombie()) {
      cerr << "Error: failed to open file " << InSIMCRootname << endl;
      return;
   }
  
  tree = (TTree*)f->Get("h10");
  if (!tree) {
    cerr << "Error: failed to get tree h10 from file " << InSIMCRootname << endl;
    f->Close();
    return;
   }

  newTree = new TTree("h10", "A modified version of the original tree");
  
  //tree->GetListOfBranches()->Print();

  nentries = tree->GetEntries();
  
  tree->SetBranchAddress("hsdelta", &hsdelta);
  tree->SetBranchAddress("hsyptar", &hsyptar);
  tree->SetBranchAddress("hsxptar", &hsxptar);
  tree->SetBranchAddress("hsytar", &hsytar);
  tree->SetBranchAddress("hsxfp", &hsxfp);
  tree->SetBranchAddress("hsxpfp", &hsxpfp);
  tree->SetBranchAddress("hsyfp", &hsyfp);
  tree->SetBranchAddress("hsypfp", &hsypfp);
  tree->SetBranchAddress("hsdeltai", &hsdeltai);
  tree->SetBranchAddress("hsyptari", &hsyptari);
  tree->SetBranchAddress("hsxptari", &hsxptari);
  tree->SetBranchAddress("hsytari", &hsytari);
  tree->SetBranchAddress("ssdelta", &ssdelta);
  tree->SetBranchAddress("ssyptar", &ssyptar);
  tree->SetBranchAddress("ssxptar", &ssxptar);
  tree->SetBranchAddress("ssytar", &ssytar);
  tree->SetBranchAddress("ssxfp", &ssxfp);
  tree->SetBranchAddress("ssxpfp", &ssxpfp);
  tree->SetBranchAddress("ssyfp", &ssyfp);
  tree->SetBranchAddress("ssypfp", &ssypfp);
  tree->SetBranchAddress("ssdeltai", &ssdeltai);
  tree->SetBranchAddress("ssyptari", &ssyptari);
  tree->SetBranchAddress("ssxptari", &ssxptari);
  tree->SetBranchAddress("ssytari", &ssytari);
  tree->SetBranchAddress("q", &q);
  tree->SetBranchAddress("nu", &nu);
  tree->SetBranchAddress("Q2", &Q2);
  tree->SetBranchAddress("W", &W);
  tree->SetBranchAddress("epsilon", &epsilon);
  tree->SetBranchAddress("epscm", &epscm);  
  tree->SetBranchAddress("Em", &Em);
  tree->SetBranchAddress("Pm", &Pm);
  tree->SetBranchAddress("thetapq", &thetapq);
  tree->SetBranchAddress("thetacm", &thetacm);
  tree->SetBranchAddress("phipq", &phipq);
  tree->SetBranchAddress("corrsing", &corrsing);
  tree->SetBranchAddress("Pmx", &Pmx);
  tree->SetBranchAddress("Pmy", &Pmy);
  tree->SetBranchAddress("Pmz", &Pmz);
  tree->SetBranchAddress("PmPar", &PmPar);
  tree->SetBranchAddress("PmPer", &PmPer);
  tree->SetBranchAddress("PmOop", &PmOop);
  tree->SetBranchAddress("fry", &fry);
  tree->SetBranchAddress("radphot", &radphot);
  tree->SetBranchAddress("sigcc", &sigcc);
  tree->SetBranchAddress("Weight", &Weight);  

  newTree->Branch("hsdelta", &hsdelta, "hsdelta/F");
  newTree->Branch("hsyptar", &hsyptar, "hsyptar/F");
  newTree->Branch("hsxptar", &hsxptar, "hsxptar/F");
  newTree->Branch("hsytar", &hsytar, "hsytar/F");
  newTree->Branch("hsxfp", &hsxfp, "hsxfp/F");
  newTree->Branch("hsxpfp", &hsxpfp, "hsxpfp/F");
  newTree->Branch("hsyfp", &hsyfp, "hsyfp/F");
  newTree->Branch("hsypfp", &hsypfp, "hsypfp/F");
  newTree->Branch("hsdeltai", &hsdeltai, "hsdeltai/F");
  newTree->Branch("hsyptari", &hsyptari, "hsyptari/F");
  newTree->Branch("hsxptari", &hsxptari, "hsxptari/F");
  newTree->Branch("hsytari", &hsytari, "hsytari/F");
  newTree->Branch("ssdelta", &ssdelta, "ssdelta/F");
  newTree->Branch("ssyptar", &ssyptar, "ssyptar/F");
  newTree->Branch("ssxptar", &ssxptar, "ssxptar/F");
  newTree->Branch("ssytar", &ssytar, "ssytar/F");
  newTree->Branch("ssxfp", &ssxfp, "ssxfp/F");
  newTree->Branch("ssxpfp", &ssxpfp, "ssxpfp/F");
  newTree->Branch("ssyfp", &ssyfp, "ssyfp/F");
  newTree->Branch("ssypfp", &ssypfp, "ssypfp/F");
  newTree->Branch("ssdeltai", &ssdeltai, "ssdeltai/F");
  newTree->Branch("ssyptari", &ssyptari, "ssyptari/F");
  newTree->Branch("ssxptari", &ssxptari, "ssxptari/F");
  newTree->Branch("ssytari", &ssytari, "ssytari/F");
  newTree->Branch("q", &q, "q/F");
  newTree->Branch("nu", &nu, "nu/F");
  newTree->Branch("Q2", &Q2, "Q2/F");
  newTree->Branch("W", &W, "W/F");
  newTree->Branch("epsilon", &epsilon, "epsilon/F");
  newTree->Branch("epscm", &epscm, "epscm/F");  
  newTree->Branch("Em", &Em, "Em/F");
  newTree->Branch("Pm", &Pm, "Pm/F");
  newTree->Branch("thetapq", &thetapq, "thetapq/F");
  newTree->Branch("thetacm", &thetacm, "thetacm/F");
  newTree->Branch("phipq", &phipq, "phipq/F");
  newTree->Branch("corrsing", &corrsing, "corrsing/F");
  newTree->Branch("Pmx", &Pmx, "Pmx/F");
  newTree->Branch("Pmy", &Pmy, "Pmy/F");
  newTree->Branch("Pmz", &Pmz, "Pmz/F");
  newTree->Branch("PmPar", &PmPar, "PmPar/F");
  newTree->Branch("PmPer", &PmPer, "PmPer/F");
  newTree->Branch("PmOop", &PmOop, "PmOop/F");
  newTree->Branch("fry", &fry, "fry/F");
  newTree->Branch("radphot", &radphot, "radphot/F");
  newTree->Branch("sigcc", &sigcc, "sigcc/F");
  newTree->Branch("Weight", &Weight, "Weight/F");
  
  //newTree = tree->CloneTree();
  
  cout << "Ending HeepReadTree() . . . " << endl;

}

void recon_hcana::EventLoop(){

  cout << "Calling EventLoop() . . . " << endl;
  
  //Convert MeV to GeV
  Ein = Ein / 1000.;     //incident beam energy
  kf0 = kf0 / 1000.;       //final electron momentum
  Pf0 = Pf0 / 1000.;       //final proton momentum

  for (Int_t i=0;i<nentries;i++) {
    
    // Progress bar
    if(i%1000==0) {	    
      int barWidth = 25;
      progress = ((float)i/(float)nentries);	    
      // cout<<i<<"/"<<nentries<<endl;
      // cout << progress << endl;
      cout << "[";
      float pos = barWidth * progress;
      for (float i = 0.; i < barWidth; ++i) {
	if (i < pos) cout << "=";
	else if (i == pos) cout << ">";
	else cout << " ";
      }
      cout << "] " << int(progress * 100.0) << " %\r";
      cout.flush();
    }	 

    tree->GetEntry(i);

    //--------Calculated Kinematic Varibales----------------

    kf = kf0*(1.0+hsdelta/100); // Corrected final electron momentum 
    Pf = Pf0*(1.0+ssdelta/100); // Corrected final proton momentum 
    
    ki = sqrt(Ein*Ein - me*me);        //initial electron momentum
    
    //redefine
    Ep = sqrt(MP*MP + Pf*Pf);
    En = sqrt(MN*MN + Pm*Pm);

    // cout << "Em: " << Em << endl;
    // cout << "Pm: " << Pm << endl;
    // cout << "Ep: " << Ep << endl;
    // cout << "En: " << En << endl;
    
    Kp = Ep - MP;
    Kn = En - MN;

    // cout << "Kp: " << Kp << endl;
    // cout << "Kn: " << Kn << endl;
    
    //Em = nu - Kp - Kn;

    // cout << "Em: " << Em << endl;
    
    // cout << "MM2: " << MM2 << endl;
    
    //W2 = W*W;

    /*
    //Use hcana formula to re-define HMS/SHMS Ztarget
    htar_z = ((h_ytar + h_yMisPoint)-xBPM*(cos(h_th*dtr)-h_yptar*sin(h_th*dtr)))/(-sin(h_th*dtr)-h_yptar*cos(h_th*dtr));
    etar_z = ((e_ytar - e_yMisPoint)-xBPM*(cos(e_th*dtr)-e_yptar*sin(e_th*dtr)))/(-sin(e_th*dtr)-e_yptar*cos(e_th*dtr));
    
    ztar_diff = htar_z - etar_z;
	  
    X = Q2 / (2.*MP*nu);                           
    th_q = acos( (ki - kf*cos(theta_e))/q );       

    //Define Dipole Exit
    xdip_hms = h_xfp - 147.48*h_xpfp;
    ydip_hms = h_yfp - 147.48*h_ypfp;
	  
    xdip_shms = e_xfp - 307.*e_xpfp;
    ydip_shms = e_yfp - 307.*e_ypfp;
    */
    
    //---------------------------------------------------

    //---------Calculate Pmx, Pmy, Pmz in the Lab, and in the q-system----------------

    //Calculate electron final momentum 3-vector
    SetCentralAngles(e_th, e_ph);
    TransportToLab(kf, hsxptar, hsyptar, kf_vec);

    // Testing with ss(x)yptar
    // Apply OOP offsets (p_oopcentral_offset = -0.00011)
    //ssxptar = ssxptar - 0.00011;
    //TransportToLab(kf, ssxptar, ssyptar, kf_vec);

    // cout << "kf_vec.X(): " << kf_vec.X() << endl;
    // cout << "kf_vec.Y(): " << kf_vec.Y() << endl;
    // cout << "kf_vec.Z(): " << kf_vec.Z() << endl;
    
    //Calculate 4-Vectors
    fP0.SetXYZM(0.0, 0.0, ki, me);  //set initial e- 4-momentum
    fP1.SetXYZM(kf_vec.X(), kf_vec.Y(), kf_vec.Z(), me);  //set final e- 4-momentum
    fA.SetXYZM(0.0, 0.0, 0.0, tgt_mass );  //Set initial target at rest
    fQ = fP0 - fP1;
    
    fScatAngle = fP0.Angle( fP1.Vect());    

    fMp.SetXYZM(0.0, 0.0, 0.0, MP);

    fMp1 = fMp + fQ;

    // Redefine higher order variables
    W2 = fMp1.M2();
    if (W2>0)  W = TMath::Sqrt(W2);
    q = fQ.P();
    Q2 = -fQ.M2();
    epsilon = 1.0 / ( 1.0 + 2.0*q*q/Q2*TMath::Power( TMath::Tan(fScatAngle/2.0), 2.0 ));
    
    fA1 = fA + fQ;   //final target (sum of final hadron four momenta)

    //Get Detected Particle 4-momentum
    SetCentralAngles(h_th, h_ph);
    // Apply OOP offsets (p_oopcentral_offset = -0.00011)
    ssxptar = ssxptar - 0.00011;
    TransportToLab(Pf, ssxptar, ssyptar, Pf_vec);

    if(reaction=="heep"){
      fX.SetVectM(Pf_vec, MP);       //SET FOUR VECTOR OF detected particle
    }else{
      // For KaonLT, this is mk (kaon)
      // ACP: u-channel needs MP (proton)
      fX.SetVectM(Pf_vec, mpi);       //SET FOUR VECTOR OF detected particle      
    }
      
    fB = fA1 - fX;                 //4-MOMENTUM OF UNDETECTED PARTICLE 

    Pmx_lab = fB.X();
    Pmy_lab = fB.Y(); 
    Pmz_lab = fB.Z();

    nu = fQ.E();
  
    Em = nu + fA.M() - fX.E();   
    
    // cout << "Pmx_lab: " << Pmx_lab << endl;
    // cout << "Pmy_lab: " << Pmy_lab << endl;
    // cout << "Pmz_lab: " << Pmz_lab << endl;

    //--------Rotate the recoil system from +z to +q-------
    // Angles of X and B wrt q-vector 
    // xq and bq are the 3-momentum vectors of X and B expressed in
    // the coordinate system where q is the z-axis and the x-axis
    // lies in the scattering plane (defined by q and e') and points
    // in the direction of e', so the out-of-plane angle lies within
    // -90<phi_xq<90deg if X is detected on the downstream/forward side of q
    rot_to_q.SetZAxis( fQ.Vect(), fP1.Vect()).Invert();

    xq = fX.Vect();
    bq = fB.Vect();

    xq *= rot_to_q;
    bq *= rot_to_q;

    //Calculate Angles of q relative to x(detected proton) and b(recoil neutron)
    th_pq = xq.Theta();   //"theta_pq"                                       
    ph_pq   = xq.Phi();   //"out-of-plane angle", "phi_pq"
    th_nq = bq.Theta();   // theta_nq                                                                                                     
    ph_nq   = bq.Phi();   //phi_nq

    thetapq = th_pq;
    phipq = ph_pq;
    
    p_miss = -bq;

    //Missing Momentum Components in the q-frame
    Pm = p_miss.Mag(); //=fB.P()
    
    // Redefine variables
    Pmx = p_miss.X();   //in-plane perpendicular component to +z
    Pmy = p_miss.Y();   //out-of-plane component (Oop)
    Pmz = p_miss.Z();   //parallel component to +z

    M_recoil = fB.M(); //recoil mass (missing mass)
    
    // Missing Mass
    missmass = sqrt(abs((Em*Em)-(Pm*Pm)));

    MM2 = missmass * missmass;

    s = (fQ+fA).M2();
    t = (fQ-fX).M2();
    u = (fQ-fB).M2();

    /************************
     ------------------------
     ---- Geometric cuts ----
     ------------------------
     ************************/
    
    /*********************
     **** End of SHMS ****
     *********************/
    // Variable to see geometric cuts at end of spectrometer
    pend_z_det = 300.0; // Approx. end of SHMS (units of cm)
    pend_x_det = ssxfp + pend_z_det*ssxpfp;
    pend_y_det = ssyfp + pend_z_det*ssypfp;
    
    /**********************
     **** SHMS AEROGEL ****
     **********************/
    paero_z_det = 231.0; // Front? of SHMS aerogel (units of cm), see PARAM/SHMS/AERO/KaonLT_PARAM/paero_geom.param
    paero_x_det = ssxfp + paero_z_det*ssxpfp;
    paero_y_det = ssyfp + paero_z_det*ssypfp;

    if (
	(InSIMCFilename.Contains("Q4p4W2p74")) || // High and low epsilon
	(InSIMCFilename.Contains("Q3p0W3p14")) || // High and low epsilon
	(InSIMCFilename.Contains("Q5p5W3p02"))    // High and low epsilon
	){
      
      // SHMS Aero Geom for n = 1.011 (DEF-files/PRODUCTION/KaonLT_DEF/Paero_1p011/Offline_Physics_Coin_Cuts.def)
      // shmsAeroxposalln    P.aero.xAtAero > -45 && P.aero.xAtAero < 45
      // shmsAeroyposalln	   P.aero.yAtAero > -30 && P.aero.yAtAero < 30
      paero_tray_cut = (paero_x_det > -45.0) & (paero_x_det < 45.0) & (paero_y_det > -30) & (paero_y_det < 30);
      
    }else{

      // SHMS Aero Geom for n = All except 1.011 (see DEF-files/PRODUCTION/KaonLT_DEF/Offline_Physics_Coin_Cuts.def)
      // shmsAeroxposalln    P.aero.xAtAero > -55 && P.aero.xAtAero < 55
      // shmsAeroyposalln	   P.aero.yAtAero > -50 && P.aero.yAtAero < 50
      paero_tray_cut = (paero_x_det > -55.0) & (paero_x_det < 55.0) & (paero_y_det > -50) & (paero_y_det < 50);

    }
    
    /**********************
     **** SHMS HGCer ****
     **********************/
    // HGCer Hole cut is now defined in lt_analysis script to be consistent with data procedure.
    // These variables are used to apply such cut.
    phgcer_z_det = 156.27; // Front? of SHMS HGcer (units of cm), see PARAM/SHMS/HGCER/KaonLT_PARAM/phgcer_geom.param
    phgcer_x_det = ssxfp + phgcer_z_det*ssxpfp;
    phgcer_y_det = ssyfp + phgcer_z_det*ssypfp;
    
    if (!(paero_tray_cut)){
      //cout << "Event outside geometric acceptance..." << endl;
      continue; // Skip events outside geometric acceptance
    }
    
    //----------
    
    // cout << "Pmx: " << Pmx << endl;
    // cout << "Pmy: " << Pmy << endl;
    // cout << "Pmz: " << Pmz << endl;
    // cout << "Pm: " << Pm << endl;
    
    newTree->Fill();  
  }
  
  cout << "Ending EventLoop() . . . " << endl;
}

void recon_hcana::WriteHist(){

  cout << "Calling WriteHist() . . . " << endl;
  
  newTree->Write("h10",TObject::kOverwrite);
  f->Close();

  cout << "Ending WriteHist() . . . " << endl;
}
//---------------AUXILIARY FUNCTIONS TO CALCULATE Pmx, Pmy, Pmz in SIMC (same as HCANA) -------------------

//_____________________________________________________
void recon_hcana::GeoToSph( Double_t  th_geo, Double_t  ph_geo, Double_t& th_sph, Double_t& ph_sph){
  
  // Convert geographical to spherical angles. Units are rad.
  
  static const Double_t twopi = 2.0*TMath::Pi();
  Double_t ct = cos(th_geo), cp = cos(ph_geo);
  Double_t tmp = ct*cp;
  th_sph = acos( tmp );
  tmp = sqrt(1.0 - tmp*tmp);
  ph_sph = (fabs(tmp) < 1e-6 ) ? 0.0 : acos( sqrt(1.0-ct*ct)*cp/tmp );
  if( th_geo/twopi-floor(th_geo/twopi) > 0.5 ) ph_sph = TMath::Pi() - ph_sph;
  if( ph_geo/twopi-floor(ph_geo/twopi) > 0.5 ) ph_sph = -ph_sph;
  
  // cout << "th_geo: " << th_geo << endl;
  // cout << "ph_geo: " << ph_geo << endl;
  // cout << "th_sph: " << th_sph << endl;
  // cout << "ph_sph: " << ph_sph << endl;
}

//_______________________________________________________________
void recon_hcana::SetCentralAngles(Double_t th_cent=0, Double_t ph_cent=0){
  
  fThetaGeo = TMath::DegToRad()*th_cent; fPhiGeo = TMath::DegToRad()*ph_cent;
  
  // cout << "th_cent: " << th_cent << endl;
  // cout << "ph_cent: " << ph_cent << endl;
    
  GeoToSph( fThetaGeo, fPhiGeo, fThetaSph, fPhiSph );
  fSinThGeo = TMath::Sin( fThetaGeo ); fCosThGeo = TMath::Cos( fThetaGeo );
  fSinPhGeo = TMath::Sin( fPhiGeo );   fCosPhGeo = TMath::Cos( fPhiGeo );
  Double_t st, ct, sp, cp;
  st = fSinThSph = TMath::Sin( fThetaSph ); ct = fCosThSph = TMath::Cos( fThetaSph );
  sp = fSinPhSph = TMath::Sin( fPhiSph );   cp = fCosPhSph = TMath::Cos( fPhiSph );
  
  Double_t norm = TMath::Sqrt(ct*ct + st*st*cp*cp);
  
  // cout << "norm: " << norm << endl;
  
  TVector3 nx( st*st*sp*cp/norm, -norm, st*ct*sp/norm );
  TVector3 ny( ct/norm,          0.0,   -st*cp/norm   );
  TVector3 nz( st*cp,            st*sp, ct            );

  // cout << "nx.X(): " << nx.X() << endl;
  // cout << "nx.Y(): " << nx.Y() << endl;
  // cout << "nx.Z(): " << nx.Z() << endl;
  
  fToLabRot.SetToIdentity().RotateAxes( nx, ny, nz );
}

//____________________________________________________________________________________
void recon_hcana::TransportToLab( Double_t p, Double_t xptar, Double_t yptar, TVector3& pvect) {
  
  TVector3 v( xptar, yptar, 1.0 );
  v *= p/TMath::Sqrt( 1.0+xptar*xptar+yptar*yptar );

  // cout << "v.X(): " << v.X() << endl;
  // cout << "v.Y(): " << v.Y() << endl;
  // cout << "v.Z(): " << v.Z() << endl;
  
  pvect = fToLabRot * v;

  // cout << "pvect.X(): " << pvect.X() << endl;
  // cout << "pvect.Y(): " << pvect.Y() << endl;
  // cout << "pvect.Z(): " << pvect.Z() << endl;  
}

//------------------------------------------------------------------------------------------


//----------------------------------UTILITIES FUNCTIONS--------------------------------------

vector <string> recon_hcana::FindString(TString keyword, TString fname){

  //Method: Finds string keyword in a given txt file. 
  //Returns the lines (stored in a vector) in which the keyword was found. Lines are counted from 0. 
  
  ifstream ifile(fname);

  vector <string> line_found; //vector to store in which lines was the keyword found
  
  int line_cnt = 0;
  string line;                  //line string to read
  
  int found = -1; //position of found keyword

  while(getline(ifile, line))
    {
      //Check 1st character of found string
      TString cmt = line[0];
      
      found = line.find(keyword);
      
      if(found<0||found>1000){found=-1;} //not found
      //if(cmt==";" || cmt=="#" || cmt=="!") {found=-1;}  //Found commented line. So Skip

      if(found!=-1){
	
	line_found.push_back(line);
	

      } //end if statement
    
      line_cnt++;
    } //end readlines

  return line_found;

}

vector <string> recon_hcana::split(string str, char del=':'){

  //method to split a string into a vetor of strings separated by a delimiter del
  //Returns a vector of strings, whose elements are separated by the delimiter del.

  string temp1, temp2;

  vector <string> parse_word;
  int del_pos;  //delimiter position
    
    for (int i=0; i < str.size(); i++){

      //Get position of delimiter
      if(str[i]==del){
	del_pos = i;
      }

    }

    for (int i=0; i < str.size(); i++){

      //append char to a string for char left of the delimiter
      if(i<del_pos){
	temp1.append(getString(str[i]));
      }      

      //append char to a string for char right of the delimiter
      else if(i>del_pos){
	temp2.append(getString(str[i]));
      }

    }
    
    parse_word.push_back(temp1);
    parse_word.push_back(temp2);
    
    return parse_word;
}

vector <float> recon_hcana::num_split(string str){

  istringstream stream(str);
  vector<float> values;
  float value;
  while (stream >> value) {
    values.push_back(value);
    if (stream.fail()) break;
  }

  return values;
}

string recon_hcana::getString(char x){
  //method to convert a character to a string
  string s(1,x);
  return s;
}

recon_hcana::~recon_hcana(){
  //Destructor

  //delete File; File = NULL;
}  
