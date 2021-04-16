//-------------------------------------------------------
// This is rewritten program of Shota Nagumo-san (U-Tokyo) 
// for the estimation of the effect of the polarized 
// electron beam at the Super Charm-Tau factory 
// 
// Generator of the events with particular e- beam polarisation: 
//  e+ e- --> (tau- --> lep- nu nu; tau+ --> pi+ pi0 nu) 
// 
//               Denis Epifanov 
//
//                 01.04.2019 
//
//------------------------------------------------------- 

#include<iostream>
#include<iomanip> 
#include<fstream> 
#include<stdlib.h> 
#include<stdio.h> 
#include<math.h> 

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TMatrixT.h" 
#include "TMatrixD.h"
#include "TMatrixTUtils.h"
#include "TRandom3.h" 
#include "TFile.h" 
#include "TTree.h" 
#include "TLorentzVector.h" 
#include "TVector3.h" 
#include "TComplex.h" 

void gen(int, int, int); 
void gen(int flr, int iomega, int Nevents){
  gROOT->ProcessLine(".! date");

  if(!(flr == 0 || flr == 1))               return; 
  if(!(iomega >= 0 && iomega <= 10))        return; 
  if(!(Nevents > 0 && Nevents <= 20000000)) return; 

  //  const int Nevents     = 1000000;   // number of generated events 
  //  const double omega    = 0.5;       // beam electron polarisation 
  //  const int flr         = 0;
  //  const double prob_max = 23.0;

  //  const double prob_max = 45.0;  // for E = 5.29 GeV
  const double prob_max = 49.0;      // for E = 1.78, 1.84, 1.89, 2.08, 2.33 GeV 
  const double Ebeam    = 1.84;      // beam energy 

  double omega = 0.1*iomega; 

  //------------------------------------------------------- 

  TString outfile = "out/gen_"; 
  outfile +=     flr; 
  outfile +=     "_"; 
  outfile +=  iomega; 
  outfile += ".root";  
  
  cout<<"---------------------------------------"<<endl; 
  cout<<"Beam energy Ebeam = "<<Ebeam<<endl; 
  cout<<"e- beam polarization: "<<omega<<endl; 
  cout<<"Flag random: "<<flr<<endl; 
  cout<<"Number of events to generate: "<<Nevents<<endl; 
  cout<<"---------------------------------------"<<endl; 

  //-------------------------------------------------------
  
  const double pi = TMath::Pi();

  // e+ e- --> tau+ tau- 

  const double M_tau = 1.77699;     // mass of tau 
  const double M_e   = 0.000511;    // mass of electron 
  const double M_mu  = 0.105658;    // mass of muon 

  double Gamma_tau = Ebeam/M_tau; 
  double P_tau     = sqrt(Ebeam*Ebeam - M_tau*M_tau);
  double Beta_tau  = P_tau/Ebeam;                       // velocity of tau in the LAB frame 
  double P_e       = sqrt(Ebeam*Ebeam - M_e*M_e);       // momentum of beam electron in the LAB frame 
  double beta      = P_e/Ebeam;                         // velocity of e- in the LAB frame 

  // tau --> lep nu nu, lep = e, mu 

  const double rho    = 0.75; 
  const double eta    = 0.00; 
  const double xi     = 1.00; 
  const double delta  = 0.75; 
  const double xi_rho = 1.00; 

  double M_lep = M_mu;                                      // mass of lepton 
  double E_max = 0.5*M_tau*(1.0 + M_lep*M_lep/M_tau/M_tau); // maximal lepton energy in tau rest frame 
  double x0    = M_lep/E_max;

  // tau --> pi pi0 nu 

  const double mpi  = 0.13957;       // mass of pi- 
  const double mpi0 = 0.13498;       // mass of pi0 

  // rho 
  const double M_rho      = 0.7743;  //  mass of rho meson 
  const double Gamma_rho  = 0.1491;  // width of rho meson 
  const double M_rho2     = M_rho*M_rho; 

  // rho' 
  const double M_rhop     = 1.370;   //  mass of rho' 
  const double Gamma_rhop = 0.386;   // width of rho' 
  const double M_rhop2    = M_rhop*M_rhop; 

  // rho-rho' mixing 
  const double betaBW     =-0.091;   // relative coupling to rho' 

  //------------------------------------------------------- 

  TLorentzVector     Ptaum(0,0,1,0); 
  TLorentzVector    Plep_0(0,0,1,0); 
  TLorentzVector      Plep(0,0,1,0); 
  TLorentzVector Ppip_lab0(0,0,1,0);
  TLorentzVector  Ppip_lab(0,0,1,0);  
  TLorentzVector Ppi0_lab0(0,0,1,0);
  TLorentzVector  Ppi0_lab(0,0,1,0);
  TVector3              Nlep(0,0,1); 
  TVector3               RTx(1,0,0);
  TVector3               RTy(0,1,0);
  TVector3               RTz(0,0,1); 

  double tht,pht;  // angle of tau- direction in the LAB frame 

  double x;        // El/Emax 
  double thl,phl;  // angle of lep direction in the tau- rest frame 

  double mpp,mpp2;                                                       
  double thn,phn;  // angle of nu_tau_bar direction in the tau+ rest frame 
  double thp,php;  // angle of pi+ direction in the rho+ rest frame 
        
  double y;
  double prob;

  int seed[2][11] = {{1,  2,  3,  5,  7, 11, 13, 17, 19, 23, 29}, {31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73}}; 
 
  TRandom3 rtht(seed[flr][0]); 
  TRandom3 rpht(seed[flr][1]); 

  TRandom3   rx(seed[flr][2]);
  TRandom3 rthl(seed[flr][3]); 
  TRandom3 rphl(seed[flr][4]); 

  TRandom3 rmpp(seed[flr][5]); 
  TRandom3 rthn(seed[flr][6]); 
  TRandom3 rphn(seed[flr][7]);
  TRandom3 rthp(seed[flr][8]);
  TRandom3 rphp(seed[flr][9]);

  TRandom3   ry(seed[flr][10]);

  TFile *fout  = new TFile(outfile,"recreate");
  TTree *trout = new TTree("tr","total");  

  double  e_tau; // energy, momentum, theta, phi of tau- in the LAB frame 
  double  p_tau;
  double th_tau;
  double ph_tau;
 
  double  e_lep; // energy, momentum, theta, phi of lep in the LAB frame 
  double  p_lep;
  double th_lep;
  double ph_lep;
  
  double  e_pip; // energy, momentum, theta, phi of pi+ in the LAB frame 
  double  p_pip;
  double th_pip;
  double ph_pip;
  
  double  e_pi0; // energy, momentum, theta, phi of pi0 in the LAB frame 
  double  p_pi0;
  double th_pi0;
  double ph_pi0;

  trout->Branch("e_tau", &e_tau, "e_tau/D");  
  trout->Branch("p_tau", &p_tau, "p_tau/D");
  trout->Branch("th_tau",&th_tau,"th_tau/D");
  trout->Branch("ph_tau",&ph_tau,"ph_tau/D");  
  trout->Branch( "e_lep", &e_lep, "e_lep/D");  
  trout->Branch( "p_lep", &p_lep, "p_lep/D");
  trout->Branch("th_lep",&th_lep,"th_lep/D");
  trout->Branch("ph_lep",&ph_lep,"ph_lep/D"); 
  trout->Branch( "e_pip", &e_pip, "e_pip/D");
  trout->Branch( "p_pip", &p_pip, "p_pip/D");
  trout->Branch("th_pip",&th_pip,"th_pip/D");
  trout->Branch("ph_pip",&ph_pip,"ph_pip/D");
  trout->Branch( "e_pi0", &e_pi0, "e_pi0/D");    
  trout->Branch( "p_pi0", &p_pi0, "p_pi0/D");
  trout->Branch("th_pi0",&th_pi0,"th_pi0/D");
  trout->Branch("ph_pi0",&ph_pi0,"ph_pi0/D");  

  int ncycle = 0;
  int ncycle_ave = 0; 

  for(int i = 0; i < Nevents; i++) {
    if((i+1)%10000 == 0){
      cout<<(i+1);
      cout<<" average = "<<ncycle_ave/10000.0<<endl;
      ncycle_ave = 0;
    }
    ncycle_ave += ncycle; 
    ncycle = 0;
    do{
      ncycle++;
      
      tht = acos(2*rtht.Rndm()-1); 
      pht = 2*pi*rpht.Rndm();

      x   = (1-x0)*rx.Rndm()+x0;
      thl = acos(2*rthl.Rndm()-1);
      phl = 2*pi*rphl.Rndm();

      mpp2= rmpp.Rndm()*(M_tau*M_tau - (mpi+mpi0)*(mpi+mpi0)) + (mpi+mpi0)*(mpi+mpi0); 
      mpp = sqrt(mpp2);

      thn = acos(2*rthn.Rndm()-1);
      phn = 2*pi*rphn.Rndm();

      thp = acos(2*rthp.Rndm()-1);
      php = 2*pi*rphp.Rndm();

      y   = ry.Rndm();

      //---------------------------------------------------------------------------------------------------- 
      // e+ e- --> tau+ tau- 

      double D0 = 1 + cos(tht)*cos(tht) + sin(tht)*sin(tht)/Gamma_tau/Gamma_tau; 

      // spin-spin correlation tensor 
      //      TMatrixD d(3,3);
      double d[3][3];
      d[0][0] = (1.0 + 1.0/Gamma_tau/Gamma_tau)*sin(tht)*sin(tht); 
      d[0][1] = 0.0; 
      d[0][2] = sin(2*tht)/Gamma_tau;                                                          
      d[1][0] = 0.0; 
      d[1][1] = -Beta_tau*Beta_tau*sin(tht)*sin(tht); 
      d[1][2] = 0.0; 
      d[2][0] = sin(2*tht)/Gamma_tau; 
      d[2][1] = 0.0; 
      d[2][2] = 1.0 + cos(tht)*cos(tht) - sin(tht)*sin(tht)/Gamma_tau/Gamma_tau;                                           

      RTx.SetXYZ(-cos(tht)*cos(pht),  sin(pht), sin(tht)*cos(pht)); 
      RTy.SetXYZ(-cos(tht)*sin(pht), -cos(pht), sin(tht)*sin(pht)); 
      RTz.SetXYZ(          sin(tht),         0,          cos(tht)); 

      // e- momentum in the tau- rest frame 
      TLorentzVector p1_taum(Ebeam*beta*sin(tht), 0.0,  Gamma_tau*Ebeam*(beta*cos(tht)-Beta_tau), Gamma_tau*Ebeam*(1-Beta_tau*beta*cos(tht)));
      // tau+ momentum in the tau- rest frame
      TLorentzVector ptaup_taum(0.0, 0.0, -2*Ebeam*Gamma_tau*Beta_tau, Gamma_tau*Ebeam*(1+Beta_tau*Beta_tau)); 
      // tau- polarisation factor in the tau- rest frame 
      TLorentzVector G = (1.0/Gamma_tau/Ebeam)*(2*p1_taum-ptaup_taum); 

      // e- momentum in the tau+ rest frame 
      TLorentzVector p1_taup(Ebeam*beta*sin(tht), 0.0,  Gamma_tau*Ebeam*(beta*cos(tht)+Beta_tau), Gamma_tau*Ebeam*(1+Beta_tau*beta*cos(tht))); 
      // tau- momentum in the tau+ rest frame
      TLorentzVector ptaum_taup(0.0, 0.0, 2*Ebeam*Gamma_tau*Beta_tau, Gamma_tau*Ebeam*(1+Beta_tau*Beta_tau)); 
      // tau+ polarisation factor in the tau+ rest frame
      TLorentzVector F = (1.0/Gamma_tau/Ebeam)*(2*p1_taup-ptaum_taup); 

      //----------------------------------------------------------------------------------------------------
      // tau- --> lep- nu nu

      Ptaum.SetPxPyPzE(0,0,1,0);
      Ptaum.SetRho(P_tau);
      Ptaum.SetTheta(tht); 
      Ptaum.SetPhi(pht); 
      Ptaum.SetE(Ebeam);
      //      TVector3 Boost_taum = Ptaum.BoostVector(); 
      TVector3 Boost_taum(0, 0, Beta_tau); 

      // lep- in the LAB frame 
      Plep_0.SetPxPyPzE(0,0,1,0); 
      Plep_0.SetRho(sqrt(x*x*E_max*E_max-M_lep*M_lep));
      Plep_0.SetTheta(thl); 
      Plep_0.SetPhi(phl); 
      Plep_0.SetE(x*E_max);
      Nlep = Plep_0.Vect().Unit(); 
      Plep_0.Boost(Boost_taum);   // boost to lab frame
      Plep.SetPxPyPzE(RTx*Plep_0.Vect(), RTy*Plep_0.Vect(), RTz*Plep_0.Vect(), Plep_0.E());
      
      double A0  = 3*x*(1 - x);
      double A1  = (2.0/3.0)*(4*x*x - 3*x - x0*x0); 
      double A2  = 3*x0*(1 - x);
      double A   = sqrt(x*x - x0*x0)*(A0 + rho*A1 + eta*A2);
     
      double B1  = 1 - x;
      double B2  = (2.0/3.0)*(4*x - 4 + sqrt(1-x0*x0)); 
      TVector3 B = xi*(x*x - x0*x0)*(B1 + delta*B2)*Nlep; 

      //----------------------------------------------------------------------------------------------------- 
      // tau+ --> pi+ pi0 nu 
     
      TLorentzVector Ptaup(0,0,1,0);
      Ptaup.SetRho(P_tau);
      Ptaup.SetTheta(pi-tht); 
      Ptaup.SetPhi(pi+pht); 
      Ptaup.SetE(Ebeam);
      //      TVector3 Boost_taup = Ptaup.BoostVector(); 
      TVector3 Boost_taup(0, 0, -Beta_tau); 

      // rho+ in the tau+ rest frame 
      TLorentzVector Prho(0,0,1,0);
      Prho.SetRho((M_tau*M_tau - mpp2)/2.0/M_tau);
      Prho.SetTheta(pi-thn); 
      Prho.SetPhi(pi+phn); 
      Prho.SetE((M_tau*M_tau + mpp2)/2.0/M_tau);
      TVector3 Boost_rho = Prho.BoostVector(); 
     
      // nu_tau_bar in the tau+ rest frame 
      TLorentzVector Pnu(0,0,1,0);
      Pnu.SetRho((M_tau*M_tau - mpp2)/2.0/M_tau);
      Pnu.SetTheta(thn); 
      Pnu.SetPhi(phn); 
      Pnu.SetE((M_tau*M_tau - mpp2)/2.0/M_tau);

      // pi+ momentum in the rho+ rest frame 
      double ppi_mpp   = sqrt((mpp2-(mpi+mpi0)*(mpi+mpi0))*(mpp2-(mpi-mpi0)*(mpi-mpi0)))/2.0/mpp; 
      double ppi_mrho  = sqrt((M_rho2-(mpi+mpi0)*(mpi+mpi0))*(M_rho2-(mpi-mpi0)*(mpi-mpi0)))/2.0/M_rho; 
      double ppi_mrhop = sqrt((M_rhop2-(mpi+mpi0)*(mpi+mpi0))*(M_rhop2-(mpi-mpi0)*(mpi-mpi0)))/2.0/M_rhop; 

      TLorentzVector Ppip(0,0,1,0); 
      Ppip.SetRho(ppi_mpp);
      Ppip.SetTheta(thp); 
      Ppip.SetPhi(php); 
      Ppip.SetE((mpp2 + mpi*mpi - mpi0*mpi0)/2.0/mpp);

      // boost to the tau+ rest frame
      Ppip.Boost(Boost_rho);
      Ppip_lab0 = Ppip;
      // boost to the LAB frame
      Ppip_lab0.Boost(Boost_taup); 
      Ppip_lab.SetPxPyPzE(RTx*Ppip_lab0.Vect(), RTy*Ppip_lab0.Vect(), RTz*Ppip_lab0.Vect(), Ppip_lab0.E());

      // pi0 momentum in the rho+ rest frame 
      TLorentzVector Ppi0(0,0,1,0); 
      Ppi0.SetRho(ppi_mpp);
      Ppi0.SetTheta(pi-thp); 
      Ppi0.SetPhi(pi+php); 
      Ppi0.SetE((mpp2 - mpi*mpi + mpi0*mpi0)/2.0/mpp);

      // boost to the tau+ rest frame
      Ppi0.Boost(Boost_rho);
      Ppi0_lab0 = Ppi0;
      // boost to the LAB frame
      Ppi0_lab0.Boost(Boost_taup); 
      Ppi0_lab.SetPxPyPzE(RTx*Ppi0_lab0.Vect(), RTy*Ppi0_lab0.Vect(), RTz*Ppi0_lab0.Vect(), Ppi0_lab0.E()); 

      double G_rho   = Gamma_rho*(M_rho/mpp)*(ppi_mpp/ppi_mrho)*(ppi_mpp/ppi_mrho)*(ppi_mpp/ppi_mrho);      // total width of rho 
      double G_rhop  = Gamma_rhop*(M_rhop/mpp)*(ppi_mpp/ppi_mrhop)*(ppi_mpp/ppi_mrhop)*(ppi_mpp/ppi_mrhop); // total width of rho' 
      TComplex BW1   = TComplex(M_rho2, 0)/TComplex(M_rho2-mpp2,-M_rho*G_rho);                              // rho  Breit-Wigner 
      TComplex BW2   = TComplex(M_rhop2,0)/TComplex(M_rhop2-mpp2,-M_rhop*G_rhop);                           // rho' Breit-Wigner 
      TComplex BWtot = (BW1 + betaBW*BW2)/(1.0 + betaBW);
      double BPS     = BWtot.Rho2()*(Prho.P()/M_tau)*(ppi_mpp/mpp);
      
      TLorentzVector Q  = Ppip - Ppi0;
      double B1_        = Q.M2() + 2*(Pnu*Q); 
      double B2_        = Q.M2() - 2*(Pnu*Q); 
      double A_         = 2*(Pnu*Q)*Q.E() - Q.M2()*Pnu.E();  
      TLorentzVector B_ = xi_rho*(B1_*Ppip + B2_*Ppi0); 

      //-----------------------------------------------------------------------------------------------------
      // total PDF 

      //      TMatrixD bb(3,3);
      //      bb[0][0] = B.X()*B_.X(); bb[0][1] = B.X()*B_.Y(); bb[0][2] = B.X()*B_.Z(); 
      //      bb[1][0] = B.Y()*B_.X(); bb[1][1] = B.Y()*B_.Y(); bb[1][2] = B.Y()*B_.Z(); 
      //      bb[2][0] = B.Z()*B_.X(); bb[2][1] = B.Z()*B_.Y(); bb[2][2] = B.Z()*B_.Z(); 
      //      TMatrixD dbb = d*bb; 
      //      double  dbb0 = dbb[0][0]+dbb[1][1]+dbb[2][2];
     
      double dbb0 = B.X()*(B_.X()*d[0][0] + B_.Y()*d[0][1] + B_.Z()*d[0][2]) + 
	            B.Y()*(B_.X()*d[1][0] + B_.Y()*d[1][1] + B_.Z()*d[1][2]) + 
	            B.Z()*(B_.X()*d[2][0] + B_.Y()*d[2][1] + B_.Z()*d[2][2]);

      prob = (D0*A*A_ + dbb0 + omega*(A*F.Vect()*B_.Vect() + A_*(G.Vect()*B)))*BPS;

      if(prob > prob_max) cout<<"Error: prob = "<<prob<<" > max = "<<prob_max<<endl;
    }while(prob_max*y > prob);           

    // tau- in the LAB frame  
     e_tau = Ptaum.E();                                       
     p_tau = Ptaum.P();       
    th_tau = Ptaum.Theta(); 
    ph_tau = Ptaum.Phi(); 
    
    // lep- in the LAB frame 
     e_lep = Plep.E();                                       
     p_lep = Plep.P();       
    th_lep = Plep.Theta(); 
    ph_lep = Plep.Phi(); 
    
    // pi+ and pi0 in the LAB frame 
     e_pip = Ppip_lab.E();                                       
     p_pip = Ppip_lab.P();       
    th_pip = Ppip_lab.Theta(); 
    ph_pip = Ppip_lab.Phi(); 
    
     e_pi0 = Ppi0_lab.E();                                       
     p_pi0 = Ppi0_lab.P();       
    th_pi0 = Ppi0_lab.Theta(); 
    ph_pi0 = Ppi0_lab.Phi(); 

    trout->Fill();
  }

  cout<<"Saving tree to the "<<outfile<<endl; 

  fout->cd();
  trout->Write("",TObject::kOverwrite);
  fout->Close(); 

  cout<<"done"<<endl; 
  cout<<"---------------------------------------"<<endl; 

  gROOT->ProcessLine(".! date"); 
  return; 
}
