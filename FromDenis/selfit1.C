#define selfit1_cxx
// The class definition in selfit1.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("selfit1.C")
// Root > T->Process("selfit1.C","some options")
// Root > T->Process("selfit1.C+")
//

#include "selfit1.h"

void selfit1::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   //----------------------------------------
   
   Ntotal = 0;
   
   //----------------------------------------
   
   L0 = new TArrayD(Nevents); 
   L1 = new TArrayD(Nevents); 
   L2 = new TArrayD(Nevents); 
   L3 = new TArrayD(Nevents); 
   L4 = new TArrayD(Nevents); 

   cout<<"---------------------------------"<<endl;
   
   if(flr == 1){
     cout<<"Calculating normalisation"<<endl; 
     AA0 = 0.0;
     AA1 = 0.0; 
     AA2 = 0.0; 
     AA3 = 0.0; 
     AA4 = 0.0; 
   }else{
     
     TString fnorm = "norm_"; 
     fnorm +=  iomega; 
     fnorm += ".dat";  

     cout<<"Reading normalisation from "<<fnorm<<endl; 

     ifstream indatafile;
     indatafile.open(fnorm,ios::in); 
     indatafile>>AA0;
     indatafile>>AA1; 
     indatafile>>AA2;
     indatafile>>AA3;
     indatafile>>AA4;
     indatafile.close();

     cout<<"AA0 = "<<AA0<<endl;
     cout<<"AA1 = "<<AA1<<endl;
     cout<<"AA2 = "<<AA2<<endl; 
     cout<<"AA3 = "<<AA3<<endl;
     cout<<"AA4 = "<<AA4<<endl; 
     
     /*
     if(ml < 0.05){
       // electron
       cout<<"tau -> e nu nu"<<endl;
       AA0 = 1.0;
       AA1 = 0.0; 
       AA2 = 4*ml/mt; 
       AA3 = 0.0; 
       AA4 = 0.0; 
     }else{
       // muon 
       cout<<"tau -> mu nu nu"<<endl;
       AA0 = 1.0;
       AA1 = 0.0; 
       AA2 = 4*ml/mt; 	  
       AA3 = 0.0; 
       AA4 = 0.0;
     }
     */
     cout<<"              done               "<<endl;
   }

   cout<<"---------------------------------"<<endl;
}

void selfit1::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t selfit1::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either selfit1::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   GetEntry(entry); 

   Ptaum.SetPxPyPzE( p_tau*sin(th_tau)*cos(ph_tau),  p_tau*sin(th_tau)*sin(ph_tau),  p_tau*cos(th_tau), e_tau); 
   Ptaup.SetPxPyPzE(-p_tau*sin(th_tau)*cos(ph_tau), -p_tau*sin(th_tau)*sin(ph_tau), -p_tau*cos(th_tau), e_tau); 

   Rx.SetXYZ(-cos(th_tau)*cos(ph_tau), -cos(th_tau)*sin(ph_tau), sin(th_tau)); 
   Ry.SetXYZ(             sin(ph_tau),             -cos(ph_tau),           0); 
   Rz.SetXYZ( sin(th_tau)*cos(ph_tau),  sin(th_tau)*sin(ph_tau), cos(th_tau)); 

   Boost_taum.SetXYZ(0, 0, Beta_tau); 
   Boost_taup.SetXYZ(0, 0,-Beta_tau); 

   // lep- 
   Plep_0.SetPxPyPzE(0,0,1,0); 
   Plep_0.SetRho(p_lep);
   Plep_0.SetTheta(th_lep); 
   Plep_0.SetPhi(ph_lep);
   Plep_0.SetE(e_lep);
   Plep.SetPxPyPzE(Rx*Plep_0.Vect(), Ry*Plep_0.Vect(), Rz*Plep_0.Vect(), Plep_0.E()); 
   Plep.Boost(-Boost_taum);    // boost to the tau- rest frame
   Nlep = Plep.Vect().Unit(); 

   // pi+
   Ppip_0.SetPxPyPzE(0,0,1,0); 
   Ppip_0.SetE(e_pip);
   Ppip_0.SetRho(p_pip);
   Ppip_0.SetTheta(th_pip); 
   Ppip_0.SetPhi(ph_pip);
   Ppip.SetPxPyPzE(Rx*Ppip_0.Vect(), Ry*Ppip_0.Vect(), Rz*Ppip_0.Vect(), Ppip_0.E()); 
   Ppip.Boost(-Boost_taup);    // boost to the tau+ rest frame 

   // pi0 
   Ppi0_0.SetPxPyPzE(0,0,1,0); 
   Ppi0_0.SetE(e_pi0);
   Ppi0_0.SetRho(p_pi0);
   Ppi0_0.SetTheta(th_pi0); 
   Ppi0_0.SetPhi(ph_pi0);
   Ppi0.SetPxPyPzE(Rx*Ppi0_0.Vect(), Ry*Ppi0_0.Vect(), Rz*Ppi0_0.Vect(), Ppi0_0.E()); 
   Ppi0.Boost(-Boost_taup);    // boost to the tau+ rest frame
  
   Prho = Ppip + Ppi0;         // in the tau+ rest frame

   // nu_tau_bar in the tau+ rest frame
   Pnu.SetPxPyPzE(-Prho.Px(), -Prho.Py(), -Prho.Pz(), M_tau-Prho.E()); 

   //-----------------------------------------------

   double d0 = 1 + cos(th_tau)*cos(th_tau) + sin(th_tau)*sin(th_tau)/Gamma_tau/Gamma_tau;

   //   TMatrixD d(3,3);
   /*
   double d[3][3]; 
   d[0][0] = (1.0 + 1.0/Gamma_tau/Gamma_tau)*sin(th_tau)*sin(th_tau); 
   d[0][1] = 0.0;
   d[0][2] = sin(2*th_tau)/Gamma_tau;                                                          
   d[1][0] = 0.0; 
   d[1][1] = -Beta_tau*Beta_tau*sin(th_tau)*sin(th_tau); 
   d[1][2] = 0.0;
   d[2][0] = sin(2*th_tau)/Gamma_tau;
   d[2][1] = 0.0;
   d[2][2] = 1 + cos(th_tau)*cos(th_tau) - sin(th_tau)*sin(th_tau)/Gamma_tau/Gamma_tau; 
   */ 

   TLorentzVector p1_taum(Ebeam*beta*sin(th_tau), 0.0, Gamma_tau*Ebeam*(beta*cos(th_tau) - Beta_tau), Gamma_tau*Ebeam*(1 - Beta_tau*beta*cos(th_tau)));
   TLorentzVector ptaup_taum(0.0, 0.0,-2*Ebeam*Gamma_tau*Beta_tau, Gamma_tau*Ebeam*(1 + Beta_tau*Beta_tau)); 
   TLorentzVector G = (1.0/Gamma_tau/Ebeam)*(2*p1_taum - ptaup_taum); 

   /*
   TLorentzVector p1_taup(Ebeam*beta*sin(th_tau), 0.0, Gamma_tau*Ebeam*(beta*cos(th_tau) + Beta_tau), Gamma_tau*Ebeam*(1 + Beta_tau*beta*cos(th_tau))); 
   TLorentzVector ptaum_taup(0.0, 0.0, 2*Ebeam*Gamma_tau*Beta_tau, Gamma_tau*Ebeam*(1 + Beta_tau*Beta_tau)); 
   TLorentzVector F = (1.0/Gamma_tau/Ebeam)*(2*p1_taup - ptaum_taup); 
   */ 

   // -------------------------------------------
   // lep-  

   double x1 = Plep.E()/E_max;
   if(x1 < x0 || x1 > 1) cout<<"Error: x1 = "<<x1<<endl;

   double A0 = sqrt(x1*x1 - x0*x0)*3.0*x1*(1 - x1);
   double A1 = sqrt(x1*x1 - x0*x0)*(2.0/3.0)*(4*x1*x1 - 3*x1 - x0*x0);
   double A2 = sqrt(x1*x1 - x0*x0)*3.0*x0*(1 - x1);
  
   double B1 = (x1*x1 - x0*x0)*(1 - x1); 
   double B2 = (x1*x1 - x0*x0)*(2.0/3.0)*(4*x1 - 4 + sqrt(1 - x0*x0)); 

   //-------------------------------------------------
   // rho+ 

   /*
   double mpp2      = Prho.M2(); 
   double mpp       = sqrt(mpp2); 
   double M_rho2    = M_rho*M_rho; 
   double M_rhop2   = M_rhop*M_rhop;
   double ppi_mpp   = sqrt((mpp2-(mpi+mpi0)*(mpi+mpi0))*(mpp2-(mpi-mpi0)*(mpi-mpi0)))/2.0/mpp; 
   double ppi_mrho  = sqrt((M_rho2-(mpi+mpi0)*(mpi+mpi0))*(M_rho2-(mpi-mpi0)*(mpi-mpi0)))/2.0/M_rho; 
   double ppi_mrhop = sqrt((M_rhop2-(mpi+mpi0)*(mpi+mpi0))*(M_rhop2-(mpi-mpi0)*(mpi-mpi0)))/2.0/M_rhop; 
   
   double G_rho     = Gamma_rho*(M_rho/mpp)*(ppi_mpp/ppi_mrho)*(ppi_mpp/ppi_mrho)*(ppi_mpp/ppi_mrho); 
   double G_rhop    = Gamma_rhop*(M_rhop/mpp)*(ppi_mpp/ppi_mrhop)*(ppi_mpp/ppi_mrhop)*(ppi_mpp/ppi_mrhop); 
   TComplex BW1     = TComplex(M_rho2, 0)/TComplex(M_rho2-mpp2,-M_rho*G_rho); 
   TComplex BW2     = TComplex(M_rhop2,0)/TComplex(M_rhop2-mpp2,-M_rhop*G_rhop); 
   TComplex BWtot   = (BW1 + betaBW*BW2)/(1.0 + betaBW); 
   double BPS       = BWtot.Rho2()*(Prho.P()/M_tau)*(ppi_mpp/mpp);

   Q = Ppip - Ppi0; 
   double B1_ = (Q*Q) + 2*(Pnu*Q);
   double B2_ = (Q*Q) - 2*(Pnu*Q);
   double A_  = 2*(Pnu*Q)*Q.E() - (Q*Q)*Pnu.E(); 
   B_ = B1_*Ppip.Vect() + B2_*Ppi0.Vect(); 
   */ 

   //------------------------------------------
   
   /*
   TMatrixD nlb(3,3);
   nlb[0][0] = Nlep.X()*B_.X();  nlb[0][1] = Nlep.X()*B_.Y(); nlb[0][2] = Nlep.X()*B_.Z(); 
   nlb[1][0] = Nlep.Y()*B_.X();  nlb[1][1] = Nlep.Y()*B_.Y(); nlb[1][2] = Nlep.Y()*B_.Z(); 
   nlb[2][0] = Nlep.Z()*B_.X();  nlb[2][1] = Nlep.Z()*B_.Y(); nlb[2][2] = Nlep.Z()*B_.Z(); 
   TMatrixD dnlb = d*nlb; 
   double  dnlb0 = dnlb[0][0]+dnlb[1][1]+dnlb[2][2]; 
   */ 
   /*
   double  dnlb0 = Nlep.X()*(B_.X()*d[0][0] + B_.Y()*d[0][1] + B_.Z()*d[0][2]) + 
                   Nlep.Y()*(B_.X()*d[1][0] + B_.Y()*d[1][1] + B_.Z()*d[1][2]) + 
                   Nlep.Z()*(B_.X()*d[2][0] + B_.Y()*d[2][1] + B_.Z()*d[2][2]); 
   */ 

   //-------------------------------------------------------------------------   

   double D0 = d0*A0;                     // 1 
   double D1 = d0*A1;                     // rho 
   double D2 = d0*A2;                     // eta 
   double D3 = omega*(G.Vect()*Nlep)*B1;  // ksi 
   double D4 = omega*(G.Vect()*Nlep)*B2;  // ksi*del 

   //------------------------------------------------------------------------- 

   if(D0<0) cout<<"Error: D0 = "<<D0<<endl; 

   L0->AddAt(D0,Ntotal); 
   L1->AddAt(D1,Ntotal); 
   L2->AddAt(D2,Ntotal); 
   L3->AddAt(D3,Ntotal);
   L4->AddAt(D4,Ntotal);

   Ntotal++; 

   if(Ntotal%100000 == 0) cout<<Ntotal<<endl;
   
   return kTRUE;
}

void selfit1::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void selfit1::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  cout<<"Total number of analyzed events is: "<<Ntotal<<endl; 

  double dAA0 = 0; 
  double dAA1 = 0; 
  double dAA2 = 0; 
  double dAA3 = 0; 
  double dAA4 = 0; 
   
  if(flr == 1){
    for(int i = 0; i < Ntotal; i++){
      double DD0 = (*L0)[i];
      double DD1 = (*L1)[i]; 
      double DD2 = (*L2)[i]; 
      double DD3 = (*L3)[i]; 
      double DD4 = (*L4)[i]; 
      double DD  = DD0 + DD1*rho_sm + DD2*eta_sm + DD3*xi_sm + DD4*xidel_sm;
      if(DD > 0){
	AA0 += DD0/DD;
	AA1 += DD1/DD;
	AA2 += DD2/DD;
	AA3 += DD3/DD;
	AA4 += DD4/DD;      
      }else{
	cout<<i<<" event DD = "<<DD<<endl;
	cout<<DD0<<" : "<<DD1<<" : "<<DD2<<" : "<<DD3<<" : "<<DD4<<endl; 
      }
    }

    AA0 = AA0/Ntotal;
    AA1 = AA1/Ntotal; 
    AA2 = AA2/Ntotal; 
    AA3 = AA3/Ntotal; 
    AA4 = AA4/Ntotal; 

    for(int i = 0; i < Ntotal; i++){
      double DD0 = (*L0)[i];
      double DD1 = (*L1)[i]; 
      double DD2 = (*L2)[i]; 
      double DD3 = (*L3)[i]; 
      double DD4 = (*L4)[i]; 
      double DD  = DD0 + DD1*rho_sm + DD2*eta_sm + DD3*xi_sm + DD4*xidel_sm;
      if(DD > 0){
	dAA0 += (DD0/DD - AA0)*(DD0/DD - AA0);
	dAA1 += (DD1/DD - AA1)*(DD1/DD - AA1);
	dAA2 += (DD2/DD - AA2)*(DD2/DD - AA2);
	dAA3 += (DD3/DD - AA3)*(DD3/DD - AA3);
	dAA4 += (DD4/DD - AA4)*(DD4/DD - AA4);      
      }else{
	cout<<i<<" event DD = "<<DD<<endl;
	cout<<DD0<<" : "<<DD1<<" : "<<DD2<<" : "<<DD3<<" : "<<DD4<<endl; 
      }
    }      

    dAA0 = sqrt(dAA0/(Ntotal-1)/Ntotal);
    dAA1 = sqrt(dAA1/(Ntotal-1)/Ntotal);
    dAA2 = sqrt(dAA2/(Ntotal-1)/Ntotal);
    dAA3 = sqrt(dAA3/(Ntotal-1)/Ntotal);
    dAA4 = sqrt(dAA4/(Ntotal-1)/Ntotal);
    
    cout<<"Calculated normalisation: "<<endl; 

    cout<<"AA0 = "<<AA0<<" +- "<<dAA0<<endl;
    cout<<"AA1 = "<<AA1<<" +- "<<dAA1<<endl;
    cout<<"AA2 = "<<AA2<<" +- "<<dAA2<<endl; 
    cout<<"AA3 = "<<AA3<<" +- "<<dAA3<<endl;
    cout<<"AA4 = "<<AA4<<" +- "<<dAA4<<endl; 

    TString fnorm = "norm_"; 
    fnorm +=  iomega; 
    fnorm += ".dat";  
    
    ofstream outdatafile;
    outdatafile.open(fnorm,ios::out); 
    outdatafile<<AA0<<endl;
    outdatafile<<AA1<<endl; 
    outdatafile<<AA2<<endl;
    outdatafile<<AA3<<endl;
    outdatafile<<AA4<<endl;
    outdatafile.close(); 

    cout<<"Write normalisation to "<<fnorm<<endl; 
  }
}
