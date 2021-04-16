//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 21 05:23:20 2014 by ROOT version 5.16/00
// from TTree tr/total
// found on file: combine.root
//////////////////////////////////////////////////////////

#ifndef selfit_h
#define selfit_h

#include<iostream> 
#include<iomanip> 
#include<fstream> 
#include<stdlib.h> 
#include<stdio.h> 
#include<math.h> 

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "TString.h"
#include "TMath.h"
#include "TTree.h"
#include "TComplex.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TMinuit.h" 



class selfit : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
 
   //----------------------
   // Your variables 

   int Ntotal;
   int nparam;

   int flr; 
   int iomega; 
   int Nevents; 

   double pi;
   double param[4];

   double Like;

   double AA0;
   double AA1;  
   double AA2;
   double AA3; 
   double AA4; 

   // Particle's parameters

   double Ebeam;  
   double omega;
   double M_e; 
   double M_mu; 
   double M_tau;
   double M_lep;
   double mpi;
   double mpi0;
   double Gamma_tau; 
   double P_tau; 
   double Beta_tau; 
   double P_e; 
   double beta;  
   double E_max;
   double x0;

   double    rho_sm; 
   double    eta_sm;
   double     xi_sm;  
   double  xidel_sm; 
   double  xirho_sm;
   
   // rho, rho' parameters

   double M_rho;
   double M_rhop;
   double Gamma_rho;
   double Gamma_rhop;
   double betaBW;

   TArrayD *L0;   
   TArrayD *L1;  // rho 
   TArrayD *L2;  // eta
   TArrayD *L3;  // ksi 
   TArrayD *L4;  // ksidel
 
   TLorentzVector Ptaum; 
   TLorentzVector Ptaup;
   TLorentzVector Plep; 
   TLorentzVector Plep_0;
   TLorentzVector Prho; 
   TLorentzVector Pnu; 
   TLorentzVector Ppip; 
   TLorentzVector Ppi0; 
   TLorentzVector Ppip_0; 
   TLorentzVector Ppi0_0; 
   TLorentzVector Q;
   
   TVector3 Boost_taum;
   TVector3 Boost_taup;
   TVector3 Nlep;
   TVector3 Rx;
   TVector3 Ry;
   TVector3 Rz;
   TVector3 B_;

   // Declaration of leave types
   Double_t         e_tau;
   Double_t         p_tau;
   Double_t        th_tau;
   Double_t        ph_tau;
   Double_t         e_lep;
   Double_t         p_lep;
   Double_t        th_lep;
   Double_t        ph_lep;
   Double_t         e_pip;
   Double_t         p_pip;
   Double_t        th_pip;
   Double_t        ph_pip;
   Double_t         e_pi0;
   Double_t         p_pi0;
   Double_t        th_pi0;
   Double_t        ph_pi0;

   // List of branches
   TBranch         *b_e_tau;   //!
   TBranch         *b_p_tau;   //!
   TBranch        *b_th_tau;   //!
   TBranch        *b_ph_tau;   //!
   TBranch         *b_e_lep;   //!
   TBranch         *b_p_lep;   //!
   TBranch        *b_th_lep;   //!
   TBranch        *b_ph_lep;   //!
   TBranch         *b_e_pip;   //!
   TBranch         *b_p_pip;   //!
   TBranch        *b_th_pip;   //!
   TBranch        *b_ph_pip;   //!
   TBranch         *b_e_pi0;   //!
   TBranch         *b_p_pi0;   //!
   TBranch        *b_th_pi0;   //!
   TBranch        *b_ph_pi0;   //!

   selfit(TTree * /*tree*/ =0) : 
   fChain(0),
   L0(0),
   L1(0),
   L2(0),
   L3(0), 
   L4(0)
   { 
     Ntotal     =  0;
     pi         =  TMath::Pi(); 
     
     nparam     =  4;
     param[0]   =  0.0;
     param[1]   =  0.0;
     param[2]   =  0.0;
     param[3]   =  0.0; 
     Like       =  0.0;
     
     Ebeam      =  2.33; 

     M_e        =  0.000511;
     M_mu       =  0.105658;
     M_tau      =  1.77699;  
     M_lep      =  M_mu; 
     mpi        =  0.13957;
     mpi0       =  0.13498;

     Gamma_tau  =  Ebeam/M_tau; 
     P_tau      =  sqrt(Ebeam*Ebeam - M_tau*M_tau);
     Beta_tau   =  P_tau/Ebeam; 
     P_e        =  sqrt(Ebeam*Ebeam - M_e*M_e);
     beta       =  P_e/Ebeam; 
     E_max      =  0.5*M_tau*(1.0 + M_lep*M_lep/M_tau/M_tau);
     x0         =  M_lep/E_max;

     M_rho      =  0.7743; 
     Gamma_rho  =  0.1491; 
     M_rhop     =  1.370; 
     Gamma_rhop =  0.386; 
     betaBW     = -0.091; 

     rho_sm     =  0.75; 
     eta_sm     =  0.00;
     xi_sm      =  1.00;  
     xidel_sm   =  0.75; 
     xirho_sm   =  1.00; 
   }
   virtual ~selfit(){ 
     if (L0) delete L0;
     if (L1) delete L1;
     if (L2) delete L2;
     if (L3) delete L3;
     if (L4) delete L4; 
   }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   //--------------------------------------------- 

   void SetParameters(int npar, double xpar[]){
     nparam = 4; 
     for(int j=0; j<nparam; j++){ 
       param[j]=xpar[j]; 
     }  
   }

   void SetFLG(int xflr, int xiomega, int xevents){ 
     flr = xflr; 
     iomega = xiomega; 
     omega  = iomega*0.1; 
     Nevents = xevents; 
   } 

   double GetLikelihood(){
     Like = 0.0;
     for(int i=0; i<Ntotal; i++){
       double DD0  = (*L0)[i]; 
       double DD1  = (*L1)[i]; 
       double DD2  = (*L2)[i]; 
       double DD3  = (*L3)[i]; 
       double DD4  = (*L4)[i]; 
       double DD   = DD0 + param[0]*DD1 + param[1]*DD2 + param[2]*DD3 + param[3]*DD4; 
       double Norm = AA0 + param[0]*AA1 + param[1]*AA2 + param[2]*AA3 + param[3]*AA4;

       if(DD > 0.0 && Norm > 0.0){
	 Like += -log(DD/Norm);
       }else{
	 cout<<"-------------------"<<endl; 
	 cout<<"DD0 = "<<DD0<<endl; 
	 cout<<"DD1 = "<<DD1<<endl; 
	 cout<<"DD2 = "<<DD2<<endl; 
	 cout<<"DD3 = "<<DD3<<endl; 
	 cout<<"DD4 = "<<DD4<<endl; 
	 cout<<"p0 = "<<param[0]<<" p1 = "<<param[1]<<" p2 = "<<param[2]<<" p3 = "<<param[3]<<endl; 
	 cout<<i<<": DD = "<<DD<<" Norm = "<<Norm<<endl;
	 Like += 0.0;
       }
     }
     return Like;
   }
   
   //-------------------------------------

   ClassDef(selfit,0);
};

#endif

#ifdef selfit_cxx
void selfit::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("e_tau", &e_tau, &b_e_tau);
   fChain->SetBranchAddress("p_tau", &p_tau, &b_p_tau);
   fChain->SetBranchAddress("th_tau", &th_tau, &b_th_tau);
   fChain->SetBranchAddress("ph_tau", &ph_tau, &b_ph_tau);
   fChain->SetBranchAddress("e_lep", &e_lep, &b_e_lep);
   fChain->SetBranchAddress("p_lep", &p_lep, &b_p_lep);
   fChain->SetBranchAddress("th_lep", &th_lep, &b_th_lep);
   fChain->SetBranchAddress("ph_lep", &ph_lep, &b_ph_lep);
   fChain->SetBranchAddress("e_pip", &e_pip, &b_e_pip);
   fChain->SetBranchAddress("p_pip", &p_pip, &b_p_pip);
   fChain->SetBranchAddress("th_pip", &th_pip, &b_th_pip);
   fChain->SetBranchAddress("ph_pip", &ph_pip, &b_ph_pip);
   fChain->SetBranchAddress("e_pi0", &e_pi0, &b_e_pi0);
   fChain->SetBranchAddress("p_pi0", &p_pi0, &b_p_pi0);
   fChain->SetBranchAddress("th_pi0", &th_pi0, &b_th_pi0);
   fChain->SetBranchAddress("ph_pi0", &ph_pi0, &b_ph_pi0);
}

Bool_t selfit::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef selfit_cxx
