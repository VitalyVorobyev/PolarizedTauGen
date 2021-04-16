//--------------------------------------------------------- 
// 
// flr = 0 : ordinary fit 
// 
// flr = 1 : calculation of the normalization coefficients, 
//           statistically independent sample should be 
//           used to calculate normalization correctly  
// 
//--------------------------------------------------------- 

#include "selfit.h" 

void fit(int, int); 
void fcn(int &, double *, double &, double *, int); 

selfit *sel = (selfit *)TSelector::GetSelector("selfit.C+");
const int NPAR = 4; 
int Nminuit; 
int flag_print_like; // print likelihood values during the fit: 1:YES, 0:NO 

void fit(int flr, int iomega){
  gROOT->ProcessLine(".! date");

  if(!(flr == 0 || flr == 1))        return; 
  if(!(iomega >= 0 && iomega <= 10)) return; 
  
  flag_print_like = 1; 

  //-------------------------------
  // open data file with tree 
  
  TString infile = "out/gen_"; 
  infile +=     flr; 
  infile +=     "_"; 
  infile +=  iomega; 
  infile += ".root";  

  cout<<"Data file: "<<infile<<endl; 

  TFile *f = new TFile(infile); 
  TTree *treein = (TTree *)f->Get("tr");

  int nentries = treein->GetEntries(); 

  sel->SetFLG(flr,iomega,nentries); 

  treein->Process(sel);

  if(flr == 1) return; // calculation of normalization ONLY 

  cout<<"         FITTING            "<<endl; 
  
  //-------------------------------------------------------
  // Use TMinuit to make fit

  TMinuit *Minu = new TMinuit(NPAR);  //initialize TMinuit with a maximum of 4 params
  Minu->SetFCN(fcn);
  double arglist[10];

  int ierflg = 0;
  arglist[0] = 0.5;         // For -Log(P) likelihood fit  
  Minu->mnexcm("SET ERR", arglist, 1, ierflg);
  
  //-------------------------------------------------------
  // Define fit parameters 
  
  static double michel[NPAR] = {0.75, 0.00, 1.00, 0.75};
  static double   step[NPAR] = {0.02, 0.02, 0.02, 0.02};
  
  Minu->mnparm(0,   "rho", michel[0], step[0], 0.25, 1.25, ierflg);
  Minu->mnparm(1,   "eta", michel[1], step[1], -0.5, 0.5, ierflg);
  Minu->mnparm(2,   "ksi", michel[2], step[2],  0.5, 1.5, ierflg);
  Minu->mnparm(3,"ksidel", michel[3], step[3], 0.25, 1.25, ierflg);

  //*****************************************************************
  // Minimization
  
  Nminuit    = 0;
  arglist[0] = 500;                            // MAX number of minimization calls
  arglist[1] = 0.1;                            // tolerance
  Minu->mnexcm("MIGRAD", arglist, 2, ierflg);  // MINIMIZATION
  Minu->mnexcm("HESSE" ,       0, 0, ierflg);  // 
  Minu->mnexcm("MINOS" ,       0, 0, ierflg);  // correct evaluation of errors 

  //*****************************************************************
  // Print results
  
  TString chnam;
  int iuint;
  double val[NPAR],err[NPAR],xlow,xup; 
  double eplus[NPAR],eminus[NPAR],eparab[NPAR],corrco;
  double fmin,fedm,errdef;
  int nvpar,nparx,icstat;

  //-----------------------------------------------------------------
  // Also save results in the file 
  
  TString outfile = "par_";
  outfile +=  iomega; 
  outfile +=  ".dat";  

  cout<<"Saving data to the "<<outfile<<" file"<<endl; 
  ofstream datafit;
  datafit.open(outfile,ios::out); 
  
  for(int i=0; i<NPAR; i++){
    Minu->mnpout(i,chnam,val[i],err[i],xlow,xup,iuint);
    Minu->mnerrs(i,eplus[i],eminus[i],eparab[i],corrco);
    cout<<chnam<<" : "<<val[i]<<" + "<<eplus[i]<<" - "<<eminus[i]<<" +- "<<eparab[i]<<" corr: "<<corrco<<endl;
    datafit<<setw(12)<<setprecision(8)<<setiosflags(ios::fixed | ios::showpoint | ios::right)<<val[i]
	   <<setw(12)<<setprecision(8)<<setiosflags(ios::fixed | ios::showpoint | ios::right)<<eplus[i]
	   <<setw(12)<<setprecision(8)<<setiosflags(ios::fixed | ios::showpoint | ios::right)<<eminus[i] 
	   <<setw(12)<<setprecision(8)<<setiosflags(ios::fixed | ios::showpoint | ios::right)<<eparab[i]
	   <<setw(12)<<setprecision(8)<<setiosflags(ios::fixed | ios::showpoint | ios::right)<<corrco<<endl; 
  }
  Minu->mnstat(fmin,fedm,errdef,nvpar,nparx,icstat);
  cout<<"After "<<Nminuit<<" steps L min: "<<fmin<<endl;
  datafit.close(); 

  //-------------------------------------
  
  gROOT->ProcessLine(".! date");
  return; 
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag){
  Nminuit++;
  sel->SetParameters(npar,par);
  f = sel->GetLikelihood();
  if(flag_print_like == 1) cout<<Nminuit<<" : like "<<f<<endl;
}
