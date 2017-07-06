
#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <fstream>





#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <cmath>
#include <cstring>

#include <unistd.h>

#include "NRGclasses.hpp"
#include "NRGfunctions.hpp"
#include "NRGOpMatRules.hpp"
#include "NRG_main.hpp"
#include "TwoChQS.hpp"



#ifndef pi
#define pi 3.141592653589793238462643383279502884197169
#endif

#ifndef _CHIN_
#define _CHIN_


double chiN(int Nsites, double Lambda)
{

  double daux[3];
  daux[0]=1.0-pow( Lambda,-((double)(Nsites)+1.0) );
  daux[1]=sqrt( 1.0-pow(Lambda,-(2.0*(double)(Nsites)+1.0)) );
  daux[2]=sqrt( 1.0-pow(Lambda,-(2.0*(double)(Nsites)+3.0)) );  

  return(daux[0]/(daux[1]*daux[2]));

}

#endif


int main (){


// Parameters for command-line passing (GetOpt)
  
  CNRGCodeHandler ThisCode;

  //#include"ModelOptMain.cpp"


  ThisCode.SaveData=false;

  // NRG objects
  
  CNRGarray Aeig;
  CNRGarray* pAeig;
  
  CNRGbasisarray AeigCut;
  CNRGbasisarray Abasis;
  CNRGbasisarray SingleSite;

  // STL vector

  CNRGmatrix HN;
  CNRGmatrix Qm1fNQ;

  //CNRGmatrix MQQp1;

  CNRGmatrix* MatArray;
  int NumNRGarrays=4;
  // Jul 09: Will this work???
  vector<CNRGmatrix> STLMatArray;
  CNRGmatrix auxNRGMat;
  // MatArray 0 is f_ch1
  // MatArray 1 is f_ch2
  // MatArray 2 is Sz
  // MatArray 3 is Sz2



  // STL vectors
  //vector <double> Params;
  vector <double> ParamsHN;
  vector <double> ParamsBetabar;
  vector<int> CommonQNs; 
  vector<int> totSpos;

  double  x[4] = {2, .5, 0.5,1};
  std::vector<double> Params(x, x + sizeof x / sizeof x[0]);
  cout << "firstvalue is " << Params[0] << '\n';
  OneChQ_SetAnderson_Hm1_old( Params, pAeig, MatArray);
  

  // Thermodynamics

  CNRGthermo Suscep;
  CNRGthermo Entropy;

  CNRGthermo *ThermoArray;
  int NumThermoArrays=2;

  double TM=0.0;
  //double DN=0.0;
  //double betabar=0.727;
  ThisCode.Nsites=0;
  ThisCode.betabar=0.727;
  // Add more than one betabar in the code... Done!
  //ThisCode.betabar=0.6; 

  // instream
  ifstream InFile;
  // outstream
  ofstream OutFile;
  
  //cout << "firstvalue is " << Params << '\n'

  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  ///                                              ///
  ///                Main code                     ///
  ///                                              ///
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////
  //OneChQ_SetAnderson_Hm1_old( Params, pAeig, MatArray);

  

}
// end main


//
// Trash
//
//////////////////////////






void OneChQ_SetAnderson_Hm1_old(vector<double> Params,
				CNRGarray* pAeig, 
				CNRGmatrix* NRGMats){
 // Set initial CNRG array (N=-1)

  double U=Params[0];
  double ed=Params[1];
  double Lambda=Params[2];
  double HalfLambdaFactor=Params[3];

  
  // First Test: the usual Anderson model

  pAeig->NQNumbers=1;
  pAeig->Nshell=-1;

  // new thing
  pAeig->totalS=false;
  /*

  pAeig->QNumbers.push_back(-1.0);
  pAeig->iDegen.push_back(0); // Sz=0  


  // Two states in this one
  pAeig->QNumbers.push_back(0.0);
  pAeig->iDegen.push_back(5); // Sz=1/2
  pAeig->iDegen.push_back(-5); // Sz=-1/2


  pAeig->QNumbers.push_back(1.0);
  pAeig->iDegen.push_back(0); // Sz=0

  pAeig->BlockBegEnd.push_back(0);pAeig->BlockBegEnd.push_back(0);
  pAeig->BlockBegEnd.push_back(1);pAeig->BlockBegEnd.push_back(2);
  pAeig->BlockBegEnd.push_back(3);pAeig->BlockBegEnd.push_back(3);

  // using version in Wilson's: the former pluse U/2 divided by Factor
  // Assuming D=1
  pAeig->dEn.push_back(0.5*U/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed+0.5*U)/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((ed+0.5*U)/(Lambda*HalfLambdaFactor));
  pAeig->dEn.push_back((2.0*ed+1.5*U)/(Lambda*HalfLambdaFactor));


  // Eigenvectors
  pAeig->dEigVec.push_back(1.0);

  pAeig->dEigVec.push_back(1.0);
  pAeig->dEigVec.push_back(0.0);
  pAeig->dEigVec.push_back(0.0);
  pAeig->dEigVec.push_back(1.0);


  pAeig->dEigVec.push_back(1.0);

  //pAeig->PrintEn();

  // Set Matrix elements

  // fN_up

  NRGMats[0].SyncNRGarray(*pAeig);
  NRGMats[0].UpperTriangular=false;

 
  NRGMats[0].MatEl.push_back(1.0);
  NRGMats[0].MatEl.push_back(0.0);
  NRGMats[0].MatBlockMap.push_back(0);
  NRGMats[0].MatBlockMap.push_back(1);
  NRGMats[0].MatBlockBegEnd.push_back(0);
  NRGMats[0].MatBlockBegEnd.push_back(1);

  NRGMats[0].MatEl.push_back(0.0);
  NRGMats[0].MatEl.push_back(1.0);
  NRGMats[0].MatBlockMap.push_back(1);
  NRGMats[0].MatBlockMap.push_back(2);
  NRGMats[0].MatBlockBegEnd.push_back(2);
  NRGMats[0].MatBlockBegEnd.push_back(3);


  // fN_dn
  NRGMats[1].SyncNRGarray(*pAeig);

  NRGMats[1].UpperTriangular=false;
 
  NRGMats[1].MatEl.push_back(0.0);
  NRGMats[1].MatEl.push_back(1.0);
  NRGMats[1].MatBlockMap.push_back(0);
  NRGMats[1].MatBlockMap.push_back(1);
  NRGMats[1].MatBlockBegEnd.push_back(0);
  NRGMats[1].MatBlockBegEnd.push_back(1);

  NRGMats[1].MatEl.push_back(-1.0);
  NRGMats[1].MatEl.push_back(0.0);
  NRGMats[1].MatBlockMap.push_back(1);
  NRGMats[1].MatBlockMap.push_back(2);
  NRGMats[1].MatBlockBegEnd.push_back(2);
  NRGMats[1].MatBlockBegEnd.push_back(3);


  cout << " OLD: Setting Sz... " << endl;

  NRGMats[2].SyncNRGarray(*pAeig);
  NRGMats[3].SyncNRGarray(*pAeig);
  /*
  int i1=0;
  for (int ibl=0; ibl<pAeig->NumBlocks(); ibl++)
    {
      int ist0=pAeig->GetBlockLimit(ibl,0);
      int ist1=pAeig->GetBlockLimit(ibl,1);

      double Qi=pAeig->GetQNumber(ibl,0);

      cout << " Q = " << Qi
	   << " ist0 = " << ist0
	   << " ist1 = " << ist1 << endl;

      // Block diagonal
      NRGMats[2].MatBlockMap.push_back(ibl);
      NRGMats[2].MatBlockMap.push_back(ibl);
      NRGMats[3].MatBlockMap.push_back(ibl);
      NRGMats[3].MatBlockMap.push_back(ibl);
      // Position in MatBlock follows ist
      NRGMats[2].MatBlockBegEnd.push_back(i1);
      NRGMats[3].MatBlockBegEnd.push_back(i1);
      for (int ist=ist0;ist<=ist1;ist++)
	{
	  double Szi=(double)(pAeig->iDegen[ist])/10.0;
	  for (int jst=ist0;jst<=ist1;jst++)
	    {
	      // Sz and Sz2 diagonal only (waste of space but...)
	      if (ist==jst)
		{
		  NRGMats[2].MatEl.push_back(Szi);
		  NRGMats[3].MatEl.push_back(Szi*Szi);
		  // Test in Sep 08: calc ndot instead WORKS!
// 		  NRGMats[2].MatEl.push_back(Qi+1.0);
// 		  NRGMats[3].MatEl.push_back(Qi+1.0);
		}
	      else
		{
		  NRGMats[2].MatEl.push_back(0.0);
		  NRGMats[3].MatEl.push_back(0.0);
		}
	      i1++;
	    } // loop in jst
	} //ist
      NRGMats[2].MatBlockBegEnd.push_back(i1-1);
      NRGMats[3].MatBlockBegEnd.push_back(i1-1);

    }
  // end Set-up of Sz,Sz2

*/
  //NRGMats[2].PrintAllBlocks();
  //NRGMats[3].PrintAllBlocks();

  // Suscep test

  /*
  CNRGbasisarray AeigCut=CutStates(pAeig, 1000);

  vector<double> ParamsBetabar;
  double betabar=0.727;
  ParamsBetabar.clear();
  ParamsBetabar.push_back(betabar);
  
  double dSz0=CalcOpAvg(ParamsBetabar,&AeigCut,&NRGMats[2],false,0);
  double dSz20=CalcOpAvg(ParamsBetabar,&AeigCut,&NRGMats[3],false,0);
  
  cout << " N=-1: Sz = " << dSz0 << " Sz2 = " << dSz20 
       << "    T_M chi = " << dSz20-dSz0*dSz0 << endl;
  */
}
/////////////////////////////////////////
/////////////////////////////////////////
/////////////////////////////////////////

///////
