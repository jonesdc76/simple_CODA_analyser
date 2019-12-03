// Analysis of data in crate "mycrate" and bank "mybank".
// See the parameters "mycrate" and "mybank" below.
// R. Michaels, Oct 2016
// Don Jones-- Dec 2019 Added argument to analyze starting at event number
#define MAXLEN  200000
#define MAXROC    50
#define MAXBANK   20
#define MAXRAW  6000
#define MAXDAT   100

#include <iostream>
#include <string>
#include <vector>
#include "THaCodaFile.h"
#include "THaEtClient.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom.h"
#include <vector>
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
//#include "./UTEvent_minimal.h"

using namespace std;

void usage();
void decode(int* data);
void clear();

// Global data 

Int_t *irn, *rocpos, *roclen;
Int_t *bankpos, *banklen;
Int_t myroc,mybank,mychan;
Int_t evlen, evtype, evnum;
Int_t rdata[MAXDAT];

int main(int argc, char* argv[])
{
  //  Int_t load_status = gSystem->Load("./UTEvent_minimal_h.so");
  //  cout << "Status: " << load_status << endl;


  THaCodaData *coda;      

  Int_t maxevent = 1000, start_event = 0;

  Int_t choice1 = 1;  /* source of data */

// Choice of ROC and BANK

// myroc=9, mybank=6  for grinchadc data
// myroc=12, mybank=3  for MPD data

  myroc=12;
  mybank=3;

  irn = new Int_t[MAXROC];
  rocpos = new Int_t[MAXROC];
  roclen = new Int_t[MAXROC];
  bankpos = new Int_t[MAXROC*MAXBANK];
  banklen = new Int_t[MAXROC*MAXBANK];

  if (argc > 1) maxevent = atoi(argv[1]);

  if (argc > 2) start_event = atoi(argv[2]);

  if (argc > 3) mychan = atoi(argv[3]);
 
  cout << "Events to process "<<maxevent<<endl;

  // Initialize root and output
  TROOT evanana("anaroot","Hall A data analysi");
  TFile hfile("ana.root","RECREATE","HA data");

  TH1F *h1 = new TH1F("h1","Example histogram",100,0,10000);

  if (choice1 == 1) 
    {  // CODA File
      
      // CODA file "run.dat" may be a link to CODA datafile on disk
      TString filename("run.dat");
      
      coda = new THaCodaFile();
      if (coda->codaOpen(filename) != 0) 
	{
	  cout << "ERROR:  Cannot open CODA data" << endl;
	  goto end1;
	}
    } 
  else 
    {         // Online ET connection
      
      int mymode = 1;
      TString mycomputer("adaq1");
      TString mysession("Compton");
      
      coda = new THaEtClient();
      if (coda->codaOpen(mycomputer, mysession, mymode) != 0) 
	{
	  cout << "ERROR:  Cannot open ET connection" << endl;
	  goto end1;
	}
      
    }

  // Loop over events
  int status, evlo, evhi;
  
  evlo = start_event;
  evhi = evlo + maxevent;

  for (int iev = evlo; iev < evhi; iev++) 
    {//the loop over the event bounds
      
      if (iev > 0 && ((iev%1000)==0) ) printf("%d events\n",iev);

      clear();

      status = coda->codaRead();  

      if (status != 0) {  // EOF or no data
	  
	  if ( status == -1)  {
	      if (choice1 == 1) {
		  cout << "End of CODA file. Bye bye." << endl;
		  evhi=iev;
	      }
	      if (choice1 == 2) cout << "CODA/ET not running. Bye bye." << endl;
	  } else {
	      cout << "ERROR: codaRread status = " << hex << status << endl;
	  }
	  goto end1;
  

      }  else  {   // have data ...

	  decode( coda->getEvBuffer() );

          h1->Fill(rdata[4]);  // semi-meaningless histogram, as an example


	}
    }
  
 end1:
  
  coda->codaClose();

  hfile.Write();
  hfile.Close();

  return 0;

}; //end of main function

void usage() {  
  cout << "Make a link to run.dat. That is your data "<<endl;
  cout << "Usage:  'xana [maxevents] [chan]' " << endl;
}; 

void clear() {
  memset(rdata,0,MAXDAT*sizeof(Int_t));
}

void decode (int* data) {

  evlen = data[0] + 1;
  evtype = data[1]>>16;
  evnum = data[4];
  static int dodump = 1;  // dump the raw data
  static int debug = 1;   // debug the decoding
  static int verbose = 1;  // some printout

  if (dodump) 
    { // Event dump
      cout << "\n\n Event number " << dec << evnum;
      cout << " length " << evlen << " type " << evtype << endl;
      int ipt = 0;
      for (int j = 0; j < (evlen/5); j++) 
	{
	  cout << dec << "\n evbuffer[" << ipt << "] = ";
	  for (int k=j; k<j+5; k++) 
	    {
	      cout << hex << data[ipt++] << " ";
	    }
	  cout << endl;
	}
      if (ipt < evlen) 
	{
	  cout << dec << "\n evbuffer[" << ipt << "] = ";
	  for (int k=ipt; k<evlen; k++) 
	    {
	      cout << hex << data[ipt++] << " ";
	    }
	  cout << endl;
	}
    }

  // Decoding for "normal" events

  if (evtype < 15) {

      // First find pointers to ROCs.  Useful for multi-crate analysis.
      
      // n1 = pointer to first word of ROC
      int pos = data[2]+3;
      int nroc = 0;
      while( pos+1 < evlen && nroc < MAXROC ) {
	  int len  = data[pos];
	  int iroc = (data[pos+1]&0xff0000)>>16;
	  if(iroc>=MAXROC || nroc >= MAXROC) 
	    {
	      cout << "Decoder error:  ";
	      cout << "  ROC num " <<dec<<iroc;
	      cout << "  nroc "<<nroc<<endl;
	      return;
	    }
	  // Save position and length of each found ROC data block
	  rocpos[iroc]  = pos;
	  roclen[iroc]  = len;
	  irn[nroc++] = iroc;
	  pos += len+1;
	}
      Int_t found = 0;
      for (int j=0; j<nroc; j++) if(myroc == irn[j]) found=1;
      if (!found) {
	  cout << "ERROR:  ROC "<<dec<<myroc<<" not in datastream !!"<<endl;
	  return;
      }

      if (debug) {
	  cout << "Roc info "<<nroc<<endl;
	  for (int j=0; j<nroc; j++) 
	    {
	      Int_t iroc = irn[j];
	      cout << "irn "<<dec<<iroc<<"   pos "<<rocpos[iroc];
	      cout <<"   len "<<roclen[iroc]<<endl;
	    }
      }

      // Go over the data in myroc
      // Find "mybank"

      pos = rocpos[myroc]+2;
      Int_t bankmax = -1;
 
      while (pos < rocpos[myroc]+roclen[myroc]) {

       Int_t blen = data[pos];  // length of bank, this includes the header
       Int_t bankhead = data[pos+1]; // header of form 0x70100 for bank 7
       Int_t banknum = (bankhead>>16)&0xffff; // bank number
       if (banknum > bankmax) bankmax=banknum;
       if (debug) cout << "bank:  pos "<<dec<<pos<<"   word 0x"<<hex<<data[pos]<<dec<<endl;
       if (debug) cout << "roc "<<dec<<myroc<<"   bank "<<banknum<<"  head 0x"<<hex<<bankhead<<"    len "<<dec<<blen<<endl;
       if (banknum >= 0 && banknum < MAXBANK) {
 	  Int_t idx = myroc*MAXBANK + banknum;
          if (idx < MAXBANK*MAXROC) {
	      bankpos[idx] = pos+2;
              banklen[idx] = blen-1;

	  } else {
	      cout << "Decoder failure.  Fire the author."<<endl;
	  }
       }
       pos += blen+1;  
  
      }

      for (int j=0; j<nroc; j++) {

	if(verbose) cout << "Roc number  "<<irn[j]<< "====================== "<<bankmax << endl;
     
        for (int bank=0; bank <= bankmax; bank++) {
  
    	  Int_t idx = irn[j]*MAXBANK + bank;

          if (banklen[idx] > 0 && banklen[idx] < MAXLEN) {

	    if (verbose) cout << " %%%%  bank number "<<bank<<"     %%%%%%%%%%%%%%%%%%% "<<endl;

           if (bank == mybank) {

// LOOP over the Data of mybank in myroc
  	      for (pos = bankpos[idx]; pos<bankpos[idx]+banklen[idx]; pos++) {
		if (verbose) cout << " data if mybank   data["<<dec<<pos<<"] = 0x"<<hex<<data[pos]<<dec<<endl;
                Int_t ix = pos-bankpos[idx];
                
                if (ix >= 0 && ix < MAXDAT) {
		    rdata[ix] = (data[pos]&0xfff);
                    if (verbose) cout << "rdata[ "<<ix<< " ] =  "<<rdata[ix]<<endl;
		}

	      }
	   }
	  }
	}
      }


  }
}


