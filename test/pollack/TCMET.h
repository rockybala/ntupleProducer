//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Sep  2 09:05:57 2013 by ROOT version 5.32/00)
//      from the StreamerInfo in file nuTuple.root
//////////////////////////////////////////////////////////


#ifndef TCMET_h
#define TCMET_h
class TCMET;

#include "TVector2.h"

class TCMET : public TVector2 {

public:
// Nested classes declaration.

public:
// Data Members.
   TVector2    _genMET;     //
   float       _sumEt;      //
   float       _muonFraction;    //
   float       _neutralHadronFraction;    //
   float       _neutralEMFraction;        //
   float       _chargedHadronFraction;    //
   float       _chargedEMFraction;        //
   float       _unCorPhi;                 //
   float       _Significance;             //
   float       _SigmaX2;                  //

   TCMET();
   TCMET(const TCMET & );
   virtual ~TCMET();

   ClassDef(TCMET,2); // Generated by MakeProject.
};
#endif
