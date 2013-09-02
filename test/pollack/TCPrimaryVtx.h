//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Sep  2 09:05:57 2013 by ROOT version 5.32/00)
//      from the StreamerInfo in file nuTuple.root
//////////////////////////////////////////////////////////


#ifndef TCPrimaryVtx_h
#define TCPrimaryVtx_h
class TCPrimaryVtx;

#include "TVector3.h"

class TCPrimaryVtx : public TVector3 {

public:
// Nested classes declaration.

public:
// Data Members.
   float       _nDof;       //
   float       _chi2;       //
   bool        _isFake;     //
   int         _nTracks;    //
   float       _sumPt2Trks;    //

   TCPrimaryVtx();
   TCPrimaryVtx(const TCPrimaryVtx & );
   virtual ~TCPrimaryVtx();

   ClassDef(TCPrimaryVtx,2); // Generated by MakeProject.
};
#endif
