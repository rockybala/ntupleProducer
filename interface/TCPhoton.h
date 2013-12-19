#ifndef _TCPHOTON_H
#define	_TCPHOTON_H

#include <memory>
#include "TObject.h"
#include "TLorentzVector.h"
#include "TArrayF.h"
#include "TCEGamma.h"
#include <vector>

using namespace std;

class TCPhoton : public TCEGamma {
private:

  // ID variables
  //float _e2OverE9;
  bool  _trackVeto;
  
  bool    _convVeto;

  //mip stuff                                                                                          
  float _mipchi2;
  float _miptoten;
  float _mipslope;
  float _mipintercept;
  float _mipnhitcone;
  float _mipishalo;

  float _roundness;
  float _angle;
  float _smin;
  float _smaj;


 public:
    TCPhoton();
    virtual ~TCPhoton();

    // "get" methods -----------

    //float E2OverE9() const; 
    bool  TrackVeto() const;

    bool  ConversionVeto() const;
    float MipChi2() const;                                                                               
    float MipTotEn() const;
    float MipSlope() const;
    float MipIntercept() const;
    float MipNHitCone() const;
    float MipIsHalo() const;

    float Roundness() const;
    float Angle() const;
    float SMin() const;
    float SMaj() const;


    // "set" methods ---------
    //void SetE2OverE9(float);
    void SetTrackVeto(bool);
    void SetConversionVeto(bool);

    void SetMipChi2(float);
    void SetMipTotEn(float);
    void SetMipSlope(float);
    void SetMipIntercept(float);
    void SetMipNHitCone(float);
    void SetMipIsHalo(float);

    void SetRoundness(float);
    void SetAngle(float);
    void SetSMin(float);
    void SetSMaj(float);

    ClassDef(TCPhoton, 1);
};

#endif


