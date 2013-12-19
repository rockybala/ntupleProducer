#include "../interface/TCPhoton.h"
#include "TCPhotonLinkDef.h"
#include <iostream>

TCPhoton::TCPhoton() {
  _roundness = -99;
  _angle = -99;
  _smin = -99;
  _smaj = -99;
  _mipchi2 = -99;
  _miptoten = -99;
  _mipslope = -99;
  _mipintercept = -99;
  _mipnhitcone = -99;
  _mipishalo = -99;

}

TCPhoton::~TCPhoton() {}

// "get" methods -------------------------------------
//float TCPhoton::E2OverE9() const { return _e2OverE9; } 
bool  TCPhoton::TrackVeto() const { return _trackVeto; }

bool  TCPhoton::ConversionVeto() const { return _convVeto; }


//mip stuff
float TCPhoton::MipChi2() const {
  return _mipchi2;
}

float TCPhoton::MipTotEn() const {
  return _miptoten;
}

float TCPhoton::MipSlope() const {
  return _mipslope;
}

float TCPhoton::MipIntercept() const{
  return _mipintercept;
}

float TCPhoton::MipNHitCone() const{
  return _mipnhitcone;
}

float TCPhoton::MipIsHalo() const{
  return  _mipishalo;
}

float TCPhoton::Roundness() const{
  return _roundness;
}

float TCPhoton::Angle() const{
  return _angle;
}

float TCPhoton::SMin() const{
  return _smin;
}

float TCPhoton::SMaj() const{
  return _smaj;
}


// "set" methods ---------------------------------------------

//void TCPhoton::SetE2OverE9(float e) { _e2OverE9 = e; } 
void TCPhoton::SetTrackVeto(bool t) { _trackVeto = t; } 
void TCPhoton::SetConversionVeto(bool v) { _convVeto = v; }

//mip stuff                                                                                         
void TCPhoton::SetMipChi2(float mchi){
  _mipchi2 = mchi;
}

void TCPhoton::SetMipTotEn(float toten){
  _miptoten = toten;
}

void TCPhoton::SetMipSlope(float slope){
  _mipslope = slope;
}

void TCPhoton::SetMipIntercept(float intercept){
  _mipintercept = intercept;
}

void TCPhoton::SetMipNHitCone(float cone){
  _mipnhitcone = cone;
}

void TCPhoton::SetMipIsHalo(float halo){
  _mipishalo = halo;
}

void TCPhoton::SetRoundness(float round){
  _roundness = round;
}

void TCPhoton::SetAngle(float ang){
  _angle = ang;
}

void TCPhoton::SetSMin(float sm){
  _smin = sm;
}

void TCPhoton::SetSMaj(float sm2){
  _smaj = sm2;
}
