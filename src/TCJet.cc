/* 
 * File:   TCJet.cc
 * Author: Anton A.
 * 
 * Created on April 30, 2010, 2:49 PM
 */

#include "../interface/TCJet.h"
#include "TCJetLinkDef.h"


TCJet::TCJet() {
}

TCJet::~TCJet() {
}

// "get" methods -------------------------------------

using namespace std;

float TCJet::ChHadFrac() const {
   return _chHadFrac;
}

float TCJet::NeuHadFrac() const {
   return _neuHadFrac;
}

float TCJet::ChEmFrac() const {
   return _chEmFrac;
}

float TCJet::NeuEmFrac() const {
   return _neuEmFrac;
}

unsigned int TCJet::NumConstit() const {
   return _numConstit;
}

unsigned int TCJet::NumChPart() const {
   return _numChPart;
}

float TCJet::VtxSumPtFrac() const {
   return _vtxSumPtFrac;
}

float TCJet::VtxSumPt() const {
   return _vtxSumPt;
}  

float TCJet::VtxTrackFrac() const {
   return _vtxTrackFrac;
}  

int TCJet::VtxNTracks() const {
   return _vtxNTracks;
}  

unsigned int TCJet::VtxSumPtIndex() const {
   return _vtxSumPtIndex;
}

unsigned int TCJet::VtxCountIndex() const {
   return _vtxCountIndex;
}

float TCJet::UncertaintyJES() const {
	return _jesUncertainty;
}

// b tagging discriminators
float TCJet::BDiscriminatorMap(string key) {
    return _bDiscrMap[key];
}

// jet flavor
int TCJet::JetFlavor() const {
    return _jetFlavor;
}

// Hgg style Jet Id vars
float TCJet::BetaStarClassic() const {
  return _betaStarClassic;
}
float TCJet::DR2Mean() const {
  return _dR2Mean;
}

int TCJet::PuJetIdFlag_cut_loose() const{
  return _id_cut_loose;
}

int TCJet::PuJetIdFlag_cut_medium() const{
  return _id_cut_medium;
}

int TCJet::PuJetIdFlag_cut_tight() const{
  return _id_cut_tight;
}

int TCJet::PuJetIdFlag_mva_loose() const{
  return _id_mva_loose;
}

int TCJet::PuJetIdFlag_mva_medium() const{
  return _id_mva_medium;
}

int TCJet::PuJetIdFlag_mva_tight() const{
  return _id_mva_tight;
}

// "set" methods ---------------------------------------------


void TCJet::SetVtxSumPtFrac(float f){
   _vtxSumPtFrac = f;
}  

void TCJet::SetVtxSumPt(float p){
   _vtxSumPt = p;
}  

void TCJet::SetVtxTrackFrac(float f){
   _vtxTrackFrac = f;
}  

void TCJet::SetVtxNTracks(int n){
   _vtxNTracks = n;
}  

void TCJet::SetVtxSumPtIndex(unsigned int i){
   _vtxSumPtIndex = i;
} 

void TCJet::SetVtxCountIndex(unsigned int i){
   _vtxCountIndex = i;
} 

void TCJet::SetChHadFrac(float c) {
   _chHadFrac = c;
}

void TCJet::SetNeuHadFrac(float n) {
   _neuHadFrac = n;
}

void TCJet::SetChEmFrac(float c) {
   _chEmFrac = c;
}

void TCJet::SetNeuEmFrac(float n) {
   _neuEmFrac = n;
}

void TCJet::SetNumConstit(unsigned int n) {
   _numConstit = n;
}

void TCJet::SetNumChPart(unsigned int n) {
   _numChPart = n;
}

void TCJet::SetUncertaintyJES(float u) {
	_jesUncertainty = u;
}
// b tagging discriminators

void TCJet::SetBDiscriminatorMap(string key, float value) {
   _bDiscrMap[key] = value;
}

void TCJet::SetJetFlavor(float f) {
    _jetFlavor = f;
}

// Hgg style Jet Id vars
void TCJet::SetBetaStarClassic(float b) {
  _betaStarClassic = b;
}
void TCJet::SetDR2Mean(float d) {
  _dR2Mean = d;
}

void TCJet::SetPuJetIdFlag_cut_loose(int i){
  _id_cut_loose = i;
}

void TCJet::SetPuJetIdFlag_cut_medium(int i){
  _id_cut_medium= i;
}

void TCJet::SetPuJetIdFlag_cut_tight(int i){
  _id_cut_tight= i;
}

void TCJet::SetPuJetIdFlag_MVA_loose(int i){
  _id_mva_loose= i;
}

void TCJet::SetPuJetIdFlag_MVA_medium(int i){
  _id_mva_medium= i;
}

void TCJet::SetPuJetIdFlag_MVA_tight(int i){
  _id_mva_tight= i;
}
