/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   PAMClass.cpp
 * Author: Hiromi WL Koh
 * 
 * Last modified on 13 June, 2018, 4:00 PM
 */

#include "globals.hpp"

#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

PAMClass::PAMClass() {
}

PAMClass::PAMClass(const PAMClass& orig) {
}

PAMClass::~PAMClass() {
}

PAMClass::PAMClass(string nodeA, string nodeB, string datatype){
    NodeA = nodeA;
    NodeB = nodeB;
    dataT = datatype;
}

vector<string> PAMClass::getNodes(){
   vector<string> ret;
   ret.push_back(NodeA);
   ret.push_back(NodeB);
   return ret;
}

void PAMClass::addDik(double d){
    d_ik.push_back(d);
}

void PAMClass::addxik(double mu){
    xbar_ik.push_back(mu);
}
void PAMClass::adds0(double s0){
    s_o = s0;
}
void PAMClass::addMk(double mk_k){
   mk.push_back(mk_k);
}
void PAMClass::addDikstar(double dik){
   dik_star.push_back(dik);
}

void PAMClass::addsi(vector<double> vec_x, vector<int> grp_size){
    double sigma=0;
    int startpos =0;
    int endpos =0;
    for(int i=0; i<grp_size.size(); i++){
	//GRPsize.push_back(grp_size.at(i));
	if(i==0) endpos =0;
	else startpos = endpos;
	endpos = startpos+grp_size.at(i);
	for(int j=startpos; j<endpos; j++){
	    sigma+= pow((vec_x.at(j) - xbar_ik[i]),2);        
	}
    }
    s_i = sqrt(sigma/ ((double)(vec_x.size() - grp_size.size())));

}

int PAMClass::getGrpsize(int K){
     int ret;
     ret = GRPsize.at(K);
     return ret;
}

void PAMClass::addgrpsize(vector<int> grp_size){
   for(int i=0; i<grp_size.size(); i++) GRPsize.push_back(grp_size.at(i));
}

void PAMClass::addxbar(double x){
    xbar_i = x;
}

void PAMClass::addLabel(string lab){
   classLab.push_back(lab);
}

double PAMClass::calNewXbarik(double thres, int group){
   double sign_dij =1;
   double dd;
   double d = dik_star[group];
   if(d<0) sign_dij= -1;
   dd = abs(d) - thres;
   if(dd<0) dd=0;
   double dij_new = sign_dij *dd;
   double b= mk[group]*(s_o+s_i)*dij_new;
   double xik_new = xbar_i + b;
  
   //if(b==0) cerr<<"0!"<<endl;
   return xik_new;
}

double PAMClass::calNewDik(double thres, int group){
   double sign_dij =1.0;
   double dd, dij_new;
   double d = dik_star[group];
   if(d<0) sign_dij= -1.0;
   dd = abs(d)- thres;
   if(dd<=0) dij_new =0;
   else dij_new = sign_dij *dd;
   
   return dij_new;
}
