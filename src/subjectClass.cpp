/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   subjectClass.cpp
 * Author: Hiromi WL Koh
 * 
 * Last modified on 13 June, 2018, 4:00 PM
 *
 */


#include "globals.hpp"
#include "geneClass.hpp"

#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

subjectClass::subjectClass() {
}

subjectClass::subjectClass(const subjectClass& orig) {
}

subjectClass::~subjectClass() {
}

subjectClass::subjectClass(string id){
    subject = id;
}

set<string> subjectClass::getCommonG(int datatype1, int datatype2){
     set<string> setA,setB;
     if(datatype1==1) setA = dna;
     if(datatype1==2) setA = rna;
     if(datatype1==3) setA = prot;
     if(datatype2==1) setB = dna;
     if(datatype2==2) setB = rna;
     if(datatype2==3) setB = prot;
     set<string> Iset = IntersectSet(setA,setB);
     return Iset;
}

void subjectClass::removeNAgenes(string g, int datatype){
   set<string>::iterator it;
   if(datatype==1){
	for(it =dna_missing.begin(); it!=dna_missing.end(); it++) if(*it==g) dna_missing.erase(it);
   }
   if(datatype==2){
	for(it = rna_missing.begin(); it!=rna_missing.end(); it++) if(*it==g) rna_missing.erase(it);
   }
   if(datatype==3){
	for(it = prot_missing.begin(); it!=prot_missing.end(); it++) if(*it==g) prot_missing.erase(it);
   }
}

void subjectClass::insertNAgenes(string g, int datatype){
  
     if(datatype==1) dna_missing.insert(g);
     if(datatype==2) rna_missing.insert(g);
     if(datatype==3) prot_missing.insert(g);
}

set<string> subjectClass::getMissingG(int datatype){
     set<string> ret;
     if(datatype==1) ret =dna_missing;
     if(datatype==2) ret = rna_missing;
     if(datatype==3) ret = prot_missing;
     return ret;
}


void subjectClass::insertMap(string g, double val ,int dataType){
    
    geneClass *curG = NULL;
    map<string,geneClass>::iterator m_iter, m_iter2;
    m_iter = geneMap.find(g);
    
    if(m_iter==geneMap.end()){
        curG = new geneClass(g);
        geneMap[g]= *curG;
        delete curG;
        m_iter2 = geneMap.find(g);
        if(dataType==1 & val!= NA_VAL) {
            m_iter2->second.addDNA(val);
            dna.insert(g);
	    dna_exist = true;
        }
        if(dataType==2 & val!= NA_VAL) {
            m_iter2->second.addRNA(val);
            rna.insert(g);
	    rna_exist = true;
        }
        if(dataType==3 & val!= NA_VAL) {
            m_iter2->second.addProt(val);
            prot.insert(g);
	    prot_exist = true;
        }
    }
    else{
        if(dataType==1 & val!= NA_VAL) {
            m_iter->second.addDNA(val);
            dna.insert(g);
	    dna_exist = true;
        }
        if(dataType==2 & val!= NA_VAL) {
            m_iter->second.addRNA(val);
            rna.insert(g);
	    rna_exist = true;
        }
        if(dataType==3 & val!= NA_VAL) {
            m_iter->second.addProt(val);
            prot.insert(g);
	    prot_exist = true;
        }
    }
}


set<string> subjectClass::getSet(int dataType){
    set<string> ret;
    if(dataType==1) ret = dna;
    if(dataType==2) ret = rna;
    if(dataType==3) ret = prot;
    return ret;
}

bool subjectClass::existSet(int dataType){
    bool ret;
    if(dataType==1) ret = dna_exist;
    if(dataType==2) ret = rna_exist;
    if(dataType==3) ret = prot_exist;
    return ret;
}



double subjectClass::extractRatio(string gene_num, string gene_den, string type){
    
    double ret;
    map<string, geneClass>::iterator g_iter, g_iter2;
    int count =0;
    g_iter = geneMap.find(gene_num);
    g_iter2 = geneMap.find(gene_den);
    if(g_iter!=geneMap.end() & g_iter2!=geneMap.end()) {
  	  double num, den, ret;
	    if(type=="PPI"){
	      num = g_iter->second.getReadings(3);
	      den = g_iter2->second.getReadings(3);
	      ret = num+den;
	      if(num==NA_VAL|den==NA_VAL) ret = NA_VAL;
	    }
	    if(type=="TF"){
	      num = g_iter->second.getReadings(3);
	      double d1 = g_iter2->second.getReadings(2);
	      double d2 = g_iter2->second.getReadings(1);
	      den = d1 - d2;
	      ret = num+den;
	      if(num==NA_VAL|den==NA_VAL) ret = NA_VAL;
	     }
	    count++;
    }
    else ret = NA_VAL;
    return ret;
}

