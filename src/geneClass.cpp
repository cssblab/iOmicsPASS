/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   geneClass.cpp
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

geneClass::geneClass() {
}

geneClass::geneClass(const geneClass& orig) {
}

geneClass::~geneClass() {
}

geneClass::geneClass(string id){
    gene = id;
}

void geneClass::addProt(double p){
    prot = p;
    prot_bool = true;
}

void geneClass::addRNA(double r){
    rna = r;
    rna_bool = true;
}

void geneClass::addDNA(double d){
    dna = d;
    dna_bool = true;
}

bool geneClass::ratioAvail(int datatype){
  bool ret = false;
	if(datatype==1) ret = dna_bool;
	if(datatype==2) ret = rna_bool;
	if(datatype==3) ret =prot_bool;
	return ret;
}

double geneClass::getReadings(int datatype){
    double ret;
    if(datatype==1) {
	if(dna_bool) ret = dna;
	else ret=NA_VAL;
    } 
    if(datatype==2) {
	if(rna_bool) ret = rna;
	else ret = NA_VAL;
    }
    if(datatype==3) {
	if(prot_bool) ret = prot;
	else ret = NA_VAL;
    }
    return ret;
}

