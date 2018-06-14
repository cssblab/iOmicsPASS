 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   geneClass.hpp
 * Author: Hiromi WL Koh
 *
 * Last modified on 13 June, 2018, 4:00 PM
 */

#ifndef GENECLASS_HPP
#define GENECLASS_HPP


#include "globals.hpp"

#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;


class geneClass {
private:
    string gene;
    bool dna_bool = false;
    bool rna_bool = false;
    bool prot_bool = false;
    double dna, rna, prot;

public:
    geneClass();
    geneClass(const geneClass& orig);
    virtual ~geneClass();
    
    geneClass(string id);
    void addDNA(double d);
    void addRNA(double r);
    void addProt(double p);
    bool ratioAvail(int datatype);
    double getReadings(int datatype);
    double getProb(int type);

};

#endif /* GENECLASS_HPP */

