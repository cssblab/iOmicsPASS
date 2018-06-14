/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   subjectClass.hpp
 * Author: Hiromi WL Koh
 *
 * Last modified on 13 June, 2018, 4:00 PM
 */

#ifndef SUBJECTCLASS_HPP
#define SUBJECTCLASS_HPP

#include "globals.hpp"
#include "geneClass.hpp"

#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;


class subjectClass {
    
private:
    string subject;
    bool dna_exist=false;
    bool rna_exist=false;
    bool prot_exist=false;
    set<string> dna;
    set<string> rna;
    set<string> prot;
    set<string> dna_missing, rna_missing, prot_missing;

public:
    subjectClass();
    subjectClass(const subjectClass& orig);
    virtual ~subjectClass();
    
    subjectClass(string id);
    void removeNAgenes(string g, int datatype);
    double extractRatio(string gene_num, string gene_den, string type);
    void insertMap(string g, double val ,int dataType);
    void insertNAgenes(string g, int datatype);
    set<string> getMissingG(int datatype);
    set<string> getSet(int dataType);
    set<string> getCommonG(int datatype1, int datatype2);
    bool existSet(int type);
    map<string, geneClass> geneMap;
};

#endif /* SUBJECTCLASS_HPP */

