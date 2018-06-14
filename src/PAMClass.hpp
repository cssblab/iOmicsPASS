 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   PAMClass.hpp
 * Author: Hiromi WL Koh
 *
 * Last modified on 13 June, 2018, 4:00 PM
 */

#ifndef PAMCLASS_HPP
#define PAMCLASS_HPP


#include "globals.hpp"

#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

class PAMClass {
private:
    string NodeA, NodeB;
    string dataT;
    vector<double> d_ik, xbar_ik, mk;
    vector<double> dik_star;
    vector<string> classLab;
    double s_o,s_i, xbar_i;
    vector<int> GRPsize;       

public:
    PAMClass();
    PAMClass(const PAMClass& orig);
    virtual ~PAMClass();
    
    PAMClass(string nodeA, string nodeB, string datatype);
    void addDik(double d);
    void addxik(double mu);
    void addMk(double mk_k);
    void adds0(double s0);
    void addDikstar(double dik);
    void addsi(vector<double> vec_x, vector<int> grp_size);
    void addgrpsize(vector<int> grp_size);
    void addxbar(double x);
    void addLabel(string lab);
    int getGrpsize(int K);
    vector<string> getNodes();
    string getDataType(){return dataT;}
    double getXbari(){return xbar_i;}
    double getSo(){return s_o;}
    double getSi(){return s_i;}
    double getMk(int K){return mk[K];}
    double getDik(int K){return d_ik[K];}
    double getXbarik(int K){return xbar_ik[K];}
    double calNewXbarik(double thres, int group);     
    double calNewDik(double thres, int group);
    vector<double> getDIKstar() {return dik_star;}
    vector<double> getDIK() {return d_ik;}
    vector<double> getXBARIK(){return xbar_ik;}
    vector<double> getMK(){return mk;}
    vector<int> getGRPSIZE(){return GRPsize;} 
    vector<string> getCLASSLAB(){return classLab;}
};

#endif /* PAMCLASS_HPP */

