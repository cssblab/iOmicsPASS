/*
 *          Copyright Joe Coder 2004 - 2006.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/* 
 * File:   globals.cpp
 * Author: Hiromi WL Koh
 * 
 * Created on 13 June, 2018, 4:00 PM
 */

#include "globals.hpp"
#include "geneClass.hpp"
#include "subjectClass.hpp"
#include "PAMClass.hpp"

#include <boost/math/distributions/hypergeometric.hpp>
#include <algorithm>
#include <random>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <regex>
#include <limits>
#include <iomanip>

using namespace std;
using namespace boost::math;

// Function splits a string on character
vector<string> stringSplit(string str, char delim) {
    vector<string> ret;
    string buf = "";
    int i = 0;
    while((unsigned) i < str.length()) {
        if(str[i] != delim) buf += str[i];
        else if(buf.length() > 0) {
            ret.push_back(buf);
            buf = "";
        }
        i++;
    }
    if(!buf.empty()) ret.push_back(buf);
    return(ret);
}

string concatenate(vector<string> vec, char delim){
    string ret=vec.at(0);
    for(int j=1; j<vec.size(); j++){
        ret += delim + vec.at(j);
    }
    return ret;
}

bool GREP(string STRING, string PATTERN){
   bool ret = false;
   string start = STRING.substr(0,PATTERN.length());
   regex rx(PATTERN);
   if(regex_match( start,rx)) ret =true;
   return ret;   
}


//trim from start 
static inline string &ltrim(string &s){
	s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int,int>(isspace))));
        return s;
}
//trim from end
static inline string &rtrim(string &s){
	s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(),s.end());
        return s;
}
//trim from both ends
static inline string &trim(string &s){
	return ltrim(rtrim(s));
}


bool allNA(vector<double> vec){
    bool ret = true;
    for(int i=0; i<vec.size(); i++){
	if(vec.at(i)!=NA_VAL) ret=false;
    }
    return ret;
}

bool NAexist(vector<double> vec){
   bool ret = false;
   for(int i=0; i <vec.size(); i++){
     if(vec.at(i)==NA_VAL){ ret =true; continue;}
   }
   return ret;
}

vector<string> IntersectVec(vector<string> v1, vector<string> v2){

      vector<string> v3;
      sort(v1.begin(), v1.end());
      sort(v2.begin(), v2.end());
    
      set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v3));
      return v3;
}

set<string> IntersectSet(set<string> v1, set<string> v2){
    
    set<string> ret;
    set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), 
	             inserter(ret,ret.begin()));
    //if(v1.empty()) ret = v2;
    //if(v2.empty()) ret = v1;
    return ret;
}
set<string> UnionSet(set<string> v1, set<string> v2){
    set<string> ret;
    set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), 
	      inserter(ret,ret.begin()));
    return ret;
}

void UniqueVec(vector<string> &v){
    sort(v.begin(), v.end());
    v.erase(unique(v.begin(), v.end()),v.end());
}

bool contains(string g, set<string> *myset){

   bool ret = false;
   set<string>::iterator it;
   it = myset->find( g);
   if(it!= myset->end()) ret = true;
   return ret;
}
 
bool containsVec(string g, vector<string> *vec){
   bool ret = false;
   for(int i =0; i< vec->size();i++){
	if(vec->at(i)==g) {
	   ret=true;
	   break;
	}
   }
   return ret;
}



bool readUserInput(string file ,Input_t &UserInput){
    bool ret = true;
    ifstream File;
    File.open(file.c_str(), ios::in);

    if(!File.is_open()){
    	cerr<<"Error! Unable to open file...\n";
	    exit(1);
    }
    string line, d;
    vector<string> v;
 
    while(!File.eof()){
	    getline(File,line);
	    string start = line.substr(0,1);
      regex expr("#");
	    if(regex_match(start,expr)) continue;

      if(GREP(line,"MAX_BLOCKSIZE")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.blocksize = stoi(d);
	    }
	    if(GREP(line,"DNA_FILE")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.dnafile = d;
      }
	    if(GREP(line, "RNA_FILE")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.rnafile = d;
	    }
	    if(GREP(line,"PROTEIN_FILE")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.proteinfile = d;
      }
	    if(GREP(line,"PPI_NETWORK")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.ppinet = d;
  	  }
	    if(GREP(line,"TF_NETWORK")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.tfnet = d;
	    }
	    if(GREP(line,"SUBTYPE_FILE")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.subtypefile = d;
	    }
      if(GREP(line,"MODULE_FILE")){
	      v = stringSplit(line,'=');
        d = trim(v.at(1));
	      UserInput.pathwayfile = d;
	    }
      if(GREP(line,"MIN_OBS")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.minobs = stoi(d);
   	  }
	    if(GREP(line,"MIN_PROP")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.minprop = stod(d);
      }
	    if(GREP(line, "INTERACT_SIGN")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      bool d_new;
	      if(d=="false") d_new = false;
	      if(d=="true") d_new = true;
	      UserInput.Interact = d_new;
	    }	

	    if(GREP(line,"LOG_TRANSFORM_DNA")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      bool d_new;
	      if(d == "false") d_new = false;
	      if(d =="true") d_new = true;
	      UserInput.log_DNA = d_new;
	      cerr<<"log DNA: "<<d<<endl;
	    }
	    if(GREP(line,"LOG_TRANSFORM_RNA")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      bool d_new;
	      if(d == "false") d_new = false;
	      if(d =="true") d_new = true;
	      UserInput.log_RNA = d_new;
	      cerr<<"log RNA: "<<d<<endl;
	    }
    	if(GREP(line,"LOG_TRANSFORM_PROT")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   bool d_new;
    	   if(d == "false") d_new = false;
    	   if(d =="true") d_new = true;
    	   UserInput.log_PROT = d_new;
    	   cerr<<"log prot: "<<d<<endl;
    	}
    	if(GREP(line,"KNN_IMPUTE")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   bool d_new;
    	   if(d == "false") d_new = false;
    	   if(d =="true") d_new = true;
    	   UserInput.knnimpute = d_new;
    	}
    	if(GREP(line,"KNN_K")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   UserInput.knnK = stoi(d);
    	}
    	if(GREP(line,"CROSS_VALIDATION")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   bool d_new;
    	   if(d =="false") d_new = false;
    	   if(d =="true") d_new = true;
    	   UserInput.crossV = d_new;
    	}
            if(GREP(line,"CV_FOLD")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   UserInput.numFold = stoi(d);
    	}
    	if(GREP(line,"MIN_THRES")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   UserInput.minthres = stod(d);
    	}
       	if(GREP(line,"BACKGROUND_PROP")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   UserInput.bgProp = stof(d);
    	}
    	if(GREP(line,"MINBG_SIZE")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   UserInput.minbgSize = stoi(d);
    	}
    	if(GREP(line,"MINSIG_SIZE")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   UserInput.minsigSize = stoi(d);
    	}
      if(GREP(line,"ANALYSE_DNA")){
	       v = stringSplit(line,'=');
	       d = trim(v.at(1));
	       bool d_new;
	       if(d =="false") d_new = false;
	       if(d=="true") d_new = true;
	       UserInput.analyseDNA = d_new;
	    }
      if(GREP(line,"Z_TRANSFORM_DNA")){
	       v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   bool d_new;
    	   if(d =="false") d_new = false;
    	   if(d=="true") d_new = true;
    	   UserInput.ztrans_dna = d_new;
    	}
      if(GREP(line,"Z_TRANSFORM_RNA")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   bool d_new;
    	   if(d =="false") d_new = false;
    	   if(d=="true") d_new = true;
    	   UserInput.ztrans_rna = d_new;
    	}
      if(GREP(line,"Z_TRANSFORM_PROT")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   bool d_new;
    	   if(d =="false") d_new = false;
    	   if(d=="true") d_new = true;
    	   UserInput.ztrans_prot = d_new;
    	}
    }
    return ret;
}

bool checkNAs(map<string, subjectClass> &subjectMap, int datatype){
    bool ret = false;
    map<string, subjectClass>::iterator iter;
    set<string> missingG;
    
    for(iter = subjectMap.begin(); iter!=subjectMap.end(); iter++){
	      missingG = iter->second.getMissingG(datatype);
        if(missingG.size()>0){
            ret = true;
            continue;
        }
    }
    return ret;
}


vector<vector<double> > creatematrix(vector<string> genes, vector<string> subjects, map<string, subjectClass>&subjectMap, int datatype){

    vector<vector<double> > ret;
    vector<double> tmp_vec;
    ret.resize(genes.size());  

    for(int i = 0; i<genes.size(); i++){
	      string currg = genes.at(i);
        tmp_vec = extractGeneVec(datatype, currg, subjects, subjectMap);
	      ret[i] = tmp_vec;
    }
    return ret;
}

bool comp(int i, int j) { return i > j;}

vector<vector<double> > getCandidateG(int pos,vector<double> x_j, vector<vector<double> >MAT){
   
     vector<vector<double> >ret;
     for(int j=0; j<MAT.size(); j++){
	      if(MAT[j][pos] != NA_VAL) ret.push_back(MAT[j]);
     }
     return ret;
}

double calEdist(vector<double> v1, vector<double> v2){

    double ret=0;
    double d;
    for(int i=0; i<v1.size(); i++){
	      if(v1.at(i)!=NA_VAL & v2.at(i)!=NA_VAL) d = pow((v1.at(i) - v2.at(i)), 2);
        ret +=d;
    } 
    return sqrt(ret);
}

set<string> getNonOverlap(vector<string> g1, vector<string> g2){

    set<string> SET;
    for(int i=0; i<g2.size(); i++){
	      if(!containsVec(g2.at(i), &g1)) SET.insert(g2.at(i));
    }
    return SET;
}



bool vecContains(int index, vector<int> vec){

    bool ret = false;
    if(vec.size()>0){
      for(int i=0; i<vec.size(); i++){
	      if(vec.at(i)==index) {
	        ret = true;
	        continue;
	      }
      }
    }
    return ret;
}

vector<double> getColvec(int pos, vector<vector<double> >MAT){
    vector<double> ret;
    for(int i=0; i <MAT.size(); i++) ret.push_back(MAT[i][pos]);
    
    return ret;
}


vector<string> SetToVec(set<string> SET){
   vector<string> ret;
   set<string>::iterator it;
   for(it = SET.begin(); it!=SET.end(); it++) ret.push_back(*it);
   return ret;
}

void copyMap(map<string, subjectClass> &MapA, map<string, subjectClass>&MapB){
    map<string, subjectClass>::iterator iter;
    for(iter = MapA.begin(); iter!=MapA.end(); iter++){
	      string sub = iter->first;
  	    MapB[sub] = iter->second;
    }
}

vector<vector<double> > removeRow(vector<double> x, vector<vector<double> >&MAT_new){

    vector<vector<double> > MAT_out;
    vector<int> keep;
    bool same;
    MAT_new.erase(unique(MAT_new.begin(), MAT_new.end()), MAT_new.end());
    for(int i=0; i<MAT_new.size(); i++){
  	    vector<double> vec = MAT_new[i];
  	    same = true;
  	    for(int j=0; j<vec.size(); j++) {
  	        if((vec.at(j)!=NA_VAL) & (x.at(j)!=NA_VAL)) if(vec.at(j)!= x.at(j)) same = false;	
        }
  	    if(!same) keep.push_back(i);
    }
    for(int k=0; k<keep.size(); k++){
	      int index = keep.at(k);
	      MAT_out.push_back(MAT_new[index]);
    }
    return MAT_out;
}


void KNNimpute(vector<int> indices,vector<string> subjects, vector<string> genes, vector<vector<double> > MAT,map<set<string>,vector<vector<double> > >&MAP, map<string, subjectClass> &subjectMap, int knnK, int datatype){

     map<string, subjectClass>::iterator s_iterI;
     map<set<string>, vector<vector<double> > >::iterator iter;
     vector<double> wdist_vec, x_j;
     vector<vector<double> > Cv, MAT_new;
     set<string> g_missing, allg; 
     double yhat;
     set<string>::iterator it;
     
     int nrow = MAT.size();
     int ncol = (MAT.at(0)).size();
     int countNA =0;
	 
     for(int h=0; h<indices.size(); h++){
         
          int i = indices.at(h);
          string curr_g = genes.at(i);
          x_j = MAT[i];
          bool search = false;
 
	      for(iter = MAP.begin(); iter!=MAP.end(); iter++){
	        set<string> curr_block = iter->first;
	        if(contains(curr_g, &curr_block)){
		        MAT_new =iter->second;
		        search=true;
		        continue;
	        }
	      }
              for(int k=0; k<x_j.size(); k++){
	         if(x_j.at(k)==NA_VAL){
		        Cv = getCandidateG(k,x_j, MAT_new);
	          	string curr_sub = subjects.at(k);
 	          	if(Cv.size()>=knnK){
              		  wdist_vec = calweightedEdist(x_j, Cv);
	            	  yhat = imputeMV(k,wdist_vec,knnK, Cv);
	          	}
	          	if(Cv.size()<knnK)  yhat = calMean(getColvec(k, MAT));
	          	s_iterI = subjectMap.find(curr_sub);
			//if(s_iterI == subjectMap.end()) cerr<<"Couldn't find subject!"<<endl;
	          	s_iterI->second.insertMap(curr_g, yhat, datatype);
	          	//s_iterI->second.removeNAgenes(curr_g,datatype);
		        countNA++;
      	 	  }
	       Cv.clear();
	       wdist_vec.clear();
	     }
             x_j.clear();
	     MAT_new.clear();
   }
   cerr<<"\n==========================================================="<<endl;
   cerr<<"Total number of values imputed: "<<countNA<<endl;
}

int whichpos(string g, vector<string> vec){
   int pos;
   for(int i=0; i<vec.size(); i++) if(vec.at(i)==g) {pos = i;continue;}
   return pos;
}

vector<double> extractGeneVec(int datatype, string g, vector<string> sub, map<string, subjectClass>&subjectMap){

     map<string, geneClass>::iterator g_iter;
     map<string, subjectClass>::iterator s_iter;
     vector<double> vec;

     for(int i = 0; i<sub.size(); i++){
	      string curr_s = sub.at(i);
        s_iter = subjectMap.find(curr_s);
	      g_iter = s_iter->second.geneMap.find(g);
        double dd = g_iter->second.getReadings(datatype);
        vec.push_back(dd);
     }
     return vec;
}

void determineCandidate(int max_p, vector<string> genes, vector<vector<double> > MAT, map<set<string>, vector<vector<double> > >&MAP){

    vector<vector<double> >::iterator iter;
    vector<string> genesSET;
    vector<int> selected;
    bool iterate = true;
    int counter=0;
    vector<vector<double> > MAT_new;
    set<string>::iterator it;

    for(int i=0; i<genes.size(); i++){
	      set<string> SETg;
	      vector<string> collect_g;
	      vector<int> curr_select, curr_select2;
	      vector<double> tmp_vec, curr_x, pw_x, dist_vec;
              if(!vecContains(i,selected)){
    	      	  string curr_g = genes.at(i);
	          curr_x = MAT.at(i);
	          for(int j=0; j<genes.size(); j++){
	             if( (j!=i) & (!vecContains(j,selected))){
	               string pw_g = genes.at(j);
	               pw_x = MAT.at(j);
	               collect_g.push_back(pw_g);
	               curr_select.push_back(j);
		             double dd = calEdist(curr_x, pw_x);
	               dist_vec.push_back(dd);
	             }
    	       } 
    	        
    	       vector<int> block = findkminPos(dist_vec,max_p-1);
    	       selected.push_back(i); 
    	       curr_select2.push_back(i);
    	       SETg.insert(curr_g);
    	       genesSET.push_back(curr_g);
    	       MAT_new.push_back(MAT[whichpos(curr_g, genes)]);
    	       for(int k=0; k<block.size(); k++){
    	          string tmpG = collect_g.at(block.at(k));
    	          int pos = whichpos(tmpG, genes);
    	          MAT_new.push_back(MAT[pos]);
    	          SETg.insert(tmpG);
    	          genesSET.push_back(tmpG);
    	          curr_select2.push_back(curr_select.at(block.at(k)));
    	          selected.push_back(curr_select.at(block.at(k)));
    	        }
    	        MAP[SETg] = MAT_new;		
    	        set<string> genesleft = getNonOverlap(genesSET, genes);
    	        //if(genesleft.size() != (genes.size()-selected.size())) cerr<<"Something wrong with size!"<<endl;
    	        if(genesleft.size()<max_p) iterate = false;
                  if (!iterate){
    	              MAT_new.clear();
    	              for(it = genesleft.begin(); it!=genesleft.end(); it++){
    		                int pos = whichpos(*it, genes);
    		                MAT_new.push_back(MAT[pos]);
    	              }
    	              MAP[genesleft] = MAT_new;
    	              break;
    	           }
    	           curr_select.clear();
    	           curr_select2.clear();
    	           collect_g.clear();
    	           SETg.clear();
    	           curr_x.clear(); 
    	           pw_x.clear(); 
    	           dist_vec.clear();
    	           MAT_new.clear();
    	    	}
    }
    //cerr<< "There are a total of " << MAP.size()<< " blocks of genes of size "<<max_p<<" created."<<endl;
   
}



vector<double> calweightedEdist(vector<double> vecA, vector<vector<double> > candidateVec){

    int size_n = vecA.size();
    int counter;
    double tmpA,tmpB,dist;
    vector<double> vecB, vec_dist;
    bool same;
    
    for(int i=0; i<candidateVec.size(); i++){
  	  vecB = candidateVec.at(i);
  	  dist=0;  
  	  counter=0;
  	  same = true;
  	  for(int j=0; j<size_n; j++){
  	    tmpA = vecA.at(j);
  	    tmpB = vecB.at(j);
        if((tmpA!=NA_VAL) & (tmpB!=NA_VAL)){
	           counter++;
	           dist += pow(tmpA-tmpB,2);	
        }
	    }
      dist = sqrt(dist/(double) counter);
	    vec_dist.push_back(dist);
    }     
    return vec_dist;
}


vector<int> findkmaxPos(vector<double> v, int k){
   int count=0;
   vector<int> ret;
   vector<double> vec_ord =v;
   sort(vec_ord.begin(), vec_ord.end());
   vec_ord.erase(unique(vec_ord.begin(), vec_ord.end()), vec_ord.end());
   for(int j=vec_ord.size(); j>0; j--){	
	    double curr = vec_ord.at(j-1);
	    if(count<k){
	      for(int h=0; h<v.size(); h++) if(v.at(h)==curr){ ret.push_back(h); count++;}
	      if(count>=k) break;
	    }
    }
    return ret;
}

vector<int> findkminPos(vector<double> v, int k){

  int count =0;
  vector<int> ret;
  vector<double> vec_ord = v;
  sort( vec_ord.begin(), vec_ord.end());
  vec_ord.erase(unique(vec_ord.begin(), vec_ord.end()), vec_ord.end());
  for(int j=0; j<vec_ord.size(); j++){
      double curr = vec_ord.at(j);
      if(count<k){
          for(int h=0; h<v.size();h++) if(v.at(h)==curr) {ret.push_back(h); count++;}
          if(count>=k) break;
      }
  }
  return ret;
}

double imputeMV(int indexNA,vector<double> weightedD, int K, vector<vector<double> > candidateVec){

    double ret=0;
    set<double> min_Kwt;
    vector<double> tmp_vec;
    double normalizingC=0;

    if(candidateVec.size()>=K){
        vector<int> cid_minK = findkminPos(weightedD,K);
        for(int k=0; k<K; k++){
  	      int index = cid_minK.at(k);
  	      double djk = weightedD.at(index);
  	      if(djk==0) cerr<<"weight D is zero"<<endl;
          normalizingC += 1/djk;
        }
        for(int i=0; i<K;i++){
  	      int index = cid_minK.at(i);
  	      tmp_vec = candidateVec.at(index);
          double djk = weightedD.at(index);
          double wk = 1/(djk * normalizingC);
          ret += wk*tmp_vec.at(indexNA);
        }
    }
    return ret;
}



vector<string> readFile(string myFile,map<string, set<string> > &subtypeMap, map<string, subjectClass> &subjectMap,Input_t &UserInput,set<string> sub,  int dataType){

    ifstream file;
    file.open(myFile.c_str(), ios::in);
    
    if(!file.is_open()){
        cerr <<"Error! Unable to open file...\n";
        exit(1);
    }
    
    bool logtrans = false;
    string datatype;
    if(dataType==1) {datatype ="DNA"; if(UserInput.log_DNA) logtrans=true; }
    if(dataType==2) {datatype ="RNA"; if(UserInput.log_RNA) logtrans=true;}
    if(dataType==3) {datatype ="Proteome"; if(UserInput.log_PROT) logtrans=true;}

    string line;
    bool skipHeader = true;
    int ctr = 1;
    vector<string> SubjectIds;
    set<string> set_sub;
    vector<string> v, v_subtypes,v_subjects,geneset;
    vector<int> v_size, v_index,cid_missing;
    map<string, subjectClass>::iterator m_iter;
    map<string, set<string> >::iterator iter;
    int counter=0;
    while(!file.eof()){
        getline(file, line);  //reading in whole line chopping off '\n'
        if(line=="") {break;}
        if(skipHeader){
            skipHeader = false;
            v = stringSplit(line, '\t');
            SubjectIds = v;
	          SubjectIds.erase(SubjectIds.begin());
            subjectClass *curSub = NULL;
            
            for(int i=0; i !=SubjectIds.size(); i++){
                string subject = SubjectIds.at(i);
            		//only insert into subjectMap if subject has a subtype
            		if(contains(subject, &sub)){
	                 set_sub.insert(subject);
                   m_iter = subjectMap.find(subject);
                   if(m_iter==subjectMap.end()){
                       //create new key
                       curSub = new subjectClass(subject);
                       subjectMap[subject] = *curSub;
                       delete curSub;  
                   }
		            }
		            if(!(contains(subject, &sub))) cerr<<"error with "<<subject<<endl;
            }
      	    set<string>::iterator it;
      	    vector<string>::iterator it2;
      	    for(iter=subtypeMap.begin(); iter!=subtypeMap.end(); iter++){
      		      v_subtypes.push_back(iter->first);   // collection of strings of subtypes
      		      set<string> tmp_mem = iter->second;
	              int n_type =0;
	              for(it = tmp_mem.begin(); it!=tmp_mem.end(); it++){
		                if(contains(*it, &set_sub)){
			                  v_subjects.push_back(*it);  // collection of members ordered by subtype
			                  it2 = find(SubjectIds.begin(), SubjectIds.end(), *it);
			                  int cid = distance(SubjectIds.begin(), it2);   //index of subject in original data
			                  v_index.push_back(cid);
		                    n_type++;  
		                  }
		            }
		            v_size.push_back(n_type);	 // number of subjects in each subtype
            }
            continue;   
        }// close header loop

        v.clear();
        v = stringSplit(line, '\t');
        string curr_g = v.at(0);
        bool keep = true;       
        if(dataType==1) cerr << "Reading in DNA..........."<<ctr<<"\r";
        if(dataType==2) cerr << "Reading in RNA..........."<<ctr<<"\r";
        if(dataType==3) cerr << "Reading in Protein........."<<ctr<<"\r";     

      	int start_pos=0;
      	int end_pos=0;
      	for(int k=0; k<v_size.size(); k++){
  	    if(k==0) end_pos=0;
  	    else start_pos = end_pos;
  	    end_pos = start_pos+v_size.at(k);
  	    int countNA=0;
        for(int j=start_pos; j<end_pos; j++){
            int h = v_index.at(j);
            string dat_str = v.at(h+1);
            if(!logtrans) if(dat_str=="NA"|dat_str=="") countNA++;
  		      if(logtrans) if(dat_str=="NA"|dat_str==""|dat_str=="0") countNA++;
  	    }
	      double tmp_p = (v_size.at(k) - countNA)/(double) v_size.at(k);
	      int diff = v_size.at(k) - countNA;
	     if((diff<UserInput.minobs)| (tmp_p<UserInput.minprop)) keep = false;	
	  }   
	
	  if(keep){
  	   geneset.push_back(curr_g);
  	   int start_pos=0;
  	   int end_pos=0;
       bool nainsert = true;
  	   for(int k=0; k<v_size.size(); k++){
  	      if(k==0) end_pos=0;
  	      else start_pos=end_pos;
  	      end_pos = start_pos +v_size.at(k);	
  	      for(int j=start_pos; j<end_pos; j++){
  		      int h= v_index.at(j);
  		      m_iter = subjectMap.find(SubjectIds.at(h));
  	        double d=0;
  		      double d_new =0;
  		      if(m_iter!=subjectMap.end()){
  		        string dat_str = v.at(h+1);
  		        if(dat_str=="NA"|dat_str==""){
  			          m_iter->second.insertNAgenes(curr_g, dataType);
  			          d_new = NA_VAL;
  			          if(nainsert) {cid_missing.push_back(counter); nainsert=false;}
  		        }
              else{
                  d=atof(dat_str.c_str()); 
  			          if(logtrans){
                      if(d==0){
  				                d_new=NA_VAL;
  				                m_iter->second.insertNAgenes(curr_g,dataType);
  			 	                 if(nainsert){cid_missing.push_back(counter); nainsert=false;}
  			              }
  			              if(d!=0) d_new =(log(d)/log(2));  
    		          }
  			          if(!logtrans) d_new = d;   
              }
              m_iter->second.insertMap(curr_g,d_new,dataType);                           
            } 
  	      }
  	    } //close for loop
        counter++;
      }// close keep loop   
      ctr++;
    }//close while loop
    if(dataType==1) UserInput.cid_dna = cid_missing;
    if(dataType==2) UserInput.cid_rna = cid_missing;
    if(dataType==3) UserInput.cid_prot = cid_missing;
    return geneset;
}



set<string> readPPINetwork(string file, map<string, set<string> >&PPIMap, vector<string> prot_g){

    set<string> ret;

    ifstream myFile;
    myFile.open(file.c_str(), ios::in);

    if(!myFile.is_open()){
	      cerr<<"Error! Unable to open file...\n";
	      exit(1);
    }

    set<string> members;
    vector<string> v;
    string line;
    bool skipHeader=true;
    map<string, set<string> >::iterator m_iter;

    while(!myFile.eof()){
    	getline(myFile, line);
    	if(line=="") break;
    	if(skipHeader){
    	    skipHeader=false;
    	    continue;
      }
  	  v = stringSplit(line,'\t');
      string geneA = v.at(0);
      string geneB = v.at(1);
  	  if(containsVec(geneA, &prot_g) & containsVec(geneB, &prot_g)){
         m_iter = PPIMap.find(geneA);
     	   if(m_iter == PPIMap.end()){
  	        members.insert(geneB);
  	        PPIMap[geneA] = members;
  	        members.clear();
  	      }
          else m_iter->second.insert(geneB);
          m_iter = PPIMap.find(geneB);
          if(m_iter == PPIMap.end()){
              members.insert(geneA);
              PPIMap[geneB] = members;
      	      members.clear();
          }
          else m_iter->second.insert(geneA);
  	      string edge = geneA+"_prot_"+geneB+"_prot";
  	      ret.insert(edge);
       }
    }

    string maxG = maximumG(PPIMap);
    //cerr<<"Gene with the most interactors is "<<maxG<<endl;
    //cerr<<"The size of the PPI-map with interactions in the data is " <<PPIMap.size()<<".\n"<<endl;
    //cerr<<"Number of edges formed is: " << ret.size()<<endl;
    return ret;
}

vector<string> getSubjects(int datatype, vector<string> genes, map<string, subjectClass> &subjectMap){

    vector<string> ret;
    map<string, subjectClass>::iterator m_iter;
    bool check;
    for(m_iter = subjectMap.begin(); m_iter!= subjectMap.end(); m_iter++){
	    string curr_sub = m_iter->first;
	    check = m_iter->second.existSet(datatype);
	    if(check) ret.push_back(curr_sub);
    }
    return ret;
}

void checkDuplicate(vector<vector<double> > &MAT){

    bool ret = false;
    vector<double> vec1,vec2;
    int count =0;
    bool samevec;
    for(int i=0; i<MAT.size(); i++){
	    vec1 = MAT[i];
	    for(int j=i; j<(MAT.size()-1); j++){
	      vec2 = MAT[j+1];
	      samevec = true;
	      for(int k=0; k<vec2.size(); k++) if(vec2.at(k)!=vec1.at(k)) samevec=false;
		    if(!samevec) continue;
	    }
	    if(samevec) {ret = true; count++;}
    }
    if(count>0) cerr<<count<<" duplicates found!"<<endl;      
}


set<string> readPATHWAY(string file, map<string, set<string> >&PATHmap, map<string, string> &PATHwayAnnot){
    set<string> ret;
    ifstream myFile;
    myFile.open(file.c_str(), ios::in);
    
    map<string, set<string> >::iterator iter;

    if(!myFile.is_open()){
      	cerr<<"Error!! Unable to open file...\n";
      	exit(1);
    }
	
    set<string> members;
    vector<string> v;
    string line;
    bool skipHeader=true;

    while(!myFile.eof()){	
      	getline(myFile, line);
      	if(line=="") break;
      	if(skipHeader){
      	   skipHeader=false;
      	   continue;
      	}
      	v = stringSplit(line,'\t');
      	string mem = v.at(0);
      	string path = v.at(1);
	string annot = v.at(2);
      	iter= PATHmap.find(path);
      	if(iter == PATHmap.end()){
      	    members.insert(mem);
      	    ret.insert(path);
      	    PATHmap[path] = members;
	    PATHwayAnnot[path] = annot;
      	}
        else iter->second.insert(mem);
      	members.clear();
    }
    string maxPath = maximumG(PATHmap);
    cerr<<"There are a total of "<< PATHmap.size() <<" number of pathways."<<endl;
    return ret;
}


set<string> readTFNetwork(bool interactType, string file, map<string, set<string> > &TFMap, map<string, int> &InteractMap, vector<string> target_g, vector<string> prot_g){
    
    set<string> ret;
    ifstream myFile;
    myFile.open(file.c_str(), ios::in);
    
    if(!myFile.is_open()){
        cerr<<"Error! Unable to open file...\n";
        exit(1);
    }
    set<string> members;
    vector<string> v;
    string line, type;
    bool skipHeader=true;
    map<string, set<string> >::iterator m_iter;

    while(!myFile.eof()){
        getline(myFile, line);
        if(line=="") break;
        if(skipHeader){
            skipHeader=false;
            continue;
        }
        v = stringSplit(line,'\t');
        string key = v.at(0);
        string mem = v.at(1);
      	if(interactType) type = v.at(2); 
      	if(containsVec(key, &prot_g) & containsVec(mem, &target_g)){
              m_iter = TFMap.find(key);
              if(m_iter == TFMap.end()){
                  members.insert(mem);
                  TFMap[key] = members;
              }
              if(m_iter!=TFMap.end())m_iter->second.insert(mem);
      	      string edge = key+"_prot_"+mem+"_mrna";
      	      if(interactType)InteractMap[edge] = stoi(type);
      	      ret.insert(edge);
              members.clear();
      	}
    }

    string maxGene = maximumG(TFMap);
    //cerr<<"Transcription factor with the largest number of targets in data is "<<maxGene<<endl;
    //cerr<<"The size of the TF-map with interactions in the data is " <<TFMap.size()<<".\n"<<endl;
    //cerr<<"Number of edges formed is: "<<ret.size()<<endl;
    return ret;
    
}

vector<string> StatusReport(map<string, subjectClass>&subjectMap, int datatype1, int datatype2){

    string type1, type2;
    if(datatype1 ==1) type1 = "DNA";
    if(datatype1 ==2) type1 = "RNA";
    if(datatype1 ==3) type1 = "Protein";
    if(datatype1 ==4) type1 = "Methylation";
    if(datatype2 ==1) type2 = "DNA";
    if(datatype2 ==2) type2 = "RNA";
    if(datatype2 ==3) type2 = "Protein";
    if(datatype2 ==4) type2 = "Methylation";

    cerr<<"reporting first 5 subjects with both "<< type1 << " and "<< type2<<"  data for checking...\n";
    cerr<<"----------------------------------------------------------"<<endl;
    int ctr=0;
    vector<string> ret;
    map<string, subjectClass>::iterator m_iter;
    for(m_iter=subjectMap.begin(); m_iter!=subjectMap.end(); m_iter++){
        set<string> setA = m_iter->second.getSet(datatype1);
        set<string> setB = m_iter->second.getSet(datatype2);
        set<string> commonG = IntersectSet(setA,setB);
        if((!setA.empty())&(!setB.empty())){
	        ret.push_back(m_iter->first);
	        if(ctr<5){
            	cerr<< "Subject: "<< m_iter->first <<endl;
            	cerr<< setA.size()<<" number of genes in "<< type1<<" file\n";
            	cerr<< setB.size()<<" number of genes in "<<type2<<" file\n";
            	cerr<< commonG.size()<<" number of common genes in both data types.\n";
            	cerr<<"----------------------------------------------------------"<<endl;	
	        }
	        ctr++;
        }
    }
    cerr<<"\nThere are "<<ctr<<" number of subjects with both "<<type1<<" and "<<type2<<" data."<<endl; 
    return ret;
}


vector<string> getKey(map<string, set<string> > *Map){
     vector<string> ret;
     map<string, set<string> >::iterator iter;
     for(iter = Map->begin(); iter!=Map->end(); iter++) ret.push_back(iter->first);
     return ret;
}


set<string> retrieveKeys(map<string, subjectClass> &MAP){
    set<string> ret;
    map<string, subjectClass>::iterator m_iter;
    for(m_iter=MAP.begin(); m_iter!=MAP.end(); m_iter++) ret.insert(m_iter->first);
    return ret;
}

string maximumG(map<string, set<string> > &MAP){

    string ret;
    map<string, set<string> >::iterator iter;
    int maxSize=0;
    for(iter=MAP.begin(); iter!=MAP.end(); iter++){
	    int tmp_size = iter->second.size();
	    string tmp_g = iter->first;
	    if(tmp_size>maxSize){
	      maxSize = tmp_size;
	      ret = tmp_g;
	    }
    }
    return ret;
}



set<string> insertSubtypeMap(string subtypeFile, map<string, set<string> > &subtypeMap){

     vector<string> ret;
     ifstream file;
     file.open(subtypeFile.c_str(), ios::in);

     if(!file.is_open()){
	      cerr<<"Error reading in subtype information...\n";
        exit(1);
     }

     bool skipHeader=true;
     set<string> sub;
     vector<string> v;
     string line;
     set<string> all_sub;
     map<string, set<string> >::iterator iter;
     int counter = 0;

     while(!file.eof()){
        	getline(file,line);
        	if(line=="") break;
        	if(skipHeader){
        	    skipHeader=false;
           	    continue;
          }
         	v = stringSplit(line,'\t');
        	string curr_s = v.at(0);
        	string curr_t = v.at(1);
        	iter=subtypeMap.find(curr_t);
        	if(iter==subtypeMap.end()){
        	    sub.insert(curr_s);
        	    all_sub.insert(curr_s);
        	    subtypeMap[curr_t] = sub;
        	    ret.push_back(curr_t);
        	}
        	if(iter!=subtypeMap.end()){
        	    iter->second.insert(curr_s);
        	    all_sub.insert(curr_s);
        	}
        	counter++;
          sub.clear();
      }
      cerr<<"Number of subtypes found is "<< ret.size()<<" with " <<all_sub.size()<<" number of subjects and the subtypes are the following:"<<endl;
      int grp_size;
      for(int i=0; i<ret.size(); i++) {
        	iter = subtypeMap.find(ret.at(i));
        	grp_size = iter->second.size();
        	cerr<<ret.at(i)<<":\t" <<grp_size<<endl;
      }	 
      iter = subtypeMap.begin();
      set<string> CommonMem = iter->second;
      for(iter = subtypeMap.begin(); iter!=subtypeMap.end(); iter++){
      	 set<string> mem = iter->second;
      	 CommonMem = IntersectSet(CommonMem, mem);
      }
      if(CommonMem.size()>0){
          cerr<<"There are common members found in the subtype map. Please remove duplicates and try again. They are: "<<endl;
          set<string>::iterator it; 
          for(it = CommonMem.begin(); it!=CommonMem.end(); it++) cerr<<*it<<"\t";
          cerr<<endl;
	        exit(1);
     }
     return all_sub;
}

double calMedian(vector<double> vec){
    double ret;
    sort(vec.begin(), vec.end());
    if((vec.size()%2)==0){ret = (vec.at((vec.size()/2)-1)+vec.at(vec.size()/2))/2;}
    else{
        div_t res = div(vec.size(), 2);
        ret = vec.at(res.quot);
    }
   return ret;
}

double calMeanDik(vector<double> vec){
    double ret =0;
    int n = 0;
    double sum = 0;
    for(int i=0; i<vec.size(); i++){
        if(vec.at(i)!=NA_VAL) {
            sum += abs(vec.at(i));
            n++;
        }
    }
    ret = sum/n;
    return ret; 
}

double calMean(vector<double> vec){
   double ret =0;
   int n = 0;
   double sum = 0;
   for(int i=0; i<vec.size(); i++){
     if(vec.at(i)!=NA_VAL){
    	    sum += vec.at(i);
    	    n++;
    	}
   }
   ret = sum/n;
   return ret; 
}

double calSdn(vector<double> vec){
   double ret=0;
   unsigned n = vec.size();
   double mu = calMean(vec);
   double num =0;
   for(int i=0; i<vec.size(); i++) num+= pow((vec.at(i)-mu),2);
   ret = num/n;
   ret = sqrt(ret);
   return ret;
}

double calSd(vector<double> vec){
   double ret=0;
   unsigned n = vec.size();
   double mu = calMean(vec);
   double num =0;
   for(int i=0; i<vec.size(); i++) num += ((vec.at(i)-mu)*(vec.at(i)-mu));
   ret = num/(n-1);
   ret = sqrt(ret);
   return ret;
}

double calCov(vector<double> x, vector<double> y){
   double ret=0;
   int size = x.size();
   double mean_x = calMean(x);
   double mean_y = calMean(y);
   for(int i =0; i<size; i++){
	    ret+= (x.at(i) - mean_x)*(y.at(i)- mean_y);
   }
   ret = ret/(size-1);
   return ret;
}



vector<double> Benjamini_Hochberg(vector<double> pvalue){
   vector<double> p_adjust;
   vector<double> ord_p;
   int rank;
   double m;
   map<double, double> Pmap;   
   m = pvalue.size();
   
   ord_p = pvalue;
   sort(ord_p.begin(), ord_p.end());
   //ord_p.erase(unique(ord_p.begin(), ord_p.end()),ord_p.end());
   map<double,double>::iterator iter;
   vector<double>::iterator it;
   vector<double> tmp_adjp;

   //check each pvalue infront.
   for(int k=0; k<ord_p.size(); k++){
    	double tmp = ord_p.at(k);
      rank = (k+1);
    	bool check =true;
    	int c=1;
    	//check if next pvalue is equal, if yes, take the index to be new rank
    	if(k<(ord_p.size()-1)){
    	    while(check){
    	    	double tmp_check = ord_p.at(k+c);
    	    	if(tmp_check == tmp) {
    		        rank = (k+c+1);
    	          if(rank==ord_p.size()) check=false;
    	      }
    	    	else check = false;
    	    	c++;
    	    }
    	}
    	double padj = (tmp*m)/rank;
    	padj = min(padj, 1.0);
    	tmp_adjp.push_back(padj);
    	Pmap[tmp] = padj;	
   }

   // to ensure monotonic increasing BH values
   for(int j=0; j<tmp_adjp.size(); j++){
    	double tmp = tmp_adjp.at(j);
    	if(j<(tmp_adjp.size()-1)){
    	    vector<double> tmp_vec = tmp_adjp;
    	    tmp_vec.erase(tmp_vec.begin(), tmp_vec.begin()+(j+1));
    	    it = min_element(tmp_vec.begin(), tmp_vec.end());
    	    double new_padj = min(tmp, *it);
    	    double key = ord_p.at(j);
    	    Pmap[key] = new_padj;
    	}
    	else{
    	    double key = ord_p.at(j);
    	    Pmap[key] = tmp;
    	}
    }
    for(int i=0; i<pvalue.size(); i++){
     	  double tmp = pvalue.at(i);
        iter = Pmap.find(tmp);
  	    double value = iter->second;
     	  p_adjust.push_back(value);
        if(value<0) cerr<<"Neg Pval, pls check!"<<endl;
    }
    return p_adjust; 
}


int countOcurrence(string str, vector<string> vec){
    int ret=0;
    for(int i=0; i<vec.size(); i++) if(vec.at(i)==str) ret++;
    return ret;
}


set<string> findParent(string child, map<string, set<string> >&Map){
     set<string> ret;
     map<string, set<string> >::iterator iter;
      
     for(iter = Map.begin(); iter!=Map.end(); iter++){
      	 string curr_key = iter->first;
      	 set<string> val = iter->second;
      	 if(contains(child,&val)) ret.insert(curr_key);
     }
     return ret; 
}





set<string> identifyGenes(int mink, map<string, set<string> >&subtypeMap, map<string, subjectClass>&subjectMap, int datatype){

     set<string> ret;
     map<string, set<string> >::iterator s_iter;
     map<string, subjectClass>::iterator sub_iter;
     set<string>::iterator it, it2;

     string curr_s,gene;
     set<string> curr_g, ALLg;
     set<string> members;
     map<string, geneClass>::iterator g_iter;

     string data;
     if(datatype==1) data = "DNA";
     if(datatype==2) data = "RNA";
     if(datatype==3) data = "Protein";

     for(s_iter=subtypeMap.begin(); s_iter!=subtypeMap.end(); s_iter++){
      	  members = s_iter->second;
      	  for(it = members.begin(); it!=members.end(); it++){
      		    sub_iter = subjectMap.find(*it);
      	      if(sub_iter->second.existSet(datatype)){
      		       curr_g = sub_iter->second.getSet(datatype);
      		       ALLg = UnionSet(curr_g, ALLg);
      		    }
	        }
     }
   
     for(it2 = ALLg.begin(); it2!=ALLg.end(); it2++){
        	gene = *it2;
        	bool satisfy = true;
        	for(s_iter = subtypeMap.begin(); s_iter!=subtypeMap.end(); s_iter++){
        	   members = s_iter->second;
        	   int counter =0;
             for(it = members.begin(); it!=members.end(); it++){
        		    sub_iter = subjectMap.find(*it);
        		    g_iter = sub_iter->second.geneMap.find(gene);
        		    if(g_iter->second.ratioAvail(datatype)) counter++;
        	   }
        	   if(counter<mink) satisfy=false;
        	 }
	         if(satisfy) ret.insert(gene);
     }
     cerr<<"In the "<<data<<" data, there are "<<ret.size()<<" number of genes with at least "<<mink <<" genes across each of the subtypes."<<endl;
     return ret;
}




void outputData(vector<string> subjects,vector<string> genes, map<string, subjectClass> &subjectMap, int datatype){

    map<string, geneClass>::iterator g_iter;
    map<string, subjectClass>::iterator iter;
    string LAB;
    if(datatype==1)LAB = "DNA";
    if(datatype==2)LAB = "RNA";
    if(datatype==3)LAB = "PROT";
   
    ofstream outF("imputedData_"+LAB+".txt");
    
    outF<<"Gene";
    for(int i=0; i<subjects.size(); i++) outF<<"\t"<<subjects.at(i);
    outF<<"\n";
    
    for(int i=0; i<genes.size(); i++){	
      	string curr_g = genes.at(i);
      	outF<<curr_g;
      	for(int j=0; j<subjects.size(); j++){
      	   iter = subjectMap.find(subjects.at(j));
      	   g_iter = iter->second.geneMap.find(curr_g);
      	   if(g_iter==(iter->second).geneMap.end()) cerr<<"no genemap"<<endl;
      	   double dd = g_iter->second.getReadings(datatype);
      	   if(dd==NA_VAL) outF<<"\t"<<"";
      	   if(dd!=NA_VAL) outF<<"\t"<<dd;
      	}
	      outF<<"\n";
    }
    outF.close();
} 


vector<double> fillPAMmap2(bool useDNA, bool interact,set<string> *Alledges,vector<string> *sub,map<string,int> *InteractMap, map<string, set<string> > *subtypeMap, map<string, subjectClass> *subjectMap, map<string, PAM_t> &PAMmap, map<string, vector<set<string> > > *FDneighborMap){
     map<string, subjectClass>::iterator s_iter;
     map<string, PAM_t>::iterator p_iter, p_iter2;
     map<string, set<string> >::iterator iter;
     map<string, geneClass>::iterator g_iterA, g_iterB;
     map<string, vector<set<string> > >::iterator iter1, iter2;
     map<string, int>::iterator i_iter;

     int count =0;
     vector<double> x_vec;
     double xbar, sigma;
     set<string> curr_grp, rna_set, prot_set, SETin;
     string curr_key,featureA, featureB, curr_e,curr_g, curr_t, str_key;
     string datatype;
     vector<string> lab = getKey(subtypeMap);
     set<string>::iterator it,it2;
     double maxDij;
     vector<string> allkeys, strVec, members,features;
     int interactINT, countPos=0, countNeg=0;

     ofstream outFile("Expressiondata_edges.txt");
     outFile<<"Edge";
    
     count =0;
     vector<set<string> > SubGroups;
     vector<int> GrpCount;
     set<string> grpSub;
     for(int j=0; j<lab.size(); j++){
    	  iter = subtypeMap->find(lab.at(j));
    	  curr_grp = iter->second;
    	  int grp_count=0;
    	  for(it = curr_grp.begin(); it!=curr_grp.end(); it++){
        		if(containsVec(*it, sub)){
        		    grpSub.insert(*it);
        		    outFile<<"\t"<<*it;
        		    grp_count++;
        		}
    	   }
      	 GrpCount.push_back(grp_count);
      	 SubGroups.push_back(grpSub);
      	 grpSub.clear();
     }
     outFile<<endl;

     vector<double> SS;
     vector<double> GRPmean;
     string geneA, geneB, typeA,typeB;
     int datatypeA, datatypeB;
     vector<string> vec_g;
     double xA,xB,xC, rr;

     for(it = Alledges->begin(); it != Alledges->end(); it++){

      	bool first = true;
      	curr_e = *it;
      	PAM_t *pam_struct = new PAM_t();
      	PAMmap[curr_e] = *pam_struct;	
      	p_iter = PAMmap.find(curr_e);
      	delete pam_struct;
      
      	outFile<<curr_e;
      	i_iter = InteractMap->find(curr_e);
      	interactINT = i_iter->second;
      
      	vec_g = stringSplit(curr_e,'_');
      	geneA = vec_g.at(0);
      	typeA = vec_g.at(1);
      	geneB = vec_g.at(2);
      	typeB = vec_g.at(3);
      	p_iter->second.identifier = curr_e;
      	p_iter->second.geneA = geneA;
      	p_iter->second.geneB = geneB;
      	if(typeB=="prot") datatype = "ppi";
      	if(typeB=="mrna") datatype = "tf";
      	p_iter->second.dataT = datatype;
      	if(datatype=="ppi"){datatypeA = 3; datatypeB=3;}
      	if(datatype=="tf"){datatypeA = 3; datatypeB=2;}	
      	for(int j=0; j<lab.size(); j++){
      	     curr_grp = SubGroups.at(j);
      	     vector<double> grp_x;
      	     for(it2 = curr_grp.begin(); it2!=curr_grp.end(); it2++){
            		 string curr_sub = *it2;
            		 s_iter = subjectMap->find(curr_sub);
            		 g_iterA = s_iter->second.geneMap.find(geneA);
            		 g_iterB = s_iter->second.geneMap.find(geneB);
            		 xA = g_iterA->second.getReadings(datatypeA);
            		 xB = g_iterB->second.getReadings(datatypeB);
            		 if(datatype=="tf"){
            		   if(useDNA) xC = g_iterB->second.getReadings(1);
            		   if(interact){
            		      if(interactINT== 1){
                  			if(first){countPos++; first=false;}
                  			if(!useDNA) rr = xA+xB;
                  			if(useDNA) rr = xA+xB-xC;
            		      }
            		      if(interactINT== -1){
            		        if(first){countNeg++; first = false;}
            		        if(!useDNA) rr = xA - xB;
            		        if(useDNA) rr = xA+ (-1*(xB-xC));
            		      }
            		   }
            		   if(!interact){
            		     if(!useDNA) rr =xA+xB;
            		     if(useDNA) rr= xA+xB-xC;
            		   }
            		 }
            		 if(datatype=="ppi") {
				if(first){countPos++; first=false;}
				rr = xA+xB;
			 } 
            		 outFile<<"\t"<<rr;
            		 grp_x.push_back(rr);
            		 x_vec.push_back(rr);
      	     }
      	     double grp_mu = calMean(grp_x);
      	     GRPmean.push_back(grp_mu);
      	     (p_iter->second.xbar_ik).push_back(grp_mu);	
      	     (p_iter->second.classLab).push_back(lab[j]);
      	     grp_x.clear();
     	    }
        	outFile<<endl;
        
        	p_iter->second.GRPsize=GrpCount;
        	p_iter->second.xbar_i = calMean(GRPmean);

          double sigma=0;
          int startpos=0;
          int endpos =0;
          for(int r=0; r<GrpCount.size(); r++){
        	   endpos = startpos+GrpCount.at(r);
        	   for(int r2= startpos;r2<endpos; r2++){
        		 sigma+=pow((x_vec.at(r2) - GRPmean.at(r)),2);
	        }
	        startpos = endpos;
      	}
      	double s_i = sqrt(sigma/((double)(x_vec.size() - GrpCount.size())));
      	
      	p_iter->second.s_i = s_i;
      	SS.push_back(s_i);

        count++;
	      x_vec.clear();
        GRPmean.clear();
     }

     outFile.close();
     if(interact) cerr<<"\nThere are a total of "<<countPos<<" positive-signed edges and "<<countNeg<<" negative-signed edges created.\n";

     double So = calMedian(SS);
     int sizeN = 0;
     vector<vector<string> > neighborG;
     for(int h=0; h<GrpCount.size(); h++) sizeN+=GrpCount.at(h);
     for(p_iter = PAMmap.begin(); p_iter!=PAMmap.end(); p_iter++){
      	p_iter->second.s_o =So;
      	double Si = p_iter->second.s_i;
      	double num, denom;
       	for(int j =0; j<lab.size(); j++){
    	    num = (p_iter->second.xbar_ik[j]) - (p_iter->second.xbar_i);
    	    double mk = (1/(double)GrpCount.at(j))+(1/(double)sizeN);
    	    mk = sqrt(mk);
    	    (p_iter->second.mk).push_back(mk);
    	    denom = mk *(Si + So); 
    	    double dik = num/denom;
    	    (p_iter->second.d_ik).push_back(dik);
  	    }
     }
      // fill up Neighbor Info 
     double curr_dik;
     vector<double> maxTHRES;
     for(int j=0; j<lab.size(); j++) maxTHRES.push_back(0.0);
     string NodeA,NodeB;
     set<string> tf_N, ppi_N, Neigh_rna, Neigh_prot;
     vector<double> collect_dik;
     vector<double> collect_dik_tf;
     vector<double> collect_dik_ppi;
     for(p_iter = PAMmap.begin(); p_iter!=PAMmap.end(); p_iter++){
          curr_e = p_iter->first;
    	    vec_g = stringSplit(curr_e,'_');
    	    geneA = vec_g.at(0);
    	    typeA = vec_g.at(1);
    	    geneB = vec_g.at(2);
    	    typeB = vec_g.at(3);
    	    datatype = p_iter->second.dataT;
    	    if(datatype=="ppi"){datatypeA = 3; datatypeB=3;}
    	    if(datatype=="tf"){datatypeA = 3; datatypeB=2;}
        
           NodeA = geneA +"_"+ typeA;
           NodeB = geneB+"_"+ typeB;
           iter1 = FDneighborMap->find(NodeA);
           iter2 = FDneighborMap->find(NodeB);
           Neigh_rna = iter1->second.at(1);
           Neigh_prot = iter1->second.at(2);
           if(Neigh_rna.size()>0){
               for(it = Neigh_rna.begin(); it!=Neigh_rna.end(); it++){
                   string tmpstr = NodeA+"_"+*it+"_mrna";
                   if(contains(tmpstr,Alledges)){
                       tf_N.insert(tmpstr);
                   }
               }
           }
           if(Neigh_prot.size()>0){
             for(it2 = Neigh_prot.begin(); it2!=Neigh_prot.end(); it2++){
                 string tmpstr = NodeA+"_"+*it2+"_prot";
                 string newstr = getEdge(tmpstr,Alledges);
                 if(!newstr.empty()) ppi_N.insert(newstr);
             }
         }
         Neigh_rna = iter2->second.at(1);
         Neigh_prot = iter2->second.at(2);
         if(Neigh_rna.size()>0){
             for(it = Neigh_rna.begin(); it!=Neigh_rna.end(); it++){
                 string tmpstr = NodeB+"_"+*it+"_mrna";
                 if(contains(tmpstr,Alledges)) tf_N.insert(tmpstr);
             }
         }
         if(Neigh_prot.size()>0){
             for(it2 = Neigh_prot.begin(); it2!=Neigh_prot.end(); it2++){
                 string tmpstr = NodeB+"_"+*it2+"_prot";
                 string newstr = getEdge(tmpstr,Alledges);
                 if(!newstr.empty()) ppi_N.insert(newstr);
             }
         }

         int NumTF = tf_N.size();
         int NumPPI = ppi_N.size();
         p_iter->second.NumFD_tf = NumTF;
         p_iter->second.NumFD_ppi = NumPPI;
         int totalNeighbors = tf_N.size() + ppi_N.size();

         for(int ss =0; ss<lab.size(); ss++){
              curr_dik = p_iter->second.d_ik.at(ss);
              for(it = tf_N.begin(); it!=tf_N.end(); it++){
                  p_iter2 = PAMmap.find(*it);
                  collect_dik.push_back((p_iter2->second.d_ik)[ss]);
                  collect_dik_tf.push_back((p_iter2->second.d_ik)[ss]);
              }
              for(it2 = ppi_N.begin(); it2!=ppi_N.end(); it2++){
                  p_iter2 = PAMmap.find(*it2);
                  //if(p_iter2==PAMmap.end()) cerr<<"Edge not found in map!"<<endl;
                  collect_dik.push_back((p_iter2->second.d_ik)[ss]);
                  collect_dik_ppi.push_back((p_iter2->second.d_ik)[ss]);
               }
               int NumAgree =0;
               double TF_dikbar=0, PPI_dikbar=0;
               if(collect_dik_tf.size()>0){
                    TF_dikbar = calMean(collect_dik_tf);
                    p_iter->second.FDtf_dik_bar.push_back(TF_dikbar);
                    for(int h=0; h<collect_dik_tf.size(); h++){
                        if(collect_dik_tf.at(h)<0 & curr_dik<0) NumAgree++;
                        if(collect_dik_tf.at(h)>0 & curr_dik>0) NumAgree++;
                    }
                }
                if(collect_dik_ppi.size()>0){
                    PPI_dikbar = calMean(collect_dik_ppi);
                    p_iter->second.FDppi_dik_bar.push_back(PPI_dikbar);
                    for(int h=0; h<collect_dik_ppi.size();h++){
                        if(collect_dik_ppi.at(h)<0 & curr_dik<0) NumAgree++;
                        if(collect_dik_ppi.at(h)>0 & curr_dik>0) NumAgree++;
                    }
                }

                double tmp_prop = (double)NumAgree/(double)(totalNeighbors);
	              double wtDikbar = (NumTF*TF_dikbar)+(NumPPI*PPI_dikbar);
            		wtDikbar = wtDikbar/(NumTF+NumPPI);
            		if((NumTF+NumPPI)==0) wtDikbar =0;
            		double phi_f = exp((tmp_prop-0.5)/0.2);
            		double phi = (2*phi_f)/(1+phi_f);
            		if(tmp_prop<0.5|(curr_dik*wtDikbar)<0) phi=0.0;
            		double dik_star = curr_dik + phi*wtDikbar;
            		if(NumTF+NumPPI==0) dik_star = curr_dik;
            		if(dik_star!=dik_star) cerr<<"dikstar problem "<< dik_star<<"\t"<<NumTF<<"\t"<<NumPPI<<endl;
            		if(abs(dik_star)>maxTHRES.at(ss)) maxTHRES.at(ss) = abs(dik_star);
                p_iter->second.prop_agree.push_back(tmp_prop);
		            p_iter->second.dikstar_new.push_back(dik_star);
                collect_dik_tf.clear();
                collect_dik_ppi.clear();
                collect_dik.clear();
            }//close subtypemap
      	    tf_N.clear();
      	    ppi_N.clear();
      	    Neigh_rna.clear();
      	    Neigh_prot.clear();
      }//close tmp_pammap
      return maxTHRES;
}


vector<double> fillPAMmap(bool useDNA, bool interact ,set<string> *Alledges,vector<string> *sub,map<string, int> *InteractMap, map<string, set<string> > *subtypeMap, map<string, subjectClass> *subjectMap, map<string, PAM_t> &PAMmap, map<string, vector<set<string> > > *FDneighborMap){

    map<string, subjectClass>::iterator s_iter;
    map<string, PAM_t>::iterator p_iter, p_iter2;
    map<string, set<string> >::iterator iter;
    map<string, geneClass>::iterator g_iterA, g_iterB;
    map<string, vector<set<string> > >::iterator iter1, iter2;
    map<string, int>::iterator i_iter;

    int count =0;
    vector<double> x_vec;
    double xbar, sigma;
    set<string> curr_grp, rna_set, prot_set, SETin;
    string curr_key,featureA, featureB, curr_e,curr_g, curr_t, str_key;
    string datatype;
    vector<string> lab = getKey(subtypeMap);
    set<string>::iterator it,it2;
    double maxDij;
    vector<string> allkeys, strVec, members,features;
    int interactINT;
    
    count =0;
    vector<set<string> > SubGroups;
    vector<int> GrpCount;
    set<string> grpSub;
    for(int j=0; j<lab.size(); j++){
        iter = subtypeMap->find(lab.at(j));
        curr_grp = iter->second;
        int grp_count=0;
        for(it = curr_grp.begin(); it!=curr_grp.end(); it++){
            if(containsVec(*it, sub)){
                grpSub.insert(*it);
                grp_count++;
            }
        }
        GrpCount.push_back(grp_count);
        SubGroups.push_back(grpSub);
        grpSub.clear();
    }
    
    vector<double> SS;
    vector<double> GRPmean;
    string geneA, geneB, typeA,typeB;
    int datatypeA, datatypeB;
    vector<string> vec_g;
    double xA,xB,xC,rr;
    
    for(it = Alledges->begin(); it != Alledges->end(); it++){
        curr_e = *it;
        PAM_t *pam_struct = new PAM_t();
        PAMmap[curr_e] = *pam_struct;
        p_iter = PAMmap.find(curr_e);
        delete pam_struct;
        
      	i_iter = InteractMap->find(curr_e);
      	interactINT = i_iter->second;
        
        vec_g = stringSplit(curr_e,'_');
        geneA = vec_g.at(0);
        typeA = vec_g.at(1);
        geneB = vec_g.at(2);
        typeB = vec_g.at(3);
        p_iter->second.identifier = curr_e;
        p_iter->second.geneA = geneA;
        p_iter->second.geneB = geneB;
        if(typeA!="prot") cerr<<"Something wrong!"<<endl;
        if(typeB=="prot") datatype = "ppi";
        if(typeB=="mrna") datatype = "tf";
        p_iter->second.dataT = datatype;
        if(datatype=="ppi"){datatypeA = 3; datatypeB=3;}
        if(datatype=="tf"){datatypeA = 3; datatypeB=2;}
        for(int j=0; j<lab.size(); j++){
            curr_grp = SubGroups.at(j);
            vector<double> grp_x;
            for(it2 = curr_grp.begin(); it2!=curr_grp.end(); it2++){
                string curr_sub = *it2;
                s_iter = subjectMap->find(curr_sub);
                g_iterA = s_iter->second.geneMap.find(geneA);
                g_iterB = s_iter->second.geneMap.find(geneB);
                xA = g_iterA->second.getReadings(datatypeA);
                xB = g_iterB->second.getReadings(datatypeB);
		            if(datatype=="tf"){
              		   if(useDNA) xC = g_iterB->second.getReadings(1);
              		   if(interact){
              		      if(interactINT ==1){
                      			if(!useDNA) rr = xA+xB;	
                      			if(useDNA) rr = xA+xB-xC;
              		      }
              		      if(interactINT == -1){
                      			if(!useDNA) rr=xA-xB;
                      			if(useDNA) rr = xA + (-1*(xB-xC));
              		      }
              		   }
              		   if(!interact){
              		      if(!useDNA) rr = xA+xB;
              		      if(useDNA) rr = xA+xB-xC;
              		   }
              	}
              	if(datatype=="ppi") rr = xA+xB;
                grp_x.push_back(rr);
                x_vec.push_back(rr);
            }
            double grp_mu = calMean(grp_x);
            GRPmean.push_back(grp_mu);
            (p_iter->second.xbar_ik).push_back(grp_mu);
            (p_iter->second.classLab).push_back(lab[j]);
            grp_x.clear();
        }
        
        p_iter->second.GRPsize=GrpCount;
        p_iter->second.xbar_i = calMean(GRPmean);
        
        double sigma=0;
        int startpos=0;
        int endpos =0;
        for(int r=0; r<GrpCount.size(); r++){
            endpos = startpos+GrpCount.at(r);
            for(int r2= startpos;r2<endpos; r2++) sigma += (pow((x_vec.at(r2) - GRPmean.at(r)),2));
            startpos = endpos;
        }
        double s_i = sqrt(sigma/((double)(x_vec.size() - GrpCount.size())));
        p_iter->second.s_i = s_i;
        SS.push_back(s_i);
        count++;
        x_vec.clear();
        GRPmean.clear();
    }
    double So = calMedian(SS);
    int counter2=0;
    int sizeN = 0;
    vector<vector<string> > neighborG;
    for(int h=0; h<GrpCount.size(); h++) sizeN+=GrpCount.at(h);
    for(p_iter = PAMmap.begin(); p_iter!=PAMmap.end(); p_iter++){
        p_iter->second.s_o =So;
        double Si = p_iter->second.s_i;
        double num, denom;
        for(int j =0; j<lab.size(); j++){
            num = (p_iter->second.xbar_ik[j]) - (p_iter->second.xbar_i);
            double mk = (1/(double)GrpCount.at(j))+(1/(double)sizeN);
            mk = sqrt(mk);
            (p_iter->second.mk).push_back(mk);
            denom = mk *(Si + So);
            double dik = num/denom;
            (p_iter->second.d_ik).push_back(dik);
        }
        counter2++;
    }
    
    counter2=0;
    // fill up Neighbor Info
    double curr_dik;
    vector<double> maxTHRES; 
    for(int j=0; j<lab.size(); j++) maxTHRES.push_back(0.0);
    string NodeA,NodeB;
    set<string> tf_N, ppi_N, Neigh_rna, Neigh_prot;
    vector<double> collect_dik;
    vector<double> collect_dik_tf;
    vector<double> collect_dik_ppi;
    for(p_iter = PAMmap.begin(); p_iter!=PAMmap.end(); p_iter++){
        curr_e = p_iter->first;
        vec_g = stringSplit(curr_e,'_');
        geneA = vec_g.at(0);
        typeA = vec_g.at(1);
        geneB = vec_g.at(2);
        typeB = vec_g.at(3);
        datatype = p_iter->second.dataT;
        if(datatype=="ppi"){datatypeA = 3; datatypeB=3;}
        if(datatype=="tf"){datatypeA = 3; datatypeB=2;}
                
        NodeA = geneA +"_"+ typeA;
        NodeB = geneB+"_"+ typeB;
        iter1 = FDneighborMap->find(NodeA);
        iter2 = FDneighborMap->find(NodeB);
        Neigh_rna = iter1->second.at(1);
        Neigh_prot = iter1->second.at(2);
        if(Neigh_rna.size()>0){
            for(it = Neigh_rna.begin(); it!=Neigh_rna.end(); it++){
                string tmpstr = NodeA+"_"+*it+"_mrna";
                if(contains(tmpstr,Alledges)){
                    tf_N.insert(tmpstr);  
                }
            }
        }
        if(Neigh_prot.size()>0){
            for(it2 = Neigh_prot.begin(); it2!=Neigh_prot.end(); it2++){
                string tmpstr = NodeA+"_"+*it2+"_prot";
                string newstr = getEdge(tmpstr,Alledges);
                if(!newstr.empty()){
                    ppi_N.insert(newstr);
                }
            }
        }
        Neigh_rna = iter2->second.at(1);
        Neigh_prot = iter2->second.at(2);
        if(Neigh_rna.size()>0){
            for(it = Neigh_rna.begin(); it!=Neigh_rna.end(); it++){
                string tmpstr = NodeB+"_"+*it+"_mrna";
                if(contains(tmpstr,Alledges)) tf_N.insert(tmpstr);
            }
        }
        if(Neigh_prot.size()>0){
            for(it2 = Neigh_prot.begin(); it2!=Neigh_prot.end(); it2++){
                string tmpstr = NodeB+"_"+*it2+"_prot";
                string newstr = getEdge(tmpstr,Alledges);
                if(!newstr.empty()) ppi_N.insert(newstr);
            }
        }
        int NumTF = tf_N.size();
        int NumPPI = ppi_N.size();
        p_iter->second.NumFD_tf = NumTF;
        p_iter->second.NumFD_ppi = NumPPI;
        int totalNeighbors = tf_N.size() + ppi_N.size();
        
        for(int ss =0; ss<lab.size(); ss++){
            curr_dik = p_iter->second.d_ik.at(ss);
            for(it = tf_N.begin(); it!=tf_N.end(); it++){
                p_iter2 = PAMmap.find(*it);
                collect_dik.push_back((p_iter2->second.d_ik)[ss]);
                collect_dik_tf.push_back((p_iter2->second.d_ik)[ss]);
            }
            for(it2 = ppi_N.begin(); it2!=ppi_N.end(); it2++){
                p_iter2 = PAMmap.find(*it2);
                collect_dik.push_back((p_iter2->second.d_ik)[ss]);
                collect_dik_ppi.push_back((p_iter2->second.d_ik)[ss]);
            }
            
            int NumAgree =0;
            double TF_dikbar=0, PPI_dikbar=0;
            if(collect_dik_tf.size()>0){
                TF_dikbar = calMean(collect_dik_tf);
                p_iter->second.FDtf_dik_bar.push_back(TF_dikbar);
                for(int h=0; h<collect_dik_tf.size(); h++){
                    if(collect_dik_tf.at(h)<0 & curr_dik<0) NumAgree++;
                    if(collect_dik_tf.at(h)>0 & curr_dik>0) NumAgree++;
                }
            }
            if(collect_dik_ppi.size()>0){
                PPI_dikbar = calMean(collect_dik_ppi);
                p_iter->second.FDppi_dik_bar.push_back(PPI_dikbar);
                for(int h=0; h<collect_dik_ppi.size();h++){
                    if(collect_dik_ppi.at(h)<0 & curr_dik<0) NumAgree++;
                    if(collect_dik_ppi.at(h)>0 & curr_dik>0) NumAgree++;
                }
            }
            
            double tmp_prop = (double)NumAgree/(double)(totalNeighbors);
            double wtDikbar = (NumTF*TF_dikbar)+(NumPPI*PPI_dikbar);
            wtDikbar = wtDikbar/(NumTF+NumPPI);
            if((NumTF+NumPPI)==0) wtDikbar =0;
            double phi_f = exp((tmp_prop-0.5)/0.2);
            double phi = (2*phi_f)/(1+phi_f);
            if(tmp_prop<0.5|(curr_dik*wtDikbar)<0) phi=0.0;
            double dik_star = curr_dik + phi*wtDikbar;
            if(NumTF+NumPPI==0) dik_star = curr_dik;
            //if(dik_star!=dik_star) cerr<<"dikstar problem "<< NumTF<<"\t"<<NumPPI<<endl;
            if(abs(dik_star)>maxTHRES.at(ss)) maxTHRES.at(ss) = abs(dik_star);
            p_iter->second.prop_agree.push_back(tmp_prop);
            p_iter->second.dikstar_new.push_back(dik_star);
            collect_dik_tf.clear();
            collect_dik_ppi.clear();
            collect_dik.clear();
        }//close subtypemap
        counter2++;
        tf_N.clear();
        ppi_N.clear();
      	Neigh_rna.clear();
      	Neigh_prot.clear();
    }//close tmp_pammap
    return maxTHRES;
}

int whichMax(vector<double> vec){
    int pos=0;
    double curr_max = vec.at(0);
    for(int i =0; i<vec.size(); i++){
	if(vec[i] >curr_max){
	   pos = i;
	   curr_max = vec[i];
	}
    }
    return pos;
}

int whichMin(vector<double> vec){
    int pos=0;
    double curr_min = vec.at(0);
    for(int i=0; i<vec.size(); i++){
      	if(vec.at(i)<curr_min){
       	    pos =i;
      	    curr_min = vec.at(i);
      	}
    }
    return pos;
}

string findClass(string sub, map<string, set<string> > *MAP){
 
   string ret;
   set<string>::iterator it;
   bool found = false;
   map<string, set<string> >::iterator iter;
   for(iter = MAP->begin(); iter!=MAP->end(); iter++){
    	set<string> SET = iter->second;
    	if(contains(sub, &SET)) {ret = iter->first;found=true;continue;}
   }
   return ret;
}


double calNewDik(vector<double> d_ik, double thres, int group){
   double sign_dij = 1.0;
   double dd;
   double d = d_ik[group];
   if(d<0) sign_dij = -1.0;
   dd = abs(d) - thres;
   if(dd<0) dd=0;
   double dij_new = sign_dij*dd;
   return dij_new;
}


double calNewXbarik(PAM_t &pam_struct,double thres, int group, int type){
   vector<double> d_ik;
   double  xik_new;
   if(type==1) d_ik= pam_struct.d_ik;
   if(type==2) d_ik = pam_struct.dikstar_new;
   vector<double> mk = pam_struct.mk;
   double s_o = pam_struct.s_o;
   double s_i = pam_struct.s_i;
   double xbar_i = pam_struct.xbar_i;
   double sign_dij = 1.0;
   double dd;
   double d = d_ik[group];
   if(d<0) sign_dij= -1.0;
   dd = abs(d) - thres;
   if(dd<=0){
    	dd=0;
    	xik_new = xbar_i;
   }
   if(dd>0){	
       double dij_new = ((double) sign_dij) *dd;
       double b= mk[group]*(s_o+s_i)*dij_new;
       xik_new = xbar_i + b;
   }
   return xik_new;
}



set<string> createNeighborMap(vector<string> *prot_genes, vector<string> *rna_genes, map<string, set<string> > *TFmap, map<string, set<string> > *TargetMap, map<string, set<string> > *PPImap, map<string, vector<set<string> > >&FDNeighborMap){

     map<string, set<string> >::iterator tf_iter, ppi_iter, target_iter;
     set<string> gene_set,rna_targets,prot_tf,prot_interact;
     set<string> parentNodes, singletons, RNA_genes,PROT_genes;
     string curr_g, curr_key;
     set<string>::iterator it, it2;
     int num_rna, num_prot, curr_max=0;
     int counter=0;
     string maxNode;
     map<string, vector<set<string> > >::iterator iter;
     set<string> features;
     vector<set<string> > Neighbors;

     for(int i=0; i<prot_genes->size(); i++) PROT_genes.insert(prot_genes->at(i));
     for(int i=0; i<rna_genes->size(); i++) RNA_genes.insert(rna_genes->at(i));

     //for rna, first degree neighbors are only TFs which are targeting the gene
     for(int i =0; i<rna_genes->size(); i++){
        curr_g = rna_genes->at(i);
        curr_key = curr_g +"_mrna";
        gene_set.insert(curr_g);
        Neighbors.push_back(gene_set);
        // for rna, insert empty set into second element
        Neighbors.push_back(rna_targets);
        target_iter = TargetMap->find(curr_g);
      	if(target_iter!=TargetMap->end()) parentNodes = target_iter->second;
              prot_tf = IntersectSet(parentNodes, PROT_genes);
              if(prot_tf.size()==0) singletons.insert(curr_key);
              Neighbors.push_back(prot_tf);
              if(prot_tf.size()>0) FDNeighborMap[curr_key] = Neighbors;
              if(prot_tf.size()>curr_max){
                  curr_max = prot_tf.size();
                  maxNode = curr_key;
              }

              gene_set.clear();
              prot_tf.clear();
              Neighbors.clear();
        }
 
       counter=0;
       //for protein, first degree neighbors are all possible RNA targets and PPI interactors
       for(int i=0; i<prot_genes->size(); i++){
          curr_g = prot_genes->at(i);
          curr_key = curr_g +"_prot";
          gene_set.insert(curr_g);
          Neighbors.push_back(gene_set);
          tf_iter = TFmap->find(curr_g);
          if(tf_iter!=TFmap->end()){
              set<string> tmp_children = tf_iter->second;
              rna_targets = IntersectSet(tmp_children, RNA_genes);
          }
          Neighbors.push_back(rna_targets);	
          ppi_iter = PPImap->find(curr_g);
          if(ppi_iter!=PPImap->end()){
               set<string> interact = ppi_iter->second;
               prot_tf = IntersectSet(interact, PROT_genes);
          }
          Neighbors.push_back(prot_tf);
          if(rna_targets.size()==0 & prot_tf.size()==0) singletons.insert(curr_key);
          if((rna_targets.size() + prot_tf.size())>0) FDNeighborMap[curr_key] = Neighbors;
          if((rna_targets.size()+prot_tf.size())>curr_max){
              curr_max = rna_targets.size()+prot_tf.size();
              maxNode = curr_key;
          }
          gene_set.clear();
          prot_tf.clear();
          rna_targets.clear();
          Neighbors.clear();
       }     
       iter = FDNeighborMap.find(maxNode);
       num_rna = (iter->second).at(1).size();
       num_prot = (iter->second).at(2).size();
  
  
       cerr<<"Generating output for information on the first-degree neighbours of all features in the data...\n"<<endl;
       cerr<<"\nThere is a total of "<<FDNeighborMap.size()<<" features with TF/PPI network interactions and "<<singletons.size()<<" features with no interactions and are removed."<<endl;
       cerr<<maxNode<<" is the feature with the most first-degree neighbours, "<<num_rna<<" RNA targets and "<< num_prot<<" PROT tfs/interactors.\n"<<endl;

       //output the first degree neighbor data
       ofstream outFile("Features_Neighbors.txt");
       outFile<<"Feature\tGene\tType\tNeighbor_inData\tNumNeigh_inData"<<endl;
  
       string currF, currG, currT;
       set<string> rna_set, prot_set;
       vector<string> mem;
       for(iter = FDNeighborMap.begin(); iter!=FDNeighborMap.end(); iter++){
         currF = iter->first;
	       features.insert(currF);
         vector<string> tmpVec = stringSplit(currF,'_');
         currG = tmpVec.at(0);
         currT = tmpVec.at(1);
         
	       rna_set = (iter->second).at(1);
         prot_set = (iter->second).at(2);
         if(rna_set.size()>0){
	          for(it = rna_set.begin(); it!=rna_set.end(); it++){
        		  string tmp = *it+"_mrna";
        		  mem.push_back(tmp);
        	  }
	        }
      	 if(prot_set.size()>0){
      	    for(it = prot_set.begin(); it!= prot_set.end(); it++){
      		      string tmp = *it + "_prot";
       	        mem.push_back(tmp);
      	    }
      	 }
      	 string str_new = concatenate(mem,';');
      	 outFile<<currF<<"\t"<<currG<<"\t"<<currT<<"\t"<<str_new<<"\t"<<mem.size()<<endl;

         mem.clear();
         rna_set.clear();
         prot_set.clear();
     }
     outFile.close();
     return features;
}


void ReverseMap(map<string, set<string> >&Map1, map<string, set<string> >&Map2){

      map<string, set<string> >::iterator iter1,iter2;
      set<string> all_children, children, parent_set;
      set<string>::iterator it;
      string parent, child, MaxP;
      int MaxParent=0;

      for(iter1 = Map1.begin(); iter1!=Map1.end(); iter1++){
	       parent = iter1->first;
         children = iter1->second;
      	 for(it = children.begin(); it!=children.end(); it++){
        		child = *it;
        		all_children.insert(child);
         		iter2 = Map2.find(child);
        		if(iter2==Map2.end()){
        		    parent_set.insert(parent);
        		    Map2[child] = parent_set;
        		}
        		if(iter2!=Map2.end())iter2->second.insert(parent);
        		parent_set.clear();
      	 }
      }
     for(iter2=Map2.begin(); iter2!=Map2.end(); iter2++){
      	if((iter2->second.size()) >MaxParent){
      		MaxParent = iter2->second.size();
      		MaxP = iter2->first;
      	}
     }
     cerr<<"There are a total of "<<Map2.size()<<" TF targets where "<<MaxP<<" has the most TF targetor (parent nodes) of "<<MaxParent<<".\n"<<endl;
}


double computeL2norm(double dd, vector<double>vec){
    double ret=0;
    for(int i=0; i<vec.size(); i++) ret += pow((vec.at(i) - dd),2);
    ret = ret/(vec.size());
    return ret;
}


double CVKfold(bool useDNA, bool interact,int Kfold,set<string> *Alledges,vector<string> *sub, int thres_size,map<string, int> *InteractMap, map<string, vector<set<string> > > *NeighborMap,  map<string, set<string> > *subtypeMap, map<string, subjectClass> *subjectMap, map<string, PAM_t> &PAMmap){

    double minClassError=999.0;
    double selectThres;
    map<int, vector<string> > CVmap;
    random_device rd;
    
    mt19937 g(rd());
    //shuffle(indices.begin(), indices.end(), g);
    //for(int j=0; j<indices.size(); j++) cerr<<indices.at(j)<<"\t";

    map<string, vector<set<string> > >::iterator N_iter;
    map<string, set<string> >::iterator s_iter;
    map<string, subjectClass>::iterator sub_iter;
    vector<string> grplabel = getKey(subtypeMap);

    //Equal allocation of subjects into Kfold CV
    vector<int> Subtype_size;
    vector<vector<string> > Subtype_grp;
    set<string> currGrp;
    set<string>::iterator it;
    vector<string> subGrp;
      
    cerr<<"Number of samples used for crossvalidation in each subtype:"<<endl;
    for(int i=0; i<grplabel.size(); i++){  
      	s_iter = subtypeMap->find(grplabel.at(i));
       	currGrp = s_iter->second;
      	for(it = currGrp.begin(); it!=currGrp.end(); it++) if(containsVec(*it,sub)) subGrp.push_back(*it);
      	Subtype_grp.push_back(subGrp);
      	Subtype_size.push_back(subGrp.size());
      	cerr<<grplabel.at(i)<<": "<<subGrp.size()<<endl;
      	subGrp.clear();
    }
    cerr<<endl;
    vector<vector<int> > indices;
    vector<int> tmpInd;
    int minSize;
    for(int i =0; i<grplabel.size(); i++){
      	int currSize = Subtype_size.at(i);
      	if(i==0) minSize = currSize;
        if(minSize<currSize) minSize = currSize;
      	for(int j=0; j<currSize; j++) tmpInd.push_back(j+1);
      	shuffle(tmpInd.begin(), tmpInd.end(),g);
      	indices.push_back(tmpInd);
      	tmpInd.clear();	
    }
    //Kfold is default to minimum group size or user-specified value
    if(Kfold > minSize) {
        Kfold = minSize;
        cerr<<"Number of cross-validation K is re-specified to the smallest group size: "<<minSize<<endl;
    }
    cerr<<endl;
    vector<string> sample_test;
    vector<int> start, end;
    for(int i=0; i<indices.size(); i++) {start.push_back(0); end.push_back(0);}
        int grp_size;
        div_t divresult;
        vector<string> curr_s;
        vector<div_t> DIVresult;
        for(int i=0; i<indices.size(); i++) DIVresult.push_back(div((indices.at(i)).size(), Kfold));
        for(int K=1; K<(Kfold+1); K++){
    	  for(int i=0; i<indices.size(); i++){
      	   divresult = DIVresult.at(i);
      	   int quot =divresult.quot;
      	   int rem = divresult.rem;
      	   if(K<=rem) grp_size = quot+1;
      	   if(K>rem) grp_size = quot;
      	   end.at(i) = start.at(i)+grp_size;
       	   for(int j=start.at(i); j<end.at(i); j++){
      	      vector<string> Curr = Subtype_grp.at(i);
      	      int pos = (indices[i]).at(j);
      	      string tmp_s = Curr[pos-1];
      	      curr_s.push_back(tmp_s);
      	      //temporary
      	      string CLASS = findClass(tmp_s,subtypeMap);
      	   }
      	   start.at(i) = end.at(i);
      	}
      	CVmap[K] = curr_s;
      	curr_s.clear();
    }   

    vector<double> maxthres;
    for(int j=0; j<grplabel.size(); j++) maxthres.push_back(0.0);
    map<int, vector<string> >::iterator iter;
    map<int, map<string, PAM_t> > BIGmap;
    map<int, map<string, PAM_t> >::iterator B_iter;
    map<string, PAM_t>::iterator piter,piter2;
    vector<string> curr_grp, train_grp, geneVec, features, Edges;
    map<string, PAM_t> tmp_map;

    int key=1;
    for(iter = CVmap.begin(); iter!=CVmap.end(); iter++){
        curr_grp = iter->second;
        for(int h=0; h<sub->size(); h++) if(!containsVec(sub->at(h), &curr_grp)) train_grp.push_back(sub->at(h));
        vector<double> cvthres = fillPAMmap(useDNA,interact,Alledges,&train_grp,InteractMap,subtypeMap, subjectMap, tmp_map, NeighborMap);
        BIGmap[key] = tmp_map;
        for(int j=0; j<grplabel.size(); j++) if(cvthres.at(j)>maxthres.at(j)) maxthres.at(j) = cvthres.at(j);
        tmp_map.clear();
        train_grp.clear();
        curr_grp.clear();
        key++;
    }

    double Avgthres = calMean(maxthres);
    vector<double> thresholdparam;
    double space = Avgthres/(double)(thres_size-1);

    for(int i=0; i<thres_size; i++) thresholdparam.push_back(i*space);

    ofstream outF("CVerrors.txt");
    outF<<"Threshold\tCVerror";
    for(int j=0; j<grplabel.size(); j++) outF<<"\tCVerror_"+grplabel.at(j);
    outF<<"\tEdgesSelected"<<endl;
  
    vector<double> errors, CVerrors, curr_grperror;
    vector<vector<double> > GRPCVerrors;
    double dij_new, dij, xik_new, xik, si, so, curr_thres;

    map<string, PAM_t> tmp_pammap;
    double curr_dik;
    string curr_E, newedge, curr_sub, geneA,geneB;

    map<string, int>::iterator i_iter;
    map<string, geneClass>::iterator g_iter;	
    vector<int> GRPsize;
    vector<vector<int> > Grpgenesurv;
    vector<int> Binvec;
    vector<double> geneSurv_vec, vec_dik;
    int genedied, ccc; 
    int groupsize, typeA, typeB, interactINT, gene_counter;
    double x_star,x_starA, x_starB, x_starC, tmp_score;
    //vector<double> prior;
    vector<double> TMPSCORE;
    double prior;

    for(int j=0; j<thres_size; j++){
        curr_thres = thresholdparam.at(j);
        outF<<curr_thres;
        for(int K=1; K<(Kfold+1); K++){
	    if(j==0) cerr<<"Carrying out cross-validation on fold......."<<K<<"\r";
            B_iter = BIGmap.find(K);
            tmp_pammap = B_iter->second;
            iter = CVmap.find(K);
            vector<string> testsample = iter->second;
            int thres_error=0;
            vector<int> THRES_error, grpsize;
	          for(int r=0; r<grplabel.size(); r++){ 
            		THRES_error.push_back(0); 
            		grpsize.push_back(0);
            		GRPCVerrors.push_back(curr_grperror);
        	   }
      	    for(int t=0; t<testsample.size();t++){
                curr_sub = testsample.at(t);
                sub_iter = subjectMap->find(curr_sub);
                string grp_sub = findClass(curr_sub, subtypeMap);
      	        gene_counter=0;
                vector<double> score;
      	        if(t==0)for(int gg=0; gg<subtypeMap->size(); gg++) Grpgenesurv.push_back(Binvec);
      
       	        for(int gg=0; gg<subtypeMap->size(); gg++) TMPSCORE.push_back(0); // score resets for every test sample
       	       
                for(it = Alledges->begin(); it!=Alledges->end(); it++){
                    string curr_edge = *it;
      		          i_iter = InteractMap->find(curr_edge);
      	  	        interactINT = i_iter->second;
                    piter = tmp_pammap.find(curr_edge);
                    geneA = piter->second.geneA;
      		          geneB = piter->second.geneB;
                    string datatype = piter->second.dataT;
            		    double xbar = piter->second.xbar_i;
            		    vec_dik = piter->second.d_ik;
            
            		    if(datatype=="tf") {typeA = 3;typeB=2;}
            		    if(datatype=="ppi"){typeA=3;typeB=3;}
      
      	            g_iter = sub_iter->second.geneMap.find(geneA);
                    x_starA = g_iter->second.getReadings(typeA);
      	            g_iter = sub_iter->second.geneMap.find(geneB);
      	            x_starB = g_iter->second.getReadings(typeB);
      	            if(useDNA) x_starC = g_iter->second.getReadings(1);

            		    if(t==0 & gene_counter==0) {
            		       GRPsize = piter->second.GRPsize;
                 			 prior = (double)1/(double) GRPsize.size();
                 	  }
                    double si = piter->second.s_i;
                    double s0 = piter->second.s_o;

            		   if(datatype=="tf"){
                			if(interact){
                			     if(interactINT==1){
                      				if(!useDNA) x_star = x_starA + x_starB;
                      				if(useDNA) x_star = x_starA+x_starB-x_starC;
                			     }
                			     if(interactINT== -1){
                      				if(!useDNA) x_star = x_starA-x_starB;
                      				if(useDNA) x_star = x_starA + (-1*(x_starB-x_starC));
                			     }
                			 }
                			 if(!interact){
                			     if(!useDNA) x_star = x_starA + x_starB;
                			     if(useDNA) x_star= x_starA + x_starB - x_starC;
                			  }
                	 }
            
            		   if(datatype=="ppi") x_star = x_starA+ x_starB;
            		     
            		   for(int gg=0; gg<subtypeMap->size(); gg++){
            		        double curr_dik = vec_dik[gg];
                			  //standarding to group threshold
                			  double modthres = (curr_thres/Avgthres) * maxthres.at(gg);
                        double dik_new = calNewDik(piter->second.dikstar_new,modthres,gg);
                			  double xij_new = calNewXbarik(piter->second,modthres,gg,2);
                			  if(t==0){
                			    if(dik_new==0) (Grpgenesurv.at(gg)).push_back(1);
                			    if(dik_new!=0) (Grpgenesurv.at(gg)).push_back(0);
                			  }	
                        TMPSCORE.at(gg) += pow((x_star - xij_new),2)/pow((si+s0),2);
                   }// close subtype
                   gene_counter++;
	           }//close edge

		         for(int gg=0; gg<subtypeMap->size(); gg++){
		            double ss = TMPSCORE.at(gg) - 2*log(prior);
                score.push_back(ss);
             }

      		    bool same = true;
      		    set<string> tmp_sc;
      		    string str;
      		    str = to_string(score.at(0));
      		    //cerr<<str<<endl;
      		    tmp_sc.insert(str);
      		    for(int g=1; g<score.size(); g++){
      		        str = to_string(score.at(g));
      		        if(!(contains(str, &tmp_sc))) same = false;
      		        tmp_sc.insert(str);
      		    }
      		    int index;
      		    if(same) thres_error++;
      		    if(!same){
      			    index = whichMin(score);
      			    if(grplabel[index]!=grp_sub) thres_error++;
      			    //cerr<<"Not all score at the same"<<endl;
      		    }
      		    for(int r=0; r<grplabel.size(); r++){
      	         if(grp_sub==grplabel.at(r)){
        			      (grpsize.at(r))++;
        			      if(same) (THRES_error.at(r))++;
        			      if(!same) {if(grplabel[index]!=grp_sub) (THRES_error.at(r))++;}
        	       }
      		    }
	            score.clear();
	            TMPSCORE.clear();
	         }// close testsample
	         double overall_error = (double)thres_error/(double)testsample.size();
	         double THRES_grpError;

        	 for(int r=0; r<grplabel.size(); r++){
        		    //cerr<<"Current thres: "<<curr_thres<<"\tGrp: "<<grplabel[r]<<"\tThres: "<<THRES_error.at(r)<<"\tGrpSize: "<<grpsize.at(r)<<endl; 
        		    THRES_grpError= (double) THRES_error.at(r)/ (double) grpsize.at(r);
        		    if(grpsize.at(r)!=0) GRPCVerrors.at(r).push_back(THRES_grpError);
        		}
		        CVerrors.push_back(overall_error);
		        genedied =0;

		        for(int c=0; c<Alledges->size(); c++){
      		     bool died = true;
      		     for(int b=0; b<Grpgenesurv.size(); b++){
            			vector<int> tmp_binvec = Grpgenesurv.at(b);
            			if(tmp_binvec.at(c)==0) {died = false; continue;}
      		     }
		          if(died) genedied++;
		        }

	          double genesurv = (double)Alledges->size() - (double)genedied;
        		geneSurv_vec.push_back(genesurv);
        		Grpgenesurv.clear();
            THRES_error.clear();
            grpsize.clear();
	      }//close CV loop
        double NumSurv = floor(calMean(geneSurv_vec));
	      double overallError = calMean(CVerrors);
	      errors.push_back(overallError);
	      if(overallError<minClassError){
	          minClassError = overallError;
	          selectThres = curr_thres;
	      }
	      outF<<"\t"<<overallError;
	      for(int r=0; r<grplabel.size(); r++) outF<<"\t"<<calMean(GRPCVerrors.at(r));
	      outF<<"\t"<<NumSurv<<endl;
	      CVerrors.clear();
	      GRPCVerrors.clear();
	      geneSurv_vec.clear();
    }//close threshold loop
    errors.clear();
    return selectThres;
}

vector<vector<string> > GetNeighbors(string str, map<string, vector<set<string> > >&Map){
    map<string, vector<set<string> > >::iterator iterA, iterB;
    set<string> set_rnaA,set_rnaB, set_protA, set_protB;
    set<string>::iterator it1,it2;
    vector<vector<string> > ret;
    vector<string> tmp;
    
    vector<string> str_vec = stringSplit(str,'_');
    string gA = str_vec[0] +"_"+str_vec[1];
    string gB = str_vec[2] +"_"+str_vec[3];

    iterA = Map.find(gA);
    iterB = Map.find(gB);
    set_rnaA = (iterA->second)[1];
    set_protA = (iterA->second)[2];
    set_rnaB = (iterB->second)[1];
    set_protB = (iterB->second)[2];

    if(set_rnaA.size()>0){
       for(it1 = set_rnaA.begin(); it1!=set_rnaA.end(); it1++){
           string curr_g = *it1+"_mrna";
           if(curr_g!=gB)tmp.push_back(curr_g);
       }
    }
    if(set_protA.size()>0){
        for(it2 = set_protA.begin(); it2!=set_protA.end(); it2++){
            string curr_g = *it2 +"_prot";
            if(curr_g!=gB) tmp.push_back(curr_g);
        }
    }
    ret.push_back(tmp);
    tmp.clear();
    if(set_rnaB.size()>0){
      	for(it1 = set_rnaB.begin(); it1!=set_rnaB.end(); it1++){
      	    string curr_g = *it1+"_mrna";
      	    if(curr_g!=gA) tmp.push_back(curr_g);
      	}
    }
    if(set_protB.size()>0){
      	for(it2=set_protB.begin(); it2!=set_protB.end(); it2++){
          	    string curr_g = *it2+"_prot";
       	    if(curr_g!=gA) tmp.push_back(curr_g);
      	}
    }
    ret.push_back(tmp);
    return ret;
    tmp.clear();
    ret.clear();

}

vector<vector<set<string> > > outputGeneSurv(bool useDNA,bool Interact,vector<double> minThres,set<string> *Alledges,vector<string> *subjects, map<string, set<string> > *subtypeMap, map<string, PAM_t> *PAMmap,map<string, subjectClass> *subjectMap, map<string, int> *InteractMap){

     vector<vector<set<string> > > edgeSurv;
     vector<set<string> >edgeSurvUp, edgeSurvDown;

     vector<string> grplabel = getKey(subtypeMap);
     map<string, PAM_t>::iterator piter;
     map<string, set<string> >::iterator s_iter;
     map<string, subjectClass>::iterator sub_iter;
     map<string, geneClass>::iterator g_iter;
     map<string, int>::iterator i_iter;
     set<string>::iterator it;
     vector<string> currG;
     string curr_k, curr_edge, geneA, geneB, NodeA,NodeB;
     string datatype;
     int typeA, typeB;
         
     int ccc =0;
     bool verbose=true;
     ofstream outFile("EdgesSelected_minThres.txt");
     outFile<<"Edge\tNodeA\tNodeB\tInteractionType";
     for(int j=0; j<grplabel.size(); j++) outFile<<"\t"+grplabel.at(j)+"_dikNew\t"+grplabel.at(j)+"_sigEdge\t"+grplabel.at(j)+"_dir\t"+grplabel.at(j)+"_absdiknew";
     outFile<<endl;
     int counter=0;
     int totalSurv=0;

     ofstream outF2("SampleClass_Probabilities.txt");
     outF2<<"Subject";
     for(int j=0; j<grplabel.size(); j++) outF2<<"\tProb_"+grplabel.at(j);
     outF2<<"\tTrueClass\tPredictedClass"<<endl;

     vector<int> Surv_GRP;
     set<string> tmp_setUp,tmp_setDown, allNodes;
     set<string> EdgeSurv;
     vector<set<string> > GrpNodes;
     vector<double> keepdik;
     vector<string> keeplab;
   
     for(int j=0; j<grplabel.size(); j++) {
    		Surv_GRP.push_back(Alledges->size()); 
    		edgeSurvUp.push_back(tmp_setUp);
    		edgeSurvDown.push_back(tmp_setDown); 
    		GrpNodes.push_back(tmp_setUp);
     }
  	for(it = Alledges->begin(); it!=Alledges->end(); it++){
  		set<string> Survivors;
  		curr_edge = *it;
  		piter = PAMmap->find(curr_edge);
  		vector<string> strVec = stringSplit(curr_edge,'_');
  		NodeA = strVec.at(0)+"_"+strVec.at(1);
  		NodeB = strVec.at(2)+"_"+strVec.at(3);
  		allNodes.insert(NodeA);
  		allNodes.insert(NodeB);
  		datatype = piter->second.dataT;
  	  	if(datatype=="tf") {typeA = 3;typeB=2;}
  	  	if(datatype=="ppi") {typeA = 3;typeB=3;} 
  		bool survive=false;
  		vector<double> curr_dik = piter->second.dikstar_new;
		

  		for(int gg=0; gg<grplabel.size(); gg++){
  		    string surv="died";	
  		    double Dik_new = calNewDik(curr_dik,minThres[gg],gg);
  		    if(Dik_new!=0) surv=datatype;
  		    keepdik.push_back(Dik_new);
  		    keeplab.push_back(surv);
  		    if(abs(Dik_new)>0) {
        			survive=true;
        			if(Dik_new>0){(edgeSurvUp.at(gg)).insert(curr_edge);}
        			if(Dik_new<0){(edgeSurvDown.at(gg)).insert(curr_edge);}
			        (GrpNodes.at(gg)).insert(NodeA);
				(GrpNodes.at(gg)).insert(NodeB);
        			EdgeSurv.insert(curr_edge);
  		    }
  		    if(Dik_new ==0) Surv_GRP.at(gg)--; 
  	   	}
  		if(survive){
  		     totalSurv++;
  		     outFile<<curr_edge<<"\t"<<NodeA<<"\t"<<NodeB<<"\t"<<datatype;
  	       		for(int gg=0;gg<grplabel.size(); gg++) {
  			   double tmp_d = keepdik.at(gg);
  			   string dir = "died";
  			   if(tmp_d>0) dir = "up";
  			   if(tmp_d<0) dir = "down";
  			   outFile<<"\t"<<tmp_d<<"\t"<<keeplab.at(gg)<<"\t"<<dir<<"\t"<<abs(tmp_d);
  		  }
  		  outFile<<endl;
  	 }
  	 counter++;
  	 keepdik.clear();
  	 keeplab.clear();
  	}//close Alledges
	cerr<<"\n\nAt the minimum threshold, a total of "<<totalSurv<<" edges survived across the subtypes and within each group, number of edges survived that are:"<<endl;
	for(int k=0; k<grplabel.size(); k++) cerr<<grplabel.at(k)<<":\t"<<Surv_GRP[k]<<endl;

	cerr<<"\nCalculating class probabilities on samples..."<<endl;
        
	string PredClass;
        set<string> edgeSelected;
        double si, s0,xbar,x_star, x_starA, x_starB, x_starC, tmp_score, prior;
	vector<double> score;
	//vector<int> classSize;
	int interactINT;
	double thres_error=0;

	//for(int i=0; i<grplabel.size(); i++) classSize.push_back(0);
 
	for(int t =0; t<subjects->size(); t++){
	    string curr_sub = subjects->at(t);
	    sub_iter = subjectMap->find(curr_sub);
	    string grp_sub = findClass(curr_sub, subtypeMap);
	    for(int gg=0; gg<grplabel.size(); gg++){
		//if(grp_sub == grplabel.at(gg)) classSize.at(gg)++;
		tmp_score =0;
		for(it = EdgeSurv.begin(); it!=EdgeSurv.end(); it++){
		    string curr_edge = *it;
		    if(Interact){
		    	i_iter = InteractMap->find(curr_edge);
		    	interactINT = i_iter->second;
		    }
		    piter = PAMmap->find(curr_edge);
		    geneA = piter->second.geneA;
		    geneB = piter->second.geneB;
		    datatype = piter->second.dataT;
		    xbar = piter->second.xbar_i;
		    if(datatype=="tf"){ typeA =3; typeB=2;}
		    if(datatype =="ppi"){ typeA=3; typeB=3;}
		    g_iter = sub_iter->second.geneMap.find(geneA);
		    x_starA = g_iter->second.getReadings(typeA);	
		    g_iter = sub_iter->second.geneMap.find(geneB);
		    x_starB = g_iter->second.getReadings(typeB);	
		    if(useDNA) x_starC = g_iter->second.getReadings(1);
		    si = piter->second.s_i;	
		    s0 = piter->second.s_o;

		    if(datatype=="tf"){
			if(Interact){
			    if(interactINT==1){
				if(!useDNA) x_star = x_starA+x_starB;
				if(useDNA) x_star = x_starA+x_starB-x_starC;
			    }
			   if(interactINT== -1){
				if(!useDNA) x_star = x_starA - x_starB;
				if(useDNA) x_star = x_starA + (-1*(x_starB-x_starC));
			   }
			}
			if(!Interact){
			    if(!useDNA) x_star = x_starA +x_starB;
			    if(useDNA) x_star = x_starA +x_starB - x_starC;
			}
		    }
		    if(datatype=="ppi") x_star = x_starA+x_starB;
		    double xij_new = calNewXbarik(piter->second, minThres.at(gg), gg, 2);
		    tmp_score += (pow((x_star-xij_new),2)/pow((si+s0),2));
		 }
		 score.push_back(tmp_score);
	   }// close subtype
	   prior = (double) 1/(double) grplabel.size();
	   double num=0, den=0, maxS=0;
	   vector<double> prob;
	   for(int gg=0; gg<grplabel.size(); gg++){
		score.at(gg) = score.at(gg) - 2*log(prior);
		double currMax = score.at(gg);
		if(currMax>maxS) maxS = currMax;
	   }
	   for(int gg=0; gg<grplabel.size(); gg++){
		double ss = score.at(gg) - maxS;
		num = exp((-0.5)*ss);
		den += num;
	   	prob.push_back(num);
	   }
	   for(int gg=0; gg<grplabel.size(); gg++) prob.at(gg) = prob.at(gg)/den;
	   set<string> tmp_sc;
	   bool same = true;
	   string str = to_string(score.at(0));
	   tmp_sc.insert(str);
	   for(int g=1; g<score.size(); g++){
		str = to_string(score.at(g));	
		if(!contains(str,&tmp_sc)) same = false;
		tmp_sc.insert(str);
	   }
	   if(same){
		thres_error++;
		PredClass = "Tie";
	   }
	   if(!same){
		int index = whichMax(prob);
		PredClass = grplabel[index];
		if(PredClass!= grp_sub) thres_error++;
	   }
	   outF2<<curr_sub;
	   for(int gg=0; gg<grplabel.size(); gg++) outF2<<"\t"<<prob.at(gg);
	   outF2<<"\t"<<grp_sub<<"\t"<<PredClass<<endl;
	   score.clear();
	   prob.clear();

       }
	double ERROR = (double) thres_error/(double)subjects->size();
        double ee = round((ERROR*1000))/10;

	cerr<<"The overall misclassification error based on the model trained in cross-validation is "<<ee<<"%."<<endl;

  	edgeSurv.push_back(edgeSurvUp);
  	edgeSurv.push_back(edgeSurvDown);
  
  	ofstream outFile2("AttributesTable.txt");
	outFile2<<"Node\tGene\tType";
	for(int j = 0; j<grplabel.size(); j++) outFile2<<"\t"+grplabel.at(j)+"_surv";
	outFile2<<endl;	
	
	for(it = allNodes.begin(); it!=allNodes.end(); it++){
		string currN = *it;
		vector<string> strVec = stringSplit(currN,'_');
		outFile2<<currN<<"\t"<<strVec.at(0)<<"\t"<<strVec.at(1);
		for(int j=0; j<grplabel.size(); j++){
		   string surv="died";
		   if(contains(currN, &GrpNodes.at(j)))surv = strVec.at(1);
		   outFile2<<"\t"<<surv;
		}
		outFile2<<endl;
  	}
	outFile.close();
	outFile2.close();
	outF2.close();
	return edgeSurv;
}



void createNNmap(set<string> *Alledges, map<string, vector<set<string> > > *FDneighborMap, map<string, vector<set<string> > > *NNmap){

    set<string>::iterator it,itA, itB;
    set<string> *Neigh_rna = NULL;
    set<string> *Neigh_prot = NULL;
    set<string> *tf_neigh = NULL;
    set<string> *ppi_neigh = NULL;
    string NodeA, NodeB, curr_edge;
    map<string, vector<set<string> > >::iterator iter1, iter2;
    vector<string> strVec;
    int counter=0;
    vector<set<string> > *NN = NULL;

    for(it = Alledges->begin(); it!=Alledges->end(); it++){
        Neigh_rna = new set<string>; 
	      Neigh_prot = new set<string>;
        tf_neigh = new set<string>;
        ppi_neigh = new set<string>;
	      NN = new vector<set<string> >;

      	curr_edge = *it;
      	strVec = stringSplit(curr_edge,'_');
      	NodeA = strVec.at(0) +"_"+ strVec.at(1);
      	NodeB = strVec.at(2) +"_"+ strVec.at(3);	
      	iter1 = FDneighborMap->find(NodeA);
      	iter2 = FDneighborMap->find(NodeB);
        *Neigh_rna = iter1->second.at(1);
      	*Neigh_prot = iter1->second.at(2);
      	if(Neigh_rna->size()>0){
      	   for(itA = Neigh_rna->begin(); itA!=Neigh_rna->end(); itA++){
        		  string tmpstr = NodeA+"_"+*itA+"_mrna";
        		  if(contains(tmpstr,Alledges)) tf_neigh->insert(tmpstr);
        	 }
        }
      	if(Neigh_prot->size()>0){
      	   for(itB = Neigh_prot->begin(); itB!=Neigh_prot->end(); itB++){
            		string tmpstr = NodeA+"_"+*itB+"_prot";
            		string newstr = getEdge(tmpstr,Alledges);
      	        if(!newstr.empty()) ppi_neigh->insert(newstr);
      	   }
        }
      	*Neigh_rna = iter2->second.at(1);
      	*Neigh_prot = iter2->second.at(2);
      	if(Neigh_rna->size()>0){
      	   for(itA = Neigh_rna->begin(); itA!=Neigh_rna->end(); itA++){
      		string tmpstr = NodeB+"_"+*itA+"_mrna";
      	        if(contains(tmpstr,Alledges)) tf_neigh->insert(tmpstr);
      	   }
      	}
      	if(Neigh_prot->size()>0){
      	   for(itB = Neigh_prot->begin(); itB!=Neigh_prot->end(); itB++){
      		string tmpstr = NodeB+"_"+*itB+"_prot";
      		string newstr = getEdge(tmpstr,Alledges);
      		if(!newstr.empty()) ppi_neigh->insert(newstr);
      	   }
      	}
      	NN->push_back(*tf_neigh);
      	NN->push_back(*ppi_neigh);
      	(*NNmap)[curr_edge] = *NN;
 
        delete tf_neigh;
      	delete ppi_neigh;
      	delete Neigh_rna;
      	delete Neigh_prot;
      	delete NN;
      	tf_neigh = NULL;
      	ppi_neigh = NULL;
      	Neigh_rna = NULL;
      	Neigh_prot = NULL;
      	counter++;
    }

}

string getEdge( string str, set<string> * Edges){
     string ret;
     if(contains(str, Edges)) ret = str;
     else{
      	vector<string> vecstr = stringSplit(str,'_');
      	string newstr = (vecstr.at(2)+"_"+vecstr.at(3))+"_"+(vecstr.at(0)+"_"+vecstr.at(1));
      	if(contains(newstr,Edges)) ret = newstr;
     }
     return ret;
}



vector<double> hypergeo(set<string> *vec, set<string> *DE, set<string> *bglist){
      
    vector<double> ret;
    unsigned int k, n1, n2, t;
    int count =0;
    n1 = DE->size();
    n2 = (bglist->size() - n1);
    set<string> overlap = IntersectSet(*vec, *DE);
    k = overlap.size();
    t = vec->size();
    double pval,pval2;
    //P(X>=0) = 1-P(X<0) = 1
    if(k==0) pval=1.0;
    //P(X >= x) = P(X > x-1) for x>0
    boost::math::hypergeometric_distribution<double> hg_dist(t, n1, n1+n2);
    if(k>0) {
      pval = 1 - boost::math::cdf<double>(hg_dist,(k-1));
      count++;
    }
    ret.push_back(pval);
    ret.push_back(t);
    ret.push_back(k);
    ret.push_back(bglist->size());
    ret.push_back(n1);
    return ret;
}

void CreateOutput(string lab,vector<string> PATHlist,int dir, set<string> *geneLIST, set<string> *bgLIST, map<string, set<string> > *PATHEdgemap, map<string, set<string> > *PATHmap, map<string, string> *PATHwayAnnot,int mink){

  string DIR, Annot;
  if(dir==0) DIR = "up";
  if(dir==1) DIR ="down";
  vector<string> collect;
  set<string> memEdges, memGenes, AllG, tmp_setN;
  set<string>::iterator it;
  map<string, string>::iterator a_iter;
  vector<string> tmpVec;
  ofstream outF(lab+"_Enrichment_"+DIR+".txt");
  outF <<"Pathway\tPathAnnotation\tHypergeoPval\tNumGene_inPathway\tNumGenes_inbg\tPropPathway\tNumEdgesformed_inbg\tEnriched_Edgesize"<<endl;
  map<string, set<string> >::iterator M_iter, m_iter;
  
  for(it = bgLIST->begin(); it!= bgLIST->end(); it++){
    tmpVec = stringSplit(*it,'_');
    AllG.insert(tmpVec.at(0));
    AllG.insert(tmpVec.at(2));
  }

  for(int i =0; i<PATHlist.size(); i++){
      string path = PATHlist.at(i);
      M_iter = PATHEdgemap->find(path);
      m_iter = PATHmap->find(path);
      a_iter = PATHwayAnnot->find(path);
      memEdges = M_iter->second;
      memGenes = m_iter->second;
      Annot = a_iter->second;
      tmp_setN = IntersectSet(AllG, memGenes);
      double pp = (double) tmp_setN.size()/ (double) memGenes.size();
      vector<double> hyperpval = hypergeo(&memEdges, geneLIST, bgLIST);
      if(hyperpval[2]>=mink){
            outF<<path<<"\t"<<Annot<<"\t"<<hyperpval.at(0)<<"\t"<< memGenes.size()<<"\t"<<tmp_setN.size()<<"\t"<<pp<<"\t"<<memEdges.size()<< "\t"<<hyperpval.at(2)<<endl;
      }
  }
  outF.close();
}

vector<string> createEdgeOrientedPathwayMap(set<string> *features,set<string> *Alledges,map<string, set<string> > *PATHmap,map<string, set<string> > *PATHmapEdge, double prop,int  minbg){
    
    map<string, set<string> >::iterator piter1;
    set<string>::iterator it;
    set<string> memEdge,tmp_set,tmp_setN, AllG;
    string curr_p;
    vector<string> tmpVec, ret;
    
    for(it = Alledges->begin(); it!= Alledges->end(); it++){
        tmpVec = stringSplit(*it,'_');
        AllG.insert(tmpVec.at(0));
        AllG.insert(tmpVec.at(2));
    }
    for(piter1=PATHmap->begin(); piter1!=PATHmap->end(); piter1++){
        curr_p = piter1->first;
        tmp_set = piter1->second;
        tmp_setN = IntersectSet(AllG, tmp_set);
	      double pp = (double) tmp_setN.size()/ (double) tmp_set.size();
	      if(pp>=prop){
            memEdge = createMEM(&tmp_setN,Alledges);
      	    if(memEdge.size() >= minbg){
      		      (*PATHmapEdge)[curr_p] = memEdge;
      		      ret.push_back(curr_p);
      	    }
      	}	
        memEdge.clear();
    }   
    return ret;
}
            


set<string> createMEM(set<string> *SET,set<string> *ALLedges){
    set<string> ret;
    set<string>::iterator it;
    string edge, edgenew, geneA, geneB;    
    vector<string> vecString;

    for(it = ALLedges->begin(); it!=ALLedges->end(); it++){
      	edge = *it;
      	vecString = stringSplit(edge,'_');
      	geneA = vecString.at(0);
      	geneB = vecString.at(2);
      	if(contains(geneA,SET) & contains(geneB,SET)) ret.insert(edge);
    }
    return ret;
}


void printBGlist(set<string> edges){

    set<string>::iterator it;
    ofstream outF("BGlist.txt");
    for(it=edges.begin(); it!=edges.end(); it++) outF<<*it<<endl;
    outF.close();
}

vector<double> Ztrans(vector<double> vec, double mean, double sd){
    vector<double> ret;
    double tmp;
    for(int k=0; k<vec.size(); k++){
       tmp = vec.at(k) - mean;
       tmp = tmp/sd;
       ret.push_back(tmp);
    }
    return ret;
    ret.clear();
}


void StandardizeFeatures(vector<string> sub, vector<string> genes, int datatype, map<string, subjectClass> &subjectMap){

     map<string, subjectClass>::iterator s_iter;
     map<string, geneClass>::iterator g_iter;
     vector<double> tmp_vec,z_vec;
     string curr_g,curr_sub;    
     double vec_mean, vec_sd;
     string lab;
     if(datatype==1) lab="DNA";
     if(datatype==2) lab="RNA";
     if(datatype==3) lab="PROT";
     ofstream outF("Ztransform_"+lab+".txt");
     outF<<"Gene";
     for(int j=0; j<sub.size(); j++) outF<<"\t"<<sub.at(j);     
     outF<<endl;
     for(int i=0; i<genes.size(); i++){
        outF<<genes.at(i);
        curr_g = genes.at(i);
        tmp_vec = extractGeneVec(datatype,curr_g,sub,subjectMap);
        vec_mean = calMean(tmp_vec);
        vec_sd = calSd(tmp_vec);
        z_vec = Ztrans(tmp_vec,vec_mean,vec_sd);
        for(int j=0; j<z_vec.size(); j++) {
      	    outF<<"\t"<<z_vec.at(j);
      	    curr_sub = sub.at(j);
            s_iter = subjectMap.find(curr_sub);
            g_iter = (s_iter->second.geneMap).find(curr_g);
            if(datatype==1)g_iter->second.addDNA(z_vec.at(j));
            if(datatype==2)g_iter->second.addRNA(z_vec.at(j));
            if(datatype==3)g_iter->second.addProt(z_vec.at(j));
      	}
      	outF<<endl;
        tmp_vec.clear();
        z_vec.clear();
     }
     outF.close();
}


