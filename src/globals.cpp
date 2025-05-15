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

// to deal with different EOL formats
std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case std::streambuf::traits_type::eof():
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

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


vector<string> stringSplit2(string str,string delimeter)
{
     vector<string> splittedString;
     int startIndex = 0;
     int  endIndex = 0;
     while( (endIndex = str.find(delimeter, startIndex)) < str.size() ){
       string val = str.substr(startIndex, endIndex - startIndex);
       splittedString.push_back(val);
       startIndex = endIndex + delimeter.size();
     }
     if(startIndex < str.size()){
       string val = str.substr(startIndex);
       splittedString.push_back(val);
     }
     return splittedString;
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

vector<string> IntersectVec2(vector<string> v1, vector<string> v2, int type){

      vector<string> v3;
      vector<string> v1_new;
      vector<string> v2_new;
      for(int i=0; i<v1.size(); i++) v1_new.push_back((stringSplit2(v1[i],"...")).at(0));
      for(int i=0; i<v1.size(); i++) v2_new.push_back((stringSplit2(v2[i],"...")).at(0));
      sort(v1_new.begin(), v1_new.end());
      sort(v2_new.begin(), v2_new.end());
      string DT, newv3;
      if(type==1) DT = "Z";
      if(type==2) DT = "Y";
      if(type==3) DT = "X";
    
      set_intersection(v1_new.begin(), v1_new.end(), v2_new.begin(), v2_new.end(), back_inserter(v3));
      for(int i=0; i<v3.size(); i++){
	newv3 = v3[i] +"...dat"+DT;
	v3[i] = newv3;
      }

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



bool readUserInput(string file ,string dir, Input_t &UserInput){
    bool ret = true;
    ifstream File;
    File.open(file.c_str(), ios::in);

    if(!File.is_open()){
    	cerr<<"Error! Unable to open file...\n";
	    exit(1);
    }
    string line, d;
    vector<string> v;
    UserInput.directory = dir;
 
    while(!File.eof()){
	    safeGetline(File,line);
	    string start = line.substr(0,1);
      regex expr("#");
	    if(regex_match(start,expr)) continue;

      if(GREP(line,"MAX_BLOCKSIZE")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.blocksize = stoi(d);
	    }
	    if(GREP(line,"DATA_Z")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.dataZ = d;
      }
	    if(GREP(line, "DATA_Y")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.dataY = d;
	    }
	    if(GREP(line,"DATA_X")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.dataX = d;
      }
	    if(GREP(line,"NETWORK_WITHIN")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.withinNet = d;
  	  }
	    if(GREP(line,"NETWORK_BTW")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.betweenNet = d;
	    }
	    if(GREP(line,"SUBTYPE_FILE")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.subtypefile = d;
	    }
	    if(GREP(line,"USE_PRIOR")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      bool d_new;
	      if(d == "false") d_new = false;
	      if(d =="true") d_new = true;
	      UserInput.prior_input = d_new;
	    }
	    if(GREP(line,"PRIOR_FILE")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      UserInput.priorfile = d;
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

	    if(GREP(line,"LOG_TRANSFORM_Z")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      bool d_new;
	      if(d == "false") d_new = false;
	      if(d =="true") d_new = true;
	      UserInput.log_Z = d_new;
	      cerr<<"log DNA: "<<d<<endl;
	    }
	    if(GREP(line,"LOG_TRANSFORM_Y")){
	      v = stringSplit(line,'=');
	      d = trim(v.at(1));
	      bool d_new;
	      if(d == "false") d_new = false;
	      if(d =="true") d_new = true;
	      UserInput.log_Y = d_new;
	      cerr<<"log RNA: "<<d<<endl;
	    }
    	if(GREP(line,"LOG_TRANSFORM_X")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   bool d_new;
    	   if(d == "false") d_new = false;
    	   if(d =="true") d_new = true;
    	   UserInput.log_X = d_new;
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
    	if(GREP(line,"ENRICHMENT")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   bool d_new;
    	   if(d =="false") d_new = false;
    	   if(d =="true") d_new = true;
    	   UserInput.enrichment = d_new;
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
      if(GREP(line,"ANALYZE_Z")){
	       v = stringSplit(line,'=');
	       d = trim(v.at(1));
	       bool d_new;
	       if(d =="false") d_new = false;
	       if(d=="true") d_new = true;
	       UserInput.analyzeZ = d_new;
	    }
     if(GREP(line,"NORMALIZE_ZBY")){
	       v = stringSplit(line,'=');
	       d = trim(v.at(1));
	       string d_new;
	       if(d =="Y"|d=="y") d_new = "Y";
	       if(d=="X"|d=="x") d_new = "X";
	       UserInput.normalize_Zby = d_new;
	    }

      if(GREP(line,"ANALYZE_Y")){
	       v = stringSplit(line,'=');
	       d = trim(v.at(1));
	       bool d_new;
	       if(d =="false") d_new = false;
	       if(d=="true") d_new = true;
	       UserInput.analyzeY = d_new;
	    }
      if(GREP(line,"Z_TRANSFORM_Z")){
	       v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   bool d_new;
    	   if(d =="false") d_new = false;
    	   if(d=="true") d_new = true;
    	   UserInput.ztrans_Z = d_new;
    	}
      if(GREP(line,"Z_TRANSFORM_Y")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   bool d_new;
    	   if(d =="false") d_new = false;
    	   if(d=="true") d_new = true;
    	   UserInput.ztrans_Y = d_new;
    	}
      if(GREP(line,"Z_TRANSFORM_X")){
    	   v = stringSplit(line,'=');
    	   d = trim(v.at(1));
    	   bool d_new;
    	   if(d =="false") d_new = false;
    	   if(d=="true") d_new = true;
    	   UserInput.ztrans_X = d_new;
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
     int countNA =0, counter=0;
     string DT;
     if(datatype==1) DT ="Z";
     if(datatype==2) DT ="Y";
     if(datatype==3) DT ="X";

     for(int h=0; h<indices.size(); h++){         
          int i = indices.at(h);
          string curr_g = genes.at(i);
	  string tag = "ID"+to_string(counter)+"...dat"+DT;
          x_j = MAT[i];
          bool search = false;
	  for(iter = MAP.begin(); iter!=MAP.end(); iter++){
	        set<string> curr_block = iter->first;
	        if(contains(tag, &curr_block)){
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
	          	if(Cv.size() <knnK)  yhat = calMean(getColvec(k, MAT));
	          	s_iterI = subjectMap.find(curr_sub);
	          	s_iterI->second.insertMap(tag, yhat, datatype);
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



vector<string> readFile(string myFile, map<string, set<string> > &subtypeMap, map<string, subjectClass> &subjectMap,Input_t &UserInput,set<string> sub,  int dataType){

    ifstream file;
    file.open(myFile.c_str(), ios::in);
    
    if(!file.is_open()){
        cerr <<"Error! Unable to open file...\n";
        exit(1);
    }
    
    bool logtrans = false;
    string datatype;
    if(dataType==1) {datatype ="Z"; if(UserInput.log_Z) logtrans=true;}
    if(dataType==2) {datatype ="Y"; if(UserInput.log_Y) logtrans=true;}
    if(dataType==3) {datatype ="X"; if(UserInput.log_X) logtrans=true;}

    string line, tag, DT;
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
        safeGetline(file, line);  
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
     	}// close skipheader loop
        v.clear();
        v = stringSplit(line, '\t');
        string curr_g = v.at(0);
        bool keep = true;       
        if(dataType==1) {cerr << "Reading in input data Z..........."<<ctr<<"\r";DT = "Z";}
        if(dataType==2) {cerr << "Reading in input data Y..........."<<ctr<<"\r";DT = "Y";}
        if(dataType==3) {cerr << "Reading in input data X........."<<ctr<<"\r";DT = "X";}     

      	int start_pos=0;
      	int end_pos=0;
      	for(int k=0; k<v_size.size(); k++){
  	    if(k==0) end_pos=0;
  	    else start_pos = end_pos;
  	    end_pos = start_pos + v_size.at(k);
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

	   tag = curr_g+"...dat"+DT;
  	   geneset.push_back(tag);
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
  			       m_iter->second.insertNAgenes(tag, dataType);
  			       d_new = NA_VAL;
  			       if(nainsert) {cid_missing.push_back(counter); nainsert=false;}
  		          }
                          else{
                               d=atof(dat_str.c_str()); 
  			       if(logtrans){
                      		  if(d==0){
  				    d_new=NA_VAL;
  				    m_iter->second.insertNAgenes(tag,dataType);
  			 	    if(nainsert){cid_missing.push_back(counter); nainsert=false;}
  			           }
  			           if(d!=0) d_new =(log(d)/log(2));  
    		                }
  			        if(!logtrans) d_new = d;   
            	           }
                           m_iter->second.insertMap(tag,d_new,dataType);                           
                       } 
  	          }
  	  } //close for loop
          counter++;
        }// close keep loop   
        ctr++;
    }//close while loop
    if(dataType==1) UserInput.cid_Z = cid_missing;
    if(dataType==2) UserInput.cid_Y = cid_missing;
    if(dataType==3) UserInput.cid_X = cid_missing;
    return geneset;
}



set<string> readwithinNetwork(string file, map<string, set<string> >&WITHINMap, map<string, int> &InteractMap,vector<string>X_g, set<string> sub, map<string, subjectClass>&subjectMap){

    set<string> ret;
    ifstream myFile;
    myFile.open(file.c_str(), ios::in);
    if(!myFile.is_open()){
	 cerr<<"Error! Unable to open file...\n";
	 exit(1);
    }
    set<string> members;
    vector<string> v, subjects;
    string line, type, Dir;
    string tagA, tagB;
    int counter=0, ctr=0;
    bool skipHeader=true;
    map<string, set<string> >::iterator m_iter;
    set<string>::iterator it;

    set<string> allS = retrieveKeys(subjectMap);
    set<string> commonS = IntersectSet(allS, sub);
    while(!myFile.eof()){
    	safeGetline(myFile, line);
    	if(line=="") break;
    	if(skipHeader){
    	    skipHeader=false;
    	    continue;
      	}
  	 v = stringSplit(line,'\t');
         string geneA = v.at(0);
         string geneB = v.at(1);
	 tagA = geneA+"...datX";
	 tagB = geneB+"...datX";
      	 type = v.at(2); 
	 if((type!="1") & (type!="-1"))  {cerr<<"Error! Entries in the third columns should only contain 1 or -1. Please try again\n"; exit(1);}
  	 if(containsVec(tagA, &X_g) & containsVec(tagB, &X_g)){
	    subjects = SetToVec(commonS);
            double rho = checkCor(tagA, tagB, subjects, subjectMap);

	    if(rho>=0) Dir = "1";
	    if(rho==0) {Dir="0";ctr++;}	
	    if(rho<0) Dir = "-1";
	    if(Dir==type){
               m_iter = WITHINMap.find(tagA);
     	       if(m_iter == WITHINMap.end()){
  	          members.insert(tagB);
  	          WITHINMap[tagA] = members;
  	          members.clear();
  	       }
               else m_iter->second.insert(tagB);
               m_iter = WITHINMap.find(tagB);
               if(m_iter == WITHINMap.end()){
                 members.insert(tagA);
                 WITHINMap[tagB] = members;
      	         members.clear();
               }
               else m_iter->second.insert(tagA);
  	       string edge = tagA+"___"+tagB;
     	       InteractMap[edge] = stoi(type);
  	       ret.insert(edge);
	    }
	    if(Dir!=type) counter++;
         }
    }
    if(ctr>0) cerr<<ctr<<" zero correlations were calculated.\n"<<endl;

    cerr<<counter<<" edges were removed due to incoherent direction of interaction specified in within-network and the correlation of nodes in the input data.\n"<<endl;
    return ret;
}


double checkCor(string tagA, string tagB, vector<string> sub, map<string, subjectClass>&subjectMap){

    set<string>::iterator it;
    vector<double> vecA, vecB;
    vector<string> tmpstrA, tmpstrB;

    tmpstrA = stringSplit2(tagA,"...");
    tmpstrB = stringSplit2(tagB,"...");

    string geneA = tmpstrA.at(0);
    string dataTypeA =tmpstrA.at(1);

    string geneB = tmpstrB.at(0);
    string dataTypeB = tmpstrB.at(1);

    int DTA, DTB;
    if(dataTypeA=="datX") DTA =3;
    if(dataTypeA=="datY") DTA = 2;
    if(dataTypeA =="datZ") DTA = 1;

    if(dataTypeB=="datX") DTB =3;
    if(dataTypeB=="datY") DTB =2;
    if(dataTypeB=="datZ") DTB =1;

    vecA = extractGeneVec(DTA, tagA, sub, subjectMap);
    vecB = extractGeneVec(DTB, tagB, sub, subjectMap);

    double rho = calPearsonCor(vecA, vecB);
    //if(rho<0)  ret="-1";
    //if(rho>0)  ret="1";


    return rho;
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


void insertPrior(string Priorfile, map<string, subjectClass> &subjectMap){

    set<string> ret;
    ifstream myFile;
    myFile.open(Priorfile.c_str(), ios::in);
    
    map<string, subjectClass>::iterator iter, iter2;

    if(!myFile.is_open()){
      	cerr<<"Error!! Unable to open file...\n";
      	exit(1);
    }

    vector<string> v, headv;
    string line;
    bool skipHeader=true;
    double pp;
    while(!myFile.eof()){	
      	safeGetline(myFile, line);
      	if(line=="") break;
      	if(skipHeader){
	   headv = stringSplit2(line,"\t");
	   headv.erase(headv.begin()); //remove first element
      	   skipHeader=false;
      	   continue;
      	}
      	v = stringSplit2(line,"\t");
      	string sub = v.at(0);
      	v.erase(v.begin());
      	iter= subjectMap.find(sub);
       	if(iter != subjectMap.end()){
		for(int j=0; j<v.size(); j++){
		    string tmpstr = v.at(j);
		    string tmphead = headv.at(j);
		    replace(tmphead.begin(), tmphead.end(), ' ', '_' );
		    if((tmpstr!="") &(tmpstr!="NA")) pp =atof(tmpstr.c_str()); 
		    if((tmpstr=="")|(tmpstr=="NA")) pp = NA_VAL;
		    iter->second.insertPrior(tmphead, pp);
		    //cerr<<headv.at(j)<<"\t"<<tmphead<<"\t"<<pp<<endl;
		}
	}
    }
    vector<string> to_rm;
    bool prior_exist;
    for(iter2=subjectMap.begin(); iter2!=subjectMap.end(); iter2++){
	string currK = iter2->first; 
	prior_exist = iter2->second.checkPrior();
	if(!prior_exist) to_rm.push_back(currK);
	//if(prior_exist){
	//	vector<string> GRP = iter2->second.getGrp();
	//	vector<double> PIE = iter2->second.getPriorVec();
	//	for(int i=0; i<GRP.size(); i++) cerr<<GRP[i]<<"\t";
	//	cerr<<endl;
	//	for(int i=0; i<PIE.size(); i++) cerr<<PIE[i]<<"\t";
	//	cerr<<endl;
	//}
    }
    if(to_rm.size()>0){
       for(int r=0; r<to_rm.size(); r++) subjectMap.erase(to_rm.at(r));
       cerr<<to_rm.size()<<" samples were further removed as prior probabilities were not available.\n"<<endl;
    }
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
    map<string, string>::iterator f_iter;
    while(!myFile.eof()){	
      	safeGetline(myFile, line);
      	if(line=="") break;
      	if(skipHeader){
      	   skipHeader=false;
      	   continue;
      	}
      	v = stringSplit2(line,"\t");
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
        else{iter->second.insert(mem);}
      	members.clear();
    }
    string maxPath = maximumG(PATHmap);
    cerr<<"There are a total of "<< PATHmap.size() <<" number of pathways."<<endl;
    return ret;
}


set<string>  readbetweenNetwork(string file, map<string, set<string> > &BTWMap, map<string, int> &InteractMap, vector<string> target_g, vector<string> X_g, set<string> sub, map<string, subjectClass>&subjectMap){
    
    set<string> ret;
    ifstream myFile;
    myFile.open(file.c_str(), ios::in);
    
    if(!myFile.is_open()){
        cerr<<"Error! Unable to open file...\n";
        exit(1);
    }
    set<string> members;
    vector<string> v, subjects;
    string line, type, Dir;
    bool skipHeader=true;
    map<string, set<string> >::iterator m_iter;
    int counter=0;


    set<string> allS = retrieveKeys(subjectMap);
    set<string> commonS = IntersectSet(allS, sub);
    while(!myFile.eof()){
        safeGetline(myFile, line);
        if(line=="") break;
        if(skipHeader){
            skipHeader=false;
            continue;
        }
        v = stringSplit2(line,"\t");
        string key = v.at(0);
        string mem = v.at(1);
      	type = v.at(2); 
	string tagK = key+"...datX";
	string tagM = mem+"...datY";
 
	if((type!="1") & (type!="-1"))  {cerr<<"Error! Entries in the third columns should only contain 1 or -1. Please try again\n";exit(1);}
      	if(containsVec(tagK, &X_g) & containsVec(tagM, &target_g)){
               
              subjects = SetToVec(commonS);
              double rho = checkCor(tagK, tagM, subjects, subjectMap);
	      if(rho>0) Dir="1";
	      if(rho<0) Dir="-1";

              if(Dir==type){
                m_iter = BTWMap.find(tagK);
                if(m_iter == BTWMap.end()){
                  members.insert(tagM);
                  BTWMap[tagK] = members;
                }
                if(m_iter!=BTWMap.end()) m_iter->second.insert(tagM);
      	        string edge = tagK+"___"+tagM;
      	        InteractMap[edge] = stoi(type);
      	        ret.insert(edge);
                members.clear();
	      }
	      if(Dir!=type) counter++;
	}

    }
    cerr<< counter<<" edges were removed due to incoherent direction of interation specified in between-network file and the correlation found between the nodes in the input data\n"<<endl;
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
        	safeGetline(file,line);
        	if(line=="") break;
        	if(skipHeader){
        	    skipHeader=false;
           	    continue;
          }
         	v = stringSplit2(line,"\t");
        	string curr_s = v.at(0);
		string curr_t = v.at(1);
		replace(curr_t.begin(), curr_t.end(), ' ', '_' );
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

double calPearsonCor(vector<double> x, vector<double> y){

   double ret = 0.0, num= 0.0;
   double denA=0.0, denB=0.0;
   double mu_x = calMean(x);
   double mu_y = calMean(y);
   
   for(int i=0; i<x.size(); i++){
     if((x.at(i)!=NA_VAL) & (y.at(i)!=NA_VAL)){
	num += ((x.at(i)-mu_x)*(y.at(i)-mu_y));
     	denA += pow((x.at(i)-mu_x),2);
     	denB += pow((y.at(i) -mu_y),2);
     }
   }
   double den = sqrt(denA*denB);
   ret = num/den;
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
        //if(value<0) cerr<<"Neg Pval, pls check!"<<endl;
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
     if(datatype==1) data = "Z";
     if(datatype==2) data = "Y";
     if(datatype==3) data = "X";

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
     cerr<<"In the dataframe "<<data<<", there are "<<ret.size()<<" features with at least "<<mink <<" non-missing samples across each of the subtypes."<<endl;
     return ret;
}




void outputData(vector<string> subjects,vector<string> genes,  map<string, subjectClass> &subjectMap, int datatype, string dir){

    map<string, geneClass>::iterator g_iter;
    map<string, subjectClass>::iterator iter;

    string LAB;
    if(datatype==1)LAB = "dataZ";
    if(datatype==2)LAB = "dataY";
    if(datatype==3)LAB = "dataX";
   
    ofstream outF(dir+"imputedData_"+LAB+".txt");
    vector<string> GG;
    outF<<"Gene";
    for(int i=0; i<subjects.size(); i++) outF<<"\t"<<subjects.at(i);
    outF<<"\n";

    for(int i=0; i<genes.size(); i++){	
      	string curr_g = genes.at(i);
        GG = stringSplit2(curr_g,"...");
	outF<<GG.at(0);
      	for(int j=0; j<subjects.size(); j++){
      	   iter = subjectMap.find(subjects.at(j));
      	   g_iter = iter->second.geneMap.find(curr_g);
      	   double dd = g_iter->second.getReadings(datatype);
      	   if(dd==NA_VAL) outF<<"\t"<<"";
      	   if(dd!=NA_VAL) outF<<"\t"<<dd;
      	}
	outF<<"\n";
    }
    outF.close();
} 


vector<double> fillPAMmap2(bool useZ, string normby, set<string> *Alledges,vector<string> *sub,map<string,int> *InteractMap, map<string, set<string> > *subtypeMap, map<string, subjectClass> *subjectMap, map<string, PAM_t> &PAMmap, map<string, vector<set<string> > > *FDneighborMap, string dir){
     map<string, subjectClass>::iterator s_iter;
     map<string, PAM_t>::iterator p_iter, p_iter2;
     map<string, set<string> >::iterator iter;
     map<string, geneClass>::iterator g_iterA, g_iterB;
     map<string, vector<set<string> > >::iterator iter1, iter2;
     map<string, int>::iterator i_iter;

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

     ofstream outFile(dir+"Expressiondata_edges.txt");
     outFile<<"Edge";

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
     string geneA, geneB, typeA,typeB, NodeA, NodeB;
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
      	vec_g = stringSplit2(curr_e,"___");
        NodeA = vec_g.at(0);
	    NodeB = vec_g.at(1);
      	geneA = stringSplit2(NodeA,"...").at(0);
      	typeA =  stringSplit2(NodeA,"...").at(1);
      	geneB =  stringSplit2(NodeB,"...").at(0);
      	typeB =  stringSplit2(NodeB,"...").at(1);
      	p_iter->second.identifier = curr_e;
      	p_iter->second.geneA = NodeA;
      	p_iter->second.geneB = NodeB;
      	if(typeB=="datX") datatype = "WITHIN";
      	if(typeB=="datY") datatype = "BTW";
      	p_iter->second.dataT = datatype;
      	if(datatype=="WITHIN"){datatypeA = 3; datatypeB=3;}
      	if(datatype=="BTW"){datatypeA = 3; datatypeB=2;}
         
      	for(int j=0; j<lab.size(); j++){
      	     curr_grp = SubGroups.at(j);
      	     vector<double> grp_x;
      	     for(it2 = curr_grp.begin(); it2!=curr_grp.end(); it2++){
            		 string curr_sub = *it2;
            		 s_iter = subjectMap->find(curr_sub);
            		 g_iterA = s_iter->second.geneMap.find(NodeA);
            		 g_iterB = s_iter->second.geneMap.find(NodeB);
            		 xA = g_iterA->second.getReadings(datatypeA);
            		 xB = g_iterB->second.getReadings(datatypeB);
			         xC = 0;
            		 if(useZ & normby=="Y") xC = g_iterB->second.getReadings(1);
            		 if(useZ & normby=="X") xC = g_iterA->second.getReadings(1);
            		 if(datatype=="BTW"){
            		   if(interactINT== 1){
                  		if(first){countPos++; first=false;}
                  		rr = xA+xB-xC;
            		   }
            		   if(interactINT== -1){
            		        if(first){countNeg++; first = false;}
            		        if(!useZ) rr = xA - xB;
            		        if(useZ & normby=="Y") rr = xA-xB+xC;
            		        if(useZ & normby=="X") rr = (xA-xC)-xB;
            		   }
            		 }          		 
            		 if(datatype=="WITHIN") {
            		   if(interactINT== 1){
                  		if(first){countPos++; first=false;}
                  		rr = xA+xB;
            		   }
            		   if(interactINT== -1){
            		        if(first){countNeg++; first = false;}
            		        rr = xA - xB;
            		   }
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
	    x_vec.clear();
        GRPmean.clear();
     }

     outFile.close();

     cerr<<"\nThere are a total of "<<countPos<<" positive-signed edges and "<<countNeg<<" negative-signed edges created.\n";

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
     set<string> BTW_N, WITHIN_N, Neigh_Y, Neigh_X;
     vector<double> collect_dik;
     vector<double> collect_dik_BTW;
     vector<double> collect_dik_WITHIN;
     cerr<<"size of PAM map: "<<PAMmap.size()<<endl;
     cerr<<"size of FD map: "<<FDneighborMap->size()<<endl;

     for(p_iter = PAMmap.begin(); p_iter!=PAMmap.end(); p_iter++){
            curr_e = p_iter->first;
    	    vec_g = stringSplit2(curr_e,"___");
            NodeA = vec_g.at(0);
            NodeB = vec_g.at(1);

    	    geneA = (stringSplit2(NodeA,"...")).at(0);
    	    typeA = (stringSplit2(NodeA,"...")).at(1);
    	    geneB = (stringSplit2(NodeB,"...")).at(0);
    	    typeB = (stringSplit2(NodeB,"...")).at(1);
    	    datatype = p_iter->second.dataT;
    	    if(datatype=="WITHIN"){datatypeA = 3; datatypeB=3;}
    	    if(datatype=="BTW"){datatypeA = 3; datatypeB=2;}
        
           iter1 = FDneighborMap->find(NodeA);
           iter2 = FDneighborMap->find(NodeB);
	   if(iter1 == (FDneighborMap->end())) cerr<<"No such Edge: "<<NodeA<<"\t"<<curr_e<<endl;
	   if(iter2 == (FDneighborMap->end())) cerr<<"No such Edge: "<<NodeB<<"\t"<<curr_e<<endl;

           if(iter1 != (FDneighborMap->end())){Neigh_Y = iter1->second.at(1);Neigh_X = iter1->second.at(2);}
           if(Neigh_Y.size()>0){
               for(it = Neigh_Y.begin(); it!=Neigh_Y.end(); it++){
                   string tmpstr = NodeA+"___"+*it;
                   if(contains(tmpstr,Alledges))BTW_N.insert(tmpstr);
               }
           }
           if(Neigh_X.size()>0){
             for(it2 = Neigh_X.begin(); it2!=Neigh_X.end(); it2++){
                 string tmpstr = NodeA+"___"+*it2;
                 string newstr = getEdge(tmpstr,Alledges);
                 if(!newstr.empty()) WITHIN_N.insert(newstr);
             }
         }
          if(iter2 != (FDneighborMap->end())){Neigh_Y = iter2->second.at(1);Neigh_X = iter2->second.at(2);}
         if(Neigh_Y.size()>0){
             for(it = Neigh_Y.begin(); it!=Neigh_Y.end(); it++){
                 string tmpstr = NodeB+"___"+*it;
                 if(contains(tmpstr,Alledges)) BTW_N.insert(tmpstr);
             }
         }
         if(Neigh_X.size()>0){
             for(it2 = Neigh_X.begin(); it2!=Neigh_X.end(); it2++){
                 string tmpstr = NodeB+"___"+*it2;
                 string newstr = getEdge(tmpstr,Alledges);
                 if(!newstr.empty()) WITHIN_N.insert(newstr);
             }
         }

         int NumBTW = BTW_N.size();
         int NumWITHIN = WITHIN_N.size();
         p_iter->second.NumFD_btw = NumBTW;
         p_iter->second.NumFD_within = NumWITHIN;
         int totalNeighbors = BTW_N.size() + WITHIN_N.size();
         for(int ss =0; ss<lab.size(); ss++){
              curr_dik = p_iter->second.d_ik.at(ss);
              for(it = BTW_N.begin(); it!=BTW_N.end(); it++){
                  p_iter2 = PAMmap.find(*it);
                  collect_dik.push_back((p_iter2->second.d_ik)[ss]);
                  collect_dik_BTW.push_back((p_iter2->second.d_ik)[ss]);
              }
              for(it2 = WITHIN_N.begin(); it2!=WITHIN_N.end(); it2++){
                  p_iter2 = PAMmap.find(*it2);
                  //if(p_iter2==PAMmap.end()) cerr<<"Edge not found in map!"<<endl;
                  collect_dik.push_back((p_iter2->second.d_ik)[ss]);
                  collect_dik_WITHIN.push_back((p_iter2->second.d_ik)[ss]);
               }
               int NumAgree =0;
               double BTW_dikbar=0, WITHIN_dikbar=0;
               if(collect_dik_BTW.size()>0){
                    BTW_dikbar = calMean(collect_dik_BTW);
                    p_iter->second.FDbtw_dik_bar.push_back(BTW_dikbar);
                    for(int h=0; h<collect_dik_BTW.size(); h++){
                        if(collect_dik_BTW.at(h)<0 & curr_dik<0) NumAgree++;
                        if(collect_dik_BTW.at(h)>0 & curr_dik>0) NumAgree++;
                    }
                }
                if(collect_dik_WITHIN.size()>0){
                    WITHIN_dikbar = calMean(collect_dik_WITHIN);
                    p_iter->second.FDwithin_dik_bar.push_back(WITHIN_dikbar);
                    for(int h=0; h<collect_dik_WITHIN.size();h++){
                        if(collect_dik_WITHIN.at(h)<0 & curr_dik<0) NumAgree++;
                        if(collect_dik_WITHIN.at(h)>0 & curr_dik>0) NumAgree++;
                    }
                }

                double tmp_prop = (double)NumAgree/(double)(totalNeighbors);
	              double wtDikbar = (NumBTW*BTW_dikbar)+(NumWITHIN*WITHIN_dikbar);
            		wtDikbar = wtDikbar/(NumBTW+NumWITHIN);
            		if((NumBTW+NumWITHIN)==0) wtDikbar =0;
            		double phi_f = exp((tmp_prop-0.5)/0.2);
            		double phi = (2*phi_f)/(1+phi_f);
            		if(tmp_prop<0.5|(curr_dik*wtDikbar)<0) phi=0.0;
            		double dik_star = curr_dik + phi*wtDikbar;
            		if(NumBTW+NumWITHIN==0) dik_star = curr_dik;
            		//if(dik_star!=dik_star) cerr<<"dikstar problem "<< dik_star<<"\t"<<NumBTW<<"\t"<<NumWITHIN<<endl;
            		if(abs(dik_star)>maxTHRES.at(ss)) maxTHRES.at(ss) = abs(dik_star);
                p_iter->second.prop_agree.push_back(tmp_prop);
		            p_iter->second.dikstar_new.push_back(dik_star);
                collect_dik_BTW.clear();
                collect_dik_WITHIN.clear();
                collect_dik.clear();
            }//close subtypemap
      	    BTW_N.clear();
      	    WITHIN_N.clear();
      	    Neigh_Y.clear();
      	    Neigh_X.clear();
      }//close tmp_pammap
      return maxTHRES;
}


vector<double> fillPAMmap(bool useZ, string normby,set<string> *Alledges,vector<string> *sub,map<string, int> *InteractMap, map<string, set<string> > *subtypeMap, map<string, subjectClass> *subjectMap, map<string, PAM_t> &PAMmap, map<string, vector<set<string> > > *FDneighborMap){

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
    string geneA, geneB, typeA,typeB, NodeA, NodeB;
    int datatypeA, datatypeB;
    vector<string> vec_g;
    double xA,xB,xC,rr;
    
    for(it = Alledges->begin(); it != Alledges->end(); it++){
	bool first = true;
        curr_e = *it;
        PAM_t *pam_struct = new PAM_t();
        PAMmap[curr_e] = *pam_struct;
        p_iter = PAMmap.find(curr_e);
        delete pam_struct;
        
      	i_iter = InteractMap->find(curr_e);
      	interactINT = i_iter->second;
        
        vec_g = stringSplit2(curr_e,"___");
            NodeA = vec_g.at(0);
            NodeB = vec_g.at(1);
    	    geneA = (stringSplit2(NodeA,"...")).at(0);
    	    typeA = (stringSplit2(NodeA,"...")).at(1);
    	    geneB = (stringSplit2(NodeB,"...")).at(0);
    	    typeB = (stringSplit2(NodeB,"...")).at(1);
        p_iter->second.identifier = curr_e;
        p_iter->second.geneA = NodeA;
        p_iter->second.geneB = NodeB;
        //cerr<<"typeB: "<<typeB<<endl;
        if(typeB=="datX") datatype = "WITHIN";
        if(typeB=="datY") datatype = "BTW";
        p_iter->second.dataT = datatype;
        if(datatype=="WITHIN"){datatypeA = 3; datatypeB=3;}
        if(datatype=="BTW"){datatypeA = 3; datatypeB=2;}
        for(int j=0; j<lab.size(); j++){
            curr_grp = SubGroups.at(j);
            vector<double> grp_x;
            for(it2 = curr_grp.begin(); it2!=curr_grp.end(); it2++){
                string curr_sub = *it2;
                s_iter = subjectMap->find(curr_sub);
                g_iterA = s_iter->second.geneMap.find(NodeA);
                g_iterB = s_iter->second.geneMap.find(NodeB);
                xA = g_iterA->second.getReadings(datatypeA);
                xB = g_iterB->second.getReadings(datatypeB);
		xC = 0;
            	if(useZ & normby=="Y") xC = g_iterB->second.getReadings(1);
            	if(useZ & normby=="X") xC = g_iterA->second.getReadings(1);
            	if(datatype=="BTW"){
            		   if(interactINT== 1){
                  		if(first){countPos++; first=false;}
                  		rr = xA+xB-xC;
            		   }
            		   if(interactINT== -1){
            		        if(first){countNeg++; first = false;}
            		        if(!useZ) rr = xA - xB;
            		        if(useZ & normby=="Y") rr = xA-xB+xC;
            		        if(useZ & normby=="X") rr = (xA-xC)-xB;
            		   }
            	}          		 
            	if(datatype=="WITHIN") {
            		   if(interactINT== 1){
                  		if(first){countPos++; first=false;}
                  		rr = xA+xB;
            		   }
            		   if(interactINT== -1){
            		        if(first){countNeg++; first = false;}
            		        rr = xA - xB;
            		   }
		}                
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
    set<string> BTW_N, WITHIN_N, Neigh_Y, Neigh_X;
    vector<double> collect_dik;
    vector<double> collect_dik_BTW;
    vector<double> collect_dik_WITHIN;
    for(p_iter = PAMmap.begin(); p_iter!=PAMmap.end(); p_iter++){
        curr_e = p_iter->first;
        vec_g = stringSplit2(curr_e,"___");
           NodeA = vec_g.at(0);
            NodeB = vec_g.at(1);
    	    geneA = (stringSplit2(NodeA,"...")).at(0);
    	    typeA = (stringSplit2(NodeA,"...")).at(1);
    	    geneB = (stringSplit2(NodeB,"...")).at(0);
    	    typeB = (stringSplit2(NodeB,"...")).at(1);
        datatype = p_iter->second.dataT;
        if(datatype=="WITHIN"){datatypeA = 3; datatypeB=3;}
        if(datatype=="BTW"){datatypeA = 3; datatypeB=2;}
                
        iter1 = FDneighborMap->find(NodeA);
        iter2 = FDneighborMap->find(NodeB);
        Neigh_Y = iter1->second.at(1);
        Neigh_X = iter1->second.at(2);
        if(Neigh_Y.size()>0){
            for(it = Neigh_Y.begin(); it!=Neigh_Y.end(); it++){
                string tmpstr = NodeA+"___"+*it;
                if(contains(tmpstr,Alledges)){
                    BTW_N.insert(tmpstr);  
                }
            }
        }
        if(Neigh_X.size()>0){
            for(it2 = Neigh_X.begin(); it2!=Neigh_X.end(); it2++){
                string tmpstr = NodeA+"___"+*it2;
                string newstr = getEdge(tmpstr,Alledges);
                if(!newstr.empty()){
                    WITHIN_N.insert(newstr);
                }
            }
        }
        Neigh_Y = iter2->second.at(1);
        Neigh_X = iter2->second.at(2);
        if(Neigh_Y.size()>0){
            for(it = Neigh_Y.begin(); it!=Neigh_Y.end(); it++){
                string tmpstr = NodeB+"___"+*it;
                if(contains(tmpstr,Alledges)) BTW_N.insert(tmpstr);
            }
        }
        if(Neigh_X.size()>0){
            for(it2 = Neigh_X.begin(); it2!=Neigh_X.end(); it2++){
                string tmpstr = NodeB+"___"+*it2;
                string newstr = getEdge(tmpstr,Alledges);
                if(!newstr.empty()) WITHIN_N.insert(newstr);
            }
        }
        int NumBTW = BTW_N.size();
        int NumWITHIN = WITHIN_N.size();
        p_iter->second.NumFD_btw = NumBTW;
        p_iter->second.NumFD_within = NumWITHIN;
        int totalNeighbors = BTW_N.size() + WITHIN_N.size();
        
        for(int ss =0; ss<lab.size(); ss++){
            curr_dik = p_iter->second.d_ik.at(ss);
            for(it = BTW_N.begin(); it!=BTW_N.end(); it++){
                p_iter2 = PAMmap.find(*it);
                collect_dik.push_back((p_iter2->second.d_ik)[ss]);
                collect_dik_BTW.push_back((p_iter2->second.d_ik)[ss]);
            }
            for(it2 = WITHIN_N.begin(); it2!=WITHIN_N.end(); it2++){
                p_iter2 = PAMmap.find(*it2);
                collect_dik.push_back((p_iter2->second.d_ik)[ss]);
                collect_dik_WITHIN.push_back((p_iter2->second.d_ik)[ss]);
            }
            
            int NumAgree =0;
            double BTW_dikbar=0, WITHIN_dikbar=0;
            if(collect_dik_BTW.size()>0){
                BTW_dikbar = calMean(collect_dik_BTW);
                p_iter->second.FDbtw_dik_bar.push_back(BTW_dikbar);
                for(int h=0; h<collect_dik_BTW.size(); h++){
                    if(collect_dik_BTW.at(h)<0 & curr_dik<0) NumAgree++;
                    if(collect_dik_BTW.at(h)>0 & curr_dik>0) NumAgree++;
                }
            }
            if(collect_dik_WITHIN.size()>0){
                WITHIN_dikbar = calMean(collect_dik_WITHIN);
                p_iter->second.FDwithin_dik_bar.push_back(WITHIN_dikbar);
                for(int h=0; h<collect_dik_WITHIN.size();h++){
                    if(collect_dik_WITHIN.at(h)<0 & curr_dik<0) NumAgree++;
                    if(collect_dik_WITHIN.at(h)>0 & curr_dik>0) NumAgree++;
                }
            }
            
            double tmp_prop = (double)NumAgree/(double)(totalNeighbors);
            double wtDikbar = (NumBTW*BTW_dikbar)+(NumWITHIN*WITHIN_dikbar);
            wtDikbar = wtDikbar/(NumBTW+NumWITHIN);
            if((NumBTW+NumWITHIN)==0) wtDikbar =0;
            double phi_f = exp((tmp_prop-0.5)/0.2);
            double phi = (2*phi_f)/(1+phi_f);
            if(tmp_prop<0.5|(curr_dik*wtDikbar)<0) phi=0.0;
            double dik_star = curr_dik + phi*wtDikbar;
            if(NumBTW+NumWITHIN==0) dik_star = curr_dik;
            if(abs(dik_star)>maxTHRES.at(ss)) maxTHRES.at(ss) = abs(dik_star);
            p_iter->second.prop_agree.push_back(tmp_prop);
            p_iter->second.dikstar_new.push_back(dik_star);
            collect_dik_BTW.clear();
            collect_dik_WITHIN.clear();
            collect_dik.clear();
        }//close subtypemap
        counter2++;
        BTW_N.clear();
        WITHIN_N.clear();
      	Neigh_Y.clear();
      	Neigh_X.clear();
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

vector<int> whichMin(vector<double> vec){
    int pos=0;
    vector<int> POS;
    double curr_min = vec.at(0);
    for(int i=0; i<vec.size(); i++){
      	if(vec.at(i)<curr_min){
       	    pos =i;
      	    curr_min = vec.at(i);
      	}
    }
    POS.push_back(pos);
    for(int i=0; i<vec.size(); i++){
      	if((vec.at(i)==curr_min) &(i!=pos)) POS.push_back(i);
    }
    return POS;
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

set<string> VecToSet(vector<string> VEC){
	set<string> ret;
	for(int i=0; i<VEC.size();i++) ret.insert(VEC.at(i));
        return ret;
}

void CleanupMap(vector<string> subjects,map<string, set<string> > &subtypeMap){

     
     map<string, set<string> >::iterator iter;
     set<string> members, commonMem;
     set<string> subjectsnew = VecToSet(subjects);
     string grp;
     for(iter = subtypeMap.begin(); iter!=subtypeMap.end(); iter++){
          grp = iter->first;
	  members = iter->second;
	  commonMem = IntersectSet(members, subjectsnew);
          iter->second = commonMem;
     }
}

set<string> createNeighborMap(string dir, vector<string> *X_g, vector<string> *Y_g, map<string, set<string> > *BTWmap, map<string, set<string> > *TargetMap, map<string, set<string> > *WITHINmap, map<string, vector<set<string> > >&FDNeighborMap){

     map<string, set<string> >::iterator BTW_iter, WITHIN_iter, target_iter;
     set<string> gene_set,Y_btw,X_btw,X_within;
     set<string> parentNodes, singletons, setY,setX;
     string curr_g;
     set<string>::iterator it, it2;
     int num_Y, num_X, curr_max=0;
     string maxNode, genename;
     map<string, vector<set<string> > >::iterator iter;
     set<string> features;
     vector<set<string> > Neighbors;

     for(int i=0; i< X_g->size(); i++) setX.insert(X_g->at(i));
     for(int i=0; i< Y_g->size(); i++) setY.insert(Y_g->at(i));

     // in the map, element zero is the feature tag, element one is the FDneighbors from data Y (between) and element two is the FDneighbors from data X (within)
    // cerr<<"Size of X_g: "<<X_g->size()<<"\t"<< "Size of Y_g: "<<Y_g->size()<<endl;

     //for data Y, first degree neighbors are those that are targeted by data X and second element should be empty
     for(int i =0; i< Y_g->size(); i++){
        curr_g = Y_g->at(i);
        gene_set.insert(curr_g);
        Neighbors.push_back(gene_set);
        Neighbors.push_back(Y_btw); //empty
        target_iter = TargetMap->find(curr_g);
      	if(target_iter!=TargetMap->end()) parentNodes = target_iter->second;
        X_btw = IntersectSet(parentNodes,setX);
        if(X_btw.size()==0) singletons.insert(curr_g);
        Neighbors.push_back(X_btw);
        if(X_btw.size()>0) FDNeighborMap[curr_g] = Neighbors;
        if(X_btw.size()>curr_max){
             curr_max = X_btw.size();
             maxNode = curr_g ;
        }

        gene_set.clear();
        X_btw.clear();
        Neighbors.clear();
     }

     //for data X, first degree neighbors are all possible within X or between X and Y
     for(int i=0; i< X_g->size(); i++){
          curr_g = X_g->at(i);
          gene_set.insert(curr_g);
          Neighbors.push_back(gene_set);
          BTW_iter = BTWmap->find(curr_g);
          if(BTW_iter!=BTWmap->end()){
              set<string> tmp_children = BTW_iter->second;
              Y_btw = IntersectSet(tmp_children, setY);
          }
          Neighbors.push_back(Y_btw);	
          WITHIN_iter = WITHINmap->find(curr_g);

          if(WITHIN_iter!=WITHINmap->end()){
               set<string> interact = WITHIN_iter->second;
               X_within = IntersectSet(interact,setX);
          }
          Neighbors.push_back(X_within);
          if(Y_btw.size()==0 & X_within.size()==0) singletons.insert(curr_g);
          if((Y_btw.size() + X_within.size())>0) FDNeighborMap[curr_g] = Neighbors;
          if((Y_btw.size()+X_within.size())>curr_max){
              curr_max = Y_btw.size()+X_within.size();
              maxNode = curr_g;
          }
          gene_set.clear();
          X_within.clear();
          Y_btw.clear();
          Neighbors.clear();
       }     
       iter = FDNeighborMap.find(maxNode);

       num_Y = (iter->second).at(1).size();
       num_X = (iter->second).at(2).size();
 
       cerr<<"Generating output for information on the first-degree neighbours of all features in the data...\n"<<endl;
       cerr<<"\nThere is a total of "<<FDNeighborMap.size()<<" features with between and/or within network interactions and "<<singletons.size()<<" features with no interactions and were removed."<<endl;
       cerr<<maxNode<<" is the feature with the most first-degree neighbours, "<<num_Y<<" between- and "<< num_X<<" within-data interactions.\n"<<endl;

       //output the first degree neighbor data
       ofstream outFile(dir+"Features_Neighbors.txt");
       outFile<<"Feature\tGene\tType\tNeighbor_inData\tNumNeigh_inData"<<endl;
  
       string currF, currG, currT;
       set<string> Y_set, X_set;
       vector<string> mem, tt;

       for(iter = FDNeighborMap.begin(); iter!=FDNeighborMap.end(); iter++){
         currF = iter->first;
	 features.insert(currF);
         vector<string> tmpVec = stringSplit2(currF,"...");
         genename = tmpVec.at(0);
         currT = tmpVec.at(1);        
	 Y_set = (iter->second).at(1);
         X_set = (iter->second).at(2);
         if(Y_set.size()>0){
	     for(it = Y_set.begin(); it!=Y_set.end(); it++){
		tt = stringSplit2(*it,"...");
		string tmp = (tt.at(0))+"...datY";
        	mem.push_back(tmp);
             }
	 }
      	 if(X_set.size()>0){
      	    for(it = X_set.begin(); it!= X_set.end(); it++){
		tt = stringSplit2(*it,"...");
		string tmp = (tt.at(0))+"...datX";
       	        mem.push_back(tmp);
      	    }
      	 }
      	 string str_new = concatenate(mem,';');
      	 outFile<<currF<<"\t"<<genename<<"\t"<<currT<<"\t"<<str_new<<"\t"<<mem.size()<<endl;
         mem.clear();
         Y_set.clear();
         X_set.clear();
     }
     outFile.close();
     return features;
}


void ReverseMap(map<string, set<string> >&Map1, map<string, set<string> >&Map2){

      map<string, set<string> >::iterator iter1, iter2;
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
}

double computeL2norm(double dd, vector<double>vec){
    double ret=0;
    for(int i=0; i<vec.size(); i++) ret += pow((vec.at(i) - dd),2);
    ret = ret/(vec.size());
    return ret;
}


double CVKfold(string dir,bool useZ, bool usePrior, string normby, int Kfold,set<string> *Alledges,vector<string> *sub, int thres_size,map<string, int> *InteractMap, map<string, vector<set<string> > > *NeighborMap,  map<string, set<string> > *subtypeMap, map<string, subjectClass> *subjectMap, map<string, PAM_t> *PAMmap){

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
        vector<double> cvthres = fillPAMmap(useZ,normby, Alledges,&train_grp,InteractMap, subtypeMap, subjectMap, tmp_map, NeighborMap);
        BIGmap[key] = tmp_map;
        
        for(int j=0; j<grplabel.size(); j++){
            if(cvthres.at(j)>maxthres.at(j)) maxthres.at(j) = cvthres.at(j);
        }
        tmp_map.clear();
        train_grp.clear();
        curr_grp.clear();
        key++;
    }

    double Avgthres = calMean(maxthres);
    vector<double> thresholdparam;
    double space = Avgthres/(double)(thres_size-1);

    for(int i=0; i<thres_size; i++) thresholdparam.push_back(i*space);

    ofstream outF(dir+ "CVerrors.txt");
    outF<<"Threshold\tCVerror";
    for(int j=0; j<grplabel.size(); j++) outF<<"\tCVerror_"+grplabel.at(j);
    outF<<"\tMeanCV_NumEdgesSelected\tNumEdgesSelected"<<endl;
  
    vector<double> errors, CVerrors, curr_grperror;
    vector<vector<double> > GRPCVerrors;
    double dij_new, dij, xik_new, xik, si, so, curr_thres;

    map<string, PAM_t> tmp_pammap;
    double curr_dik;
    string curr_E, newedge, curr_sub, geneA,geneB, grpClass;

    map<string, int>::iterator i_iter;
    map<string, geneClass>::iterator g_iter, g_iter2;	
    vector<vector<int> > Grpgenesurv;
    vector<int> Binvec, index;
    vector<double> geneSurv_vec, vec_dik;
    int genedied, ccc; 
    int groupsize, typeA, typeB, interactINT, gene_counter;
    double x_star,x_starA, x_starB, x_starC, tmp_score;
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
            
            		    	if(datatype=="BTW") {typeA = 3;typeB=2;}
            		    	if(datatype=="WITHIN"){typeA=3;typeB=3;}
      
      	           		g_iter = sub_iter->second.geneMap.find(geneA);
                    		x_starA = g_iter->second.getReadings(typeA);
      	            		g_iter2 = sub_iter->second.geneMap.find(geneB);
      	            		x_starB = g_iter2->second.getReadings(typeB);
				x_starC = 0;
      	            		if(useZ & normby=="Y") x_starC = g_iter2->second.getReadings(1);
      	            		if(useZ & normby=="X") x_starC = g_iter->second.getReadings(1);
				//if(t==0 & gene_counter==0) {
            		       	 //     GRPsize = piter->second.GRPsize;
                 		//      prior = (double)1/(double) GRPsize.size();
                 	  	//}
                   		double si = piter->second.s_i;
                   		double s0 = piter->second.s_o;

            		   	if(datatype=="BTW"){
                		      if(interactINT==1) x_star = x_starA+x_starB-x_starC;
                		      if(interactINT== -1){
                      			   if(!useZ) x_star = x_starA-x_starB;
                      			   if(useZ & normby=="Y") x_star = x_starA + (-1*(x_starB-x_starC));
                      			   if(useZ & normby=="X") x_star = (x_starA-x_starC) -x_starB;
                		      }
                	 	}
            
            		   	if(datatype=="WITHIN"){
            		   		if(interactINT== 1) x_star = x_starA+ x_starB;
          			 	if(interactINT== -1) x_star = x_starA- x_starB;
            		   	}
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
			
		        double ss,prior_set;
			int gg=0;
		        for(s_iter = subtypeMap->begin(); s_iter!=subtypeMap->end(); s_iter++){     
			    grpClass = s_iter->first;          	    
		            if(!usePrior) {
				currGrp = s_iter->second;
				prior = (double)1/(double) currGrp.size();
				ss = TMPSCORE.at(gg) - 2*log(prior);
			    }
			    if(usePrior){
				prior_set = (sub_iter->second).getPrior(grpClass);
				ss = TMPSCORE.at(gg) - 2*log(prior_set);
			    }
                	    score.push_back(ss);
			    gg++;
             		}
     			index = whichMin(score);
			 if(index.size() !=1) {thres_error++;cerr<<"Tie in the minimum scores."<<endl;}
      			 if(index.size() ==1) if(grplabel[index.at(0)]!=grp_sub) thres_error++;
     		    	for(int r=0; r<grplabel.size(); r++){
      	         		if(grp_sub==grplabel.at(r)){
        			      (grpsize.at(r))++;
				      if(index.size() !=1) (THRES_error.at(r))++;
        			      if(index.size() ==1) {if(grplabel[index.at(0)]!=grp_sub) (THRES_error.at(r))++;}
        	       		}
      		    	 }
	            	score.clear();
			index.clear();
	            	TMPSCORE.clear();
	         }// close testsample
	         double overall_error = (double)thres_error/(double)testsample.size();
	         double THRES_grpError;

        	 for(int r=0; r<grplabel.size(); r++){
        		//cerr<<"Current thres: "<<curr_thres<<"\tGrp: "<<grplabel[r]<<"\tThres: "<<THRES_error.at(r)<<"\tGrpSize: "<<grpsize.at(r)<<endl;
        		THRES_grpError= (double) THRES_error.at(r)/ (double) grpsize.at(r);
        		if(grpsize.at(r)!=0) (GRPCVerrors.at(r)).push_back(THRES_grpError);
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
            double NumSurvCV = floor(calMean(geneSurv_vec));
	    double overallError = calMean(CVerrors);
	    errors.push_back(overallError);
	    if(overallError<minClassError){
	          minClassError = overallError;
	          selectThres = curr_thres;
	    }
	    outF<<"\t"<<overallError;
	    for(int r=0; r<grplabel.size(); r++) outF<<"\t"<<calMean(GRPCVerrors.at(r));
            
	    double NumSurv = calNumSurvived(curr_thres,maxthres, Alledges, subtypeMap, PAMmap);

	    outF<<"\t"<<NumSurvCV<<"\t"<<NumSurv<<endl;
	    CVerrors.clear();
	    GRPCVerrors.clear();
	    geneSurv_vec.clear();
    }//close threshold loop
    errors.clear();
    return selectThres;
}

double calNumSurvived(double thres,vector<double> maxthres, set<string>*Alledges, map<string, set<string> > *subtypeMap,map<string, PAM_t>*PAMmap){

	string curr_edge, datatype;
	map<string, PAM_t>::iterator piter;
        map<string, int>::iterator i_iter;
        set<string>::iterator it;
        vector<vector<int> > Grpgenesurv;
	vector<int> Binvec;

        double curr_dik, modthres, dik_new, xbar;
	vector<double> vec_dik;

        double Avgthres = calMean(maxthres);
        for(int gg=0; gg<subtypeMap->size(); gg++) Grpgenesurv.push_back(Binvec);

        for(it = Alledges->begin(); it!=Alledges->end(); it++){
                curr_edge = *it;
                piter = PAMmap->find(curr_edge);
                datatype = piter->second.dataT;
            	xbar = piter->second.xbar_i;
            	vec_dik = piter->second.d_ik;
                for(int gg=0; gg<subtypeMap->size(); gg++){
            		curr_dik = vec_dik[gg];
                	//standarding to group threshold
                	modthres = (thres/Avgthres) * maxthres.at(gg);
                        dik_new = calNewDik(piter->second.dikstar_new,modthres,gg);
                	if(dik_new==0) (Grpgenesurv.at(gg)).push_back(1);
                	if(dik_new!=0) (Grpgenesurv.at(gg)).push_back(0);	
                 }// close subtype
	 }//close edge

	 int genedied =0;
         for(int c=0; c<Alledges->size(); c++){
      	      bool died = true;
      	      for(int b=0; b<Grpgenesurv.size(); b++){
            	    vector<int> tmp_binvec = Grpgenesurv.at(b);
            	    if(tmp_binvec.at(c)==0) {died = false; continue;}
      	      }
	      if(died) genedied++;
	 }

          Grpgenesurv.clear();
	  double genesurv = (double)Alledges->size() - (double)genedied;
	  return genesurv;
}


vector<vector<string> > GetNeighbors(string str, map<string, vector<set<string> > >&Map){
    map<string, vector<set<string> > >::iterator iterA, iterB;
    set<string> set_rnaA,set_rnaB, set_protA, set_protB;
    set<string>::iterator it1,it2;
    vector<vector<string> > ret;
    vector<string> tmp;
    
    vector<string> str_vec = stringSplit2(str,"___");
    string gA = str_vec[0] ;
    string gB = str_vec[1];

    iterA = Map.find(gA);
    iterB = Map.find(gB);
    set_rnaA = (iterA->second)[1];
    set_protA = (iterA->second)[2];
    set_rnaB = (iterB->second)[1];
    set_protB = (iterB->second)[2];

    if(set_rnaA.size()>0){
       for(it1 = set_rnaA.begin(); it1!=set_rnaA.end(); it1++){
           string curr_g = *it1;
           if(curr_g!=gB)tmp.push_back(curr_g);
       }
    }
    if(set_protA.size()>0){
        for(it2 = set_protA.begin(); it2!=set_protA.end(); it2++){
            string curr_g = *it2;
            if(curr_g!=gB) tmp.push_back(curr_g);
        }
    }
    ret.push_back(tmp);
    tmp.clear();
    if(set_rnaB.size()>0){
      	for(it1 = set_rnaB.begin(); it1!=set_rnaB.end(); it1++){
      	    string curr_g = *it1;
      	    if(curr_g!=gA) tmp.push_back(curr_g);
      	}
    }
    if(set_protB.size()>0){
      	for(it2=set_protB.begin(); it2!=set_protB.end(); it2++){
          	    string curr_g = *it2;
       	    if(curr_g!=gA) tmp.push_back(curr_g);
      	}
    }
    ret.push_back(tmp);
    return ret;
    tmp.clear();
    ret.clear();

}

vector<vector<set<string> > > outputGeneSurv(string dir,bool useZ, bool usePrior,string normby,vector<double> minThres,set<string> *Alledges,vector<string> *subjects, map<string, set<string> > *subtypeMap, map<string, PAM_t> *PAMmap,map<string, subjectClass> *subjectMap, map<string, int> *InteractMap){

     vector<vector<set<string> > > edgeSurv;
     vector<set<string> >edgeSurvUp, edgeSurvDown;

     vector<string> grplabel = getKey(subtypeMap);
     map<string, PAM_t>::iterator piter;
     map<string, set<string> >::iterator s_iter;
     map<string, subjectClass>::iterator sub_iter;
     map<string, geneClass>::iterator g_iter, g_iter2;
     map<string, int>::iterator i_iter;
     set<string>::iterator it;
     vector<string> currG;
     string curr_k, curr_edge, geneA, geneB, NodeA,NodeB, gnA, gnB;
     string datatype;
     int typeA, typeB, interactINT;
         
     int ccc =0;
     bool verbose=true;
     ofstream outFile(dir+ "EdgesSelected_minThres.txt");
     outFile<<"Edge\tNodeA\tNodeB\tGeneA\tGeneB\tInteractionType";
     for(int j=0; j<grplabel.size(); j++) outFile<<"\t"+grplabel.at(j)+"_dikNew\t"+grplabel.at(j)+"_sigEdge\t"+grplabel.at(j)+"_dir\t"+grplabel.at(j)+"_absdiknew";
     outFile<<endl;
     int counter=0;
     int totalSurv=0;

     ofstream outF2(dir+"SampleClass_Probabilities.txt");
     outF2<<"Subject";
     for(int j=0; j<grplabel.size(); j++) outF2<<"\tDiscriminantScore_"+grplabel.at(j);
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
		i_iter = InteractMap->find(curr_edge);
		interactINT = i_iter->second;
  		piter = PAMmap->find(curr_edge);
  		vector<string> strVec = stringSplit2(curr_edge,"___");
  		NodeA = strVec.at(0);
  		NodeB = strVec.at(1);
                gnA = (stringSplit2(NodeA,"...")).at(0);
                gnB = (stringSplit2(NodeB,"...")).at(0);

  		allNodes.insert(NodeA);
  		allNodes.insert(NodeB);
  		datatype = piter->second.dataT;
  	  	if(datatype=="BTW") {typeA = 3;typeB=2;}
  	  	if(datatype=="WITHIN") {typeA = 3;typeB=3;} 
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
  		     outFile<<curr_edge<<"\t"<<NodeA<<"\t"<<NodeB<<"\t"<<gnA<<"\t"<<gnB<<"\t"<<datatype;

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
	cerr<<"\n\nAt the specified threshold, a total of "<<totalSurv<<" edges were selected and the number of edges predictive for each group/subtype is:"<<endl;
	for(int k=0; k<grplabel.size(); k++) cerr<<grplabel.at(k)<<":\t"<<Surv_GRP[k]<<endl;

	cerr<<"\nCalculating class probabilities on samples..."<<endl;
        
	string PredClass, currGrp;
        set<string> edgeSelected;
        double si, s0,xbar,x_star, x_starA, x_starB, x_starC, tmp_score, prior, pp,xij_new;
	vector<double> score, prior_vec;
	//vector<int> classSize;
	double thres_error=0;
        int SubjectNotfound =0;

	for(int t =0; t<subjects->size(); t++){
	    string curr_sub = subjects->at(t);
	    sub_iter = subjectMap->find(curr_sub);
	    string grp_sub = findClass(curr_sub, subtypeMap);
	    if(grp_sub.empty()) SubjectNotfound++;
	    for(int gg=0; gg<grplabel.size(); gg++){
		//if(grp_sub == grplabel.at(gg)) classSize.at(gg)++;
		tmp_score =0;
		if(usePrior){
		   currGrp = grplabel.at(gg);
		   pp = sub_iter->second.getPrior(currGrp);
		   prior_vec.push_back(pp);
		}
		for(it = EdgeSurv.begin(); it!=EdgeSurv.end(); it++){
		    curr_edge = *it;
		    i_iter = InteractMap->find(curr_edge);
		    interactINT = i_iter->second;
		    piter = PAMmap->find(curr_edge);
		    geneA = piter->second.geneA;
		    geneB = piter->second.geneB;
		    datatype = piter->second.dataT;
		    xbar = piter->second.xbar_i;
		    if(datatype=="BTW"){ typeA =3; typeB=2;}
		    if(datatype =="WITHIN"){ typeA=3; typeB=3;}
		    g_iter = sub_iter->second.geneMap.find(geneA);
		    x_starA = g_iter->second.getReadings(typeA);	
      	            g_iter2 = sub_iter->second.geneMap.find(geneB);
      	            x_starB = g_iter2->second.getReadings(typeB);
		    x_starC = 0;
      	            if(useZ & normby=="Y") x_starC = g_iter2->second.getReadings(1);
      	            if(useZ & normby=="X") x_starC = g_iter->second.getReadings(1);
		    si = piter->second.s_i;	
		    s0 = piter->second.s_o;

            	    if(datatype=="BTW"){
                	if(interactINT==1) x_star = x_starA+x_starB-x_starC;
                	if(interactINT== -1){
                      		if(!useZ) x_star = x_starA-x_starB;
                      		if(useZ & normby=="Y") x_star = x_starA + (-1*(x_starB-x_starC));
                      		if(useZ & normby=="X") x_star = (x_starA-x_starC) -x_starB;
                	}
                    }       
            	    if(datatype=="WITHIN"){
            		   if(interactINT== 1) x_star = x_starA+ x_starB;
          		   if(interactINT== -1) x_star = x_starA- x_starB;
            	    }
		    xij_new = calNewXbarik(piter->second, minThres.at(gg), gg, 2);
		    tmp_score += (pow((x_star-xij_new),2)/pow((si+s0),2));
  	
		 }// close Edges loop
		 score.push_back(tmp_score);
	   }// close subtype

	   if(!usePrior) prior = (double) 1/(double) grplabel.size();
	   double num=0, den=0, maxS=0, sss;
	   vector<double> prob;

	   for(int gg=0; gg<grplabel.size(); gg++){
		if(!usePrior) sss = score.at(gg) - 2*log(prior);
		if(usePrior) sss= score.at(gg) - 2*log(prior_vec.at(gg));
		if(sss!=sss){
		    cerr<<"Current subject: "<<curr_sub<<endl;
 		    cerr<<"grp label: "<<grplabel.at(gg)<<endl;
		    cerr<<"grp score: "<<sss<<endl;	
		}
		 score.at(gg) = sss;
		double currMax = score.at(gg);
		if(currMax>maxS) maxS = currMax;
	   }
	   //cerr<<"Max S: "<<maxS<<endl;
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
       for(int gg=0; gg<grplabel.size(); gg++) outF2<<"\t"<<score.at(gg);
	   for(int gg=0; gg<grplabel.size(); gg++) outF2<<"\t"<<prob.at(gg);
	   outF2<<"\t"<<grp_sub<<"\t"<<PredClass<<endl;
	   score.clear();
	   prob.clear();
	   prior_vec.clear();
       }

	//cerr<<"Number of subjects not found is "<<SubjectNotfound<<endl;
	double ERROR = (double) thres_error/(double)subjects->size();
        double ee = round((ERROR*1000))/10;

	cerr<<"The overall misclassification error based on the model trained in cross-validation is "<<ee<<"%."<<endl;

  	edgeSurv.push_back(edgeSurvUp);
  	edgeSurv.push_back(edgeSurvDown);

     	ofstream outFile2(dir+ "PredictiveEdges_Parameters.txt");
     	outFile2<<"GeneA\tGeneB\tInteractionType\tSign\tSi\tSo";
     	for(int j=0; j<grplabel.size(); j++) outFile2<<"\t"+grplabel.at(j)+"_Xikbar";
    	outFile2<<endl;
	for(it = EdgeSurv.begin(); it!=EdgeSurv.end(); it++){
	        curr_edge = *it;
  		i_iter = InteractMap->find(curr_edge);
		interactINT = i_iter->second;
		piter = PAMmap->find(curr_edge);
		geneA = piter->second.geneA;
		geneB = piter->second.geneB;
                gnA = (stringSplit2(geneA,"...")).at(0);
                gnB = (stringSplit2(geneB,"...")).at(0);
		datatype = piter->second.dataT;
		si = piter->second.s_i;
		s0 = piter->second.s_o;
		outFile2<<gnA<<"\t"<<gnB<<"\t"<<datatype<<"\t"<<interactINT<<"\t"<<si<<"\t"<<s0;
		for(int gg=0; gg<grplabel.size(); gg++){
			xij_new = calNewXbarik(piter->second, minThres.at(gg), gg, 2);
			outFile2<<"\t"<<xij_new;
		}
		outFile2<<endl;
	}
  
  	ofstream outFile3(dir+"AttributesTable.txt");
	outFile3<<"Node\tGene\tType";
	for(int j = 0; j<grplabel.size(); j++) outFile3<<"\t"+grplabel.at(j)+"_surv";
	outFile3<<endl;	
	
	for(it = allNodes.begin(); it!=allNodes.end(); it++){
		string currN = *it;
		vector<string> strVec = stringSplit2(currN,"...");
		string genename = strVec.at(0);
		outFile3<<currN<<"\t"<<genename<<"\t"<<strVec.at(1);
		for(int j=0; j<grplabel.size(); j++){
		   string surv="died";
		   if(contains(currN, &GrpNodes.at(j))) surv = strVec.at(1);
		   outFile3<<"\t"<<surv;
		}
		outFile3<<endl;
  	}
	outFile.close();
	outFile2.close();
	outFile3.close();
	outF2.close();
	return edgeSurv;
}



void createNNmap(set<string> *Alledges, map<string, vector<set<string> > > *FDneighborMap, map<string, vector<set<string> > > *NNmap){

    set<string>::iterator it,itA, itB;
    set<string> *Neigh_Y = NULL;
    set<string> *Neigh_X = NULL;
    set<string> *BTW_neigh = NULL;
    set<string> *WITHIN_neigh = NULL;
    string NodeA, NodeB, curr_edge;
    map<string, vector<set<string> > >::iterator iter1, iter2;
    vector<string> strVec;
    int counter=0;
    vector<set<string> > *NN = NULL;

    for(it = Alledges->begin(); it!=Alledges->end(); it++){
        Neigh_Y = new set<string>; 
	Neigh_X = new set<string>;
        BTW_neigh = new set<string>;
        WITHIN_neigh = new set<string>;
	NN = new vector<set<string> >;

      	curr_edge = *it;
      	strVec = stringSplit2(curr_edge,"___");
      	NodeA = strVec.at(0);
      	NodeB = strVec.at(1);	
      	iter1 = FDneighborMap->find(NodeA);
      	iter2 = FDneighborMap->find(NodeB);
        *Neigh_Y = iter1->second.at(1);
      	*Neigh_X = iter1->second.at(2);
      	if(Neigh_Y->size()>0){
      	   for(itA = Neigh_Y->begin(); itA!=Neigh_Y->end(); itA++){
        		  string tmpstr = NodeA+"___"+*itA;
        		  if(contains(tmpstr,Alledges)) BTW_neigh->insert(tmpstr);
        	 }
        }
      	if(Neigh_X->size()>0){
      	   for(itB = Neigh_X->begin(); itB!=Neigh_X->end(); itB++){
            		string tmpstr = NodeA+"___"+*itB;
            		string newstr = getEdge(tmpstr,Alledges);
      	        if(!newstr.empty()) WITHIN_neigh->insert(newstr);
      	   }
        }
      	*Neigh_Y = iter2->second.at(1);
      	*Neigh_X = iter2->second.at(2);
      	if(Neigh_Y->size()>0){
      	   for(itA = Neigh_Y->begin(); itA!=Neigh_Y->end(); itA++){
      		string tmpstr = NodeB+"___"+*itA;
      	        if(contains(tmpstr,Alledges)) BTW_neigh->insert(tmpstr);
      	   }
      	}
      	if(Neigh_X->size()>0){
      	   for(itB = Neigh_X->begin(); itB!=Neigh_X->end(); itB++){
      		string tmpstr = NodeB+"___"+*itB;
      		string newstr = getEdge(tmpstr,Alledges);
      		if(!newstr.empty()) WITHIN_neigh->insert(newstr);
      	   }
      	}
      	NN->push_back(*BTW_neigh);
      	NN->push_back(*WITHIN_neigh);
      	(*NNmap)[curr_edge] = *NN;
 
        delete BTW_neigh;
      	delete WITHIN_neigh;
      	delete Neigh_Y;
      	delete Neigh_X;
      	delete NN;
      	BTW_neigh = NULL;
      	WITHIN_neigh = NULL;
      	Neigh_Y = NULL;
      	Neigh_X = NULL;
      	counter++;
    }

}

string getEdge( string str, set<string> * Edges){
     string ret;
     if(contains(str, Edges)) ret = str;
     else{
      	vector<string> vecstr = stringSplit2(str,"___");
      	string newstr = vecstr.at(1)+"___"+ vecstr.at(0);
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

string getTag(string str, map<string, string > *Map){
   string ret="NA";
   map<string, string>::iterator iter;
   for(iter = Map->begin(); iter!=Map->end(); iter++){
      if((iter->second)== str) ret = (iter->first);
      continue;
   }
   return ret;
}

void CreateOutput(string directory, string lab,vector<string> PATHlist,int dir, set<string> *geneLIST, set<string> *bgLIST, map<string, set<string> > *PATHEdgemap, map<string, set<string> > *PATHmap, map<string, string> *PATHwayAnnot,int mink){

  string DIR, Annot;
    string geneSET="";
  if(dir==0) DIR = "up";
  if(dir==1) DIR ="down";
  vector<string> collect;
  set<string> memEdges, memGenes, AllG,All_genesignature,common_set, tmp_setN;
  set<string>::iterator it;
  map<string, string>::iterator a_iter;
  vector<string> tmpVec,common_setV;
  ofstream outF(directory+lab+"_Enrichment_"+DIR+".txt");
  outF<<"Pathway\tPathAnnotation\tHypergeoPval\tNumGene_inPathway\tNumGenes_inbg\tPropPathway\tNumEdgesformed_inbg\tEnriched_Edgesize\tGeneset_Enriched"<<endl;
  map<string, set<string> >::iterator M_iter, m_iter;
  
  for(it = bgLIST->begin(); it!= bgLIST->end(); it++){
    tmpVec = stringSplit2(*it,"___");
    AllG.insert((stringSplit2(tmpVec.at(0),"...")).at(0));
    AllG.insert((stringSplit2(tmpVec.at(1),"...")).at(0));
  }
    for(it = geneLIST->begin(); it!= geneLIST->end(); it++){
      tmpVec = stringSplit2(*it,"___");
      All_genesignature.insert((stringSplit2(tmpVec.at(0),"...")).at(0));
      All_genesignature.insert((stringSplit2(tmpVec.at(1),"...")).at(0));
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
      common_set = IntersectSet(All_genesignature, tmp_setN);
      if(common_set.size()>0){
          common_setV = SetToVec(common_set);
          geneSET = concatenate(common_setV,';');
      }
      double pp = (double) tmp_setN.size()/ (double) memGenes.size();
      vector<double> hyperpval = hypergeo(&memEdges, geneLIST, bgLIST);
      if(hyperpval[2]>=mink){
            outF<<path<<"\t"<<Annot<<"\t"<<hyperpval.at(0)<<"\t"<<memGenes.size()<<"\t"<<tmp_setN.size() <<"\t"<<pp<<"\t"<<memEdges.size()<< "\t"<<hyperpval.at(2)<<"\t"<< geneSET<<endl;
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
    	tmpVec = stringSplit2(*it,"___");
    	AllG.insert((stringSplit2(tmpVec.at(0),"...")).at(0));
    	AllG.insert((stringSplit2(tmpVec.at(1),"...")).at(0));
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
    vector<string> vecString, tmpA, tmpB;

    for(it = ALLedges->begin(); it!=ALLedges->end(); it++){
      	edge = *it;
      	vecString = stringSplit2(edge,"___");
      	geneA =(stringSplit2(vecString.at(0),"...")).at(0);
      	geneB =(stringSplit2(vecString.at(1),"...")).at(0);
      	if(contains(geneA,SET) & contains(geneB,SET)) ret.insert(edge);
    }
    return ret;
}


void printBGlist(string dir, set<string> edges){

    set<string>::iterator it;
    ofstream outF(dir+"BGlist.txt");
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


void StandardizeFeatures(string dir, vector<string> sub, vector<string> genes, int datatype, map<string, subjectClass> &subjectMap){

     map<string, subjectClass>::iterator s_iter;
     map<string, geneClass>::iterator g_iter;
     vector<double> tmp_vec,z_vec;
     string curr_g,curr_sub;    
     double vec_mean, vec_sd;
     string lab;
     if(datatype==1) lab="dataZ";
     if(datatype==2) lab="dataY";
     if(datatype==3) lab="dataX";
     ofstream outF(dir+"Ztransform_"+lab+".txt");
     outF<<"Feature";
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
            if(datatype==1) g_iter->second.addZ(z_vec.at(j));
            if(datatype==2) g_iter->second.addY(z_vec.at(j));
            if(datatype==3) g_iter->second.addX(z_vec.at(j));
      	}
      	outF<<endl;
        tmp_vec.clear();
        z_vec.clear();
     }
     outF.close();
}


