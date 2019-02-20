/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: Hiromi WL Koh
 *
 * Last modified on 13 June, 2018, 4:00 PM
 */

#include "globals.hpp"
#include "geneClass.hpp"
#include "subjectClass.hpp"

#include <algorithm>
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <set>
#include <ctime>
#include <cmath>

using namespace std;
/*
 * 
 */
int main(int argc, char** argv) {
    
    if(argc <2){
        cerr <<"\nUSAGE: ./iOmicsPASS <<inputFile>\n";
        return 0;
    }

    clock_t begin = clock();
    clock_t time;
    double elapse, elapse2;

    string inputF = argv[1];	
    Input_t *UserInput = new Input_t();    

    map<string, subjectClass> subjectMap;
    map<string, PAM_t> PAMmap;
    map<string, set<string> > PPIMap;
    map<string, set<string> > TFMap, PATHmap;
    map<string, set<string> > subtypeMap;
    map<string, vector<set<string> > > FirstDNeighborMap;
    map<string, int> InteractMapTF;
    map<string, string> PATHwayAnnot;

    vector<string> dna_g,rna_g, prot_g;
    set<string> TF_edge, PPI_edge, pathways;
    vector<string> GENE, gene_interest;
    
    cerr << "\n\nInitiating iOmicsPASS...\n" <<endl; 
    readUserInput(inputF, *UserInput);  
    
    cerr<<"\nApplying the following filter (specified by user) for integrity of genes across subtypes:"<<endl;
    cerr<<"Minimum number of non-missing observations for each gene within each subtype: "<< UserInput->minobs<<endl;
    cerr<<"Minimum proportion of non-missing observations for each gene within each subtype: "<<UserInput->minprop<<endl; 

    cerr << "\nReading in the subtype information file....." <<endl;
    set<string> subjects = insertSubtypeMap(UserInput->subtypefile, subtypeMap);    
    vector<string> grplab = getKey(&subtypeMap);  

    if(UserInput->analyseDNA){
       cerr << "\nReading in the DNA-level file......"<<endl;
       dna_g = readFile(UserInput->dnafile, subtypeMap,subjectMap,*UserInput,subjects, 1);
       cerr<<"\nThere are "<<dna_g.size()<<" genes in the DNA-level file."<<endl;
    }

    cerr << "\nReading in the RNA-level file......"<<endl;
    rna_g = readFile(UserInput->rnafile,subtypeMap, subjectMap,*UserInput,subjects, 2);
    cerr<<"\nThere are "<<rna_g.size()<<" genes in the RNA-level file."<<endl;

    cerr<<"\nReading in the Protein-level file......"<<endl;
    prot_g = readFile(UserInput->proteinfile, subtypeMap,subjectMap,*UserInput,subjects, 3);
    cerr<<"\nThere are "<<prot_g.size()<<" genes in the Protein-level file."<<endl;

    cerr << "\nReading in the TF network file......"<<endl;
    //if DNA file is provided, we require the genes to be present in both RNA and DNA.
    vector<string> target_g;
    if(!UserInput->analyseDNA) target_g = rna_g;
    if(UserInput->analyseDNA) target_g = IntersectVec(rna_g, dna_g);
    TF_edge = readTFNetwork(UserInput->Interact, UserInput->tfnet, TFMap,InteractMapTF, target_g, prot_g);

    cerr << "\nReading in the PPI network file....." <<endl;
    PPI_edge = readPPINetwork(UserInput->ppinet, PPIMap, prot_g); 

    set<string> Alledges = UnionSet(PPI_edge,TF_edge);
    cerr<<"\nThere are a total of "<<Alledges.size()<<" edges constructed from both TF and PPI network file."<<endl;
    printBGlist(Alledges);

    cerr<<"\nReading in the Pathway module file....."<<endl;
    pathways = readPATHWAY(UserInput->pathwayfile, PATHmap, PATHwayAnnot);

    set<string> allS = retrieveKeys(subjectMap);

    cerr <<"\nAll data files have been successfully read into iOmicsPASS!\n\nProceeding to next step...\n";
    vector<string> sub_dna, sub_rna, sub_prot;
    if(UserInput->analyseDNA) sub_dna = getSubjects(1, dna_g, subjectMap);
    sub_rna = getSubjects(2, rna_g, subjectMap);
    sub_prot = getSubjects(3, prot_g, subjectMap);

    if(UserInput->analyseDNA) cerr<<"\nThere are "<< sub_dna.size()<< " subjects with subtype information and "<< dna_g.size()<<" genes after filtering in the DNA-level file"<<endl;
    cerr<<"There are "<<sub_rna.size()<<" subjects with subtype information and "<<rna_g.size()<< " genes after filtering in the RNA-level file"<<endl;
    cerr<<"There are "<<sub_prot.size()<<" subjects with subtype information and "<< prot_g.size()<<" genes after filtering in the Protein-level file"<<endl;

    //time = clock();
    //elapse = double(time - begin)/(double)CLOCKS_PER_SEC;
    //cerr<<"\nTime : "<<elapse<<" seconds.\n"<<endl; 

    bool impute=UserInput->knnimpute;
    bool checkD, checkR, checkP;
    if(impute){
       if(!UserInput->analyseDNA) checkD = false;
       if(UserInput->analyseDNA) checkD=checkNAs(subjectMap,1);
       checkR=checkNAs(subjectMap,2);
       checkP=checkNAs(subjectMap,3);
       string new_lab;
       vector<string> lab;
       if(checkD) lab.push_back("DNA");
       if(checkR) lab.push_back("RNA");
       if(checkP) lab.push_back("Protein");
       if(lab.size()==0) cerr<<"\nNo missing values detected, skipping KNN-imputation step.\n"<<endl;
       else{
         for(int j=0; j<lab.size(); j++){
  	    if(j==0) new_lab = lab.at(j);
 	    else new_lab+= " and " + lab.at(j);
         }
         cerr<<"\nKNN imputation will be carried out on "<<new_lab<<" datasets.\n"<<endl;
         int knnK = UserInput->knnK;
	 int blocksize = UserInput->blocksize;
	 map<set<string>, vector<vector<double> > > CandidateSET_d, CandidateSET_r, CandidateSET_p;

         if(checkD){
		cerr<<"Imputing DNA-level data...\n";
		vector<vector<double> > MAT_D = creatematrix(dna_g, sub_dna, subjectMap,1);
		determineCandidate(blocksize,dna_g, MAT_D, CandidateSET_d);
		time = clock();
    		elapse2 = double(time - begin)/(double)CLOCKS_PER_SEC;
    		//cerr<<"\nTook : "<<(elapse2-elapse)<<" seconds.\n"<<endl; 
		KNNimpute(UserInput->cid_dna,sub_dna, dna_g, MAT_D,CandidateSET_d, subjectMap,knnK,1 );
		outputData(sub_dna, dna_g, subjectMap,1);
		time = clock();
    		elapse = double(time - begin)/(double)CLOCKS_PER_SEC;
    		cerr<<"\nTook : "<<(elapse-elapse2)<<" seconds.\n"<<endl; 

	 }
         if(checkR){
		cerr<<"Imputing RNA-level data...\n";
	  	vector<vector<double> > MAT_R = creatematrix( rna_g,sub_rna, subjectMap, 2);
		determineCandidate(blocksize,rna_g, MAT_R, CandidateSET_r);
		time = clock();
    		elapse2 = double(time - begin)/(double)CLOCKS_PER_SEC;
    		//cerr<<"\nTook : "<<(elapse2-elapse)<<" seconds.\n"<<endl; 
		KNNimpute(UserInput->cid_rna,sub_rna, rna_g, MAT_R, CandidateSET_r,subjectMap, knnK,2);
		outputData(sub_rna, rna_g, subjectMap,2);
		time = clock();
     		elapse = double(time - begin)/(double)CLOCKS_PER_SEC;
    		cerr<<"\nTook : "<<(elapse-elapse2)<<" seconds.\n"<<endl; 

	 }
         if(checkP){
		cerr<<"Imputing Protein-level data...\n";
		vector<vector<double> > MAT_P = creatematrix( prot_g,sub_prot,subjectMap,3);
		determineCandidate(blocksize, prot_g, MAT_P, CandidateSET_p);
		time = clock();
    		elapse2 = double(time - begin)/(double)CLOCKS_PER_SEC;
    		//cerr<<"\nTook : "<<(elapse2-elapse)<<" seconds.\n"<<endl; 
		KNNimpute(UserInput->cid_prot,sub_prot, prot_g, MAT_P, CandidateSET_p, subjectMap, knnK,3);
		outputData(sub_prot, prot_g, subjectMap,3);	
		time = clock();
    		elapse = double(time - begin)/(double)CLOCKS_PER_SEC;
    		cerr<<"\nTook : "<<(elapse-elapse2)<<" seconds.\n"<<endl; 

         }
       }
    }   

    //run PAMR to identify differential genes in each subtype 
    //combining RNA and PROT level together  
    vector<string> commonSub = IntersectVec(sub_rna, sub_prot);
    if(UserInput->analyseDNA) commonSub = IntersectVec(commonSub,sub_dna);

    // standardize measurements in the map
    if(UserInput->analyseDNA & UserInput->ztrans_dna) StandardizeFeatures(commonSub,dna_g,1,subjectMap);
    if(UserInput->ztrans_prot) StandardizeFeatures(commonSub, prot_g,3, subjectMap);
    if(UserInput->ztrans_rna) StandardizeFeatures(commonSub, rna_g,2, subjectMap);
    if(!UserInput->analyseDNA)cerr<<"\n\nNumber of common subjects with RNA and PROT data is: " <<commonSub.size()<<endl;
    if(UserInput->analyseDNA) cerr<<"\n\nNumber of common subjects with all 3 types of data is: "<<commonSub.size()<<endl;
 
    map<string, set<string> > Targets2TFmap;
    ReverseMap(TFMap, Targets2TFmap);
    set<string> features = createNeighborMap(&prot_g, &target_g, &TFMap,&Targets2TFmap, &PPIMap, FirstDNeighborMap);

    cerr<<"Carrying out shrunken centroid module on the interaction edges..."<<endl;

    bool useDNA = UserInput->analyseDNA;
    vector<double> thresMax = fillPAMmap2(useDNA,UserInput->Interact, &Alledges,&commonSub, &InteractMapTF, &subtypeMap,&subjectMap, PAMmap, &FirstDNeighborMap);
    cerr<<"\nMaximum Dij (average) in mRNA and PROT data is : "<< calMean(thresMax) <<"\n"<<endl;
    cerr<<"Class-specific thresholds are:";
    for(int i=0; i<thresMax.size(); i++) cerr<<"\t"<<thresMax.at(i);
    cerr<<endl;

    double minTHRES, minThres;
    minTHRES = UserInput->minthres;
    if(UserInput->crossV){
	cerr<<"\nCarrying out cross-validation...(This may take a while)\n"<<endl;
    	minThres = CVKfold(useDNA,UserInput->Interact,UserInput->numFold, &Alledges,&commonSub,  30, &InteractMapTF, &FirstDNeighborMap, &subtypeMap, &subjectMap, PAMmap);
    }
    if(!UserInput->crossV) cerr<<"\nSkipping cross-validation and using user-specified threshold: "<<minTHRES<<"\n"<<endl;
    if(minTHRES==NA_VAL) minTHRES = minThres;
    vector<double> MODthres = thresMax;
    for(int i=0; i<thresMax.size(); i++) MODthres.at(i) = (minTHRES/calMean(thresMax))*thresMax.at(i);
    vector<vector<set<string> > > geneSurv = outputGeneSurv(useDNA, UserInput->Interact,MODthres,&Alledges, &commonSub, &subtypeMap,  &PAMmap, &subjectMap, &InteractMapTF);

    //Network-centric Enrichment
    
    cerr<<"\nStarting Network enrichment on surviving edges for each subtype..."<<endl;
    map<string, set<string> > PATHmapEdge;
    vector<string> PATHs = createEdgeOrientedPathwayMap(&features,&Alledges, &PATHmap, &PATHmapEdge,UserInput->bgProp, UserInput->minbgSize);
    set<string> EdgeSurv_grp;
    
    for(int j=0; j<2;j++){
      for(int i=0; i<grplab.size(); i++){
         EdgeSurv_grp = (geneSurv.at(j)).at(i);
         CreateOutput(grplab.at(i), PATHs, j, &EdgeSurv_grp, &Alledges, &PATHmapEdge, &PATHmap,&PATHwayAnnot, UserInput->minsigSize);
      }
    }

    cerr <<"\n\nThe program has finished running!\n" ;
    clock_t end = clock();
    double elapse_sec = double(end - begin)/(double)CLOCKS_PER_SEC;
    string timetype;
    double new_t, tmpt;
    tmpt = elapse_sec/60.0;
    double tt = round( tmpt*10)/10;
    if(tt>240){
	timetype = " hours";
	new_t = round((tmpt/60.0)*10)/10;
    }
    if(tt<=240){
	timetype =" minutes";
	new_t = tt;
    }
    cerr<<"\nTime to completion of iOmicsPASS is : "<<new_t<<timetype<<endl;
    cerr<<"----------------------------------------------------------------"<<endl;
    
    return 0;
}

