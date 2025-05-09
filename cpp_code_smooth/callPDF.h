/* 
 * File:   PDFMainMatlab.h
 * Author: jenny
 *
 * Created on December 3, 2018, 8:32 AM
 */


#ifndef CALLPDF_H
#define	CALLPDF_H


//#include "mexAdapter.hpp"
#include "InputParameters.h"
#include "InputData.h"
#include "Score.h"
#include "ScoreQZ.h"
#include "MinimizeScore.h"
#include "WriteResults.h"
#include "OutputControl.h"


#include <vector>
#include <string>


using namespace std;
//using namespace matlab::data;

class callPDF {
    
public:
    callPDF(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr);
    callPDF(const callPDF& orig);
    virtual ~callPDF();
    void makeCall(int sampleLength, double * sampleData);    
    
    void setHigh(double high);
    void setLow(double low);
    void setSURDtarget(double surd);
    void setSURDmin(double surd);
    void setSURDmax(double surd);
    void setScoreType(string type);
    void setPoints(int points);
    void setLagrangeMin(int lagrange);
    void setLagrangeMax(int lagrange);
    void setPartitionSize(int partition);
    void setFuzz(bool fuzz);
    void setVariance(bool minVariance) {this->minVariance = minVariance;}
    void setDebug(bool debug);
    void setOutlierCutoff(double cutoff);
    void setAdaptiveDx(bool adaptive);
    void setAdaptivePartition(int adaptivePartition);
    void setEstimationPoints(vector <double>  estimationPoints);
    void setSmooth(int smoothConstant);
    
    bool solutionFailed;        
    double thresholdScore;  
    double confidenceScore;
    
    int initPartitionSize;
    int partitionSize;
    int targetPartitionSize;
    
    vector <double> Vr;    
    vector <double> Vcdf;
    vector <double> Vpdf;
    vector <double> Vx;
    vector <double> Vlagrange; 
    vector <double> Vsqr;
    vector <string> VError;
private:
    InputParameters input;  
    OutputControl out;
    int adaptivePartition;
    int sampleLength;
    bool debug = false;
    bool minVariance = false;
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
};

#endif