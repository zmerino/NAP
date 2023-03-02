/* 
 * File:   callPDFMatlab.cpp
 * Author: jenny
 * 
 * Created on December 3, 2018, 8:36 AM
 */


#include "callPDF.h"


callPDF::callPDF(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr) {
    this->matlabPtr = matlabPtr;
    out.setPtr(matlabPtr);
    out.debug = debug;    
}

callPDF::callPDF(const callPDF& orig) {
}

callPDF::~callPDF() {    
}

void callPDF::setLow(double low) {
    input.lowerBoundSpecified = true;
    input.lowerBound = low;
}

void callPDF::setHigh(double high) {
    input.upperBoundSpecified = true;
    input.upperBound = high;
}

void callPDF::setSURDtarget(double surd) {
    input.SURDTarget = surd;
}

void callPDF::setSURDmin(double surd) {
    input.SURDMinimum = surd;
}

void callPDF::setSURDmax(double surd) {
    input.SURDMaximum = surd;
}

void callPDF::setScoreType(string type) {
    input.scoreType = type;
}

void callPDF::setPoints(int points) {
    input.integrationPoints = points - 1;
}

void callPDF::setDebug(bool debug) {
    this->debug = debug;
}

void callPDF::setPartitionSize(int partition) {
    input.initPartitionSize = partition;
}

void callPDF::setLagrangeMin(int lagrange) {
    input.minLagrange = lagrange;
}

void callPDF::setLagrangeMax(int lagrange) {
    input.maxLagrange = lagrange;
}

void callPDF::setOutlierCutoff(double cutoff) {
    input.outlierCutoff = cutoff;
}

void callPDF::setAdaptiveDx(bool adaptive) {
    input.adaptive = adaptive;
}

void callPDF::setEstimationPoints(vector <double> estimationPoints) {       
    input.estimatePoints = true;       
    input.setEstimationPoints(estimationPoints);
}

void callPDF::makeCall(int sampleLength, double* sampleData) {  
      
    out.debug = true;
           
    MinimizeScore minimumPDF = MinimizeScore();
    InputData data = InputData(input);
    
    data.out.setPtr(matlabPtr);
    data.out.debug = debug;
    minimumPDF.out.setPtr(matlabPtr);
    minimumPDF.out.debug = debug;
    input.out.setPtr(matlabPtr);
    input.out.debug = debug;
    vector <double> inputData;
    for (int i = 0; i < sampleLength; i++) {
        inputData.push_back(sampleData[i]);
    }
    input.writeHeader = true;
    input.writeFile = true;
    if (input.initPartitionSize == 0) {
        input.initPartitionSize = sampleLength;
    }
    data.setData(inputData);     
    if (data.processData()) {
        
        solutionFailed = minimumPDF.minimize(input, data);
        WriteResults write;        
        write.out.debug = debug;
        write.out.setPtr(matlabPtr);
        write.createSolution(input, data, minimumPDF);  
                      
        write.createQQ(minimumPDF.bestRandom, minimumPDF.N);  
        thresholdScore = minimumPDF.bestThreshold;
//        confidenceScore = score->getConfidence(minimumPDF.bestThreshold);
        initPartitionSize = minimumPDF.initPartitionSize;
        partitionSize = minimumPDF.partitionSize;
        targetPartitionSize = minimumPDF.targetPartition;            
        double * r = minimumPDF.bestRandom;
        for (int i = 0; i < minimumPDF.N; i++) {
            Vr.push_back(r[i]); 
        }
        Vlagrange = write.L;
        Vsqr = write.SQR;
        Vcdf = write.CDF;
        if (input.estimatePoints) {
            Vpdf = write.PDFPoints;
            Vx   = write.xPoints;
        } else {                
            Vpdf = write.PDF;
            Vx   = write.x;
        }
    } else {
        solutionFailed = true;
    }     
}

