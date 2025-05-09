/* 

 * File:   PDFMainMatlab.cpp
 * Author: jenny
 * 
 * Created on December 3, 2018, 8:32 AM
 */


#include "EstimatePDFmv.h"

    

void MexFunction::operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {   

    int nParameters = inputs.size();    
    matlabPtr = getEngine();
    out.setPtr(matlabPtr);
    out.debug = true;    
    
    InputParameters input;
    input.out.setPtr(matlabPtr);
            
    matlab::data::ArrayFactory factory;
    if (nParameters < 1) {        
        out.displayError(matlabPtr, "Data sample required");
    }
    if (nParameters > 5) {        
        out.displayError(matlabPtr, "Too many input parameters");
    }            

    Array temp = inputs[1];    
    int nSamples = temp[0];
    temp = inputs[2];
    int nVariables = temp[0];    
    temp = inputs[3];
    int resolution = temp[0];
    
    int writeMarginal = 0;
    if (nParameters == 5) {
        temp = inputs[4];
        writeMarginal = temp[0];
    }       
       
    vector <Variable> variables;
    variables.reserve(nVariables);
            
    TypedArray<double> doubleArray = std::move(inputs[0]);
    
    int iSample = 0;
    int iVariable = 0;
    vector <double> samples;
    samples.reserve(nSamples);
    for (auto& elem : doubleArray) {
        samples.push_back(elem);
        if (++iSample == nSamples) {          
            ostringstream vString; 
            vString << iVariable++ << "_" << writeMarginal; 
            Variable variable = Variable(input, samples, vString.str(), writeMarginal);
            variables.push_back(variable);
            samples.clear();
            iSample = 0;
        }
    }           
   
    JointProbability jp = JointProbability(variables, nSamples, resolution);
    jp.out.setPtr(matlabPtr);
    jp.out.debug = false;    
    jp.calculate();
    
    vector <double> pdf = jp.getJP();
    ArrayDimensions sz = {pdf.size(), 1, 1};
    TypedArray<double> Mpdf = factory.createArray(sz, pdf.begin(), pdf.end());
    outputs[0] = Mpdf;
    
    vector <double> x = jp.getRange();
    sz = {x.size(), 1, 1};
    TypedArray <double> Mx = factory.createArray(sz, x.begin(), x.end());
    outputs[1] = Mx;
      
}