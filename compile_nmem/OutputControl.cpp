/* 
 * PDF Estimator:  A non-parametric probability density estimation tool based on maximum entropy
 * File:   OutputControl.cpp
 * Copyright (C) 2018
 * Jenny Farmer jfarmer6@uncc.edu
 * Donald Jacobs djacobs1@uncc.edu
 * 
 * This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in 
 * the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with 
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "OutputControl.h"



OutputControl::OutputControl() {
#ifdef outputR
    debug = false;
#endif
}

OutputControl::OutputControl(const OutputControl& orig) {
}

OutputControl::~OutputControl() {
}

#ifdef outputMatlab
void OutputControl::displayError(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr, std::string errorMessage) {    
    matlab::data::ArrayFactory factory;
    matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"), 0, std::vector<matlab::data::Array>({
    factory.createScalar(errorMessage) }));
}

void OutputControl::displayOnMATLAB(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr, std::string message) {
    if (debug) {
        matlab::data::ArrayFactory factory;
        matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("fprintf"),0, std::vector<matlab::data::Array>
            ({ factory.createScalar(message)}));
        matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("fprintf"),0, std::vector<matlab::data::Array>
            ({ factory.createScalar("\n")}));
    }
}

void OutputControl::setLagrange(int first, vector <double> lagrange, int N, vector<int> indices, vector<double> x, vector<double> cdf, vector<double> data, vector<double> r, double max, double min) {
    matlab::data::ArrayFactory factory;
    matlab::data::Array firstVal = factory.createScalar(first);
    matlab::data::Array lagrangeVal = factory.createArray({lagrange.size(), 1}, lagrange.begin(), lagrange.end());  
    matlab::data::Array maxVal = factory.createScalar(max);   
    matlab::data::Array minVal = factory.createScalar(min);  
    std::vector<matlab::data::Array> args = {firstVal, lagrangeVal, maxVal, minVal};
//    matlabPtr->feval("plotLagrange", args);
    
    matlab::data::Array xVal = factory.createArray({x.size(), 1}, x.begin(), x.end());
    matlab::data::Array cdfVal = factory.createArray({cdf.size(), 1}, cdf.begin(), cdf.end());
    matlab::data::Array dataVal = factory.createArray({data.size(), 1}, data.begin(), data.end());
    matlab::data::Array rVal = factory.createArray({r.size(), 1}, r.begin(), r.end());
    matlab::data::Array nVal = factory.createScalar(N);
    matlab::data::Array indVal = factory.createArray({indices.size(), 1}, indices.begin(), indices.end());
    args = {firstVal, nVal, indVal, xVal, cdfVal, dataVal, rVal};
//    matlabPtr->feval("plotSQR", args);
    
    args = {firstVal, rVal, nVal, indVal, lagrangeVal};
    matlabPtr->feval("analyzeScore", args);
    
}

void OutputControl::setPtr(std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr) {
    this->matlabPtr = matlabPtr; 
}

 void OutputControl::print(string output) { 
     displayOnMATLAB(matlabPtr, output); 
 }
 
 void OutputControl::error(string output) { 
     displayError(matlabPtr, output); 
 }

void OutputControl::print(string output, int value) {
    ostringstream strOut;
    strOut << output << ": " << value;                                      
    print(strOut.str());
} 

void OutputControl::print(string output, double value) {
    ostringstream strOut;
    strOut << output << ": " << value;                                      
    print(strOut.str());
}

void OutputControl::error(string output, int value)  {
    ostringstream strOut;
    strOut << output << ": " << value;                                      
    error(strOut.str());
}

void OutputControl::error(string output, double value)  {
    ostringstream strOut;
    strOut << output << ": " << value;                                      
    error(strOut.str());
}
#endif


#ifdef outputR


void OutputControl::error(string output) {
    REprintf("%s\n", output.c_str());
}

void OutputControl::error(string output, int value) {
    REprintf("%s: %d\n", output.c_str(), value);
}

void OutputControl::error(string output, double value) {
    REprintf("%s: %f\n", output.c_str(), value);
}


void OutputControl::print(string output) {
    if (debug) {
        Rprintf("%s\n", output.c_str());
    }
}

void OutputControl::print(string output, int value) {
    if (debug) {
        Rprintf("%s: %d\n", output.c_str(), value);
    }
}

void OutputControl::print(string output, double value) {
    if (debug) {
        Rprintf("%s: %f\n", output.c_str(), value);
    }
}

#endif

#ifdef outputCommandLine

void OutputControl::print(string output) {
    if (debug) {
        cout << output << "\n";
    }
}

void OutputControl::print(string output, int value) {
    if (debug) {
        cout << output << ": " << value << "\n";
    }
}

void OutputControl::print(string output, double value) {
    if (debug) {
        cout << output << ": " << value << "\n";
    }
}

void OutputControl::error(string output) {
    cout << output << "\n";
}

void OutputControl::error(string output, int value) {
   cout << output << ": " << value << "\n";
}

void OutputControl::error(string output, double value) {
    cout << output << ": " << value << "\n";
}

#endif
    
   