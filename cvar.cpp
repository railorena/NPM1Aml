#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <boost/algorithm/string.hpp>
#include <string.h>
#include <stdio.h>
#include <cmath>
#include <sstream>

using namespace std;
using namespace boost::algorithm;

struct selection{
    std::vector<int> vecSelection;
    float cvCond = 0.0;
};

void printVector(vector<string> vectorToPrint){
    for (std::vector<string>::const_iterator i = vectorToPrint.begin(); i != vectorToPrint.end(); ++i)
        std::cout << *i << ' ';
}

void printVectorInt(vector<int> vectorToPrint){
    for (std::vector<int>::const_iterator i = vectorToPrint.begin(); i != vectorToPrint.end(); ++i)
        std::cout << *i << ' ';
}

vector<string> getCondition(vector<string> vectorHeader, vector<string> cond1, vector<string> cond2){
    //printVector(vectorHeader);
    vector<string> conditionsIndex;

    for (std::vector<string>::const_iterator i = vectorHeader.begin(); i != vectorHeader.end(); ++i)
        if ( std::find(cond1.begin(), cond1.end(), *i) != cond1.end() )
            conditionsIndex.push_back("0");
        else if ( std::find(cond2.begin(), cond2.end(), *i) != cond2.end() )
            conditionsIndex.push_back("1");
        else
            conditionsIndex.push_back("2");
    //printVector(conditionsIndex);
    return (conditionsIndex);
}

int checkZeros(std::vector<int> vecCond1Check, std::vector<int> vecCond2Check){
    int diffentZcond1 = 0,  diffentZcond2 = 0;
    int diffZero = 0;

    for (std::vector<int>::const_iterator i = vecCond1Check.begin(); i != vecCond1Check.end(); ++i){
        if(*i != 0){
            diffentZcond1++;
        }
    }

    for (std::vector<int>::const_iterator i = vecCond2Check.begin(); i != vecCond2Check.end(); ++i){
        if(*i != 0){
            diffentZcond2++;
        }
    }

    //if((diffentZcond1 >= (vecCond1Check.size()/2)) || (diffentZcond2 >= (vecCond2Check.size()/2)))
    if((diffentZcond1 > (vecCond1Check.size()/2)) || (diffentZcond2 > (vecCond2Check.size()/2)))
        diffZero = 1;

    return (diffZero);
}

float meanCond(vector<int> vecCond){
    float average = accumulate( vecCond.begin(), vecCond.end(), 0.0/ vecCond.size());
    return (average);
}

float sdCond(float mean1, float mean2, float mean){
    vector<float> vec = {mean1, mean2};
    float summation = 0.0, st = 0.0;

    for(int i = 0; i <2; ++i) {
        //The pow() function returns the result of the first argument raised to the power of the second argument
        //pow() is like (data-mean)²
        summation += pow(vec[i] - mean, 2);
    }
    st = sqrt(summation / 2);
    vec.clear();

    return (st);
}

float calculateCV(float mean1, float mean2){
    float stdev = 0.0, mean = 0.0, cv = 0.0;
    mean = (mean1 + mean2) /2;
    stdev = sdCond(mean1, mean2, mean);
    cv = stdev/mean;
    return (cv);
}

void writeRow(ofstream outfile, vector<string> splittedRow){
    for (std::vector<string>::const_iterator i = splittedRow.begin(); i != splittedRow.end(); ++i)
        outfile << *i << ' ';
}

selection processingRow(vector<string> vectorHeader, vector<string> splittedRow){
    std::vector<int> vecCond1;
    std::vector<int> vecCond2;
    float mean1 = 0.0, mean2 = 0.0, cvConditions = 0.0;
    vecCond1.clear();
    vecCond2.clear();
    int diff = 0;

    selection sample;

    for (long long i=1; i<splittedRow.size(); ++i){
        if(vectorHeader[i] == "0"){
            //cout << vectorHeader[i]<< " ";
            vecCond1.push_back(std::stoll(splittedRow[i]));
        }else if(vectorHeader[i] == "1"){
            //cout << vectorHeader[i]<< " ";
            vecCond2.push_back(std::stoll(splittedRow[i]));
        }

    }

    diff = checkZeros(vecCond1, vecCond2);

    if(diff == 1){
        mean1 = meanCond(vecCond1);
        mean2 = meanCond(vecCond2);


        cvConditions = calculateCV(mean1, mean2);
    }else{
        cvConditions = 0;
    }

    vecCond1.insert( vecCond1.end(), vecCond2.begin(), vecCond2.end() );

    sample.vecSelection = vecCond1;
    sample.cvCond = cvConditions;

    return (sample);
}


int main (int argc, char **argv){

    time_t t;
    struct tm * tt;
    time (&t);
    tt = localtime(&t);
    cout << "Started at "<< asctime(tt) << endl;


    //cout << argv[1] << endl;
    char * filename = argv[1];
    //char * filename = "input.tsv";

    //======================================//
    //---reading the config file---//
    string lineconfig;
    std::vector<std::string> cond1;
    std::vector<std::string> cond2;
    std::ifstream infileconf("samples_cond");
    int nrowconfig = 0;

    cond1.clear();
    while (getline (infileconf, lineconfig)){
            std::istringstream bylineconfig(lineconfig);
            std::string itemconfig;


            if(nrowconfig == 0){
                while(getline(bylineconfig, itemconfig, ','))
                    cond1.push_back(itemconfig);
            }else{
                while(getline(bylineconfig, itemconfig, ','))
                    cond2.push_back(itemconfig);
            }
            nrowconfig++;

    }

    //======================================//
    //---reading the couting table---//
    ifstream infile(filename);
    string line;
    char buffer[4000000]; //1gb
    infile.rdbuf()->pubsetbuf(buffer, sizeof(buffer));

    std::vector<std::string> splittedString;
    std::vector<std::string> headerInd;
    std::vector<std::string> header;
    ofstream outfile;
    selection current_sel;
    long long nrow = 0;


    while (getline (infile, line)){
        splittedString.clear();

        if(nrow == 0){
            std::istringstream byline(line);
            std::string item;

            while(getline (byline, item, ' '))
                splittedString.push_back(item);


            headerInd = getCondition(splittedString, cond1, cond2);
            //printVector(splittedString);

            //======================================//
            //---Printing the header in the file---//
            outfile.open("output.tsv");

            cond1.insert( cond1.end(), cond2.begin(), cond2.end() );
            cond1.insert(cond1.begin(), "tag");
            for (std::vector<string>::const_iterator i = cond1.begin(); i != cond1.end(); ++i)
                outfile << *i << '\t';
            outfile << endl;

        }else{

            std::istringstream byline(line);
            std::string item;

            while(getline(byline, item, ' '))
                splittedString.push_back(item);

                //printVector(splittedString);

            current_sel= processingRow(headerInd, splittedString);
            float resultCV = current_sel.cvCond;
            std::vector<int> resultVec = current_sel.vecSelection;


            if(resultCV >= 0.95){
                //cout<< resultCV <<endl;
                outfile << splittedString[0] << '\t';
                for (std::vector<int>::const_iterator i = resultVec.begin(); i != resultVec.end(); ++i)
                    outfile << *i << '\t'; // outfile
                outfile << endl;
            }
        }


        nrow++;
    }

    infile.close();
    outfile.close();

    time (&t);
    tt = localtime(&t);
    cout << "Finished at "<< asctime(tt) << endl;

    return 0;
}
