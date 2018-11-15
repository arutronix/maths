//============================================================================
// Name        : powerful.cpp
// Author      : AI
// Version     : 1.0
// Copyright   : math for all intellectuals
// Description : Real Power Calculation in C++, ANSI-style
//============================================================================

#include <iostream>
#include "powerful.h"
using namespace std;

int main() {
	float f1,f2;
	cout << "please feed the first real number: f1 \n";
	cin >> f1;
	cout << "please feed the second real number: f2 \n";
	cin >> f2;
	cout << "f1^f2=  " <<  powerful(f1,f2) << "\n"; // prints !!!power!!!
	return 0;
}

float powerful(float powerin1,float powerin2)
    {


        //z=x^y => log(z) = y*log(x) = z' => z = e^z'
        unsigned int index = 0;
        unsigned int max_index = 50;
        unsigned int powerfularr[max_index];

        float poweroutmp;
        float powerout;
        float factinv;

        while(index<max_index)
        {
        	powerfularr[index] = 2*index + 1;
        	index++;
        }

        if(powerin1==0)
            powerout=0;

        else if(powerin1==1.0)
            powerout=1.0;

        else if(powerin2==0.0)
            powerout=1.0;

        else
        {
            // Calculation of Natural Logarithm
            powerin1=(powerin1-1)/(powerin1+1);
            poweroutmp=powerin1*powerin1;
            powerout =(1.0f/powerfularr[max_index-1]);
            index=max_index-1;
            while(index>0)
            {
                index = index-1;
                powerout = (powerout*poweroutmp) +(1.0f/powerfularr[index]);
            }
            powerout = powerout*powerin1;
            powerout=powerout*2*powerin2;
            poweroutmp = powerout;

            // Calculation of exponential value e^x
            factinv =3.8003907548547435925936708927884e-36;
            powerout = factinv;
            index=32;
            while(index>0)
            {
                factinv =index*factinv;
                powerout = (powerout*poweroutmp) + factinv;
                index = index-1;
            }

        }
        return(powerout);
    }

