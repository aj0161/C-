// Forcasting.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Forcasting.h"
#include <math.h>       /* fabs */

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// The one and only application object
CWinApp theApp;

using namespace std;

//Function Prototype
double* Forecast_Exponential_Smoothing(double lambda,double Memory_Usages[], double* Future_vals, int length);
double Calculate_MSE(double Predicted_vals[], double Actual_vals[], int length, int weight);
double* Forecast_Moving_Average(double Memory_Usages[], double* Forcast_vals, int length, int weight);
void Calculate_Linear_Regression(double Predicted_vals[], double Actual_vals[], int length);
double Forcast_Linear_Regression(double x);
double Calculate_LinearRegression_MSE(double Predicted_vals[], double Actual_vals[], int length);
double getMean(double number[], int size);
double Calculate_LinearRegression_RSquared(double Predicted_vals[], double Actual_vals[], int length);

//global variables for temp !!
 double slope = 0;
 double intercept = 0;
 double coefficient = 0;

int _tmain(int argc, TCHAR* argv[], TCHAR* envp[])
{
	int nRetCode = 0;

	// initialize MFC and print and error on failure
	if (!AfxWinInit(::GetModuleHandle(NULL), NULL, ::GetCommandLine(), 0))
	{
		// TODO: change error code to suit your needs
		_tprintf(_T("Fatal Error: MFC initialization failed\n"));
		nRetCode = 1;
	}
	else
	{
		const int len = 7;
		double Msg_length[len] = {40,40652, 39582, 46777, 39136, 41209,38895 };  // x
		double Mem_Usage[len] = {8120,48732, 47664, 54860, 47216, 49292, 88204 };	 // y 

		double* Predicted_values = new double[len];
		double lambda = 0.9;

		//________________________Forecast_Exponential_Smoothing - Moving Average____________________________________________________
		Predicted_values= Forecast_Exponential_Smoothing(lambda, Mem_Usage, Predicted_values, len); 

		// calculate MSE !! Less the error number, better the model!!
		double MSE_error = Calculate_MSE(Predicted_values, Mem_Usage, len-1, 1);
		
		cout  << "Forecast_Exponential_Smoothing " << Predicted_values[len] << " MSE: " << MSE_error << "\n"; // the future val

		//________________________Forecast- weighted_Moving_Average________________________________________________________________________
		
		int weight = 3;// weight of weighted moving average
		Predicted_values = Forecast_Moving_Average(Mem_Usage, Predicted_values, len, weight);
		
		MSE_error = Calculate_MSE(Predicted_values, Mem_Usage, len-1, weight);
		
		cout << "Forecast_Moving_Average " << Predicted_values[len] << " MSE: " << MSE_error << "\n"; // the future val

		//________________________Linear Regression________________________________________________________________________


		//Find a linear relations between two variables.
		Calculate_Linear_Regression(Msg_length, Mem_Usage, len);

		//coeefficient strength must be greater than 0.75 to predict using linear regression.
		//higher the strength, the more precise the prediction will be
		if(fabs(coefficient) >= 0.5)
		{
			cout  << "Test the linear model !!!" << "\n";
			double Msg_length_TestingData[len] = {39136,41209,38895,41924,42090,41107};	 // testing data x
			double BufSize_ActualData[len] = {47216,49292, 88204, 130144, 91876, 1277488};	 // Actual data y
			double* Predicted_values = new double[len];

			for(int i = 0; i<len; i++) {
				Predicted_values[i] = Forcast_Linear_Regression(Msg_length_TestingData[i]);
				cout  << "Msg_length: " <<  Msg_length_TestingData[i] << " predict_y_MemUsage:" << Predicted_values[i] << " but actual is " << BufSize_ActualData[i] << "\n";
			}

			//MSE_error = Calculate_LinearRegression_MSE(Predicted_values, BufSize_ActualData, len);
			//cout <<  "linear regression MSE: " << MSE_error << "\n";  //values closer to zero are better.

			double R_squared = Calculate_LinearRegression_RSquared(Predicted_values, BufSize_ActualData, len);
			cout <<  "linear regression R_squared: " << R_squared << "\n";  //values closer to 1 are better.
		}
		else {
			cout << "linear regression prerequisites failed !! There's very low relationship between two variables " <<  "\n";
		}

		//Delete array after done!!
	}
	return nRetCode;
}

// calculate forcast using Exponential Smoothing 
// Lambda: the value is between 0 and 1, referred to as the smoothing parameter,
//		   which means that the forcast for the current peroid is obtained by adding
//		   the forcast	in the last peroid to a fraction of the error from the last 
//		   peroid	
// Formula: new Forcast_val = Lambda * Last Actual_val + (1 - Lambda)* Last Forcast_val
double* Forecast_Exponential_Smoothing(double lambda,double Memory_Usages[], double* Forcast_vals, int length){
	
	double Actual_val= 0.0; 

	for(int i = 0; i <= length; i++) {
		Actual_val = Memory_Usages[i];

		if(i >= 1 ) {
			Actual_val = Memory_Usages[i-1];
		}
		
		if(i == 0) {
			// if no forcast value for the first time is given, then we assume Forcast_val = Actual_val
			Forcast_vals[i] = Actual_val;
			}
		else {
			Forcast_vals[i] = lambda * Actual_val + (1 - lambda) * Forcast_vals[i-1];
		}

		//cout << Forcast_vals[i] << "\n";
	}
	return Forcast_vals;
}

//Calculate Moving Average with 3 items usages
// it is used to forcast by focusing on recent data.
double* Forecast_Moving_Average(double Memory_Usages[], double* Forcast_vals, int length, int weight){
	
	if(length < weight) {
		return 0;
	}

	for(int i = weight; i <= length; i++) {
		Forcast_vals[i] =  (Memory_Usages[i-3] + Memory_Usages[i-2] + Memory_Usages[i-1]) /weight;
		//cout << Memory_Usages[i] <<" Forecast_Moving_Average: "  << Forcast_vals[i] << "\n";
	}
	return Forcast_vals;
}

double Calculate_MSE(double Predicted_vals[], double Actual_vals[], int length, int weight) {
	
	double sum =0.0;
	int counter = 0;
	for(int i = weight; i<=length; i++) {
		
		double temp = fabs ( Predicted_vals[i] - Actual_vals[i]);
		temp = pow(temp, 2);
		
		 sum += temp;
		 counter++;
	}
	return (sum/counter);
}

double Calculate_LinearRegression_MSE(double Predicted_vals[], double Actual_vals[], int length) {
	
	double sum =0.0;
	int counter = 0;
	for(int i = 0; i< length; i++) {
		
		double temp = fabs ( Predicted_vals[i] - Actual_vals[i]);
		temp = pow(temp, 2);
		
		 sum += temp;
		 counter++;
	}
	return  (sqrt((sum/counter)));
}

// the coefficient of determination R2
//ref:http://www.hedgefund-index.com/d_rsquared.asp
double Calculate_LinearRegression_RSquared(double Predicted_vals[], double Actual_vals[], int length) {
	
	double mean = getMean(Actual_vals, length);
	double SSE = 0;
	double SST =0;

	for(int i = 0; i< length; i++) {
		SSE += pow(( Predicted_vals[i] - mean), 2);
		SST += pow(( Actual_vals[i] - mean), 2);
	}
	return  ((SSE/SST));
}


void Calculate_Linear_Regression(double x[], double y[], int length) {

	double xSum = 0, y2Sum = 0, x2Sum = 0, ySum =0, xySum = 0;

	for(int i =0; i < length; i++) {
		xSum = xSum + x[i];
		ySum = ySum + y[i];
		x2Sum = x2Sum + pow(x[i],2);
		y2Sum = y2Sum + pow(y[i],2);
		xySum = xySum + x[i]*y[i];
	}
     slope  = (( xSum * ySum ) - length * xySum) / ( pow(xSum,2) - length * x2Sum);

	 intercept =  ((xSum * xySum) - ySum * x2Sum )/ ( pow(xSum,2) - length * x2Sum);
	
	 double standardDev_x = sqrt((x2Sum - ( xSum * xSum)/length ) / (length -1 ));
	 double standardDev_y = sqrt((y2Sum - ( ySum * ySum)/length ) / (length -1 ));

	 coefficient = (standardDev_x/ standardDev_y) * slope;
}
 
 double Forcast_Linear_Regression(double x) {
	return (slope * x + intercept);
}

double getMean(double number[], int size){
	double sum = 0;
	for ( int i = 0; i < size; i++ ){
	sum += number[i];
	}
	return sum/size;
}
