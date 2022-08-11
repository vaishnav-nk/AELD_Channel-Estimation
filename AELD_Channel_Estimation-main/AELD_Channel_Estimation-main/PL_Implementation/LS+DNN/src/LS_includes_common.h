/*
 * LS_includes.h
 *
 *  Created on: 05-Mar-2022
 *      Author: Ranjith Srinivas
 */

#ifndef SRC_LS_INCLUDES_H_
#define SRC_LS_INCLUDES_H_

#include "complex.h"
#include "DNN_Data.h"

// Number of used sub-carriers in one OFDM symbol
#define NUM_USED_SUB_CARRIERS 52
#define nUSC NUM_USED_SUB_CARRIERS
#define SYMBOLS_USED 200
#define SNR_MAX 30


typedef float complex flt_cmplx;

// Training Symbols Used
flt_cmplx training_symbol[nUSC] = {1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, -1.0+0.0*I, -1.0+0.0*I, -1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, -1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I, 1.0+0.0*I};

// Array storing the Mean computed on the 1000 samples of data of the LS estimated channel.
float Mean_X[2*nUSC] = {-0.005635,-0.018720,0.006749,-0.002547,-0.041757,-0.030383,0.026833,0.001692,0.046953,-0.000887,-0.018531,-0.037771,-0.049036,-0.046315,-0.030443,0.031100,0.047549,0.010577,-0.001389,0.017486,-0.042558,-0.011444,-0.045186,-0.026807,-0.005365,0.008863,-0.004364,0.014753,-0.003585,-0.050581,-0.042702,0.016255,-0.003242,0.048805,0.006165,0.030723,0.010599,-0.003656,-0.005184,0.005698,-0.017835,-0.004493,0.039828,0.034683,-0.015642,0.014595,-0.031250,0.007837,-0.029570,-0.024683,-0.018596,-0.012413,0.011855,0.005550,-0.045617,-0.033131,0.037708,0.049927,-0.002711,0.042243,-0.030171,0.007264,-0.001019,0.004368,0.021409,0.042059,0.056081,0.001481,-0.012685,0.020394,0.001940,-0.058866,0.005203,-0.029451,0.046257,0.057958,0.054895,0.038731,0.000842,-0.057629,-0.046877,0.027956,0.047699,0.002691,0.052322,-0.019774,0.014237,-0.054016,-0.053671,-0.040317,-0.021597,-0.007415,0.050376,0.042251,-0.026152,-0.037997,0.011589,-0.041093,0.020157,-0.027116,0.033723,0.037517,0.038575,0.036365};

// Array storing the standard deviation computed on the 1000 samples of data of the LS estimated channel.
float SD_X[2*nUSC] = {1.000495,0.996002,0.986994,0.992332,0.975269,0.978243,0.990477,0.985497,0.972888,1.003057,0.998725,0.983856,0.971080,0.968604,0.971268,0.988302,0.970881,1.004150,1.008098,0.972477,0.975178,0.999454,0.971817,0.970904,0.977404,0.994537,1.007098,0.975490,0.991513,0.976032,0.970862,1.005752,0.989140,0.976829,1.010698,0.978475,0.985836,0.995088,1.003665,1.003134,0.981899,0.986681,0.985826,0.983105,0.996388,0.979819,0.991780,0.995376,0.976070,0.976149,0.984106,0.990955,1.010261,0.987365,1.029934,1.032121,0.984150,1.003996,0.974243,1.031226,0.981025,0.994925,0.977888,0.974924,0.980349,0.989176,1.004095,0.979386,0.981001,1.008968,0.991532,1.012977,0.983192,1.014882,0.985981,1.005408,1.020706,1.020219,1.003186,1.017016,1.023297,0.985781,0.997055,0.997881,1.021848,0.989668,1.011815,1.004966,1.015829,1.017818,1.009504,1.002296,1.007396,1.009616,0.997066,0.995940,1.009257,1.000504,1.000649,1.006989,0.992510,0.991754,0.994791,1.005704};

// Array storing the Mean computed on the 1000 samples of data of the actual channel.
float Mean_Y[2*nUSC] = {0.014400,0.001316,-0.013286,-0.022582,-0.021722,-0.010347,0.006797,0.021728,0.026917,0.019148,0.001504,-0.017736,-0.029000,-0.026280,-0.010408,0.011065,0.027513,0.030613,0.018646,-0.002550,-0.022522,-0.031479,-0.025150,-0.006772,0.014670,0.028899,0.015671,-0.005283,-0.023621,-0.030546,-0.022666,-0.003780,0.016793,0.028770,0.026200,0.010687,-0.009437,-0.023691,-0.025220,-0.014337,0.002200,0.015542,0.019792,0.014647,0.004394,-0.005441,-0.011215,-0.012198,-0.009535,-0.004647,0.001439,0.007622,-0.015515,-0.021819,-0.018247,-0.005761,0.010338,0.022557,0.024659,0.014873,-0.002801,-0.020106,-0.028389,-0.023002,-0.005961,0.014689,0.028711,0.028851,0.014685,-0.006976,-0.025430,-0.031496,-0.022167,-0.002081,0.018887,0.030588,0.027525,0.011361,-0.026528,-0.030260,-0.019508,0.000587,0.020329,0.030061,0.024952,0.007596,-0.013132,-0.026646,-0.026301,-0.012947,0.005773,0.019954,0.023006,0.014881,0.001218,-0.010627,-0.015780,-0.013723,-0.007213,0.000254,0.006353,0.010147,0.011205,0.008995};

// Array storing the standard deviation computed on the 1000 samples of data of the actual channel.
float SD_Y[2*nUSC] = {0.706533,0.706869,0.701273,0.698893,0.703418,0.705840,0.700929,0.696279,0.699233,0.704591,0.702797,0.695646,0.695241,0.703334,0.705251,0.695843,0.691908,0.701583,0.706451,0.696790,0.691166,0.700681,0.707376,0.700298,0.695158,0.702559,0.705792,0.702542,0.705966,0.709547,0.709746,0.708863,0.707261,0.707044,0.710348,0.711056,0.705831,0.703964,0.709093,0.710244,0.704334,0.702801,0.708075,0.708875,0.703540,0.702874,0.707603,0.707320,0.702443,0.702806,0.706871,0.705842,0.701943,0.700689,0.705239,0.706477,0.700843,0.697340,0.701342,0.705181,0.701628,0.695761,0.697194,0.704015,0.704117,0.695715,0.693419,0.702491,0.706012,0.696112,0.691034,0.700871,0.706818,0.698156,0.692547,0.701247,0.708208,0.702998,0.704305,0.709681,0.708164,0.706117,0.707037,0.708593,0.710428,0.710496,0.706736,0.705352,0.709796,0.710841,0.704956,0.703122,0.708485,0.709557,0.703926,0.702798,0.707824,0.708145,0.703031,0.702882,0.707296,0.706481,0.702005,0.702725};

#endif /* SRC_LS_INCLUDES_H_ */
