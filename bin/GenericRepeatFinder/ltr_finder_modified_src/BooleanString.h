#ifndef BOOLEANSTRING_H
#define BOOLEANSTRING_H

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <assert.h>

using namespace std;

class BooleanString
{

public:
    BooleanString();
    BooleanString(const int strLength);
    ~BooleanString();

    /**********************************************************/
    /* This funciton allocates space for myBooleanString and  */
    /* modifies myStringLength and numChar                    */
    /**********************************************************/
    void initialize(const int strLength);



    /**********************************************************/
    /* Given an valid index in the BooleanString, the getVal  */
    /* function returns a true if the bit at index is 1,      */
    /* otherwise it returns a false.                          */
    /* getValN is similar to getVal, instead returning true   */
    /* or false, it returns an integer 1 or 0, respectively.  */
    /* The other difference between getVal, and getValN is    */
    /* there is no if statement in getValN, which increases   */
    /* speed in large applications.                           */
    /**********************************************************/
    bool getVal(const int index) const;
    inline int getValN(const int index) const;



    /**********************************************************/
    /* Given an integer index, setVal set the bit at the      */
    /* index to 1 if indexVal is true, to 0 if indexVal is    */
    /* false.                                                 */
    /* setValN is similar to setVal, except it sets the bit   */
    /* at the index to 1 if indexVal is 1, and to 0 if        */
    /* indexVal is 0.                                         */
    /**********************************************************/
    void setVal(const int index, const bool indexVal);
    inline void setValN(const int index, const int indexVal);



    /**********************************************************/
    /* setAll, sets all the bits in the integer array to 1 if */
    /* indexVal is true, to 0 if indexVal is false.           */
    /**********************************************************/
    void setAll(const bool indexVal);



    /**********************************************************/
    /* printAll, prints all the the bit values of the         */
    /* BooleanString. It will also print strSpace between all */
    /* values.                                                */
    /* epl in the second function denotes elements per line.  */
    /**********************************************************/
    void printAll(const char* strSpace) const;
    void printAll(const char* strSpace, const int epl) const;



    /**********************************************************/
    /* isAllTrue, returns true if all the values in the       */
    /* BooleanString are true otherwise it returns a false.   */
    /* isAllFalse, returns true if all the values in the      */
    /* BooleanString are false otherwise it returns a false.  */
    /**********************************************************/
    bool isAllTrue() const;
    bool isAllFalse() const;



private:

    char* myBooleanString;
    int myStringLength;
    int numChar;
};

#endif

inline int BooleanString::getValN(const int index) const
{
    int result;
    int charNum;
    int bitNum;
    char temp;

    assert(index >= 0 && index < myStringLength &&
           myStringLength != -1);

    charNum = index >> 3;
    bitNum = index & 7;
    temp = myBooleanString[charNum];

    temp = temp >> (7 ^ bitNum);
    temp = temp << 7;
    temp = temp >> 7;
    temp = temp & 1;

    result = (int) temp;

    return (result);
}


inline void BooleanString::setValN(const int index, const int indexVal)
{
    int charNum;
    int bitNum;
    int shiftNum;
    char temp, original;

    assert(index >= 0 && index < myStringLength &&
           myStringLength != -1);

    charNum = index >> 3;
    bitNum = index & 7;
    original = myBooleanString[charNum];
    temp = (char) indexVal;
    shiftNum = 7 ^ bitNum;

    temp = temp << shiftNum;
    temp = temp ^ original;
    temp = temp & (1 << shiftNum);
    myBooleanString[charNum] = temp ^ original;
}
