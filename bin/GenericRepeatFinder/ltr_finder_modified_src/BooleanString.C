#include "BooleanString.h"


BooleanString::BooleanString()
{
    myStringLength = -1;
    numChar = -1;
}


BooleanString::BooleanString(const int strLength)
{
    myStringLength = strLength;
    numChar = strLength / 8;

    if (strLength % 8 != 0)
    {
        numChar++;
    }

    myBooleanString = (char*) malloc(sizeof(char) * numChar);
}


BooleanString::~BooleanString()
{
    if (myStringLength != -1)
    {
        delete[] myBooleanString;
    }
}


void BooleanString::initialize(const int strLength)
{
    myStringLength = strLength;
    numChar = strLength / 8;

    if (strLength % 8 != 0)
    {
        numChar = numChar + 1;
    }

    myBooleanString = (char*) malloc(sizeof(char) * numChar);
}


bool BooleanString::getVal(const int index) const
{
    int charNum;
    int bitNum;
    char temp;

    assert(index >= 0 && index < myStringLength &&
           myStringLength != -1);

    //charNum = index / 8;
    charNum = index >> 3; //;
    //bitNum = index % 8;
    bitNum = index - (charNum << 3);
    temp = myBooleanString[charNum];

    temp = temp >> 7 - bitNum;
    temp = temp << 7;
    temp = temp >> 7;

    if ((int)temp == 0)
    {
        return (false);
    }
    else if ((int)temp == -1)
    {
        return (true);
    }
    else
    {
        cout << "wrong result" << endl;
        return (false);
    }
}


void BooleanString::setVal(const int index, const bool indexVal)
{
    int charNum;
    int bitNum;
    char temp;

    assert(index >= 0 && index < myStringLength &&
           myStringLength != -1);

    //charNum = index / 8;
    //bitNum = index % 8;
    charNum = index >> 3; //;
    bitNum = index - (charNum << 3);

    if (indexVal == true)
    {
        temp = 0x01;
        temp = temp << 7 - bitNum;
        myBooleanString[charNum] = myBooleanString[charNum] | temp;
    }
    else
    {
        temp = 0x01;
        temp = temp << 7 - bitNum;
        temp = ~temp;
        myBooleanString[charNum] = myBooleanString[charNum] & temp;
    }
}


void BooleanString::setAll(const bool indexVal)
{
    assert(myStringLength > 0);

    if (indexVal == true)
    {
        memset((void*)myBooleanString, 255, numChar);
    }
    else
    {
        memset((void*)myBooleanString, 0, numChar);
    }
}


void BooleanString::printAll(const char* strSpace) const
{
    int i;

    assert(myStringLength > 0);

    for (i = 0; i < myStringLength; i++)
    {
        cout << getVal(i) << strSpace;
    }

    cout << endl;
}


void BooleanString::printAll(const char* strSpace,
                             const int epl) const
{
    int i;

    assert(myStringLength > 0);

    for (i = 0; i < myStringLength; i++)
    {
        cout << getVal(i) << strSpace;

        if (i % epl == (epl - 1))
        {
            cout << endl;
        }

    }
    cout << endl;
}


bool BooleanString::isAllTrue() const
{
    int i;

    assert(myStringLength > 0);

    for (i = 0; i < myStringLength / 8; i++)
    {
        if (myBooleanString[i] != (char) 255)
        {
            return (false);
        }

    }

    for (i = 0; i < (myStringLength % 8); i++)
    {
        if (getVal(myStringLength - i - 1) == false)
        {
            return (false);
        }

    }
    return (true);
}


bool BooleanString::isAllFalse() const
{
    int i;

    assert(myStringLength > 0);

    for (i = 0; i < myStringLength / 8; i++)
    {
        if (myBooleanString[i] != (char) 0)
        {
            return (false);
        }

    }

    for (i = 0; i < (myStringLength % 8); i++)
    {
        if (getVal(myStringLength - i - 1) == true)
        {
            return (false);
        }

    }
    return (true);
}

