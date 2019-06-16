#include <ctype.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "BooleanString.h"
#include "Timing.h"
#define INT_SIZE 32
#define MAX_INT_INDEX 31
#define MAX_NUM_CHAR 256
#define MAX_CHAR_VAL 255

/*********************************
Reads input from the file associated with inFile.
Pre:    inFile is associated with a file, i.e. open is called.
Post:   inputString holds the characters from the file,
        inputLength is the length of the string from the file.
Modify: inputString, inputLength.
*********************************/
void read_input(ifstream& inFile, char*& inputString, int& inputLength);



/*********************************
Similar to read_input, but this function reads the first sequence
in a fasta format file.
*********************************/
void read_fasta(ifstream& inFile, char*& inputString, int& inputLength);



/*********************************
Converts a char string to a int string.
Pre:    inputString is of length stringLength, and stringLength
        is non-zero and non-negative.
Post:   returns a pointer to a string who's length is stringLength
        the character A in the original string is given the value
 0. The value of other characters are calculated by subtracting
 the value of A, e.g. C has a value of 2. A value of 0 is
 added to the end of the string this acts like the character
 '$', and the length of string is incremented by 1.
Modify: stringLength++.
*********************************/
int* charToInt(const char* inputString, int& stringLength);



/*********************************
Same as the other charToInt function expect let user specify, the
smallest character by using baseChar. 
*********************************/
int* charToInt(const char* inputString, int& stringLength,
               const char baseChar);



/*********************************
Preforms the linear time suffix sort for string of integer alphabet.
Pre:    inputString is of length stringLength, and stringLength
        is non-zero and non-negative.
 inputString[stringLength - 1], i.e. the last integer of
 the string must be the smallest integer in the entire
 inputString, and it must be unique.
Post:   Returns an array A, where A[i] give the index of the
        i-th smallest suffix in the string.
Modify: None 
*********************************/
int* LinearSuffixSort(const int* inputString, const int stringLength);

/*********************************
Preforms the linear time suffix sort knowing minChar and maxChar
*********************************/
int* LinearSuffixSort(const int* inputString, const int stringLength,
                      const int minChar, const int maxChar);

/*********************************
Same as int* LinearSuffixSort(const int*, const int); but for ASCII 
characters.
Pre:    inputString is of length stringLength, and stringLength
        is non-zero and non-negative.
        ASCII character '0' is not a part of inputString.
Post:   Returns an array A, where A[i] give the index of the
        i-th smallest suffix in the string.
        ASCII character '0' is appended to the end of the inputString.
        stringLength increases by 1. 
Modify: inputString, stringLength.
*********************************/
int* LinearSuffixSort(char*& inputString, int& stringLength);







/*******************************************************************/
/* The functions in this section are functions called by           */
/* LinearSuffixSort, or helper functions during debugging, users   */
/* Should not use those functions in their program.                */
/*******************************************************************/





/*********************************
These functions are for debugging purpose.                      
These functions output the content of a string                  
*********************************/
void outstring(const int* inputString, const int stringLength);
void outstring3(const int* inputString, const int stringLength);
void outstringS(const int* inputString, const int stringLength);
void outstring(const char* inputString, const int stringLength,
               const char* strSpace);
void outsuffix(const int* intString, const int intStringLen,
               const char* charString, const int charStringLen,
               const int suffixLength);
void outsubstring(const int* intString, const int intStringLen,
                  const char* charString, const int charStringLen,
                  const BooleanString& suffixType);

/*********************************
Determines the type of each suffix. 
If suffixType[i] is true(1) then T[i] is a type S suffix.
Else if suffixType[i] is false(0) then T[i] is a type L suffix.
Pre:    inputString is a integer array of length inputLength;
        numStype, numLtype are integers.
 suffixType is a BooleanString of length inputLength;
Post:   Calculate the type of each suffix in string inputString.
        The type of the last suffix depends on which type of
 suffixes is less in number excluding the last suffix.
 Therefore, user should not compare numStype and numLtype,
 instead it should look at the type of the last suffix.
        numStype holds the number of type S suffixes;
        numLtype holds the number of type L suffixes;
 suffixType[i] gives the type of each suffixes.
Modify: numLtype, numStype, suffixType.
*********************************/
void suffix_type(const int* inputString, const int inputLength,
                 int& numStype, int& numLtype,
                 BooleanString& suffixType);

/*********************************
Determines the type of each suffix. Except the inputString is
a character string of length inputLength
*********************************/
void suffix_type(const char* inputString, const int inputLength,
                 int& numStype, int& numLtype,
                 BooleanString& suffixType);



/*********************************
Find the maximum and minimum value of an integer array.
Pre:    inputString is of length inputLength, Max and Min are
        integers.
Post:   Max will contain the value of the maximum element.
        Min will contain the value of the minimum element.
Modify: Max, Min.
*********************************/
inline void findMaxMin(const int* inputString, const int inputLength,
                       int& Max, int& Min)
{
    int i, temp;

    if (inputLength <= 0)
    {
        cout << "In function findMaxMin(int*, int, int&, int&):" << endl;
        cout << "Length of input string cannot be less than 0."
        << endl;
        exit(0);
    }

    Max = inputString[0];
    Min = inputString[0];

    for (i = 1; i < inputLength; i++)
    {
        temp = Max - inputString[i];
        temp = temp >> MAX_INT_INDEX;
        Max = Max + ((inputString[i] - Max) & temp);

        temp = inputString[i] - Min;
        temp = temp >> MAX_INT_INDEX;
        Min = Min - ((Min - inputString[i]) & temp);
    }
}

/*********************************
Same as findMaxMin, but for character strings.
*********************************/
inline void findMaxMin(const char* inputString, const int inputLength,
                       char& Max, char& Min)
{
    int i;
    char temp;

    if (inputLength <= 0)
    {
        cout << "In function findMaxMin(int*, int, int&, int&):" << endl;
        cout << "Length of input string cannot be less than 0."
        << endl;
        exit(0);
    }

    Max = inputString[0];
    Min = inputString[0];

    for (i = 1; i < inputLength; i++)
    {
        temp = Max - inputString[i];
        temp = temp >> MAX_INT_INDEX;
        Max = Max + ((inputString[i] - Max) & temp);

        temp = inputString[i] - Min;
        temp = temp >> MAX_INT_INDEX;
        Min = Min - ((Min - inputString[i]) & temp);
    }
}



/*********************************
Do counting sort on inputString. This is the same as sorting
each suffix according to their first character.
Pre:    inputString is of length inputLength;
        A is a pointer to an array of length inputLength;
 A must be initalized with "new".
 BuckA is a BooleanString of length inputLength.
 buffer is a int array of length inputLength.
Post:   A contains the index of the suffixes sorted
        according to the first character.
 BuckA[i] is set to true, if there is an bucket bundary
 between A[i] and A[i+1]. The last element of BuckA
 will always be 1, because it's the end of the array.
Modify: A, BuckA.
*********************************/
void counting_sort(const int* inputString, const int inputLength,
                   int* A, BooleanString& BuckA, int* buffer);

/*********************************
Same as counting_sort, but min & max char is provided.
*********************************/
void counting_sort(const int* inputString, const int inputLength,
                   int* A, BooleanString& BuckA, const int minChar,
                   const int maxChar, int* buffer);



/*********************************
Sort a bucket of type S substrings using counting sort.
Pre:    inputString is a character string of inputLength,
        A is an array containing the begining position of
        all type S substrings. 
        A is of length ALength - the number of type S substrings.
        BuckA is of length ALength.
        suffixType contains the type of all suffixes of inputString.
        maxDist is the maximum reverse S-distance of inputString.
Post:   A contains the sorted order of all type S substrings.
        BuckA contains the bucket boundary of A.
Modify: A, BuckA.
*********************************/
void sort_s_substringC(const char* inputString, const int inputLength,
                       int* A, BooleanString& BuckA, const int ALength,
                       const BooleanString& suffixType, const int maxDist);

/*********************************
Same as sort_s_substringC but for type L substrings.
*********************************/
void sort_l_substringC(const char* inputString, const int inputLength,
                       int* A, BooleanString& BuckA, const int ALength,
                       const BooleanString& suffixType, const int maxDist);



/*********************************
Compute the S-Distance of all suffixes according to algorithm.
Pre:    suffixType contains the type of each suffix. It is of
        length inputLength.
 Dist is initialized with "new", and have length 
 inputLength;
 DistCount is not initalized.
 maxDist is a integer.
Post:   Dist will contain the S-Distance of all suffixes.
        DistCount[i] give the number of suffixes that have
 S-Distance less than or equal to i, for 0 <= i < maxDist. 
        DistCount[maxDist] give the total number of non-zero
 S-Distance suffixes. DistCount have length maxDist+1.
        maxDist gives the maximum S-Distance length. 
 Note: DistCount[maxDist] = 0 if only 1 type S suffix exists
       in the string, maxDist is also 0 in this case. Thus
       in the case only 1 type S suffix exists this function
       should not be called.
Modify: Dist, DistCount, maxDist.
*********************************/
void s_distance(const BooleanString& suffixType,
                const int inputLength, int* Dist, int*& DistCount,
                int& maxDist);

/*********************************
Similar to function s_distance
but this calculates the l-distance of each suffix.
*********************************/
void l_distance(const BooleanString& suffixType,
                const int inputLength, int* Dist, int*& DistCount,
                int& maxDist);

/*********************************
Compute the reverse S-Distance for all type S suffix. 
i.e. the distance of a type S
suffix from the nearest type S suffix on the right instead of
left. The reverse S-Distance of the last suffix '$' is defined
to be 0. This function only returns the maximum S-Distance.
Pre:    suffixType contains the type of each suffix. It is of
        length inputLength.
        inputLength is the length of the suffixType. 
 maxDist is a integer.
Post:   maxDist gives the maximum S-Distance length. 
Modify: Dist, DistCount, maxDist.
*********************************/
void s_distanceR(const BooleanString& suffixType,
                 const int inputLength, int& maxDist);

/*********************************
Similar to function s_distance(BooleanString,Char,int,int,int)
but this calculates the reverse l-distance of each type L suffix.
*********************************/
void l_distanceR(const BooleanString& suffixType,
                 const int inputLength, int& maxDist);



/*********************************
Construct the ArrayB of type S from ArrayA.
Pre:    ArrayA is an array with all the suffixes sorted according to
        their first character.
 inputLength determines the length of ArrayA.
 BuckA is of length inputLength and BuckA[i] is true if 
 there is a bucket bundary between ArrayA[i] and ArrayA[i+1]
 The result will be stored in ArrayB, and BuckB marks the 
 bucket bundary of ArrayB. 
 suffixType[i] give the value of the i-th suffix.
Post:   BuckB have all the suffixes sorted according to their
        first character. BuckB marks the bucket bundary of ArrayB.
 BuckB[i] = true if there is a boundary between ArrayB[i]
 and ArrayB[i+1].
Modify: ArrayB, BuckB.
*********************************/
void construct_ArrayB_typeS(const int* ArrayA, const int inputLength,
                            const BooleanString& BuckA,
                            int* ArrayB, BooleanString& BuckB,
                            const BooleanString& suffixType);

/*********************************
Similar to construct_ArrayB_typeS, but this function calculates the 
ArrayB of type L from ArrayA. 
*********************************/
void construct_ArrayB_typeL(const int* ArrayA, const int inputLength,
                            const BooleanString& BuckA,
                            int* ArrayB, BooleanString& BuckB,
                            const BooleanString& suffixType);



/*********************************
Construct the m Lists.
Pre:    ArrayA is an integer array of size inputLength.
        ArrayA is all the suffixes bucketed according to their
 first character.
 Dist is an integer array of size inputLength.
 Dist contains all the S-distance, or L-distance of
 all the suffixes.
 DistCount is of size maxDist+1.
 DistCount[i] give the number of suffixes having
 S-distance of i or less (not counting 0).
 BuckA is of length inputLength.
 BuckA give the bucket boundaries of ArrayA. i.e.
 BuckA[i] = 1 if there is a boundary between
 ArrayA[i] and ArrayA[i+1].
 BuckA[inputLength - 1] is always 1.
 BuckList is of length listLength, listLength must 
 not be 0.
 listLength is equal to DistCount[maxDist], which is
 potentially 0 in some case, refer to s_distance().
Post:   ArrayA is not pointing to the List.
        Dist is no longer vaid, and is the reverse mapping
 array of the List.
 DistCount[i] now give the number of suffixes of
 S-distance of i+1 or less.
 BuckList give the bucket boundary of all buckets.
Modify: ArrayA, Dist, BuckList, DistCount
*********************************/
int* construct_list_typeS(int*& ArrayA, const int inputLength, int* Dist,
                          int* DistCount, const int maxDist,
                          const BooleanString& BuckA,
                          BooleanString& BuckList, const int listLength);

/*********************************
Similar to function construct_list_typeS
but this is for constructing the list for type L substrings.
*********************************/

int* construct_list_typeL(int*& ArrayA, const int inputLength, int* Dist,
                          int* DistCount, const int maxDist,
                          const BooleanString& BuckA,
                          BooleanString& BuckList, const int listLength);



/*********************************
Bucket the substrings by refering to the m Lists. This is for
sorting type S substrings.
Pre:    ArrayB is an integer array of size ArrayBLength.
        ArrayB is all the type S or type L substrings bucketed
 according to their first character.
 BuckB is the bucket boundaries of ArrayB.
 BuckB[i] = true if there is a boundary between ArrayB[i]
 and ArrayB[i+1].
 List is the m Lists, it is of size listLength.
 Each element List[i] refer to the beginning index of 
 the type S or type L substring it is in. 
 BuckList gives the boundaries of List.
 BuckList[i] = true if there is a boundary between List[i]
 and List[i+1].
 inputLength is the length of the string. 
 listLength is the length of List.
 ArrayBLength is the length of ArrayB.
Post:   ArrayB is all the type S or type L substrings bucketed.
        BuckB is the bucket boundaries of ArrayB.
 BuckB[i] = true if there is a boundary between ArrayB[i]
 and ArrayB[i+1].
Modify: ArrayB, BuckB.
*********************************/
void sort_by_list_typeS(int* ArrayB, BooleanString& BuckB,
                        int* List, const BooleanString& BuckList,
                        const int inputLength, const int listLength,
                        const int ArrayBLength);

/*********************************
Simular to function sort_by_list_typeS
but this is for sorting the type L substrings.
*********************************/
void sort_by_list_typeL(int* ArrayB, BooleanString& BuckB,
                        int* List, const BooleanString& BuckList,
                        const int inputLength, const int listLength,
                        const int ArrayBLength);



/*********************************
Construct the suffix array, by using the sorted order 
of type S suffixes.
Pre:    ArrayB is the array with the sorted order of
        type S suffixes, and it's of length
 ArrayBLength.
 stringT is the orginal string, and its length
 is inputLength.
 The smallest element of stringT is assumed to 
 be 0.
 And the largest element of stringT is assumed to
 be < inputLength.
 suffixType is the BooleanString with the type
 of each suffix. suffixType[i] = true iff
 stringT[i] is a type S suffix, otherwise
 suffixType[i] = false.
 suffixArray is an empty array of length
 inputLength, suffixArray must be initalized
 with operator "new". 
Post:   suffixArray will contain the sorted order of 
        all suffixes of stringT. 
Modify: suffixArray.
*********************************/
void construct_SA_typeS(int* ArrayB, const int ArrayBLength,
                        const int* stringT, const int inputLength,
                        const BooleanString& suffixType,
                        int* suffixArray);

/*********************************
Similar to construct_SA_typeS, except this function
is for using sorted order of type L suffixes to obtain
the order of the suffix array.
*********************************/
void construct_SA_typeL(int* ArrayB, const int ArrayBLength,
                        const int* stringT, const int inputLength,
                        const BooleanString& suffixType,
                        int* suffixArray);

/*********************************
counterpart of construct_SA_typeS for int stringT
*********************************/
void construct_SA_typeS(int* ArrayB, const int ArrayBLength,
                        const char* stringT, const int inputLength,
                        const BooleanString& suffixType,
                        int* suffixArray);

/*********************************
counterpart of construct_SA_typeS for int stringT
*********************************/
void construct_SA_typeL(int* ArrayB, const int ArrayBLength,
                        const char* stringT, const int inputLength,
                        const BooleanString& suffixType,
                        int* suffixArray);



/*********************************
Construct the new string T' from array B
Pre:    ArrayB is all the type S substrings bucketed
        lexicographically.
 ArrayB is of length ArrayBLength.
 BuckB marks the bucket boundaries of ArrayB.
 tPrime an integer array of length ArrayBLength.
 Memory space must be allocated for tPrime before
 the function call. 
 inputLength is the length of the original string.
 suffixType have all the types of all the suffixes.
Post:   tPrime is the new string, and tPrime is of length
        ArrayBLength.
Modify: tPrime.
*********************************/
void construct_TPrime_typeS(int* ArrayB, const int ArrayBLength,
                            const BooleanString& BuckB,
                            int* tPrime, const int inputLength,
                            const BooleanString& suffixType);

/*********************************
Similar to construct_TPrime_typeS, instead this function
construct the new string for all the type L suffixes.
*********************************/
void construct_TPrime_typeL(int* ArrayB, const int ArrayBLength,
                            const BooleanString& BuckB,
                            int* tPrime, const int inputLength,
                            const BooleanString& suffixType);

/*********************************
Construct the new string T' from array B
Pre:    ArrayB is all the type S substrings bucketed
        lexicographically.
 ArrayB is of length ArrayBLength.
 BuckB marks the bucket boundaries of ArrayB.
 tPrime an integer array of length ArrayBLength.
 Memory space must be allocated for tPrime before
 the function call. 
 inputLength is the length of the original string.
 suffixType have all the types of all the suffixes.
Post:   tPrime is the new string, and tPrime is of length
        ArrayBLength.
        Returns the length of tPrime, because some of the
        type S suffix may not be needed. 
Modify: tPrime.
*********************************/
int construct_TPrime_NU_S(int* ArrayB, const int ArrayBLength,
                          const BooleanString& BuckB,
                          int* tPrime, const int inputLength,
                          const BooleanString& suffixType);



/*********************************
Given the sorted order of all type S suffixes as a suffix in
T', the sorted order is indexed as suffixes of T' not as
type S suffixes of T. Convert the index back to reflect the
position in T. 
Pre:    ArrayB is of length ArrayBLength.
        ArrayB contains the sorted order of all suffixes of
 T'.
 suffixType gives the type of all suffixes of T.
 inputLength is the length of T or suffixType.
Post:   ArrayB now contains the sorted order of all 
        type S suffixes of T. i.e. the index is converted 
 from T' to T. 
Modify: ArrayB.
*********************************/
void reconstruct_B_typeS(int* ArrayB, const int ArrayBLength,
                         const BooleanString& suffixType,
                         const int inputLength);

/*********************************
Same as reconstruct_B_typeS, instead this function operates
on type L suffixes.
*********************************/
void reconstruct_B_typeL(int* ArrayB, const int ArrayBLength,
                         const BooleanString& suffixType,
                         const int inputLength);



void read_input(ifstream& inFile, char*& inputString, int& inputLength)
{
    int i, count;
    char temp;

    inputLength = 2048;
    inputString = (char*) malloc(inputLength * sizeof(char));

    i = 0;
    count = 0;
    inFile.get(temp);

    while (!inFile.eof())
    {
        if (isalpha(temp) || isdigit(temp))
        {
            inputString[i] = toupper(temp);
            i++;
        }
        if (i == inputLength - 1)
        {
            inputLength += 256;
            inputString = (char*) realloc((void*)inputString,
                                          inputLength * sizeof(char));
        }

        inFile.get(temp);
    }

    inputLength = i;
}

void read_fasta(ifstream& inFile, char*& inputString, int& inputLength)
{
    char temp;

    // skip the description of the sequence, if there is any.

    inFile.get(temp);

    while (isalpha(temp) == 0 && temp != '>' && !inFile.eof())
    {
        inFile.get(temp);
    }

    if (temp == '>')
    {
        inFile.get(temp);

        while (temp != '\n')
        {
            inFile.get(temp);
        }

        read_input(inFile, inputString, inputLength);
    }
    else if (isalpha(temp) != 0)
    {
        read_input(inFile, inputString, inputLength);
    }
    else if (inFile.eof())
    {
        exit(0);
    }
}
void suffix_type(const int* inputString, const int inputLength,
                 int& numStype, int& numLtype,
                 BooleanString& suffixType)
{
    int i, j, k;

    j = 0;
    numStype = 0;
    numLtype = 0;

    suffixType.setAll(true);

    for (i = 0; i < inputLength - 1; i++)
    {
        /****************************************/
        //If character T[i] does not equal to T[i+1] make a
        //decision and go back to mark all the previously
        //undecided suffix a certain type. Also increment the
        //counter for that type. And set undecided suffix
        //counter (j) to 0.
        /****************************************/

        if (inputString[i] < inputString[i + 1])
        {
            for (k = i - j; k <= i; k++)
            {
                numStype++;
                suffixType.setValN(k, 1);
            }

            j = 0;
        }
        else if (inputString[i] > inputString[i + 1])
        {
            for (k = i - j; k <= i; k++)
            {
                numLtype++;
                suffixType.setValN(k, 0);
            }

            j = 0;
        }
        else
        {
            /**************************************/
            //If two adjacent suffixes have the same first
            //character, move on, to the next, but remember the
            //number of undecided suffixes by increment j.
            /**************************************/
            j++;
        }

    }

    /******************************************/
    //The last suffix $ must be selected no matter we choose
    //to sort L or to sort S. So the type of the last suffix
    //is set to which ever type that is smaller in number.
    /******************************************/

    if (numStype < numLtype)
    {
        suffixType.setValN(inputLength - 1, 1);
        numStype++;
    }
    else
    {
        suffixType.setValN(inputLength - 1, 0);
        numLtype++;
    }
}




void suffix_type(const char* inputString, const int inputLength,
                 int& numStype, int& numLtype,
                 BooleanString& suffixType)
{
    int i, j, k;

    j = 0;
    numStype = 0;
    numLtype = 0;

    suffixType.setAll(true);

    for (i = 0; i < inputLength - 1; i++)
    {
        /****************************************/
        //If character T[i] does not equal to T[i+1] make a
        //decision and go back to mark all the previously
        //undecided suffix a certain type. Also increment the
        //counter for that type. And set undecided suffix
        //counter (j) to 0.
        /****************************************/

        if (inputString[i] < inputString[i + 1])
        {
            for (k = i - j; k <= i; k++)
            {
                numStype++;
                suffixType.setValN(k, 1);
            }

            j = 0;
        }
        else if (inputString[i] > inputString[i + 1])
        {
            for (k = i - j; k <= i; k++)
            {
                numLtype++;
                suffixType.setValN(k, 0);
            }

            j = 0;
        }
        else
        {
            /**************************************/
            //If two adjacent suffixes have the same first
            //character, move on, to the next, but remember the
            //number of undecided suffixes by increment j.
            /**************************************/
            j++;
        }

    }
    /******************************************/
    //The types of all characters, except the last one is set.
    //Because we are comparing inputString[i] and
    //inputString[i+1] so even if the last few characters
    //are the same, the very last character will determine
    //their type.
    /******************************************/

    /******************************************/
    //The last suffix $ must be selected no matter we choose
    //to sort L or to sort S. So the type of the last suffix
    //is set to which ever type that is smaller in number.
    /******************************************/

    if (numStype < numLtype)
    {
        suffixType.setValN(inputLength - 1, 1);
        numStype++;
    }
    else
    {
        suffixType.setValN(inputLength - 1, 0);
        numLtype++;
    }
}



int* LinearSuffixSort(char*& inputString, int& stringLength)
{
    int numS, numL, maxDistance, tPrimeLen;
    int* arrayA;
    int* arrayB;
    int* tPrime;
    BooleanString suffixType((stringLength + 1));
    BooleanString BuckB;

    //Append char 0 to the end of inputString.
    //And some initialization;

    // commented by Zhao XU, strings in C/C++ always has '\0' at end
    // make sure the caller has alloc one more byte for inputString
    //  inputString = (char*) realloc((void*)inputString,
    //    sizeof(char) * (stringLength + 1));
    //

    stringLength++;
    inputString[stringLength - 1] = (char) 0;

    //  outstring(inputString, stringLength, "  ");

    //Start of real program;
    suffix_type(inputString, stringLength, numS, numL, suffixType);

    if (suffixType.getValN(stringLength - 1) == 1 &&
            numS == 1)
    {
        arrayB = (int*) malloc(sizeof(int));
        arrayB[0] = stringLength - 1;
        arrayA = (int*) malloc(stringLength * sizeof(int));
        construct_SA_typeS(arrayB, numS, inputString, stringLength,
                           suffixType, arrayA);
        free(arrayB);
        return (arrayA);
    }

    if (suffixType.getValN(stringLength - 1) == 1)
    {
        //Less type S suffixes

        /*
        //This if is for debug.
        if (suffixType.getVal(stringLength - 1) == false)
        {
          suffixType.setVal(stringLength - 1, true);
          numS++;
        }
        */
        //    cout << "Size of arrayB: " << sizeof(int)*numS << endl;
        arrayB = (int *) malloc (sizeof(int) * numS);
        BuckB.initialize(numS);

        s_distanceR(suffixType, stringLength, maxDistance);
        //    suffixType.printAll("  ");
        //    cout << maxDistance << endl;

        //assign arrayB[0] to the index of '$'.
        //and for array[j], j >= 1, if inputString[i] is type S,
        //set array[j] = i, and j++.

        int j = 0;

        for (int i = 0; i < stringLength; i++)
        {
            assert(j <= numS);
            arrayB[j] = i;
            j = j + suffixType.getValN(i);
        }

        sort_s_substringC(inputString, stringLength, arrayB, BuckB, numS,
                          suffixType, maxDistance);
        /*
        //This part prints debug information
        outstring3(arrayB, numS);
        BuckB.printAll("  ");
        cout << endl;
        suffixType.printAll(" ");
        outsubstring(arrayB, numS, inputString, stringLength, suffixType);
        */

        if (BuckB.isAllTrue())
        {
            arrayA = (int*) malloc(stringLength * sizeof(int));
            construct_SA_typeS(arrayB, numS, inputString, stringLength,
                               suffixType, arrayA);
            free((void*)arrayB);
            return (arrayA);
        }

        tPrime = (int*) malloc(sizeof(int) * numS);
        construct_TPrime_typeS(arrayB, numS, BuckB, tPrime,
                               stringLength, suffixType);
        //    outstring3(tPrime, numS);

        free((void*)arrayB);

        arrayB = LinearSuffixSort(tPrime, numS);
        free((void*)tPrime);

        reconstruct_B_typeS(arrayB, numS, suffixType, stringLength);
        arrayA = (int*) malloc(stringLength * sizeof(int));
        construct_SA_typeS(arrayB, numS, inputString, stringLength,
                           suffixType, arrayA);
        free((void*)arrayB);
        return (arrayA);
    }
    else
    {
        //Less type L suffixes

        arrayB = (int *) malloc (sizeof(int) * numL);
        BuckB.initialize(numL);

        l_distanceR(suffixType, stringLength, maxDistance);

        //    suffixType.printAll("  ");
        //    cout << maxDistance << endl;

        //assign arrayB[0] to the index of '$'.
        //and for array[j], j >= 1, if inputString[i] is type S,
        //set array[j] = i, and j++.

        int j = 0;

        for (int i = 0; i < stringLength; i++)
        {
            assert(j <= numL);
            arrayB[j] = i;
            j = j - (suffixType.getValN(i) - 1);
        }

        sort_l_substringC(inputString, stringLength, arrayB, BuckB, numL,
                          suffixType, maxDistance);
        /*
        //This part prints debug information
        outstring3(arrayB, numL);
        BuckB.printAll("  ");
        cout << endl;
        suffixType.printAll(" ");
        outsubstring(arrayB, numL, inputString, stringLength, suffixType);
        */

        if (BuckB.isAllTrue())
        {
            arrayA = (int*) malloc(stringLength * sizeof(int));
            construct_SA_typeL(arrayB, numL, inputString, stringLength,
                               suffixType, arrayA);
            free((void*)arrayB);
            return (arrayA);
        }

        tPrime = (int*) malloc(sizeof(int) * numL);
        construct_TPrime_typeL(arrayB, numL, BuckB, tPrime,
                               stringLength, suffixType);
        //    outstring3(tPrime, numL);

        free((void*)arrayB);

        arrayB = LinearSuffixSort(tPrime, numL);
        free((void*)tPrime);

        reconstruct_B_typeL(arrayB, numL, suffixType, stringLength);
        arrayA = (int*) malloc(stringLength * sizeof(int));
        construct_SA_typeL(arrayB, numL, inputString, stringLength,
                           suffixType, arrayA);
        free((void*)arrayB);
        return (arrayA);
    }
}
void outstring(const char* inputString, const int stringLength,
               const char* strSpace)
{
    for (int i = 0; i < stringLength; i++)
    {
        cout << inputString[i] << strSpace;
    }

    cout << endl;
}

void outstring(const int* inputString, const int stringLength)
{
    cout << inputString[0];

    for (int i = 1; i < stringLength; i++)
    {
        cout << setw(3);
        cout << inputString[i];
    }

    cout << endl;
}

void outstring3(const int* inputString, const int stringLength)
{
    cout << inputString[0];

    for (int i = 1; i < stringLength; i++)
    {
        cout << setw(3);
        cout << inputString[i];
    }

    cout << endl;
}

void outstringS(const int* inputString, const int stringLength)
{
    for (int i = 0; i < stringLength; i++)
    {
        cout << inputString[i] << ' ';

        if (i % 10 == 9)
        {
            cout << endl;
        }

    }
    cout << endl;
}

void outsuffix(const int* intString, const int intStringLen,
               const char* charString, const int charStringLen,
               const int suffixLength)
{
    int j;

    for (int i = 0; i < intStringLen; i++)
    {
        cout << setw(10);
        cout << intString[i];
        cout << ' ';
        j = intString[i];

        if (j > 1) //output Left(j)
            cout << '[' << charString[j - 1] << ']';

        while ((j < suffixLength + intString[i]) && (j < charStringLen))
        {
            cout << charString[j];
            j++;
        }

        cout << endl;
    }

    cout << endl;
}

void outsubstring(const int* intString, const int intStringLen,
                  const char* charString, const int charStringLen,
                  const BooleanString& suffixType)
{
    int j;
    bool continueloop;

    for (int i = 0; i < intStringLen; i++)
    {
        cout << setw(6);
        cout << intString[i];
        cout << ' ';
        j = intString[i];

        if (suffixType.getValN(j) == 0)
        {
            cout << charString[j];
            j++;
            continueloop = true;

            while (j < charStringLen && continueloop)
            {
                if (suffixType.getValN(j) == 1)
                {
                    cout << charString[j];
                    j++;
                }
                else
                {
                    cout << charString[j];
                    continueloop = false;
                }

            }

        }
        else
        {
            cout << charString[j];
            j++;
            continueloop = true;

            while (j < charStringLen && continueloop)
            {
                if (suffixType.getValN(j) == 0)
                {
                    cout << charString[j];
                    j++;
                }
                else
                {
                    cout << charString[j];
                    continueloop = false;
                }

            }

        }
        cout << endl;
    }

    cout << endl;
}
void s_distanceR(const BooleanString& suffixType,
                 const int inputLength, int& maxDist)
{
    int i, temp, prevDist;

    //Set the reverse S-distance for the last suffix to 0.

    maxDist = 0;

    prevDist = 1;

    for (i = inputLength - 2; i >= 0; i--)
    {
        //Calculate maxDist

        temp = maxDist - prevDist;
        temp = temp >> MAX_INT_INDEX;
        maxDist = maxDist + ((prevDist - maxDist) & temp);

        //If suffix i is type L, then increment prevDist by 1.
        //else set prevDist to 1.

        temp = suffixType.getValN(i) - 1;
        prevDist = (prevDist & temp) + 1;
    }
}


//This function computes the reverse l_distance.

void l_distanceR(const BooleanString& suffixType,
                 const int inputLength, int& maxDist)
{
    int i, temp, prevDist;

    //Set the reverse l-distance for the last suffix to 0;

    maxDist = 0;

    prevDist = 1;

    for (i = inputLength - 2; i >= 0; i--)
    {
        //Calculate maxDist.

        temp = maxDist - prevDist;
        temp = temp >> MAX_INT_INDEX;
        maxDist = maxDist + ((prevDist - maxDist) & temp);

        //If suffix i is type S, then increment prevDist by 1.
        //else set prevDist to 1.

        temp = 0 - suffixType.getValN(i);
        prevDist = (prevDist & temp) + 1;
    }
}



void sort_s_substringC(const char* inputString, const int inputLength,
                       int* A, BooleanString& BuckA, const int ALength,
                       const BooleanString& suffixType, const int maxDist)
{
    int i, prevCount, temp, tmpIndex, offset;
    int start, end, prevPos, bufferLen;
    int* buffer;
    int* skipVal;
    int* tempBucket;

    //set BuckA to all false.
    //invALength is the negative value of ALength.
    //and other initializations.

    bufferLen = MAX_NUM_CHAR * 2;
    BuckA.setAll(false);
    buffer = (int*) malloc (bufferLen * sizeof(int));
    skipVal = (int*) malloc (ALength * sizeof(int));
    memset ((void*)skipVal, 0, sizeof(int) * ALength);
    tempBucket = (int*) malloc (ALength * sizeof(int));

    //Set skipVal[0] to ALength, i.e. the entire array is an
    //unsorted bucket.

    skipVal[0] = ALength;

    //while it is not all sorted. continue to sort.

    for (offset = 0; offset <= maxDist; offset++)
    {
        start = 0;
        prevPos = 0;
        //continue sort until reach the end of array A.

        while (start < ALength)
        {
            //identify the first bucket need work done.
            //also go back to previous bucket to increase
            //its skipVal.
            //start is the index where the bucket starts.
            //end is the index where the bucket ends, i.e.
            //the index of the last element in the array.

            prevPos = start;

            while (skipVal[start] < 0 && start < ALength)
            {
                start = -skipVal[start];
            }

            end = skipVal[start] - 1;
            skipVal[prevPos] = -start;
            memset((void*)buffer, 0, sizeof(int) * bufferLen);

            //Sort the bucket if start <= ALength - 1.

            if (start < ALength)
            {
                //count the occurence, also copy the elements
                //into a temporary bucket.
                //If the next character is type S. Then add 1 to 2*c+1.
                //Otherwise add 1 to 2*c.

                for (i = start; i <= end; i++)
                {
                    //copy the element into temporary bucket.
                    tempBucket[i] = A[i];

                    //count the occurence.
                    tmpIndex = A[i] + offset;
                    temp = ((int) inputString[tmpIndex]) << 1;
                    temp = temp + suffixType.getValN(tmpIndex);
                    buffer[temp]++;
                }

                //calculate the new starting value of each bucket.
                prevCount = buffer[0];

                buffer[0] = start;

                for (i = 1; i < bufferLen; i++)
                {
                    temp = buffer[i];
                    buffer[i] = buffer[i - 1] + prevCount;
                    prevCount = temp;
                }

                //put the elements in there correct positions.
                //must read from the temporary bucket, because array
                //A is not reliable.

                for (i = start; i <= end; i++)
                {
                    tmpIndex = tempBucket[i] + offset;
                    temp = ((int) inputString[tmpIndex]) << 1;
                    temp = temp + suffixType.getValN(tmpIndex);
                    A[buffer[temp]] = tempBucket[i];
                    buffer[temp]++;
                }

                //set the bucket boundary of array A, and the correct
                //skipVal. skipVal[i] is negative if no more work required for
                //the bucket anymore (all the buckets whose last character is
                //type S).
                //An empty bucket will not affect a previously non-empty
                //bucket. Because suppose the empty bucket is i, then
                //buffer[i] actually points to the beginning of the
                //next bucket. buffer[i-1] points to either the beginning
                //of the previously non-empty bucket, or it points to
                //the beginning of the next bucket. In either case, it
                //is set correctly.

                i = 1;

                if (offset > 0)
                {
                    //if the first bucket has a value greater than start then
                    //it is not empty.

                    if (buffer[1] > start)
                    {
                        BuckA.setValN(buffer[1] - 1, 1);
                        skipVal[start] = - buffer[0];
                    }

                    //for all other buckets if its not empty, its value is greater
                    //than the previous bucket.
                    for (i = 1; i < bufferLen; i++)
                    {
                        //if the bucket only have 1 element, then skipVal is negative

                        if (buffer[i] == buffer[i - 1] + 1)
                        {
                            BuckA.setValN(buffer[i] - 1, 1);
                            skipVal[buffer[i - 1]] = - buffer[i];
                        }
                        else if (buffer[i] > buffer[i - 1] + 1)
                        {
                            BuckA.setValN(buffer[i] - 1, 1);
                            temp = - (i & 1);
                            temp = (buffer[i] ^ temp) - temp;
                            skipVal[buffer[i - 1]] = temp;
                        }

                    }

                }
                else
                {
                    //if the first bucket has a value greater than start then
                    //it is not empty.
                    //Note that the first bucket is 1 because '$' is type S,
                    //and it's count is in 1 not 0.

                    if (buffer[1] > start)
                    {
                        BuckA.setValN(buffer[1] - 1, 1);
                        skipVal[start] = - buffer[0];
                    }

                    //for all other buckets if its not empty, its value is greater
                    //than the previous bucket.
                    for (i = 1; i < bufferLen; i++)
                    {
                        //if the bucket only have 1 element, then skipVal is negative

                        if (buffer[i] == buffer[i - 1] + 1)
                        {
                            BuckA.setValN(buffer[i] - 1, 1);
                            skipVal[buffer[i - 1]] = - buffer[i];
                        }
                        else if (buffer[i] > buffer[i - 1] + 1)
                        {
                            BuckA.setValN(buffer[i] - 1, 1);
                            skipVal[buffer[i - 1]] = buffer[i];
                        }

                    }

                }

                //the above while loop will not be executed if the
                //bucket contains only 1 element. In that case we set
                //the values in the if statement below.

                if (start == end)
                {
                    skipVal[start] = -(end - 1);
                    BuckA.setValN(start, 1);
                }

                //set the start to the point to the next bucket.

                start = end + 1;
            }

        }

    }

    free(skipVal);
    free(tempBucket);
    free(buffer);

}



void sort_l_substringC(const char* inputString, const int inputLength,
                       int* A, BooleanString& BuckA, const int ALength,
                       const BooleanString& suffixType, const int maxDist)
{
    int i, prevCount, temp, tmpIndex, offset;
    int start, end, prevPos, bufferLen;
    int* buffer;
    int* skipVal;
    int* tempBucket;

    //set BuckA to all false.
    //invALength is the negative value of ALength.
    //and other initializations.

    bufferLen = MAX_NUM_CHAR * 2;
    BuckA.setAll(false);
    buffer = (int*) malloc (bufferLen * sizeof(int));
    skipVal = (int*) malloc (ALength * sizeof(int));
    memset ((void*)skipVal, 0, sizeof(int) * ALength);
    tempBucket = (int*) malloc (ALength * sizeof(int));

    //set value for A[0] which is the first position in
    //array A. Which should be occupied by '$'.

    skipVal[0] = ALength;

    //while it is not all sorted. continue to sort.

    for (offset = 0; offset <= maxDist; offset++)
    {
        start = 0;
        prevPos = 0;
        //continue sort until reach the end of array A.

        while (start < ALength)
        {
            //identify the first bucket need work done.
            //also go back to previous bucket to increase
            //its skipVal.
            //start is the index where the bucket starts.
            //end is the index where the bucket ends, i.e.
            //the index of the last element in the array.

            prevPos = start;

            while (skipVal[start] < 0 && start < ALength)
            {
                start = -skipVal[start];
            }

            end = skipVal[start] - 1;
            skipVal[prevPos] = -start;
            memset((void*)buffer, 0, sizeof(int) * bufferLen);

            //Sort the bucket if start <= ALength - 1.

            if (start < ALength)
            {
                //count the occurence, also copy the elements
                //into a temporary bucket.
                //If the next character is type S. Then add 1 to 2*c+1.
                //Otherwise add 1 to 2*c.

                for (i = start; i <= end; i++)
                {
                    //copy the element into temporary bucket.
                    tempBucket[i] = A[i];

                    //count the occurence.
                    tmpIndex = A[i] + offset;
                    temp = ((int) inputString[tmpIndex]) << 1;
                    temp = temp + suffixType.getValN(tmpIndex);
                    buffer[temp]++;
                }

                //calculate the new starting value of each bucket.
                prevCount = buffer[0];

                buffer[0] = start;

                for (i = 1; i < bufferLen; i++)
                {
                    temp = buffer[i];
                    buffer[i] = buffer[i - 1] + prevCount;
                    prevCount = temp;
                }

                //put the elements in there correct positions.
                //must read from the temporary bucket, because array
                //A is not reliable.

                for (i = start; i <= end; i++)
                {
                    tmpIndex = tempBucket[i] + offset;
                    temp = ((int) inputString[tmpIndex]) << 1;
                    temp = temp + suffixType.getValN(tmpIndex);
                    A[buffer[temp]] = tempBucket[i];
                    buffer[temp]++;
                }

                //set the bucket boundary of array A, and the correct
                //skipVal. skipVal[i] is negative if no more work required for
                //the bucket anymore (all the buckets whose last character is
                //type L).
                //An empty bucket will not affect a previously non-empty
                //bucket. Because suppose the empty bucket is i, then
                //buffer[i] actually points to the beginning of the
                //next bucket. buffer[i-1] points to either the beginning
                //of the previously non-empty bucket, or it points to
                //the beginning of the next bucket. In either case, it
                //is set correctly.

                if (offset > 0)
                {
                    //if the first bucket has a value greater than start then
                    //it is not empty.

                    if (buffer[0] > start)
                    {
                        BuckA.setValN(buffer[0] - 1, 1);
                        skipVal[start] = - buffer[0];
                    }

                    //for all other buckets if its not empty, its value is greater
                    //than the previous bucket.
                    for (i = 1; i < bufferLen; i++)
                    {
                        //if the bucket only have 1 element, then skipVal is negative

                        if (buffer[i] == buffer[i - 1] + 1)
                        {
                            BuckA.setValN(buffer[i] - 1, 1);
                            skipVal[buffer[i - 1]] = - buffer[i];
                        }
                        else if (buffer[i] > buffer[i - 1] + 1)
                        {
                            BuckA.setValN(buffer[i] - 1, 1);
                            temp = - ((i & 1) ^ 1);
                            temp = (buffer[i] ^ temp) - temp;
                            skipVal[buffer[i - 1]] = temp;
                        }

                    }

                }
                else
                {
                    //if the first bucket has a value greater than start then
                    //it is not empty.

                    if (buffer[0] > start)
                    {
                        BuckA.setValN(buffer[0] - 1, 1);
                        skipVal[start] = - buffer[0];
                    }

                    //for all other buckets if its not empty, its value is greater
                    //than the previous bucket.
                    for (i = 1; i < bufferLen; i++)
                    {
                        //if the bucket only have 1 element, then skipVal is negative

                        if (buffer[i] == buffer[i - 1] + 1)
                        {
                            BuckA.setValN(buffer[i] - 1, 1);
                            skipVal[buffer[i - 1]] = - buffer[i];
                        }
                        else if (buffer[i] > buffer[i - 1] + 1)
                        {
                            BuckA.setValN(buffer[i] - 1, 1);
                            skipVal[buffer[i - 1]] = buffer[i];
                        }

                    }

                }

                //set the start to the point to the next bucket.

                start = end + 1;
            }

        }

    }
    free(skipVal);
    free(tempBucket);
    free(buffer);
}
void construct_TPrime_typeS(int* ArrayB, const int ArrayBLength,
                            const BooleanString& BuckB,
                            int* tPrime, const int inputLength,
                            const BooleanString& suffixType)
{
    int* Buckets;
    int i, j, currBuck;
    int tempVal, tempValInv;

    // initalize array Buckets, and for each type
    // S suffix calculate their bucket number

    Buckets = (int*) malloc(inputLength * sizeof(int));

    currBuck = 0;

    for (i = 0; i < ArrayBLength; i++)
    {
        Buckets[ArrayB[i]] = currBuck;
        currBuck = currBuck + BuckB.getValN(i);
    }

    // construct tPrime
    j = 0;

    for (i = 0; i < inputLength; i++)
    {
        tempVal = - suffixType.getValN(i);
        tempValInv = ~tempVal;
        tPrime[j] = (tPrime[j] & tempValInv) | (Buckets[i] & tempVal);
        j = j + (1 & tempVal);
    }

    free(Buckets);
}



void construct_TPrime_typeL(int* ArrayB, const int ArrayBLength,
                            const BooleanString& BuckB,
                            int* tPrime, const int inputLength,
                            const BooleanString& suffixType)
{
    int* Buckets;
    int i, j, currBuck;
    int tempVal, tempValInv;

    // initalize array Buckets, and for each type
    // L suffix calculate their bucket number

    Buckets = (int*) malloc(inputLength * sizeof(int));

    currBuck = 0;

    for (i = 0; i < ArrayBLength; i++)
    {
        Buckets[ArrayB[i]] = currBuck;
        currBuck = currBuck + BuckB.getValN(i);
    }

    // construct tPrime
    j = 0;

    for (i = 0; i < inputLength; i++)
    {
        tempVal = (suffixType.getValN(i) - 1);
        tempValInv = ~tempVal;
        tPrime[j] = (tPrime[j] & tempValInv) | (Buckets[i] & tempVal);
        j = j + (1 & tempVal);
    }

    free(Buckets);
}




void construct_ArrayB_typeS(const int* ArrayA, const int inputLength,
                            const BooleanString& BuckA,
                            int* ArrayB, BooleanString& BuckB,
                            const BooleanString& suffixType)
{
    int i, j;
    int temp;

    BuckB.setAll(false);

    j = 0;

    for (i = 0; i < inputLength; i++)
    {
        temp = ArrayA[i];

        if (suffixType.getValN(temp) == 1)
        {
            ArrayB[j] = ArrayA[i];
            j++;
        }
        if (BuckA.getValN(i) == 1 && j - 1 >= 0)
        {
            BuckB.setValN(j - 1, 1);
        }

    }
}

void construct_ArrayB_typeL(const int* ArrayA, const int inputLength,
                            const BooleanString& BuckA,
                            int* ArrayB, BooleanString& BuckB,
                            const BooleanString& suffixType)
{
    int i, j;
    int temp;

    BuckB.setAll(false);

    j = 0;

    for (i = 0; i < inputLength; i++)
    {
        temp = ArrayA[i];

        if (suffixType.getValN(temp) == 0)
        {
            ArrayB[j] = ArrayA[i];
            j++;
        }
        if (BuckA.getValN(i) == 1 && j - 1 >= 0)
        {
            BuckB.setValN(j - 1, 1);
        }

    }
}

void construct_SA_typeS(int* ArrayB, const int ArrayBLength,
                        const int* stringT, const int inputLength,
                        const BooleanString& suffixType,
                        int* suffixArray)
{
    int* count;
    int sigma, maxChar, minChar;
    int i, j, temp, prevCount, prevChar, charBuck;

    // count the size of the alphabet in stringT.

    findMaxMin(stringT, inputLength, maxChar, minChar);
    sigma = maxChar - minChar + 1;

    // for each character in the alphabet calculate
    // the beginning position of their bucket.

    count = (int*) malloc(sigma * sizeof(int));
    memset((void*)count, 0, sigma * sizeof(int));

    for (i = 0; i < inputLength; i++)
    {
        temp = stringT[i] - minChar;
        count[temp]++;
    }

    prevCount = count[0];
    count[0] = 0;

    for (i = 1; i < sigma; i++)
    {
        temp = count[i];
        count[i] = count[i - 1] + prevCount;
        prevCount = temp;
    }

    // initalize the suffix array

    memset((void*)suffixArray, -1, inputLength * sizeof(int));

    // initalize j, which means the at the beginning of array B

    j = 0;

    // move all the suffixes into the suffix array in
    // order

    for (i = 0; i < inputLength; i++)
    {
        // if suffixArray[i] is -1, then an element of
        // array B should be put in that place. And its
        // previous suffix moved to the approperiate
        // place.
        // Otherwise just move the previous suffix to
        // the approperiate place.

        if (suffixArray[i] == -1)
        {
            assert(j < ArrayBLength);
            suffixArray[i] = ArrayB[j];
            j++;
            prevChar = suffixArray[i] - 1;

            if (prevChar >= 0)
            {
                if (suffixType.getValN(prevChar) == 0)
                {
                    charBuck = stringT[prevChar];

                    if (count[charBuck] > i)
                    {
                        suffixArray[count[charBuck]] = prevChar;
                        count[charBuck]++;
                    }

                }

            }

        }
        else
        {
            prevChar = suffixArray[i] - 1;

            if (prevChar >= 0)
            {
                if (suffixType.getValN(prevChar) == 0)
                {
                    charBuck = stringT[prevChar];

                    if (count[charBuck] > i)
                    {
                        suffixArray[count[charBuck]] = prevChar;
                        count[charBuck]++;
                    }

                }

            }

        }

    }
    free(count);
}





void construct_SA_typeL(int* ArrayB, const int ArrayBLength,
                        const int* stringT, const int inputLength,
                        const BooleanString& suffixType,
                        int* suffixArray)
{
    int* count;
    int sigma, maxChar, minChar, prevChar, charBuck;
    int i, j, temp;

    // count the size of the alphabet in stringT.

    findMaxMin(stringT, inputLength, maxChar, minChar);
    sigma = maxChar - minChar + 1;

    // for each character in the alphabet calculate
    // the ending position of their bucket.

    count = (int*) malloc(sigma * sizeof(int));
    memset((void*)count, 0, sigma * sizeof(int));

    for (i = 0; i < inputLength; i++)
    {
        temp = stringT[i] - minChar;
        count[temp]++;
    }

    count[0] = count[0] - 1;

    for (i = 1; i < sigma; i++)
    {
        count[i] = count[i - 1] + count[i];
    }

    // initalize the suffix array

    memset((void*)suffixArray, -1, inputLength * sizeof(int));

    // initalize j to be at the end of array B.

    j = ArrayBLength - 1;

    // move all the suffixes into the suffix array in
    // order

    for (i = inputLength - 1; i >= 0; i--)
    {
        // if suffixArray[i] is -1, then an element of
        // array B should be put in that place. And its
        // previous suffix moved to the approperiate
        // place.
        // Otherwise just move the previous suffix to
        // the approperiate place.

        if (suffixArray[i] == -1)
        {
            assert(j < ArrayBLength);
            suffixArray[i] = ArrayB[j];
            j--;
            prevChar = suffixArray[i] - 1;

            if (prevChar >= 0)
            {
                if (suffixType.getValN(prevChar) == 1)
                {
                    charBuck = stringT[prevChar];

                    if (count[charBuck] < i)
                    {
                        suffixArray[count[charBuck]] = prevChar;
                        count[charBuck]--;
                    }

                }

            }

        }
        else
        {
            prevChar = suffixArray[i] - 1;

            if (prevChar >= 0)
            {
                if (suffixType.getValN(prevChar) == 1)
                {
                    charBuck = stringT[prevChar];

                    if (count[charBuck] < i)
                    {
                        suffixArray[count[charBuck]] = prevChar;
                        count[charBuck]--;
                    }

                }

            }

        }

    }
    free(count);
}
void construct_SA_typeS(int* ArrayB, const int ArrayBLength,
                        const char* stringT, const int inputLength,
                        const BooleanString& suffixType,
                        int* suffixArray)
{
    int* count;
    int i, j, temp, prevCount, prevChar, charBuck;

    // for each character in the alphabet calculate
    // the beginning position of their bucket.

    count = (int*) malloc(MAX_NUM_CHAR * sizeof(int));
    memset((void*)count, 0, MAX_NUM_CHAR * sizeof(int));

    for (i = 0; i < inputLength; i++)
    {
        temp = (int) stringT[i];
        count[temp]++;
    }

    prevCount = count[0];
    count[0] = 0;

    for (i = 1; i < MAX_NUM_CHAR; i++)
    {
        temp = count[i];
        count[i] = count[i - 1] + prevCount;
        prevCount = temp;
    }

    // initalize the suffix array

    memset((void*)suffixArray, -1, inputLength * sizeof(int));

    // initalize j, which means the at the beginning of array B

    j = 0;

    // move all the suffixes into the suffix array in
    // order

    for (i = 0; i < inputLength; i++)
    {
        // if suffixArray[i] is -1, then an element of
        // array B should be put in that place. And its
        // previous suffix moved to the approperiate
        // place.
        // Otherwise just move the previous suffix to
        // the approperiate place.

        if (suffixArray[i] == -1)
        {
            assert(j < ArrayBLength);
            suffixArray[i] = ArrayB[j];
            j++;
            prevChar = suffixArray[i] - 1;

            if (prevChar >= 0)
            {
                if (suffixType.getValN(prevChar) == 0)
                {
                    charBuck = (int) stringT[prevChar];

                    if (count[charBuck] > i)
                    {
                        suffixArray[count[charBuck]] = prevChar;
                        count[charBuck]++;
                    }

                }

            }

        }
        else
        {
            prevChar = suffixArray[i] - 1;

            if (prevChar >= 0)
            {
                if (suffixType.getValN(prevChar) == 0)
                {
                    charBuck = (int) stringT[prevChar];

                    if (count[charBuck] > i)
                    {
                        suffixArray[count[charBuck]] = prevChar;
                        count[charBuck]++;
                    }

                }

            }

        }

    }
    free(count);
}





void construct_SA_typeL(int* ArrayB, const int ArrayBLength,
                        const char* stringT, const int inputLength,
                        const BooleanString& suffixType,
                        int* suffixArray)
{
    int* count;
    int prevChar, charBuck;
    int i, j, temp;

    // for each character in the alphabet calculate
    // the ending position of their bucket.

    count = (int*) malloc(MAX_NUM_CHAR * sizeof(int));
    memset((void*)count, 0, MAX_NUM_CHAR * sizeof(int));

    // count the occurrence.

    for (i = 0; i < inputLength; i++)
    {
        temp = (int) stringT[i];
        count[temp]++;
    }

    count[0] = count[0] - 1;

    for (i = 1; i < MAX_NUM_CHAR; i++)
    {
        count[i] = count[i - 1] + count[i];
    }

    // initalize the suffix array

    memset((void*)suffixArray, -1, inputLength * sizeof(int));

    // initalize j to be at the end of array B.

    j = ArrayBLength - 1;

    // move all the suffixes into the suffix array in
    // order

    for (i = inputLength - 1; i >= 0; i--)
    {
        // if suffixArray[i] is -1, then an element of
        // array B should be put in that place. And its
        // previous suffix moved to the approperiate
        // place.
        // Otherwise just move the previous suffix to
        // the approperiate place.

        if (suffixArray[i] == -1)
        {
            assert(j < ArrayBLength);
            suffixArray[i] = ArrayB[j];
            j--;
            prevChar = suffixArray[i] - 1;

            if (prevChar >= 0)
            {
                if (suffixType.getValN(prevChar) == 1)
                {
                    charBuck = (int) stringT[prevChar];

                    if (count[charBuck] < i)
                    {
                        suffixArray[count[charBuck]] = prevChar;
                        count[charBuck]--;
                    }

                }

            }

        }
        else
        {
            prevChar = suffixArray[i] - 1;

            if (prevChar >= 0)
            {
                if (suffixType.getValN(prevChar) == 1)
                {
                    charBuck = (int) stringT[prevChar];

                    if (count[charBuck] < i)
                    {
                        suffixArray[count[charBuck]] = prevChar;
                        count[charBuck]--;
                    }

                }

            }

        }

    }
    free(count);
}
int* LinearSuffixSort(const int* inputString, const int stringLength)
{
    BooleanString suffixType(stringLength);
    BooleanString BuckA(stringLength);
    BooleanString BuckB, BuckList;
    int numS, numL, maxDist;
    int* ArrayA;
    int* ArrayB;
    int* Dist;
    int* DistCount;
    int* theList;
    int* tPrime;
    int* int_buffer;
    int listLength;

    if (stringLength <= 0)
    {
        cout << "In function LinearSuffixSort(int*, int):" << endl;
        cout << "Length of input string cannot be less than 0."
        << endl;
        exit(0);
    }
    //cout << endl;

    suffix_type(inputString, stringLength, numS, numL, suffixType);

    if (suffixType.getValN(stringLength - 1) == 1 &&
            numS == 1)
    {
        ArrayB = (int*) malloc(sizeof(int));
        ArrayB[0] = stringLength - 1;
        ArrayA = (int*) malloc(stringLength * sizeof(int));
        construct_SA_typeS(ArrayB, numS, inputString, stringLength,
                           suffixType, ArrayA);
        free(ArrayB);
        return (ArrayA);
    }

    ArrayA = (int*) malloc(stringLength * sizeof(int));
    int_buffer = (int*) malloc(stringLength * sizeof(int));
    counting_sort(inputString, stringLength, ArrayA, BuckA, int_buffer);
    free(int_buffer); //add by Zhao XU

    if (suffixType.getValN(stringLength - 1) == 1)
    {
        //cout << "Sorting type S suffixes" << endl;
        //cout << numS << " " << numL << " " << stringLength << endl;
        ArrayB = (int*) malloc (numS * sizeof(int));
        BuckB.initialize(numS);
        construct_ArrayB_typeS(ArrayA, stringLength, BuckA, ArrayB,
                               BuckB, suffixType);

        Dist = (int*) malloc(stringLength * sizeof(int));

        s_distance(suffixType, stringLength, Dist,
                   DistCount, maxDist);

        BuckList.initialize(DistCount[maxDist]);
        listLength = DistCount[maxDist];

        theList = construct_list_typeS(ArrayA, stringLength, Dist,
                                       DistCount, maxDist, BuckA, BuckList,
                                       listLength);

        free(DistCount);
        free(Dist);

        sort_by_list_typeS(ArrayB, BuckB, theList, BuckList, stringLength,
                           listLength, numS);

        free(theList);

        if (BuckB.isAllTrue())
        {
            ArrayA = (int*) malloc(stringLength * sizeof(int));
            construct_SA_typeS(ArrayB, numS, inputString, stringLength,
                               suffixType, ArrayA);
            free(ArrayB);
            return (ArrayA);
        }

        tPrime = (int*) malloc(numS * sizeof(int));

        construct_TPrime_typeS(ArrayB, numS, BuckB, tPrime,
                               stringLength, suffixType);
        free(ArrayB);

        ArrayB = LinearSuffixSort(tPrime, numS);
        free(tPrime);

        reconstruct_B_typeS(ArrayB, numS, suffixType, stringLength);
        ArrayA = (int*) malloc(stringLength * sizeof(int));
        construct_SA_typeS(ArrayB, numS, inputString, stringLength,
                           suffixType, ArrayA);
        free(ArrayB);
        return (ArrayA);
    }
    else
    {
        //cout << "Sorting type L suffixes" << endl;
        //cout << numS << " " << numL << " " << stringLength << endl;
        ArrayB = (int*) malloc(numL * sizeof(int));
        BuckB.initialize(numL);

        construct_ArrayB_typeL(ArrayA, stringLength, BuckA, ArrayB,
                               BuckB, suffixType);

        Dist = (int*) malloc(stringLength * sizeof(int));

        l_distance(suffixType, stringLength, Dist,
                   DistCount, maxDist);

        BuckList.initialize(DistCount[maxDist]);
        listLength = DistCount[maxDist];

        theList = construct_list_typeL(ArrayA, stringLength, Dist,
                                       DistCount, maxDist, BuckA, BuckList,
                                       listLength);
        free(DistCount);
        free(Dist);

        sort_by_list_typeL(ArrayB, BuckB, theList, BuckList, stringLength,
                           listLength, numL);
        free(theList);

        if (BuckB.isAllTrue())
        {
            ArrayA = (int*) malloc(stringLength * sizeof(int));
            construct_SA_typeL(ArrayB, numL, inputString, stringLength,
                               suffixType, ArrayA);
            free(ArrayB);
            return (ArrayA);
        }

        tPrime = (int*) malloc(numL * sizeof(int));

        construct_TPrime_typeL(ArrayB, numL, BuckB, tPrime,
                               stringLength, suffixType);
        free(ArrayB);

        ArrayB = LinearSuffixSort(tPrime, numL);
        free(tPrime);

        reconstruct_B_typeL(ArrayB, numL, suffixType, stringLength);
        ArrayA = (int*) malloc(stringLength * sizeof(int));
        construct_SA_typeL(ArrayB, numL, inputString, stringLength,
                           suffixType, ArrayA);
        free(ArrayB);
        return (ArrayA);
    }
}



void reconstruct_B_typeS(int* ArrayB, const int ArrayBLength,
                         const BooleanString& suffixType,
                         const int inputLength)
{
    int* convertion;
    int i, j, tempVal, tempValInv;

    // build a conversion array, such that conversion[i] = j
    // where i is the index of a suffix of T', and j is
    // the index of the corresponding type S suffix in T

    convertion = (int*) malloc(ArrayBLength * sizeof(int));

    j = 0;

    for (i = 0; i < inputLength; i++)
    {
        tempVal = suffixType.getValN(i);
        tempVal = tempVal << MAX_INT_INDEX;
        tempVal = tempVal >> MAX_INT_INDEX;
        tempValInv = ~tempVal;

        convertion[j] = (i & tempVal) | (convertion[j] & tempValInv);
        j = j + (1 & tempVal);
    }

    // use the conversion array to calculate the actual
    // index of the type S suffix of ArrayB[i]

    for (i = 0; i < ArrayBLength; i++)
    {
        tempVal = ArrayB[i];
        ArrayB[i] = convertion[tempVal];
    }

    free(convertion);
}



void reconstruct_B_typeL(int* ArrayB, const int ArrayBLength,
                         const BooleanString& suffixType,
                         const int inputLength)
{
    int* convertion;
    int i, j, tempVal, tempValInv;

    // build a conversion array, such that conversion[i] = j
    // where i is the index of a suffix of T', and j is
    // the index of the corresponding type L suffix in T

    convertion = (int*) malloc(ArrayBLength * sizeof(int));

    j = 0;

    for (i = 0; i < inputLength; i++)
    {
        tempVal = suffixType.getValN(i);
        tempVal = tempVal << MAX_INT_INDEX;
        tempVal = tempVal >> MAX_INT_INDEX;
        tempValInv = ~tempVal;

        convertion[j] = (i & tempValInv) | (convertion[j] & tempVal);
        j = j + (1 & tempValInv);
    }

    // use the conversion array to calculate the actual
    // index of the type S suffix of ArrayB[i]

    for (i = 0; i < ArrayBLength; i++)
    {
        tempVal = ArrayB[i];
        ArrayB[i] = convertion[tempVal];
    }

    free(convertion);
}
int* construct_list_typeS(int*& ArrayA, const int inputLength, int* Dist,
                          int* DistCount, const int maxDist,
                          const BooleanString& BuckA,
                          BooleanString& BuckList, const int listLength)
{
    int i, j;
    int posList;
    int temp;
    int startB, endB;

    BuckList.setAll(false);
    i = 0;

    while (i < inputLength)
    {
        startB = i;

        while (BuckA.getValN(i) != 1 && i < inputLength)
        {
            //use Dist to be the temporary reverse mapping array of
            //the List.
            temp = Dist[ArrayA[i]];

            if (temp > 0)
            {
                posList = DistCount[temp - 1];
                Dist[ArrayA[i]] = posList;
                BuckList.setValN(posList, 1);
                DistCount[temp - 1]++;
            }
            else
            {
                Dist[ArrayA[i]] = -1;
            }

            i++;
        }

        //set the values for the last member in the bucket.

        temp = Dist[ArrayA[i]];

        if (temp != 0)
        {
            posList = DistCount[temp - 1];
            Dist[ArrayA[i]] = posList;
            BuckList.setValN(posList, 1);
            DistCount[temp - 1]++;
        }
        else
        {
            Dist[ArrayA[i]] = -1;
        }

        endB = i;

        //for all the member in the middle of the bucket in list
        //reset its BuckList to false. i.e. show that there is
        //no boundary between i and i+1 in list. It is possible
        //that the boundary of the last bucket of a S-distance
        //is not drawn properly.

        for (j = startB; j < endB; j++)
        {
            posList = Dist[ArrayA[j]];

            if (posList >= 0 && posList != listLength - 1)
            {
                if (BuckList.getValN(posList + 1) == 1)
                {
                    BuckList.setValN(posList, 0);
                }

            }

        }
        i++;
    }
    //reconstruct the list from Dist

    for (i = 0; i < inputLength; i++)
    {
        if (Dist[i] >= 0)
        {
            ArrayA[Dist[i]] = i;
        }

    }

    //At this point DistCount[i] give the total number
    //of suffixes having S-distance of i+1. i.e.
    //DistCount[0] give the number of suffixes with
    //S-distance of 1.

    for (i = 0; i < maxDist; i++)
    {
        BuckList.setValN(DistCount[i] - 1, 1);
    }

    //Calculate the index of the beginning of the type S
    //substring.

    for (i = 0; i < maxDist; i++)
    {
        if (i == 0)
        {
            j = 0;
        }
        else
        {
            j = DistCount[i - 1];
        }
        while (j < DistCount[i])
        {
            ArrayA[j] = ArrayA[j] - i - 1;
            j++;
        }

    }
    return (ArrayA);
}




int* construct_list_typeL(int*& ArrayA, const int inputLength, int* Dist,
                          int* DistCount, const int maxDist,
                          const BooleanString& BuckA,
                          BooleanString& BuckList, const int listLength)
{
    int i, j;
    int posList;
    int temp;
    int startB, endB;
    bool firstElement;

    BuckList.setAll(false);
    i = inputLength - 1;

    while (i >= 0)
    {
        endB = i;

        if (i > 0)
        {
            if (BuckA.getValN(i - 1) != 1)
            {
                firstElement = false;
            }
            else
            {
                firstElement = true;
            }

        }
        else
        {
            firstElement = true;
        }

        while (!firstElement)
        {

            //use Dist to be the temporary reverse mapping array of
            //the List.

            temp = Dist[ArrayA[i]];

            if (temp > 0)
            {
                posList = DistCount[temp - 1];
                Dist[ArrayA[i]] = posList;
                BuckList.setValN(posList, 1);
                DistCount[temp - 1]++;
            }
            else
            {
                Dist[ArrayA[i]] = -1;
            }

            // Decrement i, and determine if i is the first element
            // of the bucket.

            i--;

            if (i > 0)
            {
                if (BuckA.getValN(i - 1) != 1)
                {
                    firstElement = false;
                }
                else
                {
                    firstElement = true;
                }

            }
            else
            {
                firstElement = true;
            }

        }

        //set the values for the first member in the bucket.

        temp = Dist[ArrayA[i]];

        if (temp != 0)
        {
            posList = DistCount[temp - 1];
            Dist[ArrayA[i]] = posList;
            BuckList.setValN(posList, 1);
            DistCount[temp - 1]++;
        }
        else
        {
            Dist[ArrayA[i]] = -1;
        }

        startB = i;

        //for all the member in the middle of the bucket in list
        //reset its BuckList to false. i.e. show that there is
        //no boundary between i and i+1 in list. It is possible
        //that the boundary of the last bucket of a L-distance
        //is not drawn properly.

        for (j = endB; j >= startB; j--)
        {
            posList = Dist[ArrayA[j]];

            if (posList >= 0 && posList != listLength - 1)
            {
                if (BuckList.getValN(posList + 1) == 1)
                {
                    BuckList.setValN(posList, 0);
                }

            }

        }
        i--;
    }

    //reconstruct the list from Dist

    for (i = 0; i < inputLength; i++)
    {
        if (Dist[i] >= 0)
        {
            ArrayA[Dist[i]] = i;
        }

    }

    //At this point DistCount[i] give the total number
    //of suffixes having L-distance of i+1. i.e.
    //DistCount[0] give the number of suffixes with
    //L-distance of 1.

    for (i = 0; i < maxDist; i++)
    {
        BuckList.setValN(DistCount[i] - 1, 1);
    }

    //Calculate the index of the beginning of the type L
    //substring.

    for (i = 0; i < maxDist; i++)
    {
        if (i == 0)
        {
            j = 0;
        }
        else
        {
            j = DistCount[i - 1];
        }
        while (j < DistCount[i])
        {
            ArrayA[j] = ArrayA[j] - i - 1;
            j++;
        }

    }
    return (ArrayA);
}




void sort_by_list_typeS(int* ArrayB, BooleanString& BuckB,
                        int* List, const BooleanString& BuckList,
                        const int inputLength, const int listLength,
                        const int ArrayBLength)
{
    int* Rev;
    int* Left;
    int i, j, newBuckNum, BucketNum;
    int BucketRight;

    Rev = (int*) malloc(inputLength * sizeof(int));
    Left = (int*) malloc(ArrayBLength * sizeof(int));

    memset((void*)Rev, -1, inputLength*sizeof(int));
    memset((void*)Left, -1, ArrayBLength*sizeof(int));

    // initalize Rev and Left

    BucketRight = ArrayBLength - 1;

    for (i = ArrayBLength - 1; i > 0; i--)
    {
        Rev[ArrayB[i]] = BucketRight;

        if (BuckB.getValN(i - 1) == 1)
        {
            Left[BucketRight] = i;
            BucketRight = i - 1;
        }

    }

    // initalize Rev and Left for the first element of ArrayB

    Rev[ArrayB[0]] = BucketRight;
    Left[BucketRight] = 0;

    // sort the type S substrings according to the list
    // bucket by bucket.

    i = 0;

    while (i < listLength)
    {

        // count the number of elements to move in each bucket
        // and also set the value of Left for all the elements
        // that being moved.

        j = i;

        while (BuckList.getValN(j) == 0)
        {
            Left[Rev[List[j]]]++;
            j++;
        }

        // count for the last element of the bucket.

        Left[Rev[List[j]]]++;

        // moving the elements by re-assigning the Rev
        // after this we still need to update Left

        j = i;

        while (BuckList.getValN(j) == 0)
        {
            newBuckNum = Left[Rev[List[j]]] - 1;
            Rev[List[j]] = newBuckNum;
            j++;
        }

        // update Rev for the last element of the bucket

        newBuckNum = Left[Rev[List[j]]] - 1;

        Rev[List[j]] = newBuckNum;

        // correct the values of Left for all affected
        // buckets

        j = i;

        while (BuckList.getValN(j) == 0)
        {
            newBuckNum = Rev[List[j]];

            if (Left[newBuckNum] == -1)
            {
                Left[newBuckNum] = newBuckNum;
            }
            else
            {
                Left[newBuckNum]--;
            }

            BuckB.setValN(newBuckNum, 1);
            j++;
        }

        // correct the values of Left for the last element

        newBuckNum = Rev[List[j]];

        if (Left[newBuckNum] == -1)
        {
            Left[newBuckNum] = newBuckNum;
        }
        else
        {
            Left[newBuckNum]--;
        }

        BuckB.setValN(newBuckNum, 1);

        // set i to point to the first element of the next bucket

        i = j + 1;
    }

    // Reconsturct ArrayB from Rev.

    for (i = 0; i < inputLength; i++)
    {
        BucketNum = Rev[i];

        if (BucketNum > -1)
        {
            ArrayB[Left[BucketNum]] = i;
            Left[BucketNum]++;
        }

    }
    free(Rev);
    free(Left);
}



void sort_by_list_typeL(int* ArrayB, BooleanString& BuckB,
                        int* List, const BooleanString& BuckList,
                        const int inputLength, const int listLength,
                        const int ArrayBLength)
{
    int* Rev;
    int* Right;
    int i, j, newBuckNum, BucketNum;
    int BucketLeft;

    Rev = (int*) malloc(inputLength * sizeof(int));
    Right = (int*) malloc(ArrayBLength * sizeof(int));

    memset((void*)Rev, -1, inputLength*sizeof(int));
    memset((void*)Right, -1, ArrayBLength*sizeof(int));

    // initalize Rev and Right

    BucketLeft = 0;

    for (i = 0; i < ArrayBLength; i++)
    {
        Rev[ArrayB[i]] = BucketLeft;

        if (BuckB.getValN(i) == 1)
        {
            Right[BucketLeft] = i;
            BucketLeft = i + 1;
        }

    }

    // sort the type L substrings according to the list
    // bucket by bucket.

    i = 0;

    while (i < listLength)
    {

        // count the number of elements to move in each bucket
        // and also set the value of Right for all the elements
        // that being moved.

        j = i;

        while (BuckList.getValN(j) == 0)
        {
            Right[Rev[List[j]]]--;
            j++;
        }

        // count for the last element of the bucket

        Right[Rev[List[j]]]--;

        // moving the elements by re-assigning the Rev
        // after this we still need to update Right

        j = i;

        while (BuckList.getValN(j) == 0)
        {
            newBuckNum = Right[Rev[List[j]]] + 1;
            Rev[List[j]] = newBuckNum;
            j++;
        }

        // update Rev for the last element of the bucket

        newBuckNum = Right[Rev[List[j]]] + 1;

        Rev[List[j]] = newBuckNum;

        // correct the values of Right for all affected
        // buckets

        j = i;

        while (BuckList.getValN(j) == 0)
        {
            newBuckNum = Rev[List[j]];

            if (Right[newBuckNum] == -1)
            {
                Right[newBuckNum] = newBuckNum;
            }
            else
            {
                Right[newBuckNum]++;
            }
            if (newBuckNum > 0)
            {
                BuckB.setValN(newBuckNum - 1, 1);
            }

            j++;
        }

        // correct the values of Right for the last element

        newBuckNum = Rev[List[j]];

        if (Right[newBuckNum] == -1)
        {
            Right[newBuckNum] = newBuckNum;
        }
        else
        {
            Right[newBuckNum]++;
        }
        if (newBuckNum > 0)
        {
            BuckB.setValN(newBuckNum - 1, 1);
        }

        // set i to point to the first element of the next bucket

        i = j + 1;
    }

    // Reconsturct ArrayB from Rev.

    for (i = 0; i < inputLength; i++)
    {
        BucketNum = Rev[i];

        if (BucketNum > -1)
        {
            ArrayB[Right[BucketNum]] = i;
            Right[BucketNum]--;
        }

    }
    free(Rev);
    free(Right);
}
void counting_sort(const int* inputString, const int inputLength,
                   int* A, BooleanString& BuckA, int* buffer)
{
    //The int array buffer is used as a array of counters for
    //each character in the alphabet.

    int max, min, sigmaSize;
    int Temp, i, prevCount;

    assert(inputLength > 0);

    findMaxMin(inputString, inputLength, max, min);

    sigmaSize = max - min + 1;

    memset((void*) buffer, 0, sizeof(int)*sigmaSize);

    //Count the number of occurences of each character;

    for (i = 0; i < inputLength; i++)
    {
        Temp = inputString[i] - min;
        buffer[Temp]++;
    }

    //Convert buffer from a list that have the count of each
    //characters in the alphabet to a list that points to the
    //left boundaries of A. So it can be used when making A.

    prevCount = buffer[0];

    buffer[0] = 0;

    for (i = 1; i < sigmaSize; i++)
    {
        Temp = buffer[i];
        buffer[i] = prevCount + buffer[i - 1];
        prevCount = Temp;
    }

    //Constructing A. First find out which bucket a suffix goes,
    //then using buffer to calculate where should the suffix
    //be put in A.

    for (i = 0; i < inputLength; i++)
    {
        Temp = inputString[i] - min;
        A[buffer[Temp]] = i;
        buffer[Temp]++;
    }

    //Draw the bucket boundaries. BuckA[i] is defined to be true,
    //if there is an boundary between A[i] and A[i+1].

    BuckA.setAll(false);

    for (i = 0; i < sigmaSize; i++)
    {
        BuckA.setValN(buffer[i] - 1, 1);
    }
}



void counting_sort(const char* inputString, const int inputLength,
                   int* A, BooleanString& BuckA, int* buffer)
{
    int Temp, i, prevCount, sigmaSize;

    assert(inputLength > 0);

    sigmaSize = MAX_NUM_CHAR;

    memset((void*) buffer, 0, sigmaSize * sizeof(int));

    //Count the number of occurences of each character;

    for (i = 0; i < inputLength; i++)
    {
        Temp = (int) inputString[i];
        buffer[Temp]++;
    }

    //Convert buffer from a list that have the count of each
    //characters in the alphabet to a list that points to the
    //left boundaries of A. So it can be used when making A.

    prevCount = buffer[0];

    buffer[0] = 0;

    for (i = 1; i < sigmaSize; i++)
    {
        Temp = buffer[i];
        buffer[i] = prevCount + buffer[i - 1];
        prevCount = Temp;
    }

    //Constructing A. First find out which bucket a suffix goes,
    //then using buffer to calculate where should the suffix
    //be put in A.

    for (i = 0; i < inputLength; i++)
    {
        Temp = (int) inputString[i];
        A[buffer[Temp]] = i;
        buffer[Temp]++;
    }

    //Draw the bucket boundaries. BuckA[i] is defined to be true,
    //if there is an boundary between A[i] and A[i+1].

    BuckA.setAll(false);

    for (i = 0; i < sigmaSize; i++)
    {
        BuckA.setValN(buffer[i] - 1, 1);
    }
}



void counting_sort(const int* inputString, const int inputLength,
                   int* A, BooleanString& BuckA, const int minChar,
                   const int maxChar, int* buffer)
{
    int sigmaSize;
    int Temp, i, prevCount;

    assert(inputLength > 0);

    sigmaSize = maxChar - minChar + 1;

    memset((void*) buffer, 0, sizeof(int)*sigmaSize);

    //Count the number of occurences of each character;

    for (i = 0; i < inputLength; i++)
    {
        Temp = inputString[i] - minChar;
        buffer[Temp]++;
    }

    //Convert buffer from a list that have the count of each
    //characters in the alphabet to a list that points to the
    //left boundaries of A. So it can be used when making A.

    prevCount = buffer[0];

    buffer[0] = 0;

    for (i = 1; i < sigmaSize; i++)
    {
        Temp = buffer[i];
        buffer[i] = prevCount + buffer[i - 1];
        prevCount = Temp;
    }

    //Constructing A. First find out which bucket a suffix goes,
    //then using buffer to calculate where should the suffix
    //be put in A.

    for (i = 0; i < inputLength; i++)
    {
        Temp = inputString[i] - minChar;
        A[buffer[Temp]] = i;
        buffer[Temp]++;
    }

    //Draw the bucket boundaries. BuckA[i] is defined to be true,
    //if there is an boundary between A[i] and A[i+1].

    BuckA.setAll(false);

    for (i = 0; i < sigmaSize; i++)
    {
        BuckA.setValN(buffer[i] - 1, 1);
    }
}
void s_distance(const BooleanString& suffixType,
                const int inputLength, int* Dist, int*& DistCount,
                int& maxDist)
{
    int i, j;
    int prevCount, temp, prevDist;

    Dist[0] = 0;
    maxDist = 0;

    //Find the first S type, everything before that have Dist[i]=0.

    i = 0;

    while (suffixType.getValN(i) == 0)
    {
        Dist[i] = 0;
        i++;
    }

    Dist[i] = 0;

    if (i < inputLength - 1)
    {
        maxDist = 1;
    }

    //Find Dist for the rest of the string.

    j = i + 1;

    prevDist = 1;

    for (i = j; i < inputLength; i++)
    {
        Dist[i] = prevDist;

        //If suffix i is type L, then increment prevDist by 1.
        //else set prevDist to 1.

        prevDist = prevDist - prevDist * suffixType.getValN(i) + 1;

        //Calculate maxDist.

        temp = maxDist - Dist[i];
        temp = temp >> MAX_INT_INDEX;
        maxDist = maxDist + ((Dist[i] - maxDist) & temp);
    }

    //Initalize and count the number of suffix of each s-distant.

    DistCount = (int*) malloc((maxDist + 1) * sizeof(int));

    for (i = 0; i <= maxDist; i++)
    {
        DistCount[i] = 0;
    }

    // memset((void*) DistCount, 0, sizeof(int) * (maxDist + 1));

    //Skip the zero values.

    j = 0;

    while (Dist[j] == 0)
    {
        j++;
    }

    for (i = j; i < inputLength; i++)
    {
        DistCount[Dist[i] - 1]++;
    }

    //Compute the total number of suffixes that has a lesser s-distant.

    prevCount = DistCount[0];

    DistCount[0] = 0;

    for (i = 1; i <= maxDist; i++)
    {
        temp = DistCount[i];
        DistCount[i] = prevCount + DistCount[i - 1];
        prevCount = temp;
    }
}


//This function computes the l_distance.

void l_distance(const BooleanString& suffixType,
                const int inputLength, int* Dist, int*& DistCount,
                int& maxDist)
{
    int i, j;
    int prevCount, temp, prevDist;

    Dist[0] = 0;
    maxDist = 0;

    //Find the first L type, everything before that have Dist[i]=0.

    i = 0;

    while (suffixType.getValN(i) == 1)
    {
        Dist[i] = 0;
        i++;
    }

    Dist[i] = 0;

    if (i < inputLength - 1)
    {
        maxDist = 1;
    }

    //Find Dist for the rest of the string.

    j = i + 1;

    prevDist = 1;

    for (i = j; i < inputLength; i++)
    {
        Dist[i] = prevDist;

        //If suffix i is type S, then increment prevDist by 1.
        //else set prevDist to 1.

        temp = suffixType.getValN(i) - 1;
        prevDist = prevDist - (prevDist & temp) + 1;

        //Calculate maxDist.

        temp = maxDist - Dist[i];
        temp = temp >> MAX_INT_INDEX;
        maxDist = maxDist + ((Dist[i] - maxDist) & temp);
    }

    //Initalize and count the number of suffix of each l-distant

    DistCount = (int*) malloc((maxDist + 1) * sizeof(int));

    for (i = 0; i <= maxDist; i++)
    {
        DistCount[i] = 0;
    }

    //Skip the zero values.

    j = 0;

    while (Dist[j] == 0)
    {
        j++;
    }

    for (i = j; i < inputLength; i++)
    {
        DistCount[Dist[i] - 1]++;
    }

    //Compute the total number of suffixes that has a lesser l-distant.

    prevCount = DistCount[0];

    DistCount[0] = 0;

    for (i = 1; i <= maxDist; i++)
    {
        temp = DistCount[i];
        DistCount[i] = prevCount + DistCount[i - 1];
        prevCount = temp;
    }
}



