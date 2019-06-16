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



