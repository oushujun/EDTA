#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

using namespace std;

#ifndef TIMING_H
#define TIMING_H

int gettimeofday(timeval* time);

class Timing
{

public:
    Timing();

    void outtime();
    void outtime(const char* message);
    void outtime(const char* m_before, const char* m_after);
    int markbeg();
    int markend();

private:

    timeval before;
    timeval after;
};

#endif
