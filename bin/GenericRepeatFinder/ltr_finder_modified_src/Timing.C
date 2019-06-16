#include "Timing.h"

int gettimeofday(timeval* time)
{

    struct timezone temp;

    return (gettimeofday(time, &temp));
}


Timing::Timing()
{
    markbeg();
}


int Timing::markbeg()
{
    return (gettimeofday(&before));
}


int Timing::markend()
{
    return (gettimeofday(&after));
}


void Timing::outtime()
{
    if ((after.tv_usec - before.tv_usec) / 1000 >= 0)
    {
        cout << after.tv_sec - before.tv_sec << '.';

        if ((after.tv_usec - before.tv_usec) / 1000 < 100)
        {
            cout << '0';
        }
        if ((after.tv_usec - before.tv_usec) / 1000 < 10)
        {
            cout << '0';
        }

        cout << (after.tv_usec - before.tv_usec) / 1000;
    }
    else
    {
        cout << after.tv_sec - before.tv_sec - 1 << '.';

        if ((1000 + (after.tv_usec - before.tv_usec) / 1000) < 100)
        {
            cout << '0';
        }
        if ((1000 + (after.tv_usec - before.tv_usec) / 1000) < 10)
        {
            cout << '0';
        }

        cout << 1000 + (after.tv_usec - before.tv_usec) / 1000;
    }

    cout << " second" << endl;
}


void Timing::outtime(const char* message)
{
    cout << message << " ";

    if ((after.tv_usec - before.tv_usec) / 1000 >= 0)
    {
        cout << after.tv_sec - before.tv_sec << '.';

        if ((after.tv_usec - before.tv_usec) / 1000 < 100)
        {
            cout << '0';
        }
        if ((after.tv_usec - before.tv_usec) / 1000 < 10)
        {
            cout << '0';
        }

        cout << (after.tv_usec - before.tv_usec) / 1000;
    }
    else
    {
        cout << after.tv_sec - before.tv_sec - 1 << '.';

        if ((1000 + (after.tv_usec - before.tv_usec) / 1000) < 100)
        {
            cout << '0';
        }
        if ((1000 + (after.tv_usec - before.tv_usec) / 1000) < 10)
        {
            cout << '0';
        }

        cout << 1000 + (after.tv_usec - before.tv_usec) / 1000;
    }

    cout << " second" << endl;
}


void Timing::outtime(const char* m_before, const char* m_after)
{
    cout << m_before << " ";

    if ((after.tv_usec - before.tv_usec) / 1000 >= 0)
    {
        cout << after.tv_sec - before.tv_sec << '.';

        if ((after.tv_usec - before.tv_usec) / 1000 < 100)
        {
            cout << '0';
        }
        if ((after.tv_usec - before.tv_usec) / 1000 < 10)
        {
            cout << '0';
        }

        cout << (after.tv_usec - before.tv_usec) / 1000;
    }
    else
    {
        cout << after.tv_sec - before.tv_sec - 1 << '.';

        if ((1000 + (after.tv_usec - before.tv_usec) / 1000) < 100)
        {
            cout << '0';
        }
        if ((1000 + (after.tv_usec - before.tv_usec) / 1000) < 10)
        {
            cout << '0';
        }

        cout << 1000 + (after.tv_usec - before.tv_usec) / 1000;
    }

    cout << m_after << endl;
}
