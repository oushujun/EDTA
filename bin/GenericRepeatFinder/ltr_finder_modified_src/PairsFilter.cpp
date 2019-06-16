/*
* =====================================================================================
*
*       Filename:  PairsFilter.cpp
*
*    Description:  
*
*        Version:  1.0
*        Created:  2006年10月20日 09时07分25秒 CST
*       Revision:  none
*       Compiler:  gcc
*
*         Author:   (), 
*        Company:  
*
* =====================================================================================
*/
#include "PairsFilter.h"
#include "LinearSuffixSort.h"
#include "lcp.h"
#include <algorithm>
using namespace std;
void outstr(char *str, int pos, int len)
{
    for (int i = pos; i < pos + len && str[i] != '\0'; ++i)
        cout << str[i];
}
/* Construct suffix array and LCP on myString, find all exact match sequence pairs in different position
 * Delete all pairs that length less than Lex or longer than Lmax
 * Deltet all pairs that distance between them less than Dmin or longer than Dmax
 * Default is to find Left-Max match sequence pairs, but you can set the last parameter to 'FALSE' 
 *     to find all match sequence pairs
 * pairs returned in Pair
 */
int GetPairs (char *myString, int stringLength, int Lex, int Lmax, int Dmin, int Dmax,
              vector < candidate > &Pair, bool LeftMaxAlign)
{
    int *suffixArray;
    //cerr<<"seq len:"<<stringLength<<endl;
    suffixArray = LinearSuffixSort (myString, stringLength);
    ////suffixArray,total line, string, string lenth, substr length
    //outsuffix (suffixArray, 800, myString, stringLength - 1, 21);
    int *lcpA;
    lcpA = lcp (suffixArray, myString, stringLength);
    vector < vector < int > >Bk;
    //vector< vector<int>[5] > Lsets;
    vector < int >Lset[5];  //
    vector < int >bk;
    int CI[128];

    for (int i = 0; i < 128; ++i)
        CI[i] = 4;

    CI['A'] = 0;

    CI['C'] = 1;

    CI['G'] = 2;

    CI['T'] = 3;

    CI['N'] = 4;

    int flag = 0;

    for (int i = 0; i < stringLength; ++i)
    {
        if (lcpA[i] >= Lex)
        {
            if (flag == 0)
            {
                flag = 1;
                bk.clear ();
                bk.push_back (suffixArray[i - 1]);
            }

            bk.push_back (suffixArray[i]);
        }
        else if (flag == 1)
        {
            flag = 0;
            sort (bk.begin (), bk.end ());
            Bk.push_back (bk);
        }

    }
    candidate candi;

    for (int i = 0; i < Bk.size (); ++i)
    {
        for (int j = 0; j < Bk[i].size (); ++j)
        {
            if (LeftMaxAlign)
                Lset[CI[myString[Bk[i][j] - 1]]].push_back (Bk[i][j]);
            else
            {
                Lset[0].push_back (Bk[i][j]);
                Lset[1].push_back (Bk[i][j]);
            }

        }
        //check highly repeat region
        //if repeat more than 15 times
        //and consentrate in a short region 15x3000bp
        //then, it is a highly repeat region

        for (int m = 0; m < 5; ++m)
        {
            if (checkCentriole && Lset[m].size() >= 15)
            {
                int times = 0;

                for (int n = 0;n < Lset[m].size() - 1;++n)
                {
                    if (Lset[m][n + 1] - Lset[m][n] < 3000)
                        times++;
                    else
                        times = 0;

                    if (times >= 14)
                        break;
                }
                if (times >= 14) //repeat too much, centriole
                    Lset[m].clear();
            }
            //else if(Lset[m].size()>0)
            // cout<<"reserve:("<<Lset[m].size()<<") ave space:"<<(Lset[m][Lset[m].size()-1]-Lset[m][0])/Lset[m].size()<<endl;
        }
        for (int m = 0; m < 5; ++m)
        {
            for (int n = m + 1; n < 5; ++n)
            {
                for (int mi = 0; mi < Lset[m].size (); ++mi)
                {
                    for (int ni = 0; ni < Lset[n].size (); ++ni)
                    {
                        int dist = abs (Lset[m][mi] - Lset[n][ni]);

                        if (dist < (Dmax + Lmax) && dist > (Dmin + Lex)) //changed on 11.4. old: dist<Dmax&&dist>Dmin
                        {
                            //Pair.push(Lset[m][mi],Lset[n][ni]);
                            int len = Lex;

                            bool hasN = false; //change in v0.84, we need check N in pairs

                            for (int k = 0; k < len; ++k)
                            {
                                if (myString[Lset[m][mi] + k] == 'N' ||
                                        myString[Lset[n][ni] + k] == 'N')
                                {
                                    hasN = true;
                                    break;
                                }

                            }

                            if (hasN)
                                continue;

                            while (Lset[m][mi] + len < stringLength - 1 &&
                                    Lset[n][ni] + len < stringLength - 1 &&
                                    myString[Lset[m][mi] + len] ==
                                    myString[Lset[n][ni] + len] )
                            {
                                len++;
                            }
                            //cout<<" len:"<<len<<endl;
                            //if(len<=Lmax) //I think this will be faster, need test
                            {
                                candi.pos1 = Lset[m][mi];
                                candi.pos2 = Lset[n][ni];
                                candi.len = len;
                                Pair.push_back (candi);
                            }

                        }
                        //cout<<" "<<Lset[m][mi]<<"-"<<Lset[n][ni]<<endl;
                    }

                }

            }

        }

        //Lsets.push(Lset);

        for (int j = 0; j < 5; ++j)
            Lset[j].clear ();
    }

    for (int i = 0; i < Pair.size (); ++i)
    {
        if (Pair[i].pos1 > Pair[i].pos2)
        {
            int t = Pair[i].pos1;
            Pair[i].pos1 = Pair[i].pos2;
            Pair[i].pos2 = t;
        }

    }


    //cout<<"out put pair axis for draw(before del)"<<endl;

    if (CHECK_PAIRS)
    {
        sort (Pair.begin (), Pair.end (), pos1_sort);

        for (int i = 0;i < Pair.size();++i)
        {
            cout << Pair[i].pos1 << "-"
            << Pair[i].pos1 + Pair[i].len - 1 << " * "
            << Pair[i].pos2 << "-" << Pair[i].pos2 + Pair[i].len - 1
            << " len:"
            << Pair[i].len << endl;
        }

    }

    //del repeat candi
    //cout << "before del1, Pair size:" << Pair.size () << endl;
    sort (Pair.begin (), Pair.end (), pos2_sort);

    for (int i = Pair.size () - 1; i > 0; --i)
    {
        if (Pair[i].pos2 == Pair[i - 1].pos2 &&   //one side begin pos same
                Pair[i].pos1 + Pair[i].len - 1 <= Pair[i - 1].pos1 + Pair[i - 1].len - 1) //one side end pos included
            Pair.erase (Pair.begin () + i);
    }

    //cout << "before del2, Pair size:" << Pair.size () << endl;
    sort (Pair.begin (), Pair.end (), pos1_sort);

    for (int i = Pair.size () - 1; i > 0; --i)
    {
        if (Pair[i].pos1 == Pair[i - 1].pos1 &&
                Pair[i].pos2 + Pair[i].len - 1 <=
                Pair[i - 1].pos2 + Pair[i - 1].len - 1)
        {
            //cout<<Pair[i-1].pos1<<' '<<Pair[i-1].pos2<<' '<<Pair[i-1].len<<endl;
            //cout<<"del:"<<Pair[i].pos1<<' '<<Pair[i].pos2<<' '<<Pair[i].len<<endl;
            Pair.erase (Pair.begin () + i);
        }

    }

    if (showPairNum)
        cout << "Pairs:" << Pair.size () << endl;

    /*
    cout<<"out put pair axis for draw(after del)"<<endl;
    for(int i=0;i<Pair.size();++i)
    {
          cout<<Pair[i].pos1<<"\t"
            <<Pair[i].pos2<<"\t"
            <<Pair[i].len<<endl;
    }
    */
    free (lcpA);

    free (suffixArray);

    return Pair.size ();
}

void JoinPairs(char* myString, int stringLength, vector<candidate>& Pair, vector<stick>& sticks)
{
    stick st;
    vector < bool > used (Pair.size ()); //default init to false
    //    cout << "Pair size:" << Pair.size () << endl;
    //    for(int i=0;i<Pair.size();++i)
    //    {
    //        cout<<i<<" pos1:"<<Pair[i].pos1<<" pos2"<<Pair[i].pos2<<endl;
    //    }
    //    getchar();
    int total_used = 0;
    int beginPos = 0;
    bool is_extend = false;

    while (total_used < Pair.size ())
    {
        st.end1 = -1;
        st.candi.clear ();
        //cout<<"\nnew round\n";

        for (int i = beginPos; i < Pair.size (); ++i)
        {
            if (used[i])
                continue;

            if (st.end1 < 0) //first ele
            {
                st.candi.push_back (Pair[i]);
                st.pos1 = Pair[i].pos1;
                st.pos2 = Pair[i].pos2;
                st.end1 = Pair[i].len + st.pos1 - 1;
                st.len1 = Pair[i].len;
                st.len2 = st.len1;
                st.end2 = st.pos2 + st.len2 - 1;
                //st.gap = 0;       //Indels in pair
                st.match_score = st.len1; //in fact, the count of matched base
                st.match_len = st.len1; //in fact, the total match length
                used[i] = true;
                total_used++;
                is_extend = false;
                beginPos = i + 1;

                if (CHECK_PAIRS)
                    cout << "first pair: " << st.pos1 << "-" << st.end1 << " * " << st.pos2 << "-" << st.end2 << endl;
            }
            else
            {
                //int offset = Pair[i].pos2 - Pair[i].pos1;
                //    if (Pair[i].pos1 > st.end1 + MaxGap)
                //     break;
                //    if (Pair[i].pos2 > st.end2 + MaxGap || Pair[i].pos1 == st.candi[st.candi.size () - 1].pos1) //???
                //     continue;

                if (Pair[i].pos1 - st.pos1 > Lmax) //they are too long to become a LTR
                    break;

                if (CHECK_PAIRS)
                    cout << " judge:" << Pair[i].pos1 << "-" << Pair[i].pos1 + Pair[i].len - 1 << " * "
                    << Pair[i].pos2 << "-" << Pair[i].pos2 + Pair[i].len - 1 << endl;

                if (Pair[i].pos2 + Pair[i].len > st.end2 &&
                        Pair[i].pos1 + Pair[i].len > st.end1 &&
                        Pair[i].pos2 > st.pos2)
                {
                    int gap2 = (Pair[i].pos2 - (st.pos2 + st.len2));
                    int gap1 = (Pair[i].pos1 - (st.pos1 + st.len1));
                    int max_gap = 0;
                    float similarity = 0;
                    int match_score = 0;

                    if (gap1 <= 0 && gap1 <= gap2) //if gap2<0 and gap2 is bigger, gap1 must <0
                    {
                        gap2 += ( -gap1);
                        max_gap = gap2;
                        //gap_score = -(gap_open + gap_ext * (gap2-1) );
                        //score = score_match * Pair[i].len + gap_score;
                        int ext_len = Pair[i].len + gap1; //add a -xxx
                        similarity = (float)(ext_len + st.match_score) / (ext_len + max_gap + st.match_len);
                        match_score = ext_len;

                        if (CHECK_PAIRS)
                            cout << "similarity 1:" << similarity << endl;
                    }
                    else if (gap2 <= 0 && gap2 <= gap1)
                    {
                        gap1 += ( -gap2);
                        max_gap = gap1;
                        //gap_score = -(gap_open + gap_ext * (gap1-1) );
                        //score = score_match * Pair[i].len + gap_score;
                        int ext_len = Pair[i].len + gap2; //add a -xxx
                        similarity = (float)(ext_len + st.match_score) / (ext_len + max_gap + st.match_len);
                        match_score = ext_len;

                        if (CHECK_PAIRS)
                            cout << "similarity 2:" << similarity << endl;

                    }
                    else if ( abs(gap1 - gap2) > MaxGap ) //check the MaxGap first, before alignment, for performance
                    {
                        continue;
                    }
                    else //if the first MaxGap pass, then align, then check MaxGap again
                    {
                        //cerr<<"global align:"<<st.pos1+st.len1<<" to "<<Pair[i].pos1-1
                        //    <<"\n    target:"<<st.pos2+st.len2<<" to "<<Pair[i].pos2-1<<endl;
                        int match_count =
                            global_align (myString, st.pos1 + st.len1,
                                          Pair[i].pos1 - 1, myString,
                                          st.pos2 + st.len2,
                                          Pair[i].pos2 - 1, &max_gap);
                        //cout<<"align:"<<Pair[i].pos1 - 1-(st.pos1 + st.len1)+1<<" to "<<
                        //    Pair[i].pos2 - 1-(st.pos2 + st.len2)+1<<" score:"<<gap_score<<endl;
                        //score = score_match * Pair[i].len + gap_score;
                        //cerr<<"global align end\n";
                        int ext_len = gap1 > gap2 ? gap1 : gap2;
                        ext_len += Pair[i].len;
                        similarity = (double)(Pair[i].len + match_count + st.match_score) / (ext_len + st.match_len);
                        match_score = Pair[i].len + match_count;

                        if (CHECK_PAIRS)
                        {
                            cout << "ext_len:" << ext_len << " max_gap" << max_gap << " pair len" << Pair[i].len << endl;
                            cout << "similarity 3:" << similarity << endl;
                        }
                        //max_gap=0;
                    }

                    //check MaxGap
                    if (max_gap > MaxGap) //check MaxGap again
                        continue;

                    //cerr<<"score: "<<score<<endl;
                    if (similarity > JoinThreshold) //0)
                        //if ((gap1 <= 0 || gap2 <= 0) || //only Insertion or Deletion
                        // (abs (gap2 - gap1) / (gap1 + gap2) < 0.1 && gap_score > 0)) //InDel less than 1/5
                    {

                        if (similarity < SplitThreshold ) //|| max_gap > 10 )          //-(MaxGap*gap_open*(1-Sensitive)) )//gap_score is too low ??? need a define
                        {
                            //cerr <<"might be a gap, sim:"<<similarity<<" mgap:"<<max_gap<<endl;
                            is_extend = true;
                            sticks.push_back (st);

                            if (CHECK_PAIRS)
                                cout << "might be a gap, sim:" << similarity << " mgap:" << max_gap << endl;
                        }

                        if (Pair[i].pos1 + Pair[i].len - st.pos1 > Lmax ||
                                Pair[i].pos2 + Pair[i].len - st.pos2 > Lmax ||
                                st.pos2 - (Pair[i].pos1 + Pair[i].len) < Dmin)
                            continue; // too long

                        st.candi.push_back (Pair[i]);

                        int new_end1 = Pair[i].pos1 + Pair[i].len - 1;

                        //st.gap += gap;
                        int new_end2 = Pair[i].pos2 + Pair[i].len - 1;

                        st.match_score += match_score;

                        st.match_len = new_end1 - st.pos1 > new_end2 - st.pos2 ? new_end1 - st.pos1 + 1 : new_end2 - st.pos2 + 1;

                        st.end1 = new_end1;

                        st.end2 = new_end2;

                        st.len1 = st.end1 - st.pos1 + 1;

                        st.len2 = st.end2 - st.pos2 + 1;

                        if (!is_extend) //only if this is good align, marked as used.
                        {
                            used[i] = true;
                            total_used++;

                            if (CHECK_PAIRS)
                                cout << "sim:" << similarity << " mgap:" << max_gap << " could extend" << endl;
                        }
                        //cout<<"add one: "<<i<<endl;
                    }

                }

            }

        }

        if ((st.len1 + 2 * Lex > Lmin && st.len1 < Lmax) &&
                (st.len2 + 2 * Lex > Lmin && st.len2 < Lmax) &&
                st.pos2 - st.end1 > Dmin && st.pos2 - st.end1 < Dmax)
        {
            sticks.push_back (st);

            if (CHECK_PAIRS)
                cout << " join piar:" << st.pos1 << "-" << st.end1 << " * " << st.pos2 << "-" << st.end2 << endl;
        }
        else
        {
            if (CHECK_PAIRS)
            {
                cout << "DIDN'T join piar:" << st.pos1 << "-" << st.end1 << " * " << st.pos2 << "-" << st.end2 << endl;
                cout << "Lex:" << Lex << " Lmin:" << Lmin << " Lmax:" << Lmax << " Dmin:" << Dmin << " Dmax:" << Dmax << endl;
                cout << "L:" << st.pos1 << "-" << st.end1 << "  D:" << st.pos2 - st.end1 << endl;
            }

        }

    }
}
int FindEdge(char* str, int len, int win_size, int direct)
{
    int maxDiff = 0;
    int count1 = 0;
    int count2 = 0;
    int pos = -1;

    if (len < 2*win_size)
        return pos;

    for (int j = 0;j < win_size;++j)
    {
        if (str[j] == '|')
            count1++;

        if (str[j + win_size] == '|')
            count2++;
    }

    maxDiff = direct * (count2 - count1);
    pos = win_size;

    for (int i = 1;i <= len - 2*win_size;++i)
    {
        count1 -= (str[i - 1] == '|');
        count1 += (str[i + win_size - 1] == '|');
        count2 -= (str[i + win_size - 1] == '|');
        count2 += (str[i + 2 * win_size - 1] == '|');

        if (direct == 1) //5'
        {

            if (maxDiff <= (count2 - count1)*direct) //the most right result
            {
                maxDiff = (count2 - count1) * direct;
                pos = i + win_size;
            }

        }
        else
        {
            if (maxDiff < (count2 - count1)*direct) //the most left result
            {
                maxDiff = (count2 - count1) * direct;
                pos = i + win_size;
            }

        }

    }

    if (direct == 1)
        return pos;
    else
        return pos -1;
}
void ScoreCompensate(char* myString, TG_CA_TSR& tct)
{
    tct.status = 0;

    if (myString[tct.tg_pos1] == 'T' &&
            myString[tct.tg_pos1 + 1] == 'G')
        tct.status = tct.status | M5TG;

    if (myString[tct.tg_pos2] == 'T' &&
            myString[tct.tg_pos2 + 1] == 'G')
        tct.status = tct.status | M3TG;

    if (myString[tct.ca_pos1] == 'C' &&
            myString[tct.ca_pos1 + 1] == 'A')
        tct.status = tct.status | M5CA;

    if (myString[tct.ca_pos2] == 'C' &&
            myString[tct.ca_pos2 + 1] == 'A')
        tct.status = tct.status | M3CA;

    tct.signal_score = 0;

    if (myString[tct.tg_pos1] != 'T')
        tct.score -= 5;
    else
        tct.signal_score += 0.5;

    if (myString[tct.tg_pos2] != 'T')
        tct.score -= 5;
    else
        tct.signal_score += 0.5;

    if (myString[tct.ca_pos1] != 'C')
        tct.score -= 5;
    else
        tct.signal_score += 0.5;

    if (myString[tct.ca_pos2] != 'C')
        tct.score -= 5;
    else
        tct.signal_score += 0.5;

    if (myString[tct.tg_pos1 + 1] != 'G')
        tct.score -= 5;
    else
        tct.signal_score += 0.5;

    if (myString[tct.tg_pos2 + 1] != 'G')
        tct.score -= 5;
    else
        tct.signal_score += 0.5;

    if (myString[tct.ca_pos1 + 1] != 'A')
        tct.score -= 5;
    else
        tct.signal_score += 0.5;

    if (myString[tct.ca_pos2 + 1] != 'A')
        tct.score -= 5;
    else
        tct.signal_score += 0.5;

}
float Sharpness(char* str, char* aux, int len, int pos, int win_size)
{
    if (aux != NULL) //refine pos
    {
        int i = 0;

        for (int c = 0;i < len && c <= pos;++i)
        {
            if (aux[i] != '-')
                ++c;
        }

        pos = i - 1;
    }

    int cb = 0;
    int ce = 0;
    int wb = 0;
    int we = 0;
    int i = pos - win_size;

    if (i < 0)
        i = 0;

    for (;i < pos && i < len;++i)
    {
        wb++;

        if (str[i] == '|')
            cb++;
    }
    for (i = pos + 1;i < pos + 1 + win_size && i < len;++i)
    {
        we++;

        if (str[i] == '|')
            ce++;
    }

    return fabs((float)ce / we - (float)cb / wb);
}
void ShorterAlignStr(string& str)
{
    int output_len = outAlignLen; //global variable
    vector<int> vpos;
    vpos.push_back( -1);
    int pos = 0;

    while ( (pos = str.find_first_of('\n', pos)) != string::npos )
    {
        vpos.push_back(pos);
        ++pos;
    }

    int len = vpos[2] - vpos[1] - 1;
    int p1 = str.find_first_of('|', 0) - vpos[0];
    int p2 = str.find_first_of('*', vpos[2]) - vpos[2];
    int begin = (int)((p1 + p2) * 0.5) - output_len;
    int end = (int)((p1 + p2) * 0.5) + output_len;

    if (begin < 0)
        begin = 0;

    if (end >= len)
        end = len - 1;

    string ret;

    for (int i = 0;i < 5;++i)
    {
        int this_end = end;

        if (this_end > vpos[i + 1])
            this_end = vpos[i + 1] - 1;

        ret += str.substr(vpos[i] + begin + 1, this_end - begin + 1);

        ret += "\n";
    }

    int n_count = 0;

    for (int i = ret.length() - 1;i >= 0;--i)
    {
        if (ret[i] == '\n')
            n_count++;
        else
            break;
    }

    ret = ret.substr(0, ret.length() - n_count + 1);
    ret += "\n";
    str = ret;
}
void RefineSimilarity(char* str, char* aux, int len, int pos1, int pos2, stick& st)
{
    int i = 0;

    for (int c = 0;i < len && c <= pos1;++i)
    {
        if (aux[i] != '-')
            ++c;
    }

    pos1 = i - 1;
    i = 0;

    for (int c = 0;i < len && c <= pos2;++i)
    {
        if (aux[i] != '-')
            ++c;
    }

    pos2 = i - 1;

    //cerr<<"str:"<<str<<endl;
    //cerr<<"len:"<<len<<" pos1:"<<pos1<<" pos2:"<<pos2<<" ";
    //cerr<<"'";

    if (pos1 < pos2)
    {
        for (int i = pos1;i <= pos2;++i)
        {
            if (str[i] == '|')
                ++st.match_score;

            ++st.match_len;

            //cerr<<str[i];
        }

        --st.match_len;
        --st.match_score;
    }
    else if (pos1 > pos2)
    {
        for (int i = pos2;i <= pos1;++i)
        {
            if (str[i] == '|')
                --st.match_score;

            --st.match_len;

            //cerr<<str[i];
        }

        ++st.match_len;
        ++st.match_score;
    }
    //cerr<<"'";
    //cerr<<"   match:"<<st.match_score<<"  len:"<<st.match_len<<endl;
}
bool ExtendPairs(char* myString, int stringLength, stick& st)
{
    //some const
    int in_len = 80; //4*Lex;
    int out_len = 120; //6*Lex;
    int win_size = 70;

    bool full5end = true;
    bool full3end = true;
    st.sharpness5 = 1;
    st.sharpness3 = 1;

    int aln_begin[4]; //

    //check in_len
    //if(in_len>st.candi[0].len)
    //{
    //    in_len=st.candi[0].len;
    //    win_size=in_len-(int)(Lex*0.5);
    //}

    int begin1 = st.pos1 - out_len;

    if (begin1 < 0)
    {
        full5end = false;
        begin1 = 0;
    }

    int end1 = st.pos1 + in_len - 1;

    if (end1 >= stringLength)
        end1 = stringLength - 1;

    int len1 = end1 - begin1 + 1;

    int begin2 = st.pos2 - out_len;

    if (begin2 < 0)
        begin2 = 0;

    int end2 = st.pos2 + in_len - 1;

    if (end2 >= stringLength)
        end2 = stringLength - 1;

    int len2 = end2 - begin1 + 1;

    if (len2 > len1) //length must equal
    {
        len2 = len1; //cut seq2
        begin2 = end2 - len2 + 1;
    }

    AlnAln *aln_global1, *aln_global2;
    aln_begin[0] = begin1;
    aln_begin[1] = begin2;
    aln_global1 = aln_stdaln_aux(myString + begin1, myString + begin2, &aln_param_nt2nt, 1,
                                 len1, len2); //right align first
    int pos = FindEdge(aln_global1->outm, aln_global1->path_len, win_size, 1);

    if (pos < 0)
        pos = st.pos1 - begin1;

    if (full5end) //pre calculate sharpness, del obvious false
        st.sharpness5 = Sharpness(aln_global1->outm, 0, aln_global1->path_len, pos, win_size);

    if (st.sharpness5 < LowerSharpness)
    {
        st.score = -1;
        aln_free_AlnAln(aln_global1);
        return false;
    }

    int posm1 = pos;
    int pos1 = 0;

    for (int i = 0;i < pos;++i)
        if (aln_global1->out1[i] != '-')
            ++pos1;

    //cerr<<"pre pos1:"<<st.pos1<<endl;
    st.pos1 = begin1 + pos1;

    //cerr<<"final pos1:"<<st.pos1<<endl;
    int pos2 = 0;

    for (int i = 0;i < pos;++i)
        if (aln_global1->out2[i] != '-')
            ++pos2;

    //cerr<<"pre pos2:"<<st.pos2<<endl;
    st.pos2 = begin2 + pos2;

    //cerr<<"final pos2:"<<st.pos2<<endl;


//    //check in_len
//    assert(st.candi.size() > 0);

    //if(in_len>st.candi[st.candi.size()-1].len)
    //{
    //    in_len=st.candi[st.candi.size()-1].len;
    //    win_size=in_len-(int)(Lex*0.5);
    //}

    begin1 = st.end1 - in_len + 1;

    if (begin1 < 0)
        begin1 = 0;

    end1 = st.end1 + out_len;

    if (end1 >= stringLength)
        end1 = stringLength - 1;

    len1 = end1 - begin1 + 1;

    begin2 = st.end2 - in_len + 1;

    if (begin2 < 0)
        begin2 = 0;

    end2 = st.end2 + out_len;

    if (end2 >= stringLength)
    {
        full3end = false;
        end2 = stringLength - 1;
    }

    len2 = end2 - begin2 + 1;

    if (len1 > len2)
    {
        len1 = len2;
        end1 = begin1 + len1 - 1;
    }

    aln_begin[2] = begin1;
    aln_begin[3] = begin2;
    reverse(myString + begin1, myString + (begin1 + len1)); //left align first
    reverse(myString + begin2, myString + (begin2 + len2));
    aln_global2 = aln_stdaln_aux(myString + begin1, myString + begin2, &aln_param_nt2nt, 1,
                                 len1, len2);
    reverse(myString + begin1, myString + (begin1 + len1));
    reverse(myString + begin2, myString + (begin2 + len2));
    reverse(aln_global2->out1, aln_global2->out1 + aln_global2->path_len);
    reverse(aln_global2->out2, aln_global2->out2 + aln_global2->path_len);
    reverse(aln_global2->outm, aln_global2->outm + aln_global2->path_len);
    pos = FindEdge(aln_global2->outm, aln_global2->path_len, win_size, -1);

    if (pos < 0)
        pos = in_len - 1;

    if (full3end) //pre calculate, del obvious false
        st.sharpness3 = Sharpness(aln_global2->outm, 0, aln_global2->path_len, pos, win_size);

    if (st.sharpness3 < LowerSharpness)
    {
        st.score = -1;
        aln_free_AlnAln(aln_global1);
        aln_free_AlnAln(aln_global2);
        return false;
    }
    if (st.sharpness3 < HigherSharpness && st.sharpness5 < HigherSharpness) //no one edge satisfy the score limit
    {
        st.score = -1;
        aln_free_AlnAln(aln_global1);
        aln_free_AlnAln(aln_global2);
        return false;
    }

    int posm2 = pos;
    pos1 = 0;

    for (int i = 0;i < pos;++i)
        if (aln_global2->out1[i] != '-')
            ++pos1;

    st.end1 = begin1 + pos1;

    pos2 = 0;

    for (int i = 0;i < pos;++i)
        if (aln_global2->out2[i] != '-')
            ++pos2;

    st.end2 = begin2 + pos2;

    st.len1 = st.end1 - st.pos1 + 1;

    st.len2 = st.end2 - st.pos2 + 1;

    //find TG..CA
    vector<int> tg_pos1, tg_pos2, ca_pos1, ca_pos2;

    int tsr_len = -1;

    TG_CA_TSR tg_ca_tsr;

    vector<TG_CA_TSR> tct;

    FindTwoChar(aln_global1, tg_pos1, tg_pos2, "TG");

    FindTwoChar(aln_global2, ca_pos1, ca_pos2, "CA");

    if (CHECK_PAIRS)
    {
        cout << "st.pos1:" << st.pos1 + 1 << " st.end2:" << st.end2 + 1 << endl;
    }

    for (int i = 0;i < tg_pos1.size();++i)
        for (int j = 0;j < ca_pos2.size();++j)
        {
            tsr_len = FindTSR(myString, aln_begin[0] + tg_pos1[i],
                              aln_begin[3] + ca_pos2[j]); //even if tsr_len==0
            {
                tg_ca_tsr.tg_pos1 = aln_begin[0] + tg_pos1[i];
                tg_ca_tsr.tg_pos2 = aln_begin[1] + tg_pos2[i];
                tg_ca_tsr.ca_pos1 = aln_begin[2] + ca_pos1[j];
                tg_ca_tsr.ca_pos2 = aln_begin[3] + ca_pos2[j];
                tg_ca_tsr.tsr_pos1 = tg_ca_tsr.tg_pos1 - tsr_len;
                tg_ca_tsr.tsr_pos2 = tg_ca_tsr.ca_pos2 + 2;
                tg_ca_tsr.score = -abs(tg_ca_tsr.tg_pos1 - st.pos1)
                                  - abs(tg_ca_tsr.ca_pos2 + 1 - st.end2);
                ScoreCompensate(myString, tg_ca_tsr); //set score lower, if not exact TG..CA

                if (tsr_len > 0)
                {
                    tg_ca_tsr.signal_score += 1;
                    tg_ca_tsr.status = tg_ca_tsr.status | MTSR;
                }
                if (CHECK_PAIRS)
                {
                    cout << "tg ca score:" << tg_ca_tsr.score << " tsr_len:" << tsr_len << endl;
                    cout << "  TG offset:" << tg_ca_tsr.tg_pos1 - st.pos1
                    << " TG pos:" << tg_ca_tsr.tg_pos1 + 1 << endl;
                    cout << "  CA offset:" << tg_ca_tsr.ca_pos2 + 1 - st.end2
                    << " CA pos:" << tg_ca_tsr.ca_pos2 + 2 << endl;
                }

                if (tg_ca_tsr.score > -40) //???, this score is the min limit of tg_ca_tsr
                    tct.push_back(tg_ca_tsr);
            }

        }

    if (tct.size())
    {
        sort(tct.begin(), tct.end());

        if (CHECK_PAIRS)
            cout << "choose tg ca score:" << tct[0].score << " tsr len:" << tct[0].tg_pos1 - tct[0].tsr_pos1 << endl;

        st.tg_pos1 = tct[0].tg_pos1;

        st.tg_pos2 = tct[0].tg_pos2;

        st.ca_pos1 = tct[0].ca_pos1;

        st.ca_pos2 = tct[0].ca_pos2;

        st.tsr_pos1 = tct[0].tsr_pos1;

        st.tsr_pos2 = tct[0].tsr_pos2;

        st.tsr_len = st.tg_pos1 - st.tsr_pos1;

        st.pos1 = st.tg_pos1;

        st.pos2 = st.tg_pos2;

        st.end1 = st.ca_pos1 + 1;

        st.end2 = st.ca_pos2 + 1;

        st.len1 = st.end1 - st.pos1 + 1;

        st.len2 = st.end2 - st.pos2 + 1;

        //first time to asign st.score and st.status
        st.score = tct[0].signal_score; //TG 1, CA 1

        st.status = tct[0].status;

        //        //after boundary settled, calculate Sharpness and LTR Similarity
        //        if(full5end)
        //            st.sharpness5=Sharpness(aln_global1->outm,aln_global1->out1,
        //                    aln_global1->path_len,st.pos1-aln_begin[0],win_size);
        //        if(st.sharpness5<LowerSharpness)
        //        {
        //            st.score=-1;
        //            return false;
        //        }
        //
        //        if(full3end)
        //            st.sharpness3=Sharpness(aln_global2->outm,aln_global2->out1,
        //                    aln_global2->path_len,st.end1-aln_begin[2],win_size);
        //        if(st.sharpness3<LowerSharpness)
        //        {
        //            st.score=-1;
        //            return false;
        //        }
        //        if(st.sharpness3<HigherSharpness && st.sharpness5<HigherSharpness)//no one edge satisfy the score limit
        //        {
        //            st.score=-1;
        //            return false;
        //        }


//        RefineSimilarity(aln_global1->outm, aln_global1->out1, aln_global1->path_len,
//                         st.pos1 - aln_begin[0], st.candi[0].pos1 - aln_begin[0], st);
//
//        RefineSimilarity(aln_global2->outm, aln_global2->out1, aln_global2->path_len,
//                         st.candi[(int)st.candi.size() - 1].pos1 + st.candi[(int)st.candi.size() - 1].len - 1 - aln_begin[2],
//                         st.end1 - aln_begin[2], st);
        
        RefineSimilarity(aln_global1->outm, aln_global1->out1, aln_global1->path_len,
                 st.pos1 - aln_begin[0], st.pos1 - aln_begin[0], st);

        RefineSimilarity(aln_global2->outm, aln_global2->out1, aln_global2->path_len,
                 st.end2 - aln_begin[2],
                 st.end1 - aln_begin[2], st);

        //cerr<<"--------------\n";

        char buf[1024];

        int i, c;

        for (i = 0, c = 0;aln_begin[1] + c < st.pos2 || aln_global1->out2[i] == '-';++i)
        {
            if (aln_global1->out2[i] != '-')
                ++c;

            st.LTR5 += " ";
        }

        sprintf(buf, "|%d\n", st.pos2 + 1);
        st.LTR5 += buf;
        st.LTR5 += aln_global1->out2;
        st.LTR5 += "\n";
        aln_global1->outm[posm1] = '*';
        st.LTR5 += aln_global1->outm;
        st.LTR5 += "\n";
        st.LTR5 += aln_global1->out1;
        st.LTR5 += "\n";

        for (i = 0, c = 0;aln_begin[0] + c < st.tsr_pos1 || aln_global1->out1[i] == '-';++i)
        {
            if (aln_global1->out1[i] != '-')
                ++c;

            st.LTR5 += " ";
        }
        for (c = 0;c < st.tsr_len || aln_global1->out1[i] == '-';++i)
        {
            if (aln_global1->out1[i] != '-')
            {
                ++c;
                st.LTR5 += "*";
            }
            else
                st.LTR5 += "-";
        }

        sprintf(buf, "|%d\n", st.pos1 + 1);
        st.LTR5 += buf;


        sprintf(buf, "%d|", st.end2 + 1);
        //sprintf(buf,"|%d <--",st.end2+1);
        int d_len = strlen(buf);
        //d_len=1;//no need to calculate this

        for (i = 0, c = 0;aln_begin[3] + c < st.end2 - d_len + 1 || aln_global2->out2[i] == '-';++i)
        {
            if (aln_global2->out2[i] != '-')
                ++c;

            st.LTR3 += " ";
        }
        for (int j = i;j < i + d_len;++j)
            if (aln_global2->out2[j] == '-')
                st.LTR3 += " ";

        st.LTR3 += buf;

        i += d_len;

        for (c = 0;c < st.tsr_len || (aln_global2->out2[i] == '-' && st.tsr_len);++i)
        {
            if (aln_global2->out2[i] != '-')
            {
                ++c;
                st.LTR3 += "*";
            }
            else
                st.LTR3 += "-";
        }

        st.LTR3 += "\n";
        st.LTR3 += aln_global2->out2;
        st.LTR3 += "\n";
        aln_global2->outm[posm2] = '*';
        st.LTR3 += aln_global2->outm;
        st.LTR3 += "\n";
        st.LTR3 += aln_global2->out1;
        st.LTR3 += "\n";
        sprintf(buf, "%d|", st.end1 + 1);
        //sprintf(buf,"|%d <--",st.end1+1);
        d_len = strlen(buf);
        //d_len=1;

        for (i = 0, c = 0;aln_begin[2] + c < st.end1 - d_len + 1 || aln_global2->out1[i] == '-';++i)
        {
            if (aln_global2->out1[i] != '-')
                ++c;

            st.LTR3 += " ";
        }
        for (int j = i;j < i + d_len - 1 || aln_global2->out1[j] == '-';++j)
        {
            if (aln_global2->out1[j] == '-')
                st.LTR3 += " ";
        }

        st.LTR3 += buf;
        st.LTR3 += "\n";

        ShorterAlignStr(st.LTR5);
        ShorterAlignStr(st.LTR3);


    }
    else
    {
        //try to find TSR independent
    }

    aln_free_AlnAln(aln_global1);
    aln_free_AlnAln(aln_global2);
    return true;
}
int FindTSR( char* myString, int pos1, int pos2 )
{
    int len = 0;
    int stringLen = strlen(myString);

    for ( int i = pos1 - TSRmax;i <= pos1 - TSRmin;++i )
    {
        if ( i < 0 || pos2 + 2 >= stringLen )
            continue;

        if ( myString[ i ] == myString[ pos2 + 2 ]
                && myString[ i ] != 'N' ) //N can not in TSR
        {
            int m = 3, flag = 1;

            for ( int j = i + 1;j < pos1;++m, ++j )
            {
                if ( pos2 + m >= stringLen || (myString[ j ] != myString[ pos2 + m ] ||
                                                        myString[ j ] == 'N' ||
                                                        myString[ pos2 + m ] == 'N') )
                {
                    flag = 0;
                    break;
                }

            }

            if ( flag )
            {
                len = pos1 - i;
                break;
            }

        }

    }
    return len;
}
void FindTwoChar(AlnAln *aln, vector<int>& pos1, vector<int>& pos2, char* sub) //sub = 'TG' or 'CA'
{
    int p1 = 0, p2 = 0;

    for (int i = 0;i < aln->path_len - 1;++i)
    {
        //if(aln->outm[i]=='|'
        //        && aln->outm[i+1]=='|'
        //        && aln->out1[i]==sub[0]
        //        && aln->out1[i+1]==sub[1])

        if ( (aln->out2[i] == sub[0]
                && aln->out2[i + 1] == sub[1]
                && aln->out1[i] != '-'
                && aln->out1[i + 1] != '-')
                ||
                (aln->out1[i] == sub[0]
                 && aln->out1[i + 1] == sub[1]
                 && aln->out2[i] != '-'
                 && aln->out2[i + 1] != '-'
                )
           )
        {
            pos1.push_back(p1);
            pos2.push_back(p2);
        }
        if (aln->out1[i] != '-')
            ++p1;

        if (aln->out2[i] != '-')
            ++p2;
    }

}

void FindSignal(char* myString, int maxLen, PBS& pbs, PPT& ppt, stick& st)
{
    PBS_PPT pp;
    string tmp_name;
    int tmp_begin = 0;
    int tmp_end = 0;
    int tmp_match = 0;
    int tmp_count = 0;

    //PBS
    pp.PBS = false;
    pp._PBS = false;
    pp.PPT = false;
    pp._PPT = false;
    pp._PPT_count = 0;
    pp.PPT_count = 0;

    int search_len = PBS_region;

    if (search_len + st.pos1 + st.len1 > maxLen)
    {
        search_len = maxLen - (st.pos1 + st.len1);
    }
    if ( pbs.Search(myString + (st.pos1 + st.len1), search_len, tmp_name,
                    tmp_begin, tmp_end, tmp_match, pp.pbs_str))
    {
        pp.PBS = true;
        pp.PBS_name = tmp_name;
        pp.PBS_begin = tmp_begin + st.pos1 + st.len1;
        pp.PBS_end = tmp_end + st.pos1 + st.len1;
        pp.PBS_match = tmp_match;
    }

    search_len = PBS_region;

    if (st.pos2 - search_len < 0)
    {
        search_len = st.pos2;
    }
    if ( pbs.rSearch(myString + (st.pos2 - search_len), search_len, tmp_name,
                     tmp_begin, tmp_end, tmp_match, pp._pbs_str))
    {
        pp._PBS = true;
        pp._PBS_name = tmp_name;
        pp._PBS_begin = tmp_begin + st.pos2 - search_len;
        pp._PBS_end = tmp_end + st.pos2 - search_len;
        pp._PBS_match = tmp_match;
    }

    search_len = PPT_region;

    if (search_len + st.pos1 + st.len1 > maxLen)
    {
        search_len = maxLen - (st.pos1 + st.len1);
    }
    if ( ppt.rSearch(myString + (st.pos1 + st.len1), search_len,
                     tmp_begin, tmp_count))
    {
        pp._PPT = true;
        pp._PPT_begin = tmp_begin + st.pos1 + st.len1;
        pp._PPT_count = tmp_count;
        //cout<<"_ppt: "<<pp._PPT_begin<<" "<<tmp_count<<endl;
    }

    search_len = PPT_region;

    if (st.pos2 - search_len < 0)
    {
        search_len = st.pos2;
    }
    if ( ppt.Search(myString + (st.pos2 - search_len), search_len,
                    tmp_begin, tmp_count))
    {
        pp.PPT = true;
        pp.PPT_begin = tmp_begin + st.pos2 - search_len;
        pp.PPT_count = tmp_count;
        //cout<<"ppt: "<<pp.PPT_begin<<" "<<tmp_count<<endl;
    }
    /*
        if(pp.PBS && pp.PPT)
        {
            st.PBS_name= pp.PBS_name;
            st.PBS_begin=pp.PBS_begin;
            st.PBS_end=pp.PBS_end;
            st.PBS_str=pbs_str;
            st.PPT_begin=pp.PPT_begin;
            st.PPT_count=pp.PPT_count;
            st.strand=1;
        }
        else if(pp._PBS && pp._PPT)
        {
            st.PBS_name= pp._PBS_name;
            st.PBS_begin=pp._PBS_begin;
            st.PBS_end=pp._PBS_end;
            st.PBS_str=_pbs_str;
            st.PPT_begin=pp._PPT_begin;
            st.PPT_count=pp._PPT_count;
            st.strand=-1;
        }
        else if(pp.PBS || pp.PPT)
        {
            st.PBS_name= pp.PBS_name;
            st.PBS_begin=pp.PBS_begin;
            st.PBS_end=pp.PBS_end;
            st.PBS_str=pbs_str;
            st.PPT_begin=pp.PPT_begin;
            st.PPT_count=pp.PPT_count;
            st.strand=1;
        }
        else if(pp._PBS || pp._PPT)
        {
            st.PBS_name= pp._PBS_name;
            st.PBS_begin=pp._PBS_begin;
            st.PBS_end=pp._PBS_end;
            st.PBS_str=_pbs_str;
            st.PPT_begin=pp._PPT_begin;
            st.PPT_count=pp._PPT_count;
            st.strand=-1;
        }
    */
    st.pp = pp;
}

//this function is lack of thinking
void EraseOverlap(vector<stick>& sticks)
{
    sort(sticks.begin(), sticks.end(), stick_pos1_sort);

    for (int i = 0;i < sticks.size();++i)
    {
        if (sticks[i].score < 0)
            continue;

        for (int j = i + 1;j < sticks.size();++j)
        {
            if (sticks[j].pos1 > sticks[i].end2)
                break;

            if (sticks[j].score < 0)
                continue;

            stick *st1, *st2;

            st1 = &sticks[i];

            st2 = &sticks[j];

            cout << endl;

            cout << "i:" << st1->pos1 + 1 << "-" << st1->end1 + 1 << endl;

            cout << "j:" << st2->pos1 + 1 << "-" << st2->end1 + 1 << endl;

            int oLen1 = st1->end2 - st1->pos1 + 1;

            int oLen2 = st2->end2 - st2->pos1 + 1;

            int oLen = oLen1 < oLen2 ? oLen1 : oLen2;

            int begin = st1->pos1 > st2->pos1 ? st1->pos1 : st2->pos1;

            int end = st1->end2 < st2->end2 ? st1->end2 : st2->end2;

            if ( (float)(end - begin + 1) / oLen >= 0.5 )
            {
                cout << "judge\n";

                if ( (st1->status & 7) == (st2->status & 7) ) //both has/hasn't protein domain
                {

                    if ( (st1->status & MTSR) == (st2->status & MTSR) ) //both has TSR
                    {

                        if (st1->score < st2->score)
                            st1->score = -1; //del st1
                        else if (st1->score > st2->score)
                            st2->score = -1;

                        //if the score is same, do nothing
                    }
                    else if ( (st1->status & MTSR) < (st2->status & MTSR) ) //st2 has TSR
                        st1->score = -1;
                    else
                        st2->score = -1;
                }
                else if ( (st1->status & 7) == (st2->status & 7) ) //st2 has domain
                    st1->score = -1;
                else
                    st2->score = -1;
            }

        }

    }
}

