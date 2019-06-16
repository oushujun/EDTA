// Copyright (C) 2016 shij@miamioh.edu

#include <fstream>
#include <sstream>
#include <cctype>
#include <tuple>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <ctime>
#include <omp.h>
#include "DetectIntersperse.h"

using namespace std;

DetectIntersperse::DetectIntersperse(const parameter & param) : p(param) {
}

void DetectIntersperse::run() {
    // start
    time_t rawtime;
    time(&rawtime);
    cout << ctime(&rawtime);
    cout << "start" << endl;
    // read input data
    fastaread(chroms, genomeData, p.input);
    // score, vector index
    unordered_map<string, unsigned> index;
    // calculate cum score for all seeds in each chrom
    for (size_t i = 0; i < chroms.size(); i++) {
        time(&rawtime);
        cout << ctime(&rawtime);
        cout << "calculating cumulative scores for chromosome " << chroms[i]
                << "..." << endl;
        score(i);
        time(&rawtime);
        cout << ctime(&rawtime);
        cout << "grouping sequences with same scores..." << endl;
        groupScore(index, i);
    }
    index.clear();
    cout << "found " << groupedScore.size() << " groups" << endl;
    // find repeats within groups
    time(&rawtime);
    cout << ctime(&rawtime);
    cout << "detecting candidate repeats within groups..." << endl;
    result.resize(groupedScore.size());
#pragma omp parallel for
    for (size_t i = 0; i < groupedScore.size(); i++) {
        findCopy(groupedScore[i], result[i]);
    }
    groupedScore.clear();
    groupedScore.shrink_to_fit();
    // size for result2
    unsigned n = 0;
    for (auto & v : result) {
        n += v.size();
    }
    result2.resize(n);
    unsigned i = 0;
    for (auto & v : result) {
        for (auto & g: v) {
            result2[i] = move(g);
            i++;
        }
        v.clear();
        v.shrink_to_fit();
    }    
    result.clear();
    result.shrink_to_fit();
    cout << "found " << result2.size() << " candidate repeat groups" << endl;
    if (result2.size() == 0) {
        return;
    }
    // extend seed region
    time(&rawtime);
    cout << ctime(&rawtime);
    cout << "extending candidate repeats..." << endl;    
#pragma omp parallel for
    for (size_t i = 0; i < result2.size(); i++) {        
        // extend right
        extendR(result2[i]);
        // extend left
        extendL(result2[i]);
        // remove repeats with > max mismatches
        filter(result2[i]);
    }
    // merge gruop and output
    time(&rawtime);
    cout << ctime(&rawtime);
    cout << "merging groups and printing..." << endl;
    output();
    // end
    time(&rawtime);
    cout << ctime(&rawtime);
    cout << "end" << endl;
}

void DetectIntersperse::fastaread(vector<string> & id, vector<string> & seq,
        const string & file) {
    ifstream fid(file);
    if (!fid.good()) {
        cerr << "can't open " << file << endl;
        exit(1);
    }
    string line;
    string header = "";
    while (getline(fid, line)) {
        if (line[0] == '>') {
            istringstream ss(line);
            ss >> header;
            header = header.substr(1);
            id.push_back(header);
            seq.push_back("");
        } else {
            seq.back() += line;
        }
    }
    fid.close();
}

void DetectIntersperse::score(size_t count) {
    string & gseq = genomeData[count];
    // 1. Check genome sequence.
    size_t gsize = gseq.size();
    if (gsize < p.seed) {
        cerr << "Chromosome " << chroms[count]
                << " is too short. Can't perform detection." << endl;
        exit(1);
    }
    for (char & c : gseq) {
        c = std::toupper(c);
    }
    // 2. Mapping.
    // a t c g
    vector<tuple<unsigned, unsigned, unsigned, unsigned> > cum(gsize);
    mapVector(gseq, cum);
    // 3. Calculate cumulative score.
    cumsum(cum);
    // 4. Calculate cum scores for sub seqs
    auto & v = cumScore;
    v.resize(gsize - p.seed + 1);    
    v[0] = cum[p.seed - 1];
    transform(gseq.substr(0, p.seed), v[0]);
#pragma omp parallel for
    for (size_t i = 1; i < v.size(); i++) {
        // cum score
        get<0>(v[i]) = get<0>(cum[i + p.seed - 1]) - get<0>(cum[i - 1]);
        get<1>(v[i]) = get<1>(cum[i + p.seed - 1]) - get<1>(cum[i - 1]);
        get<2>(v[i]) = get<2>(cum[i + p.seed - 1]) - get<2>(cum[i - 1]);
        get<3>(v[i]) = get<3>(cum[i + p.seed - 1]) - get<3>(cum[i - 1]);      
        transform(gseq.substr(i, p.seed), v[i]);
    }    
}

void DetectIntersperse::mapVector(const string & seq,
        vector<tuple<unsigned, unsigned, unsigned, unsigned> > & nseq) {
    for (size_t i = 0; i < seq.size(); i++) {
        switch (seq[i]) {
            case 'A':
                get<0>(nseq[i]) = 1;
                break;
            case 'T':
                get<1>(nseq[i]) = 1;
                break;
            case 'C':
                get<2>(nseq[i]) = 1;
                break;
            case 'G':
                get<3>(nseq[i]) = 1;
                break;
            default:
                break;
        }
    }
}

void DetectIntersperse::cumsum(
        vector<tuple<unsigned, unsigned, unsigned, unsigned> > & v) {
    // score of previous position
    unsigned a = 0, t = 0, c = 0, g = 0;
    for (size_t i = 0; i < v.size(); i++) {
        a += get<0>(v[i]);
        t += get<1>(v[i]);
        c += get<2>(v[i]);
        g += get<3>(v[i]);
        v[i] = make_tuple(a, t, c, g);
    }
}

void DetectIntersperse::transform(const string & s,
        tuple<unsigned, unsigned, unsigned, unsigned> & t) {
    // filter out seq with N and low complexity seq
    if (!filterSeq(s) && !lowComplex(s)) {
        // transform scores
        if (get<0>(t) < get<1>(t)) {
            swap(get<0>(t), get<1>(t));
            swap(get<2>(t), get<3>(t));
        } else if (get<0>(t) == get<1>(t) && get<2>(t) < get<3>(t)) {
            swap(get<2>(t), get<3>(t));
        }
    } else {
        // mask
        t = make_tuple(0, 0, 0, 0);
    }
}

bool DetectIntersperse::lowComplex(const string & s) {
    double gc = gcContent(s);
    if (gc < 0.2 || gc > 0.8 || seqComplexity(s) < 0.675) {
        return true;
    }    
    return false;
}

bool DetectIntersperse::filterSeq(const string & s) {
    for (auto c : s) {
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
            return true;
        }
    }
    return false;
}

double DetectIntersperse::gcContent(const string & s) {
    double r = 0;
    for (auto i : s) {
        if (i == 'G' || i == 'C') {
            r++;
        }
    }
    return r / s.size();
}

double DetectIntersperse::seqComplexity(const string & seq) {
    unsigned int c = 1;
    string S(1, seq[0]);
    string Q = "", SQ = "";
    size_t seqLen = seq.size();
    for (size_t i = 1; i < seqLen; i++) {
        Q += seq[i];
        SQ = S + Q;
        string SQv = SQ.substr(0, SQ.size() - 1);
        if (SQv.find(Q) == string::npos) {
            S = SQ;
            Q = "";
            c++;
        }
    }
    return (c / (seqLen / (log(seqLen) / log(4))));
}

void DetectIntersperse::groupScore(unordered_map<string, unsigned> & index, 
        size_t chrom) {
    for (size_t i = 0; i < cumScore.size(); i++) {
        tuple<unsigned, unsigned, unsigned, unsigned> & t = cumScore[i];
        if (t != make_tuple(0, 0, 0, 0)) {
            string key = to_string(get<0>(t)) + ":" 
                    + to_string(get<1>(t))+ ":"
                    + to_string(get<2>(t)) + ":" 
                    + to_string(get<3>(t));
            if(index.find(key) == index.end()) {
                size_t n = groupedScore.size();
                groupedScore.resize(n + 1);
                groupedScore.back().push_back(make_pair(chrom, i));
                index[key] = n;
            } else {
                groupedScore[index[key]].push_back(make_pair(chrom, i));
            }
        }        
    }
    cumScore.clear();
    cumScore.shrink_to_fit();
}

void DetectIntersperse::findCopy(vector<pair<size_t, size_t> > & v,
        vector<group> & r) {
    if (v.size() < p.copy) {
        return;
    }
    // seq, copies
    unordered_map<string, vector<repeat> > index;
    for (size_t i = 0; i < v.size(); i++) {
        string seq1 = genomeData[v[i].first].substr(v[i].second, p.seed);
        string seq2 = reverseComp(seq1);
        repeat r;
        r.chrom = v[i].first;
        r.start = v[i].second;
        r.end = r.start + p.seed - 1;
        if (index.find(seq1) != index.end()) {
            r.strand = true;
            index[seq1].push_back(r);
        } else if (index.find(seq2) != index.end()) {
            r.strand = false;
            index[seq2].push_back(r);
        } else {
            r.strand = true;
            index[seq1].push_back(r);
        }
    }
    v.clear();
    v.shrink_to_fit();
    for (auto & i : index) {
        if (i.second.size() >= p.copy) {
            group tmp;
            tmp.con = move(i.first);
            tmp.pos = move(i.second);
            r.push_back(move(tmp));
        }
    }
}

string DetectIntersperse::reverseComp(const string & s) {
    size_t last = s.size() - 1;
    string r(s.size(), ' ');
    for (size_t i = 0; i < s.size(); i++) {
        r[i] = complement(s[last - i]);
    }
    return r;
}

char DetectIntersperse::complement(char c) {
    char r;
    switch (c) {
        case 'A':
            r = 'T';
            break;
        case 'T':
            r = 'A';
            break;
        case 'C':
            r = 'G';
            break;
        case 'G':
            r = 'C';
            break;
        default:
            r = 'N';
    }
    return r;
}

void DetectIntersperse::extendR(group & g) {
    vector<repeat> & v = g.pos;
    // undetermined base in extended seq
    unsigned n = 0;
    // extended seq
    vector <char> ext;
    while (n <= p.max_n) {
        // frequency of each base: a, t, c, g, n
        int freq[] = {0, 0, 0, 0, 0};
        // extend base one by one
        for (size_t i = 0; i < v.size(); i++) {
            const string & chromSeq = genomeData[v[i].chrom];
            char tmp;
            // reach boundary
            if ((v[i].strand && v[i].end == chromSeq.size() - 1)
                    || (!v[i].strand && v[i].start == 0)) {
                // mask
                v[i].flag = false;
                continue;
            }
            if (v[i].strand) {
                v[i].end++;
                tmp = chromSeq[v[i].end];
            } else {
                v[i].start--;
                tmp = complement(chromSeq[v[i].start]);
            }
            if (tmp == 'A') {
                freq[0]++;
            } else if (tmp == 'T') {
                freq[1]++;
            } else if (tmp == 'C') {
                freq[2]++;
            } else if (tmp == 'G') {
                freq[3]++;
            } else {
                freq[4]++;
            }
        }
        // find base with largest frequency
        char conbase = findCon(base, freq, 5, v.size() * p.min_con);
        ext.push_back(conbase);
        if (conbase == 'N') {
            n++;
        }
    }
    // trim N in the end of extended seq
    unsigned len = 0;
    for (int i = ext.size() - 1; i >= 0; i--) {
        if (ext[i] == 'N') {
            len++;
        } else {
            break;
        }
    }
    g.con += string(ext.begin(), ext.end() - len);
    for (size_t i = 0; i < v.size(); i++) {
        if (v[i].strand) {
            v[i].end -= len;
        } else {
            v[i].start += len;
        }
    }
}

void DetectIntersperse::extendL(group & g) {
    vector<repeat> & v = g.pos;
    // undetermined base in extended seq
    unsigned n = 0;
    // extended seq
    vector <char> ext;
    while (n <= p.max_n) {
        // frequency of each base: a, t, c, g, n
        int freq[] = {0, 0, 0, 0, 0};
        // extend base one by one
        for (size_t i = 0; i < v.size(); i++) {
            const string & chromSeq = genomeData[v[i].chrom];
            char tmp;
            // reach boundary
            if ((v[i].strand && v[i].start == 0) || 
                    (!v[i].strand && v[i].end == chromSeq.size() - 1)) {                
                // mask
                v[i].flag = false;
                continue;
            }
            if (v[i].strand) {
                v[i].start--;
                tmp = chromSeq[v[i].start];
            } else {
                v[i].end++;
                tmp = complement(chromSeq[v[i].end]);
            }
            if (tmp == 'A') {
                freq[0]++;
            } else if (tmp == 'T') {
                freq[1]++;
            } else if (tmp == 'C') {
                freq[2]++;
            } else if (tmp == 'G') {
                freq[3]++;
            } else {
                freq[4]++;
            }
        }
        // find base with largest frequency
        char conbase = findCon(base, freq, 5, v.size() * p.min_con);
        ext.push_back(conbase);
        if (conbase == 'N') {
            n++;
        }
    }
    // trim N in the end of extended seq
    unsigned len = 0;
    for (int i = ext.size() - 1; i >= 0; i--) {
        if (ext[i] == 'N') {
            len++;
        } else {
            break;
        }
    }
    g.con = string(ext.rbegin() + len, ext.rend()) + g.con;
    for (size_t i = 0; i < v.size(); i++) {
        if (v[i].strand) {
            v[i].start += len;         
        } else {
            v[i].end -= len;
        }
    }
}

char DetectIntersperse::findCon(char base[], int freq[], int len, double t) {
    for (int i = 0; i < len; i++) {
        if ((double) freq[i] > t) {
            return base[i];
        }
    }
    return 'N';
}

void DetectIntersperse::filter(group & g) {
    vector<repeat> & v = g.pos;
    for (size_t i = 0; i < v.size(); i++) {
        if (v[i].flag) {
            unsigned count = 0;
            string s = genomeData[v[i].chrom].substr(v[i].start, g.con.size());
            if (!v[i].strand) {
                s = reverseComp(s);
            }
            for (size_t j = 0; j < g.con.size(); j++) {
                if (g.con[j] != s[j] && g.con[j] != 'N') {
                    count++;
                }
            }
            if (count > p.max_m) {
                // mask
                v[i].flag = false;
            }
        }
    }
}

void DetectIntersperse::output() {
    unordered_map<string, unordered_set<string> > index;
    ofstream out(p.output + "/interspersed_repeat.out");
    if (!out.good()) {
        cerr << "can't open file: " << p.output << endl;
        exit(1);
    }
    // merge repeats with same consensus sequence to groups
    for (group & g : result2) {
        unordered_set<string> & s = index[g.con];
        for (repeat & r : g.pos) {
            if (r.flag) {
                string key = to_string(r.strand) + '\t' + to_string(r.chrom) 
                        + '\t' + to_string(r.start) + '\t' + to_string(r.end);
                if (s.find(key) == s.end()) {
                    s.insert(key);
                }
            }
        }
        g.con = "";
        g.pos.clear();
        g.pos.shrink_to_fit();
    }
    result2.clear();
    result2.shrink_to_fit();
    // print
    unsigned n = 0;
    for (auto & i : index) {
        if (i.second.size() >= p.copy) {
            n++;
            out << "--------------------------------------------------" << endl;
            for (auto & r : i.second) {
                istringstream is(r);            
                string situation, str;
                is >> situation;
                is >> str;
                unsigned chrom = stoul(str);
                is >> str;
                unsigned start = stoul(str);
                is >> str;
                unsigned end = stoul(str);
                char strand = '+';
                if (situation == "0") {
                    strand = '-';
                }
                out << ">" << chroms[chrom] << ":"
                        << start + 1 << ":"
                        << end + 1 << ":"
                        << strand << endl;
                if (p.format == 0) {
                    if (strand == '+') {
                        out << genomeData[chrom]
                            .substr(start, end - start + 1) << endl;
                    } else {
                        out << reverseComp(genomeData[chrom]
                            .substr(start, end - start + 1)) << endl;
                    }
                }
            }
        }
    }
    out.close();
    cout << "printed " << n << " unique groups" << endl;
}
