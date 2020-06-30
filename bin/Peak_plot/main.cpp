#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include <zlib.h>
#include <pthread.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>

#include "bgzf.h"
#include "tabix.h"
#include "knetfile.h"

using namespace std;

int len_upstream   = 1000;
int len_downstream = 1000;
int bin_size       = 20;
int body_bin       = 0;
int bp_Total = (len_upstream + len_downstream)/bin_size;
int bin_Total= (len_upstream + len_downstream)/bin_size + body_bin;

string query       = "GCA.GCC.GCT,ACG.TCG";
int threadNum      = 1;

double **Array_rat;
int **Array_cnt;
map< string,vector<string> > M_chr_pos;
vector <string> V_query;

int * Array_query_peak_i;
int peak_tot = 0;
int i_query = 0;

int read_bed(string &bed_file, map< string,vector<string> > &M_chr_pos){
    ifstream infile;
    infile.open(bed_file.c_str());
    if ( ! infile ){
        cerr << "Cannot find interval file" << bed_file << endl;
        exit(0);
    }
    
    map< string,vector<string> >::iterator I_it;
    vector<string> V_pos;
    
    string lineStr;
    while (getline(infile,lineStr,'\n')){
        if (lineStr[0] == ' ' || lineStr[0] == '\n'){
            continue;
        }
        vector<string> LineVec;
        boost::split(LineVec, lineStr, boost::is_any_of("\t\n"), boost::token_compress_on);
        string chrom = LineVec[0];
        if (chrom.substr(0, 3) != "chr"){
            continue;
        }
        I_it = M_chr_pos.find(chrom);
        if (I_it==M_chr_pos.end()){
            M_chr_pos.insert(map< string,vector<string> >::value_type(chrom, V_pos));
        }
        M_chr_pos[chrom].push_back(lineStr);
        peak_tot += 1;
    }
    infile.close();
    return peak_tot;
}

int pos_bin(string &strand, int32_t p_beg, int32_t t_beg, int32_t t_end){
    // The center of t_beg and t_end
    int32_t t_mid = (t_beg + t_end) / 2;
    
    int bin = (p_beg - (t_mid-len_upstream))/(bin_size);
    if (strand == "-"){
        bin = bp_Total - bin;
    }
//    cout << p_beg << "\t"<< t_mid << "\t" << bin << "\t" <<  (p_beg - (t_mid-len_upstream)) << endl;
    return bin;
}

void load_bin(tabix_t *t, string &region, int32_t t_beg, int32_t t_end, string &strand){
    int i, len;
    ti_iter_t iter;
    const char *s;
    const ti_conf_t *idxconf;

    if (ti_lazy_index_load(t) < 0 ) {
        fprintf(stderr,"[tabix] failed to load the index file.\n");
    }
    idxconf = ti_get_conf(t->idx);

    int32_t tid, beg, end;
    
    vector<double> V_value;
    
    int ovlp_cnt = 0;
    int bin = -1;
    if (ti_parse_region(t->idx, region.c_str(), &tid, &beg, &end) == 0) {
        iter = ti_queryi(t, tid, beg, end);
        while ((s = ti_read(t, iter, &len)) != 0) {
            string lineCons = s;
            vector<string> V_Cons;
//            cout << lineCons << endl;
            boost::split(V_Cons, lineCons, boost::is_any_of(" \t\n"), boost::token_compress_on);
            int32_t p_beg = boost::lexical_cast<int32_t>(V_Cons[1]);
            int umt = boost::lexical_cast<int>(V_Cons[3]);
            int met = boost::lexical_cast<int>(V_Cons[4]);
            int bin = pos_bin(strand, p_beg, t_beg, t_end);
            
            if (umt + met >= 1){
// point coverage >=3!!!
                double ratio = met / boost::lexical_cast<double> (umt+met);
                int idx  = Array_query_peak_i[i_query] + i_query*peak_tot;
                Array_rat[idx][bin] += ratio;
                Array_cnt[idx][bin] += 1;
            }
        }
        ti_iter_destroy(iter);
    }
}

void parse_line(string &lineStr, tabix_t *t_reg){
    vector<string> lineVec;
    boost::split(lineVec, lineStr,boost::is_any_of("\t"), boost::token_compress_on);
    string  t_chr = lineVec[0];
    int32_t t_beg = boost::lexical_cast<int32_t>(lineVec[1]);
    int32_t t_end = boost::lexical_cast<int32_t>(lineVec[2]);
    int32_t t_mid = (t_beg + t_end) / 2;
    string strand = lineVec[3];
//    string region = t_chr + ":" + lineVec[1] + "-" + lineVec[2];
    ostringstream stringStream;
    stringStream << t_chr << ":" << t_mid-len_upstream << "-" << t_mid+len_downstream;
    string region = stringStream.str();
    
    load_bin(t_reg, region, t_beg, t_end, strand);
}


void interval_matrix(int i_query, string chrom, string prefix_dir){
    string reg_file = prefix_dir + "/" + chrom + "." + V_query[i_query] + ".bed.gz";
    tabix_t *t_reg;
    cout << reg_file << endl;
    cout << (t_reg = ti_open(reg_file.c_str(), 0)) << endl;
    for(int i_chrLine = 0; i_chrLine < M_chr_pos[chrom].size(); i_chrLine++){
        int idx  = Array_query_peak_i[i_query] + i_query*peak_tot;
        Array_rat[idx] = new double[bp_Total]();
        Array_cnt[idx] = new int[bp_Total]();
        string lineStr = M_chr_pos[chrom][i_chrLine];
        parse_line(lineStr, t_reg);
        Array_query_peak_i[i_query]++;
    }
    ti_close(t_reg);
}

void output_matrix(string &prefix_out_dir){
    
    string out_file_avg = prefix_out_dir + ".ratio.avg.xls";
    ofstream f_out_avg;
    f_out_avg.open(out_file_avg.c_str());
    for (i_query=0; i_query<V_query.size(); ++i_query){
        f_out_avg << V_query[i_query];
        int peak_out_i = 0;
        string out_file_mat = prefix_out_dir + "." + V_query[i_query] + ".ratio.mat";
        
        ofstream f_out_mat;
        f_out_mat.open(out_file_mat.c_str());
            
        int *array_cnt_tot = new int[bp_Total]();
        double *array_rat_tot = new double[bp_Total]();
        
        map< string,vector<string> >::iterator I_it;
        for (I_it=M_chr_pos.begin(); I_it!=M_chr_pos.end(); ++I_it){
            string chrom = I_it->first;
            if (M_chr_pos[chrom].size() == 0){
                continue;
            }
            for (int i_pos=0; i_pos<M_chr_pos[chrom].size(); ++i_pos){
                string lineStr = M_chr_pos[chrom][i_pos];
                
                vector<string> lineVec;
                boost::split(lineVec, lineStr,boost::is_any_of("\t"), boost::token_compress_on);
                string region = lineVec[0] + ":" + lineVec[1] + "-" + lineVec[2];
                f_out_mat << region;
                for (int i_bin=0; i_bin<bp_Total; ++i_bin){
                   
                    double rat_meth = Array_rat[peak_out_i+i_query*peak_tot][i_bin];
                    int    cnt_meth = Array_cnt[peak_out_i+i_query*peak_tot][i_bin];
                    
                    if(cnt_meth == 0){
                        f_out_mat << "\tnan";
                    }
                    else{
                        double avg_meth = rat_meth / boost::lexical_cast<double>(cnt_meth);
                        array_cnt_tot[i_bin] += 1;
                        array_rat_tot[i_bin] += avg_meth;
                        f_out_mat << "\t" << avg_meth;
                    }
                }
                
                delete [] Array_rat[peak_out_i+i_query*peak_tot];
                delete [] Array_cnt[peak_out_i+i_query*peak_tot];
                    
                peak_out_i += 1;
                f_out_mat << endl;
            } 
        }
        for (int i_bin=0; i_bin<bp_Total; ++i_bin){
            double rat_meth = array_rat_tot[i_bin];
            int    cnt_meth = array_cnt_tot[i_bin];
            
            if(cnt_meth == 0){
                f_out_avg << "\tnan";
            }
            else{
                double avg_meth = rat_meth / boost::lexical_cast<double>(cnt_meth);
                f_out_avg << "\t" << avg_meth;
            }
        }
        delete [] array_cnt_tot;
        delete [] array_rat_tot;
        
        f_out_mat.close();
        f_out_avg << endl;
    }
    delete [] Array_rat;
    delete [] Array_cnt;
    f_out_avg.close();
}

void usage()
{
    cout << "awk -v OFS=\"\\t\" '{print $1,$2,$3,\"+\" }' /date/huboqiang/NORM_seq_PGC/mouse/03.NDR/distal_region/mPGC_Week135_rep3_1.query.NDR.distal.cnt10.bed | ./peak_plot /date/huboqiang/NORM_seq_PGC/mouse/02.SingleC/mPGC_Week135_rep3_1/singleC  /dev/stdin test_out /date/huboqiang/NORM_seq_PGC/mouse/04.Flank_plot/mPGC_Week135_rep3_1" << endl;
    cout << "   -s <string>query strings, default=" << len_upstream << endl;
    cout << "   -U <int>   length of upstream, default=" << len_upstream << endl;
    cout << "   -D <int>   length of downstream, default=" << len_downstream << endl;
    cout << "   -B <int>   bin size, default=" << bin_size << endl;
    cout << "   -b <int>   body bin counts, default=" << body_bin << endl;
    cout << "  -h        get help information"   << endl;
    cout << " U \% b or D \% b should better to be zero to avoid mistake edge-caculation in TSS or TES" << endl;
    exit (0);
}

int main(int argc, char *argv[])
{
    int c;
    while ( (c=getopt(argc,argv,"s:U:D:B:b:p:h")) != -1 ){
        switch(c){
            case 's' : query          = optarg;break;
            case 'U' : len_upstream   = atoi(optarg);break;
            case 'D' : len_downstream = atoi(optarg);break;
            case 'B' : bin_size       = atoi(optarg);break;
            case 'b' : body_bin       = atoi(optarg);break;
            case 'p' : threadNum      = atoi(optarg);break;
            case 'h' : usage();break;
            default : usage();
        }
    }
    if (argc < 4) usage();
    string prefix_dir = argv[optind++];
    string file_bed = argv[optind++];
    string prefix_out_dir = argv[optind++];
    
    bp_Total = (len_upstream + len_downstream)/bin_size;
    bin_Total= (len_upstream + len_downstream)/bin_size + body_bin;
    
    boost::split(V_query, query, boost::is_any_of(","), boost::token_compress_on);
    peak_tot = read_bed(file_bed, M_chr_pos);
    cout << peak_tot << endl;
    Array_rat = new double *[peak_tot*V_query.size()];
    Array_cnt = new int *[peak_tot*V_query.size()];
    Array_query_peak_i = new int [V_query.size()]();
    
    map< string,vector<string> >::iterator I_it;
    for (I_it=M_chr_pos.begin();  I_it!=M_chr_pos.end(); I_it++ ){
        string chrom = I_it->first;
        for (i_query=0; i_query<V_query.size(); ++i_query){
            interval_matrix(i_query, chrom, prefix_dir);
        }
        
    }
    
    output_matrix(prefix_out_dir);
    
    delete [] Array_query_peak_i;
}
