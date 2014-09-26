//
//  data_to_plot.h
//  MAB
//
//  Created by Alex  Jennings on 7/3/14.
//  Copyright (c) 2014 ___Robotics___. All rights reserved.
//

#ifndef MAB_data_to_plot_h
#define MAB_data_to_plot_h

/// DOCUMENTATION:
/// @AJ

using namespace std;
#include <iostream>
#include<math.h> 
#include<numeric>
#include<functional>



class statistics_library{ // a class containing statistical member functions devloped in the code below
    /// Programmer-user functions
public:
    void carriage_return();
    void take_value(double);
    void run_stats_library(char* filename);
    void test_stats();
    
    /// Hidden from the programmer-user
private: 
    vector<double> single_stat_run;
    vector< vector<double> > master_table;
    
    vector<double> means;
    vector<double> medians;
    vector<double> stdevs;
    vector<double> running_means;
    vector<double> z_episode_running_means;
    void prep();
    void reset();
    
    void calc_means();
    void calc_medians();
    void calc_stdevs();
    void calc_running_means();
    void calculate_z_episode_running_means(int);
    void calculate_all_statistics();
    

    FILE* pFile;
    void push_to_file(char* filename);
};

void statistics_library::take_value(double action_value){
    single_stat_run.push_back (action_value); //creates a vector that takes the action value for current episode
    /// @AJ
}

void statistics_library::carriage_return(
){
    master_table.push_back(single_stat_run); //adds vector take_value to another vector as a table (a row)
    single_stat_run.clear(); // clears the take_value vector for the next run
    /// @AJ
}


void statistics_library::calc_means(){
    int nmbr_runs = master_table.size();
    int nmbr_iterations = master_table.at(0).size();
    
    double total;
    double mean;
    vector <double> columns;
    
    for (int k=0; k<nmbr_iterations; k++){
        mean=0;
        total=0;
        columns.clear();
        for ( int n=0; n<nmbr_runs; n++)
        {
            columns.push_back(master_table.at(n).at(k));
        }
        
        for (int n=0; n<nmbr_runs; n++){
            total+=columns.at(n);
        }
        mean=total/nmbr_runs;
        means.push_back(mean);
    }
    
      /// AJ
}

void statistics_library::calc_medians(){
    int nmbr_runs = master_table.size();
    int nmbr_iterations = master_table.at(0).size();
    int x=nmbr_runs/2;
    int w=(nmbr_runs+1)/2;
    double median;
    vector <double> columns;
    
    if(nmbr_runs%2 == 0){ /// @AJ
        
        for (int k=0; k<nmbr_iterations; k++){
            median=0;
            columns.clear();
            for ( int n=0; n<nmbr_runs; n++)
            {
                columns.push_back(master_table.at(n).at(k));
            }
            sort(columns.begin(),columns.end());
            median=(columns.at(x)+columns.at(x+1))/2;
            medians.push_back(median);
            //cout << "THIS IS A FLAG #349875" << endl;
            }
        }

    
    if(nmbr_runs%2 == 1){
        
        for (int k=0; k<nmbr_iterations; k++){
            median =0;
            columns.clear();
            for ( int n=0; n<nmbr_runs; n++)
            {
                columns.push_back(master_table.at(n).at(k));
            }
            sort(columns.begin(),columns.end());
            median=(columns.at(w));
            medians.push_back(median);
        }
    }
        
            }


void statistics_library::calc_stdevs(){
    int nmbr_runs = master_table.size();
    int nmbr_iterations = master_table.at(0).size();
    
    
    vector<double> stdevs_comp;
    double N;                   // defines variables and vectors for intermediate calculations
    double u=0;
    double v;
    double y;
    double sum;
    double stdevs_total;
    vector <double> columns;
    
    for (int k=0; k<nmbr_iterations; k++){
        stdevs_comp.clear();
        y=0;
        u=0;
        v=0;
        sum=0;
        N=0;
        stdevs_total=0;
        columns.clear();
        u=means.at(k);
        for ( int n=0; n<nmbr_runs; n++)
        {
            columns.push_back(master_table.at(n).at(k));
            v=columns.at(n);
            y=pow(v-u,2);
            stdevs_comp.push_back(y);
        }
        sum=accumulate(stdevs_comp.begin(), stdevs_comp.end(),0);
        //cout<< sum << endl;
        N=nmbr_runs;
        stdevs_total=pow(sum/N,1.0/2.0);
        //cout<< stdevs_total << " " << N << endl;
        stdevs.push_back(stdevs_total);
    }
}

    
    // This is the url for the standard deviation formula used for this calculation:
        //www.mathsisfun.com/data/standard-deviation-formulas.html
    


void statistics_library::calc_running_means(){
    int nmbr_runs = master_table.size();
    int nmbr_iterations = master_table.at(0).size();
        double c=0;
        double running_avg;
        vector <double> columns;
        double N=0;
    
    
    for (int k=0; k<nmbr_iterations; k++){
        //columns.clear();
        c=0;
        N=0;
        running_avg=0;
        for ( int n=0; n<nmbr_runs; n++)
        {
            columns.push_back(master_table.at(n).at(k));
        }
        for (int n=0; n<columns.size(); n++)
        {
            c++;
            N+=columns.at(n);
            running_avg=N/c;
        }
        running_means.push_back(running_avg);
            }
}

void statistics_library:: calculate_z_episode_running_means(int z=10){
    int nmbr_runs = master_table.size();
    int nmbr_iterations = master_table.at(0).size();
    int a =0;
    double c=0;
    double running_avg=0;
    vector <double> columns;
    double N=0;
    double M=0;
    
    // calculations based off of mean calculations.
    for(int i=0; i<means.size(); i++){
        if(i==0){
            z_episode_running_means.push_back(means.at(0));
        }
        if(0<i && i<z){
            // running mean
            double rm = accumulate(means.begin(),means.begin()+i,0);
            rm/=i;
            z_episode_running_means.push_back(rm);
        }
        if(i>=z){
            // z-ep running mean
            double rm = accumulate(means.begin()+i-z,means.begin()+i,0);
            rm/=z;
            z_episode_running_means.push_back(rm);
        }
    }
    
    }



    /// @AJ running means with a certain length of memory. The one from above does 1:nmbr_iterations. This one should do (nmbr_iterations-n):nmbr_iterations


void statistics_library:: calculate_all_statistics(){
    calc_means();
    calc_medians();
    calc_stdevs();
    calc_running_means();
    calculate_z_episode_running_means();
    
    
    /// @AJ have this function execute all the other statistics functions!
}

void statistics_library::prep(){
    
    single_stat_run.clear();
    means.clear();
    medians.clear();
    stdevs.clear();
    running_means.clear();
    
    
    /// @AJ Clears all the vectors associated with this class, but not the master table.
}

void statistics_library::reset(){
    
    single_stat_run.clear();
    means.clear();
    medians.clear();
    stdevs.clear();
    running_means.clear();
    master_table.clear();
    
    /// @AJ Clears all the vectors associated with this class, including the master table.
}


void statistics_library::push_to_file(char* filename){
    /// Run function calculate_all_statistics() first.
    int nmbr_iterations = master_table.at(0).size();
    pFile = fopen(filename,"wt");
    for(int i=0; i<nmbr_iterations; i++){
        fprintf (pFile, "%.4f\t", means.at(i));
    }
    fprintf (pFile, "\b \b\n");
    for(int i=0; i<nmbr_iterations; i++){
        fprintf (pFile, "%.4f\t", medians.at(i));
    }
    fprintf (pFile, "\b \b\n");
    for(int i=0; i<nmbr_iterations; i++){
        fprintf (pFile, "%.4f\t\t", stdevs.at(i));
    }
    fprintf (pFile, "\b \b\n");
    for(int i=0; i<nmbr_iterations; i++){
        fprintf (pFile, "%.4f\t", running_means.at(i));
    }
    fprintf (pFile, "\b \b\n");
    for(int i=0; i<nmbr_iterations; i++){
        fprintf (pFile, "%.4f\t", z_episode_running_means.at(i));
    }
    fclose (pFile);
    /// @AJ
}


void statistics_library::run_stats_library(char* filename){
    prep();
    calculate_all_statistics();
    push_to_file(filename);
    // reset();
}

    



#endif



