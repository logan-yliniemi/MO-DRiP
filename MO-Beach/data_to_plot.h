//
//  data_to_plot.h
//  MAB
//
//  Created by Alex  Jennings on 7/3/14.
//  Copyright (c) 2014 ___Robotics___. All rights reserved.

// *** Last Updated: 9/24/2014 ***

// *** Instructions: ***

//********************************************************************************************************************//
//Main:
//I have generated some value I want to use for statistical anlaysis
//I place take_value(some value), a member function, directly beneath variable generation (depending on flow control used)
//I place carriage_return() at the end of the code used for a single statistical run
//I place run_stats_library(), a member function, at the end of main
//I place push_to_file(desired filename), a member function, below run_stats_library()
//A text file with the desired filename should be created, allowing for simple plotting in MATLAB or Excel
//********************************************************************************************************************//

#ifndef MAB_data_to_plot_h
#define MAB_data_to_plot_h


using namespace std;
#include <iostream>
#include<math.h> 
#include<numeric>
#include<functional>


//A class containing statistical member functions devloped in the code below

class statistics_library{
    
//Programmer-user functions
    
public:
    void carriage_return();
    void take_value(double);
    void run_stats_library(char* filename);
    void test_stats();
    
//Hidden from the programmer-user
    
private:
    
//Vectors used in class, statistics library
    
    vector<double> single_stat_run;
    vector< vector<double> > master_table;
    
    vector<double> means;
    vector<double> medians;
    vector<double> stdevs;
    vector<double> running_means;
    vector<double> z_episode_running_means;
    
//Member functions in class, statistics library
    
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

 //Creates a vector that stores the desired value for current episode

// *** User must place this member function in the appropiate place in main ***

// *** The appropiate place to use the take_value member function in main, is after the statistically significant value is generated (The value "I" use for statistical analysis)***

void statistics_library::take_value(double action_value){
    single_stat_run.push_back (action_value);
   
}

//Creates a master table of all statistical values

//Adds vector single_stat_run to another vector as a table (a row)

void statistics_library::carriage_return(
){
    master_table.push_back(single_stat_run);
    
//Clears the take_value vector for the next run
    
    single_stat_run.clear();
}

//Calculates means

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
        
//Stores means in a vector
        
        means.push_back(mean);
    }
    
}

//Calculates medians

void statistics_library::calc_medians(){
    int nmbr_runs = master_table.size();
    int nmbr_iterations = master_table.at(0).size();
    int x=nmbr_runs/2;
    int w=(nmbr_runs+1)/2;
    double median;
    vector <double> columns;
    
//For even amount of numbers in a set
    
    if(nmbr_runs%2 == 0){
        
        for (int k=0; k<nmbr_iterations; k++){
            median=0;
            columns.clear();
            for ( int n=0; n<nmbr_runs; n++)
            {
                columns.push_back(master_table.at(n).at(k));
            }
            sort(columns.begin(),columns.end());
            median=(columns.at(x)+columns.at(x+1))/2;
            
//Stores median values in a vector
            
            medians.push_back(median);
            }
        }

//For odd amount of numbers in a set
    
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
            
//Stores median values in a vector
            
            medians.push_back(median);
            
        }
    }
        
            }

//Calculates standard deviations

// This is the url for the standard deviation formula used for this calculation:
//www.mathsisfun.com/data/standard-deviation-formulas.html


void statistics_library::calc_stdevs(){
    int nmbr_runs = master_table.size();
    int nmbr_iterations = master_table.at(0).size();
    
   
    vector<double> stdevs_comp;
    double N;
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
        
//Stores standard deviation values in a vector
        
        stdevs.push_back(stdevs_total);
    }
}

    
//Calculates running means

void statistics_library::calc_running_means(){
    int nmbr_runs = master_table.size();
    int nmbr_iterations = master_table.at(0).size();
        double c=0;
        double running_avg;
        vector <double> columns;
        double N=0;
    
    
    for (int k=0; k<nmbr_iterations; k++){
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
//Stores running mean values in a vector
        
        running_means.push_back(running_avg);
            }
}

//Calculates epsiode-number-running means

void statistics_library:: calculate_z_episode_running_means(int z=10){
    int nmbr_runs = master_table.size();
    int nmbr_iterations = master_table.at(0).size();
    int a =0;
    double c=0;
    double running_avg=0;
    vector <double> columns;
    double N=0;
    double M=0;
    
    
    for(int i=0; i<means.size(); i++){
        if(i==0){
            
//Stores running means in a vector
            
            z_episode_running_means.push_back(means.at(0));
        }
        
//Calculates running means component
        
        if(0<i && i<z){
            double rm = accumulate(means.begin(),means.begin()+i,0);
            rm/=i;
            
//Stores running means in a vector
            
            z_episode_running_means.push_back(rm);
        }
//Calculates z-ep running means component
        
        if(i>=z){
            double rm = accumulate(means.begin()+i-z,means.begin()+i,0);
            rm/=z;
            
//Stores running means in a vector
            
            z_episode_running_means.push_back(rm);
        }
    }
    
    }


//Calculates all statistical values within the class

void statistics_library:: calculate_all_statistics(){
    
    calc_means();
    calc_medians();
    calc_stdevs();
    calc_running_means();
    calculate_z_episode_running_means();
}

//Clears all vectors except for the master table

void statistics_library::prep(){
    
    single_stat_run.clear();
    means.clear();
    medians.clear();
    stdevs.clear();
    running_means.clear();
}

//Clears all vectors including the master table

void statistics_library::reset(){
    
    single_stat_run.clear();
    means.clear();
    medians.clear();
    stdevs.clear();
    running_means.clear();
    master_table.clear();
}

//Creates text files for plotting by using raw data from statisical value vectors

// *** User must imput the desired filename of the text file ***

// *** Must calculate all other statistics first ***

// *** Values from the master table are pushed to a text file ***

void statistics_library::push_to_file(char* filename){
    int nmbr_iterations = master_table.at(0).size();
    pFile = fopen(filename,"wt");
    for(int i=0; i<nmbr_iterations; i++){
        fprintf (pFile, "%.4f\t", means.at(i));
    }
    fprintf (pFile, "\n");
    for(int i=0; i<nmbr_iterations; i++){
        fprintf (pFile, "%.4f\t", medians.at(i));
    }
    fprintf (pFile, "\n");
    for(int i=0; i<nmbr_iterations; i++){
        fprintf (pFile, "%.4f\t", stdevs.at(i));
    }
    fprintf (pFile, "\n");
    for(int i=0; i<nmbr_iterations; i++){
        fprintf (pFile, "%.4f\t", running_means.at(i));
    }
    fprintf (pFile, "\n");
    for(int i=0; i<nmbr_iterations; i++){
        fprintf (pFile, "%.4f\t", z_episode_running_means.at(i));
    }
    fclose (pFile);
    
}

//Calculates all statisitcs and creates a text file of the statistical values for simple plotting

// *** User must appropiatel place take_value() function in main to perform statiscal operations and store statistical values

// *** User must imput a desired name for the text file when using this member function in main ***

void statistics_library::run_stats_library(char* filename){
    prep();
    calculate_all_statistics();
    push_to_file(filename);
    
//User can clear master table (all stored statistical values for text file creation) by uncommenting the reset member function below
    
    // reset();
}

    



#endif



