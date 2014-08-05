#include <cstdlib>
#include <math.h>
#include <time.h>

#include <stdio.h>
#include <iostream>
using namespace std;

#ifndef VECTOR_INCLUDE
#define VECTOR_INCLUDE
#include <vector>
#endif

#define STAT_RUNS 30
#define EPISODES 1000
#define STEPS 10
#define NUM_AGENTS 700
#define BEACHES 10
#define CAPACITY 15
#define FAMILY_CAPACITY 35
#define ACTIONS 3

#define LYRAND (double)rand()/RAND_MAX
#define SMALL 0.0001
#define HUNDI_ROUND(x) (double)floor(x*100)/100

/// Normalization Limits (calculated in determine_limits():
double D1MIN=999999;
double D2MIN=999999;
double D1MAX=-999999;
double D2MAX=-999999;
double G1MAX;
double G2MAX;
double L1MAX;
double L2MAX;

#include "QLearner.h"

bool pretty_print = true;

/// All ``command" booleans initialized to false because they are cycled in main()
bool command_global = false;
bool command_difference = false;
bool command_global_PBRS_gzmi = false;
bool command_global_PBRS_hand = false;
bool command_difference_PBRS_hand = false;
bool command_local = false;

using namespace std;

double MO_Combine(double obj_1, double obj_2, double beta){
    return obj_1*beta + (obj_2*(1-beta));
}

double normalize(double value, double min, double max){
    return (value - min) / (max-min);
}

class beach {
public:
    void start();
    vector<double> beach_lookup;
    vector<double> beach_family_lookup;
    void make_capacity_lookup();
    void make_family_capacity_lookup();
    void determine_limits();
    
    void get_attendance(vector<QLearner>*);
    
    void evaluate();
    void console_attendance();

    vector<int> attendance;
    vector<double> sect_local;
    vector<double> sect_difference;
    vector<double> sect_family_local;
    vector<double> sect_family_difference;
    double global;
    double family_global;
};

void beach::console_attendance() /// Put attendance values out to console.
{
    for (int i = 0; i < attendance.size(); i++) {
        cout << attendance.at(i);
        cout << "\t";
    }
    cout << endl;
}

void beach::make_capacity_lookup() /// Creates "L" Lookup table, so that we can skip exponential calculations.
{
    double max_reward = 10.0;
    beach_lookup.clear();
    //cout << "BAR LOOKUP CHART: " << endl;
    for (int x = 0; x < NUM_AGENTS + 1; x++) {
        double xx = x;
        beach_lookup.push_back(xx * exp(-xx / CAPACITY));
        if(x>0){
            double diff = beach_lookup.at(x) - beach_lookup.at(x-1);
            if(diff > D1MAX){D1MAX = diff;}
            if(diff < D1MIN){D1MIN = diff;}
        }
//        if (x < CAPACITY){
//            beach_lookup.push_back(max_reward);
//        } else {
//            beach_lookup.push_back(max_reward/pow(x + 1 - CAPACITY,2));
//        }
        //cout << "If attendance = "<< x << ", reward = " << beach_lookup.back() << "\n";
    }
}

void beach::make_family_capacity_lookup()
{
    double max_reward = 10.0;
    beach_family_lookup.clear();
    //cout << "BAR LOOKUP CHART: " << endl;
    for (int x = 0; x < NUM_AGENTS + 1; x++) {
        double xx = x;
        beach_family_lookup.push_back(xx * exp(-xx / FAMILY_CAPACITY));
//        if (x < FAMILY_CAPACITY){
//            beach_lookup.push_back(max_reward);
//        } else {
//            beach_lookup.push_back(max_reward/pow(x + 1 - CAPACITY,2));
//        }
        //cout << "If attendance = "<< x << ", reward = " << beach_family_lookup.back() << "\n";
    }
}

void beach::determine_limits(){
        for (int x = 1; x < NUM_AGENTS + 1; x++) {
            double diff = beach_family_lookup.at(x) - beach_family_lookup.at(x-1);
            if(diff > D2MAX){D2MAX = diff;}
            if(diff < D2MIN){D2MIN = diff;}
            double diff2 = beach_lookup.at(x) - beach_lookup.at(x-1);
            if(diff2 > D1MAX){D1MAX = diff;}
            if(diff2 < D1MIN){D1MIN = diff;}
        }
        
        double c = CAPACITY;
        double fc = FAMILY_CAPACITY;
        double glut = NUM_AGENTS - c*BEACHES;
        double family_glut = NUM_AGENTS - fc * BEACHES;
        
        G1MAX = c*exp(-c/c) * (BEACHES-1) + glut*exp(-glut/c);
        G2MAX = fc*exp(-fc/fc) * (BEACHES-1) + family_glut*exp(-family_glut/fc);
        
        L1MAX = c*exp(-c/c);
        L2MAX = fc*exp(-fc/fc);
}

void beach::start() { /// Initializes values.
    sect_local.clear();
    sect_difference.clear();
    attendance.clear();
    global = 0;
}

void beach::get_attendance(vector<QLearner>* pA) { /// Determines agent attendance per beach
    attendance.clear();
    attendance.resize(BEACHES, 0);
    for (int i = 0; i < NUM_AGENTS; i++) {
        attendance.at(pA->at(i).state)++;
    }
}

void beach::evaluate() { /// Evaluates local, global, difference rewards
    for (int beach = 0; beach < attendance.size(); beach++) {
        sect_local.push_back(beach_lookup.at(attendance.at(beach)));
        sect_family_local.push_back(beach_family_lookup.at(attendance.at(beach)));
        global += sect_local.back();
        family_global += sect_family_local.back();
        if (attendance.at(beach) == 0) {
            sect_difference.push_back(0);
            sect_family_difference.push_back(0);
        } else {
            sect_difference.push_back(/**/ beach_lookup.at(attendance.at(beach)) - beach_lookup.at(attendance.at(beach) - 1) /**/);
            sect_family_difference.push_back(beach_family_lookup.at(attendance.at(beach)) - beach_family_lookup.at(attendance.at(beach) - 1));
        }
    }
}

void sense(vector<QLearner>* pA) {
    for (int agent = 0; agent < NUM_AGENTS; agent++) {
        pA->at(agent).previousState = pA->at(agent).state;
    }
}

void decide(vector<QLearner>* pA) { ///agents decide which day they are attending.
    for (int agent = 0; agent < NUM_AGENTS; agent++) {
        pA->at(agent).decay_alpha();
        pA->at(agent).decay_epsilon();
        pA->at(agent).choose_egreedy_action();
    }
}

void act(vector<QLearner>* pA, beach* pE) { /// agents attend, beach grabs attendance
    int chosen_agent = LYRAND*NUM_AGENTS; //Choose a random agent
    pA->at(chosen_agent).randomise_learnability_action();

    
    for (int agent = 0; agent < NUM_AGENTS; agent++) {
        //Calculate next state for all agents
        pA->at(agent).state += pA->at(agent).action - 1;
        if (pA->at(agent).state < 0) {
            pA->at(agent).state = 0;
        } else if (pA->at(agent).state >= BEACHES) {
            pA->at(agent).state = BEACHES-1;
        }
        
        //Calculate next noise state
        if (agent == chosen_agent){
            pA->at(agent).noise_state = pA->at(agent).state;
        } else {
            pA->at(agent).noise_state += pA->at(agent).learnability_action - 1;
        }
        if (pA->at(agent).noise_state < 0) {
            pA->at(agent).noise_state = 0;
        } else if (pA->at(agent).noise_state >= BEACHES) {
            pA->at(agent).noise_state = BEACHES-1;
        }
                
        //Calculate next signal state
        if (agent == chosen_agent){
            pA->at(agent).signal_state += pA->at(agent).learnability_action - 1;
        } else {
            pA->at(agent).signal_state = pA->at(agent).state;
        }
        if (pA->at(agent).signal_state < 0) {
            pA->at(agent).signal_state = 0;
        } else if (pA->at(agent).signal_state >= BEACHES) {
            pA->at(agent).signal_state = BEACHES-1;
        }
    }
    
    pE->start();
    pE->get_attendance(pA);
}

void react(vector<QLearner>* pA, beach* pE) { /// reward calculations, Q updates.
    pE->evaluate();
    
    int bta=0;
    int num_assigned_to=0;
    
    for (int agent = 0; agent < NUM_AGENTS; agent++) {
        double L1,L2;
        double G1,G2;
        double D1,D2;
        double GZMI1,GZMI2;

        L1 = pE->sect_local.at(pA->at(agent).state);
        G1 = pE->global;
        D1 = pE->sect_difference.at(pA->at(agent).state);
        
        L2 = pE->sect_family_local.at(pA->at(agent).state);
        G2 = pE->family_global;
        D2 = pE->sect_family_difference.at(pA->at(agent).state);

        GZMI1 = G1 - D1; /// TODO THIS IS A HACK
        GZMI2 = G2 - D2;
        
        if(command_difference_PBRS_hand || command_global_PBRS_hand)
        {
            //Manual heuristic encouraging agents maximise capacity on all but one beach 
            // TODO: Generalise - this version is specific to 5 beaches
            //int num_excess_agents = NUM_AGENTS-(CAPACITY*(BEACHES-1));
            //int encouraged_state = 2; //All agents to the middle beach except... 
            //if (0 <= agent && agent < CAPACITY) {encouraged_state = 0;}
            //else if (CAPACITY <= agent && agent < CAPACITY*2) {encouraged_state = 1;}
            //else if (NUM_AGENTS/2 <= agent && agent < CAPACITY + NUM_AGENTS/2) {encouraged_state = 3;}
            //else if (CAPACITY + NUM_AGENTS/2 <= agent && agent < CAPACITY * 2 + NUM_AGENTS/2) {encouraged_state = 4;}
            
            //Manual heuristic for MOB problem, generalized to BEACHES sections.
            int num_excess_agents = NUM_AGENTS - (CAPACITY*(BEACHES-1));
            num_assigned_to++;
            if(num_assigned_to == CAPACITY -1){
                num_assigned_to=0;
                bta++;
            }
            if(bta >= BEACHES){bta = BEACHES-1;}
            int encouraged_state = bta;
            
            
            //Encourage even spread
//            int encouraged_state = 2; //All agents to the middle beach except... 
//            if (0 <= agent && agent < NUM_AGENTS/5) {encouraged_state = 0;}
//            else if (NUM_AGENTS/5 <= agent && agent <  2*NUM_AGENTS/5) {encouraged_state = 1;}
//            else if (3* NUM_AGENTS/5 <= agent && agent < 4* NUM_AGENTS/5) {encouraged_state = 3;}
//            else if (4* NUM_AGENTS/5 <= agent && agent < NUM_AGENTS) {encouraged_state = 4;}
            
            //Encourage staying still
//            int encouraged_state = 3; 
//            if (0 <= agent && agent < NUM_AGENTS/2) {encouraged_state = 1;}
            
            //Encourage everyone to super beach
            //int encouraged_state = 2; //All agents to the middle beach 
            
            //Debug message
            //cout << "Agent " << agent << " encouraged to go on day " << encouraged_state << endl; 

            //For potential gradient 
            pA->at(agent).previousPhi = (BEACHES - abs(pA->at(agent).previousState - encouraged_state))*100;
            pA->at(agent).currentPhi = (BEACHES - abs(pA->at(agent).state - encouraged_state))*100;

            //For discrete potential
//            if (pA->at(agent).previousState == encouraged_state){   //For setting initial state potential
//                pA->at(agent).previousPhi = 10;                   //Works ok in theory as this is a static potential function
//            } else {
//                pA->at(agent).previousPhi = 0;
//            }         
//            if (pA->at(agent).state == encouraged_state){
//                pA->at(agent).currentPhi = 10;        
//            } else {
//                pA->at(agent).currentPhi = 0;
//            }
            
            //Dynamic Manual Negative If attendance is between CAPACITY and 2*CAPACITY
//            int attendanceAtCurrentState = pE->attendance.at(pA->at(agent).state);
//            if (CAPACITY < attendanceAtCurrentState && attendanceAtCurrentState < 2*CAPACITY) {
//                pA->at(agent).currentPhi = -10;
//            } else {
//                pA->at(agent).currentPhi = 0;
//            }
            
            //Dynamic Manual Positive If attendance is NOT between CAPACITY and 2*CAPACITY
//            int attendanceAtCurrentState = pE->attendance.at(pA->at(agent).state);
//            if (CAPACITY < attendanceAtCurrentState && attendanceAtCurrentState < 2*CAPACITY) {
//                pA->at(agent).currentPhi = 0;
//            } else {
//                pA->at(agent).currentPhi = 10;
//            }
            
            //Dynamic Manual encourage attendance is between CAPACITY and 2*CAPACITY
            //int attendanceAtCurrentState = pE->attendance.at(pA->at(agent).state);
            //if (CAPACITY < attendanceAtCurrentState && attendanceAtCurrentState < 2*CAPACITY) {
            //    pA->at(agent).currentPhi = 10;
            //} else {
            //    pA->at(agent).currentPhi = 0;
            //}
        }
                       
        //Automated multi-agent potential function        
        if(command_global_PBRS_gzmi){
            pA->at(agent).currentPhi = GZMI1;   //Potential-based difference reward
        }
                
        //Calculate potential based reward
        double PBRS = pA->at(agent).gamma * pA->at(agent).currentPhi - pA->at(agent).previousPhi;   
        double shapedReward=0;
        if(command_global_PBRS_gzmi || command_global_PBRS_hand)
        {
            shapedReward = MO_Combine(normalize(G1,0,G1MAX),normalize(G2,0,G2MAX),0.5) + PBRS;
        }
        if(command_difference_PBRS_hand)
        {
            shapedReward = MO_Combine(normalize(D1,D1MIN,D1MAX),normalize(D2,D2MIN,D2MAX),0.5) + PBRS;
        }
        
//        //Debug output for PBRS
//        if (agent == 20 
//            //&& pA->at(agent).currentPhi != pA->at(agent).previousPhi
//            //&& pA->at(agent).state == pA->at(agent).previousState
//           ) 
//        {
//            cout << "Current state : " << pA->at(agent).state << "\t Previous state : " << pA->at(agent).previousState << endl;
//            cout << PBRS << "Equals " << pA->at(agent).gamma << " * " << pA->at(agent).currentPhi << " - " << pA->at(agent).previousPhi << endl;
//        }

        pA->at(agent).set_local(MO_Combine(normalize(L1,0,L1MAX),normalize(L2,0,L2MAX),0.5));
        pA->at(agent).set_global(MO_Combine(normalize(G1,0,G1MAX),normalize(G2,0,G2MAX),0.5));
        pA->at(agent).set_difference(MO_Combine(normalize(D1,D1MIN,D1MAX),normalize(D2,D2MIN,D2MAX),0.5));
        pA->at(agent).set_gzmi(MO_Combine(normalize(GZMI1,0,L1MAX),normalize(GZMI2,0,L2MAX),0.5));
        pA->at(agent).set_shaped_reward(shapedReward);

        if(command_local){pA->at(agent).learn_with_local();}
        if(command_global){pA->at(agent).learn_with_global();}
        if(command_difference){pA->at(agent).learn_with_difference();}
        if(command_difference_PBRS_hand || command_global_PBRS_gzmi || command_global_PBRS_hand)
        {pA->at(agent).learn_with_shaped_reward();}
        
        
        pA->at(agent).Qupdate();  
    }
}

void report(FILE* pFILE, double global) { /// report to text file
    fprintf(pFILE, "%.5f\t", global);
}

int main() {
    srand(time(NULL));    
    
    FILE* pFILE;
    FILE* pFILE2;

    for(int method=0; method < 5; method++)
    {
        if(method==0)
        {
            cout << "HERE BEGINS GLOBAL REWARDS" << endl;
            pFILE = fopen("global.txt", "w");
            pFILE2 = fopen("global2.txt", "w");
            command_global = true;
            command_difference=false;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = false;
            command_local = false;
        }
        if(method==1)
        {
            cout << "HERE BEGINS DIFFERENCE REWARDS" << endl;
            pFILE = fopen("difference.txt", "w");
            pFILE2 = fopen("difference2.txt", "w");
            command_global = false;
            command_difference=true;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = false;
            command_local = false;
        }
        if(method==2){
            cout << "HERE BEGINS GLOBAL + PBRS (AUTO) REWARDS" << endl;
            pFILE = fopen("global_PBRS_gzmi.txt", "w");
            pFILE2 = fopen("global_PBRS_gzmi2.txt", "w");
            command_global = false;
            command_difference=false;
            command_global_PBRS_gzmi = true;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = false;
            command_local = false;
        }
        if(method==3){
            cout << "HERE BEGINS GLOBAL + PBRS (HAND) REWARDS" << endl;
            pFILE = fopen("global_PBRS_hand.txt", "w");
            pFILE2 = fopen("global_PBRS_hand2.txt", "w");
            command_global = false;
            command_difference=false;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = true;
            command_difference_PBRS_hand = false;
            command_local = false;
        }
        if(method==4){
            cout << "HERE BEGINS DIFFERENCE + PBRS (HAND) REWARDS" << endl;
            pFILE = fopen("difference_pbrs_hand.txt", "w");
            pFILE2 = fopen("difference_pbrs_hand2.txt", "w");
            command_global = false;
            command_difference=false;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = true;
            command_local = false;
        }
        if(method==5){
            cout << "HERE BEGINS LOCAL REWARDS" << endl;
            pFILE = fopen("local.txt", "w");
            pFILE2 = fopen("local2.txt", "w");
            command_global = false;
            command_difference=false;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = false;
            command_local = true;
        }

        
        for(int stat_run=0; stat_run < STAT_RUNS; stat_run++) {
            beach Pebble;
            Pebble.make_capacity_lookup();
            Pebble.make_family_capacity_lookup();
            Pebble.determine_limits();
            beach* pE = &Pebble;

            vector<QLearner> Agents;
            vector<QLearner>* pA = &Agents;

            for (int i = 0; i < NUM_AGENTS; i++) {
                QLearner Q;
                Q.id = i;
                Q.start();
                pA->push_back(Q);
            }

            for (int episode = 0; episode < EPISODES; episode++){
                if (episode % (EPISODES / 10) == 0) {
                    cout << "Run No." << stat_run << " is " << (double) episode / EPISODES * 100 << " % Complete!" << endl;
                }
                //cout << "\n---\n";
                for (int time = 0; time < STEPS; time++) {
    //                if (time % (STEPS / 10) == 0) {
    //                    cout << "Episode No." << episode << " is " << (double) time / STEPS * 100 << " % Complete!" << endl;
    //                }

                    //cout << "begin time " << time << endl;
                    sense(pA);
                    //cout << "sensed time " << time << endl;
                    decide(pA);
                    //cout << "decided time " << time << endl;
                    act(pA, pE);
                    //pE->console_attendance();         //Print behaviour to console
                    //cout << "acted time " << time << endl;
                    react(pA, pE);
                    //cout << "reacted time " << time << endl;
    //                if (time % (STEPS / 1000) == 0) {
    //                    report(pFILE, pE->global);      //Report only occasionally
    //                    //cout << pE->global;
    //                }  
                    //fprintf(pFILE,"\n");
                    //cout << "end time " << time << endl;
               }

               //For beautiful graphs
               if (pretty_print) {
                   report(pFILE, pE->global); // Report every result
                   report(pFILE2, pE->family_global); // Report every result
               } else {
                   //For Coarse Results
                   if (episode % (EPISODES / 1000) == 0) {
                        report(pFILE, pE->global);      //Report only occasionally
                        report(pFILE2, pE->family_global);      //Report only occasionally
                        //cout << pE->global;
                   }
               }

               // Pick an action from the final state
               // This action will cause a transition to the absorbing state
               sense(pA);
               decide(pA);
               for (int i = 0; i < NUM_AGENTS; i++) {
                   pA->at(i).final_Qupdate();       //Updates the final state->absorbing state transition
                   pA->at(i).restart();
               }
            }        

            //Start a new line in output file for next run
            fprintf(pFILE,"\n");
            fprintf(pFILE2,"\n");

            //Debug output
    //        cout << "\n---\nAgent 1s Qtable: Should go to beach 0" << endl;
    //        pA->at(1).console_2vector(pA->at(1).Qtable); 
    //        cout << "\n---\nAgent 8s Qtable: Should stay in beach 1" << endl;
    //        pA->at(8).console_2vector(pA->at(8).Qtable); 
    //        cout << "\n---\nAgent 20s Qtable: Should go to beach 2" << endl;
    //        pA->at(1).console_2vector(pA->at(20).Qtable); 
    //        cout << "\n---\nAgent 51s Qtable: Should stay in beach 3" << endl;
    //        pA->at(8).console_2vector(pA->at(51).Qtable); 
    //        cout << "\n---\nAgent 58s Qtable: Should go to beach 4" << endl;
    //        pA->at(1).console_2vector(pA->at(58).Qtable); 

            cout << endl << "Lane attendance:" << endl;
            pE->console_attendance();       //Print final behaviour
            cout << endl << "Final performance = " << pE->global << endl << endl;     //Final global reward  
        }
        fclose(pFILE);
        fclose(pFILE2);
    }
    return 0;
}

