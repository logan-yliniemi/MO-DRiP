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

#define EPISODES 2500
#define STEPS 50
#define STAT_RUNS 30

#define NUM_AGENTS 6
#define ACTIONS 5
#define NUM_POIS 4

#define LYRAND (double)rand()/RAND_MAX
#define SMALL 0.0001
#define HUNDI_ROUND(x) (double)floor(x*100)/100


#define XMAX 10
#define YMAX 10
#define XSTATES 10
#define YSTATES 10
#define TOTAL_STATES XSTATES * YSTATES

#define VELOCITY 1.0

#include "QLearner_Grid.h"

bool pretty_print = false;

bool command_global = false;
bool command_difference=false;
bool command_global_PBRS_gzmi = false;
bool command_global_PBRS_hand = false;
bool command_difference_PBRS_hand = false;
bool command_difference_PBRS_auto = false;
bool command_local = false;
bool command_local_PBRS_gzmi = false;
bool command_local_PBRS_hand = false;

using namespace std;

class gridworld {
public:
    void start();
    int start_pois();
    void start_epsiode();
    
    vector<double> gridworld_lookup;
    void make_lookup();
    double getReward(double distanceToPOI);

    vector<double> rover_x;
    vector<double> rover_y;
    void obtain_locations(vector<QLearner>*);
    vector<double> rover_local;
    vector<double> rover_difference;
    vector<double> rover_gzmi;
    double global;

    void console_locations();
    double rover_poi_dist(int,int);
    int find_closest_rover_index(int poidex);
    int find_closest_rover_index_not(int exclude, int poidex);

    void evaluate();
    
    vector<double> POI_values;
    vector<double> POI_X;
    vector<double> POI_Y;
};

void gridworld::make_lookup() /// Creates "L" Lookup table, so that we can skip exponential calculations.
{
    for(int i=0; i<(XMAX+YMAX); i++)
    {
        gridworld_lookup.push_back(getReward(i));
    }
};

double gridworld::getReward(double distanceToPOI) 
{
        int poimin = 2;
        int poimax = XMAX; 
        if(distanceToPOI<poimin){
            return 1.0;
        }
        else if(poimin<=distanceToPOI && distanceToPOI<poimax){
            return 1.0/(distanceToPOI*distanceToPOI);
        }
        else { //if(poimax<=distanceToPOI)
            return 0.0;
        }
};  
    
int gridworld::start_pois(){
    POI_X.clear();
    POI_Y.clear();
    POI_values.clear();
    
    
    POI_X.resize(NUM_POIS,0);
    POI_Y.resize(NUM_POIS,0);
    POI_values.resize(NUM_POIS,0);
    
    /// CORNER POIS
    POI_X.at(0)=1;
    POI_Y.at(0)=1;
    POI_values.at(0)=1;
    
    POI_X.at(1)=XMAX-2;
    POI_Y.at(1)=1;
    POI_values.at(1)=1;
    
    POI_X.at(2)=XMAX-2;
    POI_Y.at(2)=YMAX-2;
    POI_values.at(2)=1;
    
    POI_X.at(3)=1;
    POI_Y.at(3)=YMAX-2;
    POI_values.at(3)=1;
    
    
    
//    POI_X.at(4)=XMAX/2;
//    POI_Y.at(4)=YMAX-2;
//    POI_values.at(4)=1;
//        
//    POI_X.at(5)=XMAX/2;
//    POI_Y.at(5)=1;
//    POI_values.at(5)=1;
//    
//    POI_X.at(6)=1;
//    POI_Y.at(6)=YMAX/2;
//    POI_values.at(6)=1;
//    
//    POI_X.at(7)=XMAX-2;
//    POI_Y.at(7)=YMAX/2;
//    POI_values.at(7)=1;
    
  
    //Set one POI randomly to a higher value
    POI_values.at(rand() % NUM_POIS) = 5;
    
}

double gridworld::rover_poi_dist(int roverdex, int poidex){
    double dx = fabs(rover_x.at(roverdex) - POI_X.at(poidex));
    double dy = fabs(rover_y.at(roverdex) - POI_Y.at(poidex));
    return dx+dy;
}

void gridworld::evaluate()
{
    
    global=0;
    rover_local.clear();
    rover_difference.clear();
    rover_gzmi.clear();
    rover_local.resize(NUM_AGENTS,0);
    rover_gzmi.resize(NUM_AGENTS,0);
    rover_difference.resize(NUM_AGENTS,0);

    for(int p=0; p<NUM_POIS; p++){
        int closest=find_closest_rover_index(p);
        int second=find_closest_rover_index_not(closest,p);
        double d1 = rover_poi_dist(closest,p);
        double d2 = rover_poi_dist(second,p);
        double best = getReward(d1)*POI_values.at(p);
        double backup = getReward(d2)*POI_values.at(p);
        //cout << "POI #" << p << ": " << closest << " " << second << "; " << d1 << " <= " << d2 << " -> " << best << " , " << backup << endl; 
        rover_local.at(closest)+=best;
        rover_local.at(second)+=backup; /// locally the second rover thinks he's doing a good job.
        
        //the loop below makes all agents receive rewards from all POIs when using local 
        // i.e. zero consideration of other agents, just how close am I to each POI
//        for(int agent=0; agent<NUM_AGENTS; agent++) {
//            if (agent != closest || agent != second) {
//                rover_local.at(agent)+=getReward(rover_poi_dist(agent,p))*POI_values.at(p);
//            }
//        }
        
        rover_gzmi.at(closest)+=backup;
        rover_difference.at(closest)+=best-backup;
        rover_difference.at(second)+=0; /// no reward for second place
        global+=best;
    }
}

int gridworld::find_closest_rover_index(int poidex){
    int roverdex=-1;
    int mindist=99999999;
    for(int i=0; i<NUM_AGENTS; i++){
        int a=rover_poi_dist(i,poidex);
        if(a<mindist){
            mindist=a;
            roverdex=i;
        }
    }
    return roverdex;
}

int gridworld::find_closest_rover_index_not(int exclude, int poidex){
    int roverdex=-1;
    int mindist=99999999;
    for(int i=0; i<NUM_AGENTS; i++){
        if(i==exclude){continue;}
        int a=rover_poi_dist(i,poidex);
        if(a<mindist){
            mindist=a;
            roverdex=i;
        }
    }
    return roverdex;
}


void gridworld::start() { /// Initializes values.
    rover_local.clear();
    rover_difference.clear();
    rover_gzmi.clear();
    rover_x.clear();
    rover_y.clear();
    global=0;
    start_pois();
      
    rover_local.resize(NUM_AGENTS,0);
    rover_gzmi.resize(NUM_AGENTS,0);
    rover_difference.resize(NUM_AGENTS,0);
}

// Same as above without reinitialising the POIs
// So that the randomly chosen higher value POI remains the same throughout a run
void gridworld::start_epsiode() { 
    rover_local.clear();
    rover_difference.clear();
    rover_gzmi.clear();
    rover_x.clear();
    rover_y.clear();
    global=0;
      
    rover_local.resize(NUM_AGENTS,0);
    rover_gzmi.resize(NUM_AGENTS,0);
    rover_difference.resize(NUM_AGENTS,0);
}

void gridworld::obtain_locations(vector<QLearner>* pA) { /// Determines agent attendance per lane
    rover_x.clear();
    rover_y.clear();

    for (int i = 0; i < NUM_AGENTS; i++) {
        rover_x.push_back(pA->at(i).x);
        rover_y.push_back(pA->at(i).y);
    }
}


void sense(vector<QLearner>* pA) {
    for (int agent = 0; agent < NUM_AGENTS; agent++) {
        pA->at(agent).previousState = pA->at(agent).state;
    }
}

void decide(vector<QLearner>* pA) { ///agents decide on an action
    for (int agent = 0; agent < NUM_AGENTS; agent++) {
        pA->at(agent).decay_alpha();
        pA->at(agent).decay_epsilon();
        pA->at(agent).choose_egreedy_action();
    }
}

void act(vector<QLearner>* pA, gridworld* pE) { //agents make their action
    
    for (int agent = 0; agent < NUM_AGENTS; agent++) {
        //Move agent's 95 % of the time
        //Other 5% movement fails
        //if (rand() % 20 != 0){
            double x = pA->at(agent).x;
            double y = pA->at(agent).y;

            //Calculate next state for all agents
            if(pA->at(agent).action==0){x=x;y=y;}
            if(pA->at(agent).action==1){x=x;y=y+VELOCITY;}
            if(pA->at(agent).action==2){x=x+VELOCITY;y=y;}
            if(pA->at(agent).action==3){x=x;y=y-VELOCITY;}
            if(pA->at(agent).action==4){x=x-VELOCITY;y=y;}
            /// Edge correction
            if(x<0){x=0.0;y=y;}
            if(x>=XMAX){x=XMAX-VELOCITY;y=y;}
            if(y<0){x=x;y=0.0;}
            if(y>=YMAX){x=x;y=YMAX-VELOCITY;}
        
            pA->at(agent).x = x;
            pA->at(agent).y = y;
        //}
        pA->at(agent).set_state();
        //cout << pA->at(agent).x << "," << pA->at(agent).y << "\t";
    }
    //cout << endl;
    pE->obtain_locations(pA);
}

void react(vector<QLearner>* pA, gridworld* pE) { /// reward calculations, Q updates.
    pE->evaluate();
    for (int agent = 0; agent < NUM_AGENTS; agent++) {
        double L;
        double G;
        double D;
        double GZMI;

        L = pE->rover_local.at(agent);
        G = pE->global;
        D = pE->rover_difference.at(agent);

        GZMI = pE->rover_gzmi.at(agent);

        //Manual heuristic encouraging agents maximise capacity on all but one lane
        if(command_difference_PBRS_hand || command_global_PBRS_hand || command_local_PBRS_hand)
        {
            pA->at(agent).currentPhi = fabs(pA->at(agent).x - 5.0) + fabs(pA->at(agent).y - 5.0);
        }

        //Automated multi-agent potential function
        if(command_global_PBRS_gzmi || command_local_PBRS_gzmi || command_difference_PBRS_auto){
            pA->at(agent).currentPhi = GZMI;   //Potential-based difference reward
        }

        //Calculate potential based reward
        double PBRS = pA->at(agent).gamma * pA->at(agent).currentPhi - pA->at(agent).previousPhi;
        double shapedReward=0;
        if(command_global_PBRS_gzmi || command_global_PBRS_hand)
        {
            shapedReward = G + PBRS;
        }
        if(command_difference_PBRS_hand || command_difference_PBRS_auto)
        {
            shapedReward = D + PBRS;
        }
        if(command_local_PBRS_hand || command_local_PBRS_gzmi)
        {
            shapedReward = L + PBRS;            
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

        pA->at(agent).set_local(L);
        pA->at(agent).set_global(G);
        pA->at(agent).set_difference(D);
        pA->at(agent).set_gzmi(GZMI);
        pA->at(agent).set_shaped_reward(shapedReward);

        if(command_local){pA->at(agent).learn_with_local();}
        if(command_global){pA->at(agent).learn_with_global();}
        if(command_difference){pA->at(agent).learn_with_difference();}
        if(command_difference_PBRS_hand || command_difference_PBRS_auto ||
           command_global_PBRS_gzmi || command_global_PBRS_hand ||
           command_local_PBRS_hand || command_local_PBRS_gzmi)
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

    for(int method=0; method < 9; method++)
    {
        if(method==0)
        {
            cout << "HERE BEGINS GLOBAL REWARDS" << endl;
            pFILE = fopen("global.txt", "w");
            command_global = true;
            command_difference=false;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = false;
            command_difference_PBRS_auto = false;
            command_local = false;
            command_local_PBRS_hand = false;
            command_local_PBRS_gzmi = false;
        }
        if(method==1)
        {
            cout << "HERE BEGINS DIFFERENCE REWARDS" << endl;
            pFILE = fopen("difference.txt", "w");
            command_global = false;
            command_difference=true;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = false;
            command_difference_PBRS_auto = false;
            command_local = false;
            command_local_PBRS_hand = false;
            command_local_PBRS_gzmi = false;
        }
        if(method==2){
            cout << "HERE BEGINS GLOBAL + PBRS (AUTO) REWARDS" << endl;
            pFILE = fopen("global_PBRS_gzmi.txt", "w");
            command_global = false;
            command_difference=false;
            command_global_PBRS_gzmi = true;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = false;
            command_difference_PBRS_auto = false;
            command_local = false;
            command_local_PBRS_hand = false;
            command_local_PBRS_gzmi = false;
        }
        if(method==3){
            cout << "HERE BEGINS GLOBAL + PBRS (HAND) REWARDS" << endl;
            pFILE = fopen("global_PBRS_hand.txt", "w");
            command_global = false;
            command_difference=false;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = true;
            command_difference_PBRS_hand = false;
            command_difference_PBRS_auto = false;
            command_local = false;
            command_local_PBRS_hand = false;
            command_local_PBRS_gzmi = false;
        }
        if(method==4){
            cout << "HERE BEGINS DIFFERENCE + PBRS (HAND) REWARDS" << endl;
            pFILE = fopen("difference_pbrs_hand.txt", "w");
            command_global = false;
            command_difference=false;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = true;
            command_difference_PBRS_auto = false;
            command_local = false;
            command_local_PBRS_hand = false;
            command_local_PBRS_gzmi = false;
        }
        if(method==5){
            cout << "HERE BEGINS LOCAL REWARDS" << endl;
            pFILE = fopen("local.txt", "w");
            command_global = false;
            command_difference=false;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = false;
            command_difference_PBRS_auto = false;
            command_local = true;
            command_local_PBRS_hand = false;
            command_local_PBRS_gzmi = false;
        }
        if(method==6){
            cout << "HERE BEGINS LOCAL + PBRS(AUTO) REWARDS" << endl;
            pFILE = fopen("local_PBRS_gzmi.txt", "w");
            command_global = false;
            command_difference=false;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = false;
            command_difference_PBRS_auto = false;
            command_local = false;
            command_local_PBRS_hand = false;
            command_local_PBRS_gzmi = true;
        }
        if(method==7){
            cout << "HERE BEGINS LOCAL + PBRS(HAND) REWARDS" << endl;
            pFILE = fopen("local_PBRS_hand.txt", "w");
            command_global = false;
            command_difference=false;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = false;
            command_difference_PBRS_auto = false;
            command_local = false;
            command_local_PBRS_hand = true;
            command_local_PBRS_gzmi = false;
        }
        if(method==8){
            cout << "HERE BEGINS DIFFERENCE + PBRS (AUTO) REWARDS" << endl;
            pFILE = fopen("difference_pbrs_auto.txt", "w");
            command_global = false;
            command_difference=false;
            command_global_PBRS_gzmi = false;
            command_global_PBRS_hand = false;
            command_difference_PBRS_hand = false;
            command_difference_PBRS_auto = true;
            command_local = false;
            command_local_PBRS_hand = false;
            command_local_PBRS_gzmi = false;
        }
        
        
        double episode_global;
        for(int stat_run=0; stat_run < STAT_RUNS; stat_run++) {
            gridworld grid;
            //grid.make_lookup(); Obsolete? //TODO: Remove?
            grid.start();
            gridworld* pE = &grid;

            vector<QLearner> Agents;
            vector<QLearner>* pA = &Agents;

            for (int i = 0; i < NUM_AGENTS; i++) {
                QLearner Q;
                Q.id = i;
                Q.start();
                pA->push_back(Q);
            }

            for (int episode = 0; episode < EPISODES; episode++){
                grid.start_epsiode();
                episode_global=0;
                
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
                    episode_global+=pE->global;
               }

               //For beautiful graphs
               if (pretty_print) {
                   report(pFILE,episode_global); // Report every result
               } else {                
                   //For Coarse Results
                   if (episode % (EPISODES / 1000) == 0) {
                       report(pFILE,episode_global);      //Report only occasionally
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
//               if(episode%100 ==0){
//               cout << "GLOBAL EPISODE: " << episode << "\t" << episode_global << endl;
//               }
            }

            //Start a new line in output file for next run
            fprintf(pFILE,"\n");

            cout << endl << "Final performance = " << episode_global << endl << endl;     //Final global reward
        }
        fclose(pFILE);
    }
    return 0;
}

