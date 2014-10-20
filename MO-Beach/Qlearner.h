#ifndef QLEARNER_H
#define	QLEARNER_H

#ifndef VECTOR_INCLUDE
#define VECTOR_INCLUDE
#include <vector.h>
#include <math.h>
#endif

class QLearner {
    double learning_signal; /// The value used to update the Q table. G, L, or D are assigned to this value.

    double global_reward;
    double local_reward;
    double difference_reward;
    double gzmi; /// G(z_{-i});
    double shaped_reward;
    
    void console_1vector(vector<double>);

    int rand_action();
    int greedy_action();
    void initial_Qtable();
    void set_initial_state();
    void show_action();
    

    double Initial_Q_Value;

public:
    /// AGENT DOMAIN SPECIFIC:
    double instant_local_moves;
    double instant_global_moves;
    double total_local_moves;
    double total_global_moves;
    
    /// AGENT GENERIC:
    void console_2vector(vector<vector<double> >);
    vector< vector<double> > Qtable;
    int id_num;
    int previousState;
    int state;
    int noise_state;
    int signal_state;
    int action;
    int learnability_action;
    int epsilon_counter;
    double alpha;
    double gamma;
    double epsilon;
    double currentPhi;
    double previousPhi;
    void start();       //Used before first run of a statistical run / repeat
    void restart();     //Used at start of episode (i.e. does not change Qtable or learning parameters)
    void timestep_restart();
    void Qupdate();
    void final_Qupdate();
    void choose_egreedy_action();
    void choose_greedy_action();
    void randomise_learnability_action();
    void set_local(double);
    void set_global(double);
    void set_difference(double);
    void set_gzmi(double);
    void set_shaped_reward(double);
    void learn_with_global();
    void learn_with_difference();
    void learn_with_local();
    void learn_with_shaped_reward();
    void decay_epsilon();
    void decay_alpha();
};

void QLearner::decay_epsilon()
{
    epsilon_counter++;
    epsilon=1.0/epsilon_counter;
}
void QLearner::decay_alpha()
{
    alpha*=0.9999;
}

void QLearner::set_initial_state()
{
    //state = id % BEACHES; //uniform initial distribution
    
    // Initialise states - NOTE: This version is specific to 5 sections
    /*
    if (id_num < NUM_AGENTS / 2) {
        state = 1;
    } else {
        state = 3;
    }
     */
    
    /// All same initial distribution
    state = BEACHES/2;
}

void QLearner::start() {
    total_global_moves = 0;
    total_local_moves = 0;
    instant_local_moves = 0;
    instant_global_moves = 0;
    
    previousState = 0;
    set_initial_state();
   
    noise_state = 0; 
    signal_state = 0;
    action = 0;
    learnability_action = 0;
    alpha = 0.05;
    epsilon_counter=0;
    epsilon = 1;
    gamma = 0.9;
    Initial_Q_Value = -1;
    learning_signal = 0;

    global_reward = 0;
    local_reward = 0;
    difference_reward = 0;
    gzmi = 0;
    shaped_reward = 0;

    currentPhi = 0.0;
    previousPhi = 0.0;
    
    //Initialise Q-table, SECTIONS(STATES)x3(ACTIONS)    
    for (int i = 0; i < BEACHES; i++) {
        vector<double> vec;
        vec.resize(ACTIONS,0); /// create inner vector   
        Qtable.push_back(vec); /// push inner vector into outer vector
    }
    initial_Qtable();
}

void QLearner::restart() {
    previousState = 0;
    set_initial_state();
    
    total_global_moves = 0;
    total_local_moves = 0;
   
    noise_state = 0; //Todo: update if calculating learnability? 
    signal_state = 0; //TODO: Same as above, should these just equal state?
    action = 0;
    learnability_action = 0;

    learning_signal = 0;
    global_reward = 0;
    local_reward = 0;
    difference_reward = 0;
    gzmi = 0;
    shaped_reward = 0;

    currentPhi = 0.0;
    previousPhi = 0.0;
}

void QLearner::timestep_restart(){
    instant_local_moves=0;
    instant_global_moves=0;
    
    learning_signal = 0;
    global_reward = 0;
    local_reward = 0;
    difference_reward = 0;
    gzmi = 0;
    shaped_reward = 0;
}

void QLearner::initial_Qtable() { 
   for (int i = 0; i < Qtable.size(); i++) {
        for (int j = 0; j < Qtable.at(i).size(); j++) {
                Qtable.at(i).at(j) = Initial_Q_Value + LYRAND*SMALL - LYRAND*SMALL;
        }
   }
}

void QLearner::console_2vector(vector< vector<double> > a) {
    for (int i = 0; i < a.size(); i++) {
        console_1vector(a.at(i));
        cout << endl;
    }
}

void QLearner::console_1vector(vector<double> a) {
    for (int i = 0; i < a.size(); i++) {
        cout << a.at(i);
        cout << "\t";
    }
}

void QLearner::choose_egreedy_action() {
    double a = (double) rand() / RAND_MAX;
    if (a < epsilon) {
        action = rand_action();
    } else {
        action = greedy_action();
    }
    
    a = (double) rand() / RAND_MAX;
    if (a < epsilon) {
        learnability_action = rand_action();
    } else {
        learnability_action = greedy_action();
    }
}

void QLearner::choose_greedy_action() {
    action = greedy_action();
}

void QLearner::randomise_learnability_action() {
    learnability_action = rand_action();
}

int QLearner::rand_action() {
    int a;
    a = rand() % ACTIONS;
    return a;
}

int QLearner::greedy_action() {
    int LL = ACTIONS;
    double best = -9999999999;
    int bestdex = -1;
    for (int i = 0; i < LL; i++) {
        if (Qtable.at(state).at(i) > best) {
            best = Qtable.at(state).at(i);
            bestdex = i;
        }
    }
    return bestdex;
}

void QLearner::show_action() {
   // cout << "Agent " << index << ": " << action << endl;
}

void QLearner::learn_with_global() {
    learning_signal = global_reward;
    currentPhi= 0;
}

void QLearner::learn_with_difference() {
    learning_signal = difference_reward;
    currentPhi= 0;
}

void QLearner::learn_with_local() {
    learning_signal = local_reward;
    currentPhi= 0;
}

void QLearner::learn_with_shaped_reward() {
    learning_signal = shaped_reward; 
}

void QLearner::Qupdate() {
    double Q = Qtable.at(previousState).at(action);
    
    double Qmax = -9999999999;
    for (int i = 0; i < ACTIONS; i++) {
        if (Qtable.at(state).at(i) > Qmax) {
            Qmax = Qtable.at(state).at(i);
        }
    }
    //cout << "Q before: " << Q << endl;
    //cout << "In state: " << state << endl;
    //cout << "For action: " << action << endl;
    Q = Q + alpha * (learning_signal + gamma * Qmax - Q);
    //cout << "Q after: "  << Q << endl;
    Qtable.at(previousState).at(action) = Q;
    previousPhi = currentPhi;   //Handy for maintaining theoretical guarantees with dynamic potential functions
}

void QLearner::final_Qupdate() {
    double Q = Qtable.at(state).at(action);
    //cout << "Q before: " << Q << endl;
    //cout << "In state: " << state << endl;
    //cout << "For action: " << action << endl;
    Q = Q + alpha * (-previousPhi - Q);
    //cout << "Q after: "  << Q << endl;
    Qtable.at(state).at(action) = Q;
}


void QLearner::set_local(double L) {
    local_reward = L;
}

void QLearner::set_global(double G) {
    global_reward = G;
}

void QLearner::set_difference(double D) {
    difference_reward = D;
}

void QLearner::set_gzmi(double GZMI) {
    gzmi = GZMI;
}

void QLearner::set_shaped_reward(double R) {
    shaped_reward = R;
}



#endif	/* QLEARNER_H */

