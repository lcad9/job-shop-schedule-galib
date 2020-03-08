/******************************************************************************
*   Researchers involved:													  *
*   Fabio Henrique Pereira													  *
*   Flavio Grassi															  *
*   Luis Carlos dos Santos Junior                                             *
*																		      *
*   Contributors:															  *
*   Pedro Henrique Triguis													  *												  *
*   																		  *
*   Affiliation:															  *
*   Post-Graduation Program of Industrial Engineering, 						  *
*   Universidade Nove de Julho, Sao Paulo, SP, Brazil.						  *
*																			  *
*   Short Description:														  *
*   Program that uses Genetic Algorithms to find the optimal schedule 		  *
*   sequence for a job shop environment (JSSP) with method 3 phase.           *
*   It uses a develop methodology named "DSGA - Dynamic Seed Genetic 		  *
*   Algorithm", and runs using GAlib - a Genetic Algorithm library developed  *
*   in C++ by Matthew Wall, from the Massachusetts Institute of Technology.   *
******************************************************************************/

//Includes
#include "fact_LA2.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <queue>
#include <iostream>
#include <cstdlib>
#include <set>
#include <list>

//#include <ga/ga.h>
#include <ga/std_stream.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <ga/GASimpleGA.h>       // we're going to use the simple GA
#include <ga/GA2DBinStrGenome.h> // and the 2D binary string genome
#include <ga/GASStateGA.h>       // and Steady State GA
#include <ga/GARealGenome.h>
#include <ga/GARealGenome.C>
#include <windows.h> // for SHGetFolderPath function
#include <shlobj.h>  //for SHGetFolderPath function

//Defines
#define cout STD_COUT
#define ostream STD_OSTREAM

#define PRINT

using namespace std;

//User-defined function declarations
int factivel(int *vPriorities, int *vSequences, int iPrioMode);
void carregadados();
void checkmachines(int timestep, int iter, int iPrioMode);
void checkprocesses(int timestep, int iter, int iPrioMode);
double geraVarNormal(double media, double desvio);
float Objective(GAGenome &);
vector<int> gera_sequencia(int rota[]);
void setup();
void localSearch();
vector<int> localSearch(vector<int> seed, int originalScore, int bestLocalScore);
void localSearch(const GAStatistics &);
// Calcula e atualiza o caminho critico
void findCriticalPath(vector<int> s);
// Apenas calcula o caminho critico
vector<int> calcCriticalPath(vector<int> s);

int contLocalSearchCall;
void setCriticalPath(vector<int> s, int theBestScore);
void setSeed(vector< vector<int> > s);
void changeSeed(const GAStatistics &);
vector<int> gera_seed_ga();
void return_seed(vector<int> &seed, GAGenome &g);
float ObjectiveSeed(GAGenome &g);

//void changeSeed(const GAStatistics &);
string desktopFolder();
int projectDuration(vector<int> S);
int contSolucao = 0;
///Global variable declarations

int menorEncontrado = 9999;
int somaTotal = 0;

/// /////////////////////////////////////////////////////////////////////////
/// ******   RETORNA SEMENTE DEFINIDA PELO ATRASO OTIMO DO AG *********** ///
/// /////////////////////////////////////////////////////////////////////////
// vetor para a semente
    vector< vector<int> > SS_rseed(1);
    vector< vector<int> > S_rseed(MACHINE);
    vector< vector<int> > SSFinal_rseed(MACHINE);
    // vetor para a semente temporaria
    vector< vector<int> > Stemp_rseed(MACHINE);
    // filas nas m�quinas
    vector< vector<int> > Queue_rseed(MACHINE);
    // tempos de processamento
    vector< vector<int> > processTime_rseed(MACHINE);
    vector< int > seedTemp(MACHINE * JOB);
    // eventos de fim de atividade (idle time)
    vector<float> idleTime_rseed(MACHINE, 0);
    // contador de alocacao das maquinas
    vector<int> seized_rseed(MACHINE, 0);
    // vetor de pr�ximas m�quinas
    vector<int> nextStation_rseed(JOB);
    //Numero de fila
    vector<int> numberInQueue_rseed;
    // vetor de eventos programados para m�quinas
    vector<int> numberOfEventsPerMachine_rseed(MACHINE, 0);
    // vetor de passos do job
    vector<int> step_rseed(JOB, 0);
    // STATUS DAS M�QUINAS
    vector<int> status_rseed(MACHINE, 0);
    // STATUS ANTERIOR DAS M�QUINAS
    vector<int> previousStatus_rseed(MACHINE, 0);
    //Job a ser visitado
    vector<int> numberOfNodes_rseed;
    // Rel�gio de simula��o
    float t_rseed = 0;
    // Lista de eventos futuros
    list<float> futureEvents_rseed(1, 0.0);
    // Contador de opera��es escalonadas
    int k_seed = 0;
    //Contem o fitness
    float fitness_rseed = 0;
    //Armazena o delay
    float delay;
    //The job
    int theJob_rseed;
    //Tempo do processo até o momento
    int time_rseed;
    //Variavel que contrala a qtde de nós a serem processados
    int indiceVet = -1;
    // Está retornando
    bool rollback = false;
    //Maquina
    int mq = 0;
    //Numero na fila
    int nf = 0;
    int contTotal = 0;
    int score_temp = MACHINE * JOB * 1000;

/// Protótipos de funções phase 3
    void init(GAGenome &g);
    void phaseA();
    void phaseB(GAGenome &g);
    void phaseC(GAGenome &g);
    void runProcess(int mqn);
    void armazenarProcesso();
    void retornarProcesso();
    void atendeJobFIFO();
    void atendeJob(int indexJob);
    void backtracking(int sizeInQueue);
    void return_seed_ga(vector<int> &seed, GAGenome &g);
    void refinaIndividuos();
    void localSearchIn3phases();

class storeProcess{
    public:
         // vetor para a semente
        vector< vector<int> > S_class;
        // filas nas m�quinas
        vector< vector<int> > Queue_class;
        // tempos de processamento
        vector< vector<int> > processTime_class;
        // eventos de fim de atividade (idle time)
        vector<float> idleTime_class;
        // contador de alocacao das maquinas
        vector<int> seized_class;
        // vetor de pr�ximas m�quinas
        vector<int> nextStation_class;
        // vetor de eventos programados para m�quinas
        vector<int> numberOfEventsPerMachine_class;
        // vetor de passos do job
        vector<int> step_class;
        // STATUS DAS M�QUINAS
        vector<int> status_class;
        // STATUS ANTERIOR DAS M�QUINAS
        vector<int> previousStatus_class;
        //Armazena a ser visitado
        vector<int> numberOfNodes_class;
        //Numero de Jobs na fila
        vector<int> numberInQueue_class;
        // Rel�gio de simula��o
        float t_class;
        // Lista de eventos futuros
        list<float> futureEvents_class;
        // Contador de opera��es escalonadas
        int k_class;
        // valor do fitness currente
        float fitness_class;
        // Numero da maquina
        int mq;
        //Numero de job na fila
        int nf;
};

vector < storeProcess > storeStack;

//ft20
//int T[JOB * MACHINE] = {29,9,49,62,44,43,75,69,46,72,91,39,90,12,45,81,71,9,85,22,14,22,26,21,72,84,52,48,47,6,46,61,32,32,30,31,46,32,19,36,76,76,85,40,26,85,61,64,47,90,78,36,11,56,21,90,11,28,46,30,85,74,10,89,33,95,99,52,98,43,6,61,69,49,53,2,95,72,65,25,37,13,21,89,55,86,74,88,48,79,69,51,11,89,74,13,7,76,52,45};
//int R[MACHINE * JOB] = {1,2,3,4,5,1,2,4,3,5,2,1,3,5,4,2,1,5,3,4,3,2,1,4,5,3,2,5,1,4,2,1,3,4,5,3,2,1,4,5,1,4,3,2,5,2,3,1,4,5,2,4,1,5,3,3,1,2,4,5,1,3,2,4,5,3,1,2,4,5,1,2,5,3,4,2,1,4,5,3,1,3,2,4,5,1,2,5,3,4,2,3,1,4,5,1,2,3,4,5};
//int MAKESPAN = 1165;

// JSSP 5x5 Mirshekarian (2016) Correlation of job-shop scheduling problem features.
//int R[MACHINE * JOB] = {3, 1, 2, 5, 4, 3, 2, 1, 4, 5, 3, 2, 5, 4, 1, 5, 2, 3, 1, 4, 5, 2, 3, 4, 1};
//int T[JOB * MACHINE] = {96, 70, 6, 58, 86, 37, 26, 23, 14, 34, 21, 83, 29, 66, 25, 18, 73, 28, 29, 89, 87, 41, 32, 77, 77};
//int MAKESPAN = 465;

//int T[MACHINE * JOB] = {29,78,9,36,49,11,62,56,44,21,43,90,75,11,69,28,46,46,72,30,91,85,39,74,90,10,12,89,45,33,81,95,71,99,9,52,85,98,22,43,14,6,22,61,26,69,21,49,72,53,84,2,52,95,48,72,47,65,6,25,46,37,61,13,32,21,32,89,30,55,31,86,46,74,32,88,19,48,36,79,76,69,76,51,85,11,40,89,26,74,85,13,61,7,64,76,47,52,90,45};
//int R[JOB * MACHINE] = {0,1,2,3,4,5,6,7,8,9,0,2,4,9,3,1,6,5,7,8,1,0,3,2,8,5,7,6,9,4,1,2,0,4,6,8,7,3,9,5,2,0,1,5,3,4,8,7,9,6,2,1,5,3,8,9,0,6,4,7,1,0,3,2,6,5,9,8,7,4,2,0,1,5,4,6,8,9,7,3,0,1,3,5,2,9,6,7,4,8,1,0,2,6,8,9,5,3,4,7};
//int MAKESPAN = 500

int R[MACHINE*JOB]= {1,4,2,5,3,5,3,1,2,4,2,3,5,1,4,3,2,5,1,4,5,1,4,3,2,2,1,5,4,3,5,2,4,1,3,2,1,3,4,5,5,1,3,2,4,5,3,2,4,1};
int T[JOB*MACHINE]={20, 87, 31, 76, 17,25, 32, 24, 18, 81,72, 23, 28, 58, 99,86, 76, 97, 45, 90,27, 42, 48, 17, 46,67, 98, 48, 27, 62,28, 12, 19, 80, 50,63, 94, 98, 50, 80,14, 75, 50, 41, 55,72, 18, 37, 79, 61};
int MAKESPAN = 655;

//LA3
//int R[MACHINE*JOB]={2,3,1,5,4,3,2,1,5,4,3,4,5,1,2,5,1,3,2,4,5,1,2,4,3,5,1,2,3,4,4,3,1,5,2,5,2,1,3,4,5,1,4,3,2,5,2,1,3,4};
//int T[JOB*MACHINE]={23, 45, 82, 84, 38,21, 29, 18, 41, 50,38, 54, 16, 52, 52,37, 54, 74, 62, 57,57, 81, 61, 68, 30,81, 79, 89, 89, 11,33, 20, 91, 20, 66,24, 84, 32, 55, 8,56, 7, 54, 64, 39,40, 83, 19, 8, 7};
//int MAKESPAN = 597;

////////////////////////////////////////////////////////
//LA4
//int R[MACHINE*JOB]={1,3,4,5,2,2,4,5,/3,1,2,1,4,5,3,3,5,1,4,2,2,4,5,1,3,4,3,1,5,2,3,2,1,4,5,2,4,1,5,3,3,5,1,2,4,3,5,4,2,1};
//int T[JOB*MACHINE]={12, 94, 92, 91, 7,19, 11, 66, 21, 87,14, 75, 13, 16, 20,95, 66,  7,  7, 77,45,  6, 89, 15, 34,77, 20, 76, 88, 53, 74, 88, 52, 27,  9, 88, 69, 62, 98, 52, 61,  9, 62, 52, 90, 54,  5, 59, 15, 88};
//
/// /////////////////////////////////////////////////////
//LA5
//int R[MACHINE*JOB]={2,1,5,3,4,5,4,1,3,2,2,4,3,1,5,1,4,5,2,3,5,3,4,2,1,4,1,5,2,3,1,4,2,5,3,5,3,4,2,1,3,4,2,1,5,3,4,1,5,2};
//int T[JOB*MACHINE]={72, 87, 95, 66, 60, 5, 35, 48, 39, 54, 46, 20, 21, 97, 55, 59, 19, 46, 34, 37, 23, 73, 25, 24, 28, 28, 45, 5, 78, 83, 53, 71, 37, 29, 12, 12, 87, 33, 55, 38, 49, 83, 40, 48, 7, 65, 17, 90, 27, 23};
/// /////////////////////////////////////////////////////
//LA6
//int R[MACHINE*JOB]={2,3,5,1,4,4,5,2,3,1,3,1,2,4,5,4,2,5,1,3,5,4,3,2,1,3,2,1,4,5,1,4,2,5,3,1,2,3,5,4,3,4,5,1,2,1,5,4,2,3,5,3,1,4,2,1,5,3,2,4,5,4,2,3,1,5,2,1,3,4,1,2,3,5,4};
//int T[JOB*MACHINE]={21, 34, 95, 53, 55, 52, 16, 71, 26, 21, 31, 12, 42, 39, 98, 77, 77, 79, 55, 66, 37, 34, 64, 19, 83, 43, 54, 92, 62, 79, 93, 69, 87, 77, 87, 60, 41, 38, 83, 24, 98, 17, 25, 44, 49, 96, 77, 79, 75, 43, 28, 35, 95, 76, 7, 61, 10, 95, 9, 35, 59, 16, 91, 59, 46, 43, 52, 28, 27, 50, 87, 45, 39, 9, 41};
/// /////////////////////////////////////////////////////
//LA7
//int R[MACHINE*JOB]={1,5,2,4,3,1,2,5,4,3,4,1,3,2,5,1,2,5,4,3,4,2,1,3,5,2,3,1,5,4,3,2,1,5,4,3,4,5,1,2,5,1,3,2,4,5,1,2,4,3,5,1,2,3,4,4,3,1,5,2,5,2,1,3,4,5,1,4,3,2,5,2,1,3,4};
//int T[JOB*MACHINE]={47, 57, 71, 96, 14, 75, 60, 22, 79, 65, 32, 33, 69, 31, 58, 44, 34, 51, 58, 47, 29, 44, 62, 17, 8, 15, 40, 97, 38, 66, 58, 39, 57, 20, 50, 57, 32, 87, 63, 21, 56, 84, 90, 85, 61, 15, 20, 67, 30, 70, 84, 82, 23, 45, 38, 50, 21, 18, 41, 29, 16, 52, 52, 38, 54, 37, 54, 57, 74, 62, 57, 61, 81, 30, 68};
/// /////////////////////////////////////////////////////
//LA8
//int R[MACHINE*JOB]={4,3,1,5,2,3,2,1,4,5,2,4,1,5,3,3,5,1,2,4,3,5,4,2,1,5,4,3,2,1,5,4,1,2,3,4,3,1,2,5,4,1,5,3,2,5,3,4,1,2,1,2,5,4,3,1,5,3,4,2,1,4,5,3,2,4,2,1,5,3,3,1,4,2,5};
//int T[JOB*MACHINE]={92, 94, 12, 91, 7, 21, 19, 87, 11, 66, 14, 13, 75, 16, 20, 95, 66, 7, 77, 7, 34, 89, 6, 45, 15, 88, 77, 20, 53, 76, 9, 27, 52, 88, 74, 69, 52, 62, 88, 98, 90, 62, 9, 61, 52, 5, 54, 59, 88, 15, 41, 50, 78, 53, 23, 38, 72, 91, 68, 71, 45, 95, 52, 25, 6, 30, 66, 23, 36, 17, 95, 71, 76, 8, 88};
/// /////////////////////////////////////////////////////
//LA9
//int R[MACHINE*JOB]={2,4,3,1,5,4,2,3,5,1,5,4,2,3,1,1,2,3,4,5,1,5,3,4,2,1,4,2,3,5,4,3,5,2,1,3,2,4,5,1,3,5,1,2,4,3,5,4,2,1,5,4,3,1,2,2,1,4,5,3,5,4,1,2,3,4,2,3,1,5,1,2,3,5,4};
//int T[JOB*MACHINE]={66, 85, 84, 62, 19, 59, 64, 46, 13, 25, 88, 80, 73, 53, 41, 14, 67, 57, 74, 47, 84, 64, 41, 84, 78, 63, 28, 46, 26, 52, 10, 17, 73, 11, 64, 67, 97, 95, 38, 85, 95, 46, 59, 65, 93, 43, 85, 32, 85, 60, 49, 41, 61, 66, 90, 17, 23, 70, 99, 49, 40, 73, 73, 98, 68, 57, 9, 7, 13, 98, 37, 85, 17, 79, 41};
/// /////////////////////////////////////////////////////
//LA10
//int R[MACHINE*JOB]={2,3,4,1,5,2,1,5,4,3,1,2,3,5,4,4,2,3,1,5,3,1,2,4,5,4,5,3,1,2,2,5,1,3,4,3,4,2,5,1,1,4,5,2,3,3,5,4,1,2,1,5,4,3,2,3,1,2,5,4,4,3,2,5,1,2,3,5,1,4,4,3,1,5,2};
//int T[JOB*MACHINE]={58, 44, 5, 9, 58, 89, 97, 96, 77, 84, 77, 87, 81, 39, 85, 57, 21, 31, 15, 73, 48, 40, 49, 70, 71, 34, 82, 80, 10, 22, 91, 75, 55, 17, 7, 62, 47, 72, 35, 11, 64, 75, 50, 90, 94, 67, 20, 15, 12, 71, 52, 93, 68, 29, 57, 70, 58, 93, 7, 77, 27, 82, 63, 6, 95, 87, 56, 36, 26, 48, 76, 36, 36, 15, 8};
/// /////////////////////////////////////////////////////
//LA11
//int R[MACHINE*JOB]={3,2,1,4,5,1,4,2,5,3,1,2,3,5,4,3,4,5,1,2,1,5,4,2,3,5,3,1,4,2,1,5,3,2,4,5,4,2,3,1,5,2,1,3,4,1,2,3,5,4,1,4,2,5,3,5,3,1,2,4,2,3,5,1,4,3,2,5,1,4,5,1,4,3,2,2,1,5,4,3,5,2,4,1,3,2,1,3,4,5,5,1,3,2,4,5,3,2,4,1};
//int T[JOB*MACHINE]={34, 21, 53, 55, 95, 21, 52, 71, 16, 26, 12, 42, 31, 98, 39, 66, 77, 79, 55, 77, 83, 37, 34, 19, 64, 79, 43, 92, 62, 54, 93, 77, 87, 87, 69, 83, 24, 41, 38, 60, 25, 49, 44, 98, 17, 96, 75, 43, 77, 79, 95, 76, 7, 28, 35, 10, 95, 61, 9, 35, 91, 59, 59, 46, 16, 27, 52, 43, 28, 50, 9, 87, 41, 39, 45, 54, 20, 43, 14, 71, 33, 28, 26, 78, 37, 89, 33, 8, 66, 42, 84, 69, 94, 74, 27, 81, 45, 78, 69, 96};
/// /////////////////////////////////////////////////////

//LA16
//int R[MACHINE*JOB]={2,7,10,9,8,3,1,5,4,6,5,3,6,10,1,8,2,9,7,4,4,3,9,2,5,10,8,7,1,6,2,4,3,8,9,10,7,1,6,5,3,1,6,7,8,2,5,10,4,9,3,4,6,10,5,7,1,9,2,8,4,3,1,2,10,9,7,6,5,8,2,1,4,5,7,10,9,6,3,8,5,3,9,6,4,8,2,7,10,1,9,10,3,5,4,1,8,7,2,6};
//int T[JOB*MACHINE]={21, 71, 16, 52, 26, 34, 53, 21, 55, 95, 55, 31, 98, 79, 12, 66, 42, 77, 77, 39, 34, 64, 62, 19, 92, 79, 43, 54, 83, 37, 87, 69, 87, 38, 24, 83, 41, 93, 77, 60, 98, 44, 25, 75, 43, 49, 96, 77, 17, 79, 35, 76, 28, 10, 61, 9, 95, 35, 7, 95, 16, 59, 46, 91, 43, 50, 52, 59, 28, 27, 45, 87, 41, 20, 54, 43, 14, 9, 39, 71, 33, 37, 66, 33, 26, 8, 28, 89, 42, 78, 69, 81, 94, 96, 27, 69, 45, 78, 74, 84};
//int MAKESPAN = 945;
//
/// /////////////////////////////////////////////////////
//LA17
//int R[MACHINE*JOB]={5,8,10,3,4,9,6,7,2,1,9,6,2,8,3,4,7,10,5,1,3,5,4,2,9,7,8,1,10,6,1,9,4,8,6,3,5,7,2,10,10,1,5,9,7,3,6,4,8,2,4,3,6,1,8,5,9,2,7,10,2,8,9,4,5,6,7,1,3,10,2,8,3,1,9,7,4,10,6,5,3,4,5,10,1,7,8,9,2,6,2,1,6,4,10,8,9,3,7,5};
//int T[JOB*MACHINE]={18, 21, 41, 45, 38, 50, 84, 29, 23, 82, 57, 16, 52, 74, 38, 54, 62, 37, 54, 52, 30, 79, 68, 61, 11, 89, 89, 81, 81, 57, 91, 8, 33, 55, 20, 20, 32, 84, 66, 24, 40, 7, 19, 7, 83, 64, 56, 54, 8, 39, 91, 64, 40, 63, 98, 74, 61, 6, 42, 15, 80, 39, 24, 75, 75, 6, 44, 26, 87, 22, 15, 43, 20, 12, 26, 61, 79, 22, 8, 80, 62, 96, 22, 5, 63, 33, 10, 18, 36, 40, 96, 89, 64, 95, 23, 18, 15, 64, 38, 8};
//int MAKESPAN = 784;
/// /////////////////////////////////////////////////////
//LA18
//int R[MACHINE*JOB]={7,1,5,4,8,9,2,6,3,10,4,10,7,6,1,9,5,3,8,2,5,2,9,1,8,7,6,4,10,3,10,2,5,4,9,3,7,1,8,6,4,3,7,10,8,1,5,6,2,9,2,5,1,3,10,7,8,9,6,4,2,4,1,3,10,8,9,5,7,6,6,4,7,2,1,8,9,10,3,5,2,1,8,5,4,6,10,9,7,3,5,9,3,4,2,7,8,10,6,1};
//int T[JOB*MACHINE]={54, 87, 48, 60, 39, 35, 72, 95, 66, 5, 20, 46, 34, 55, 97, 19, 59, 21, 37, 46, 45, 24, 28, 28, 83, 78, 23, 25, 5, 73, 12, 37, 38, 71, 33, 12, 55, 53, 87, 29, 83, 49, 23, 27, 65, 48, 90, 7, 40, 17, 66, 25, 62, 84, 13, 64, 46, 59, 19, 85, 73, 80, 41, 53, 47, 57, 74, 14, 67, 88, 64, 84, 46, 78, 84, 26, 28, 52, 41, 63, 11, 64, 67, 85, 10, 73, 38, 95, 97, 17, 60, 32, 95, 93, 65, 85, 43, 85, 46, 59};
//int MAKESPAN = 848;
/// /////////////////////////////////////////////////////
//LA19
//int R[MACHINE*JOB]={3,4,6,5,1,8,9,10,2,7,5,8,2,9,1,4,3,6,10,7,10,7,5,4,2,1,9,3,8,6,2,3,8,6,9,5,4,7,10,1,7,2,4,1,3,9,5,8,10,6,8,6,9,3,5,7,4,2,10,1,7,2,5,6,3,4,8,9,10,1,1,6,9,10,4,7,5,8,3,2,6,3,4,7,5,8,9,10,2,1,10,5,7,8,1,3,9,6,4,2};
//int T[JOB*MACHINE]={44, 5, 58, 97, 9, 84, 77, 96, 58, 89, 15, 31, 87, 57, 77, 85, 81, 39, 73, 21, 82, 22, 10, 70, 49, 40, 34, 48, 80, 71, 91, 17, 62, 75, 47, 11, 7, 72, 35, 55, 71, 90, 75, 64, 94, 15, 12, 67, 20, 50, 70, 93, 77, 29, 58, 93, 68, 57, 7, 52, 87, 63, 26, 6, 82, 27, 56, 48, 36, 95, 36, 15, 41, 78, 76, 84, 30, 76, 36, 8, 88, 81, 13, 82, 54, 13, 29, 40, 78, 75, 88, 54, 64, 32, 52, 6, 54, 82, 6, 26};
//int MAKESPAN = 842;
/// /////////////////////////////////////////////////////
//LA20
//int R[MACHINE*JOB]={7,2,5,3,9,4,1,6,10,8,8,3,10,5,2,6,9,1,4,7,3,6,1,4,2,7,5,9,8,10,5,7,8,1,3,6,4,2,10,9,1,7,5,2,3,4,10,9,6,8,3,7,4,6,2,9,1,10,5,8,5,4,2,6,7,8,9,10,1,3,2,8,4,5,7,10,9,1,3,6,4,9,1,3,2,6,5,10,8,7,1,3,4,6,7,10,9,5,8,2};
//int T[JOB*MACHINE]={9, 81, 55, 40, 32, 37, 6, 19, 81, 40, 21, 70, 65, 64, 46, 65, 25, 77, 55, 15, 85, 37, 40, 24, 44, 83, 89, 31, 84, 29, 80, 77, 56, 8, 30, 59, 38, 80, 41, 97, 91, 40, 88, 17, 71, 50, 59, 80, 56, 7, 8, 9, 58, 77, 29, 96, 45, 10, 54, 36, 70, 92, 98, 87, 99, 27, 86, 96, 28, 73, 95, 92, 85, 52, 81, 32, 39, 59, 41, 56, 60, 45, 88, 12, 7, 22, 93, 49, 69, 27, 21, 61, 68, 26, 82, 71, 44, 99, 33, 84};
//int MAKESPAN = 902;
/// ////////////////////////////////////////////////////
//LA21
/// ////////////////////////////////////////////////////
// LA22
/*
int R[MACHINE*JOB]={10,6,5,3,8,4,2,1,9,7,
                    4,3,5,2,10,1,7,6,8,9,
                    9,8,3,1,10,6,7,4,2,5,
                    4,3,7,5,8,9,6,10,1,2,
                    5,7,2,3,8,1,9,6,4,10,
                    7,1,5,4,8,9,2,6,3,10,
                    4,10,7,6,1,9,5,3,8,2,
                    5,2,9,1,8,7,6,4,10,3,
                    10,2,5,4,9,3,7,1,8,6,
                    4,3,7,10,8,1,5,6,2,9,
                    2,5,1,3,10,7,8,9,6,4,
                    2,4,1,3,10,8,9,5,7,6,
                    6,4,7,2,1,8,9,10,3,5,
                    2,1,8,5,4,6,10,9,7,3,
                    5,9,3,4,2,7,8,10,6,1};

int T[JOB*MACHINE]={66, 91, 87, 94, 21, 92, 7, 12, 11, 19,
                    13, 20, 7, 14, 66,  75, 77, 16, 95, 7,
                    77, 20, 34, 15, 88, 89, 53, 6, 45, 76,
                    27, 74, 88, 62, 52, 69, 9, 98, 52, 88,
                    88, 15, 52, 61, 54, 62, 59, 9, 90,  5,
                    71, 41, 38, 53, 91, 68, 50, 78, 23,72,
                    95, 36, 66, 52, 45, 30, 23, 25, 17, 6,
                    65, 8, 85, 71, 65, 28, 88, 76, 27, 95,
                    37, 37, 28, 51, 86, 9, 55, 73, 51, 90,
                    39, 15, 83, 44, 53, 16, 46, 24, 25,82,
                    72, 48, 87, 66, 5, 54, 39, 35, 95, 60,
                    46, 20, 97, 21, 46, 37, 19, 59, 34, 55,
                    23, 25, 78, 24, 28, 83, 28, 5, 73, 45,
                    37, 53, 87, 38, 71, 29, 12, 33, 55, 12,
                    90, 17, 49, 83, 40, 23, 65, 27, 7, 48};
int MAKESPAN = 927;
*/

///Configurando os parametros

//vector<int> vetorTempo = T;

//LA20
//float MIN_VALUE = 0.0;
//float MAX_VALUE = 12.0;
//float INC = 0.001;

float MIN_VALUE = 0.0;
float MAX_VALUE = 60.0;
float INC = 10.0;

/////////////////////////////////////////////////////////////
///Salva popula��o inicial para testes dissertacao Valdemar
/////////////////////////////////////////////////////////////
string filePop = "PreDefinedIndividual." + (string)PROBLEMA + ".txt";
bool saveInitialPop = true;
int contPop = 0;

//SALVA OS ARQUIVOS NO PROPRIO DIRETORIO
string ArqCon = "Convergence." + (string)PROBLEMA + ".txt";
string ArqBsi = "BestSequenceIdentified." + (string)PROBLEMA + ".txt";
/////////////////////////////////////////////////////////////
int score_ant = ITER;
bool primeira_geracao_s = true, primeira_geracao;
int *s;
int horas, minutos, segundos, horas_seg = 3600;
GAGenome g_temp;
int cont;
int nBestIndividuals = 1;
int newBestIndividuals = 10;
vector<int> score_rodada_anterior(nBestIndividuals);
vector<int> bestSeedsIndex(nBestIndividuals);
vector<int> scoreTopList(nBestIndividuals);
std::vector<int>::iterator it;
int bestCurrent = ITER;
int firstBest = 0;
int seedIndexCurrent = 0;
vector< vector<int> > S(nBestIndividuals);
int outIter;
vector<int> currentCriticalPath;
bool STOP = false;
bool cpBasedSeed = false;
////////// PROPORCAO DE FACTIVEIS  /////////////////
int numberOfEvaluations = 0;
int numberOfFactibles = 0;
///////////// MATRIZ PARA CONTROLE DAS PERMUTA��ES ////////////
vector<int> M(MACHINE *JOB);

/******************************************************************************
* 									Main Program							  *
*******************************************************************************/
int main(int argc, char **argv)
{
    cout << "Novo teste para o problema:  " << PROBLEMA << "." << endl;
    cout << "Ira tentar encontrar o melhor sequenciamento de jobs.\n\n";
    cout.flush();

    remove(ArqCon.c_str());
    remove(ArqBsi.c_str());

    double start = clock();

    numberOfFactibles = 0;
    numberOfEvaluations = 0;

    int height = MACHINE;
    int width = JOB;

    GA2DBinaryStringGenome genome(width, height, Objective);

    GASteadyStateGA ga(genome);

    ga.minimize();

    // Crossover operators
    //genome.crossover(GA2DBinaryStringGenome::UniformCrossover);
    genome.crossover(GA2DBinaryStringGenome::EvenOddCrossover);
    //genome.crossover(OnePointCrossover); //default

    ga.populationSize(60);
    ga.pReplacement(0.9);
    ga.nGenerations(100);
    ga.pMutation(0.05);
    ga.pCrossover(0.85);
    ga.scoreFilename(ArqCon.c_str());
    ga.scoreFrequency(10);
    ga.selectScores(GAStatistics::Minimum);
    ga.flushFrequency(10);

    ga.nBestGenomes(newBestIndividuals);

    ga.pConvergence(1);
    ga.nConvergence((int)MACHINE * JOB * 3);
    ga.terminator(GAGeneticAlgorithm::TerminateUponConvergence);

    int outIter = 100;

    for (int i = 1; i <= outIter; i++)
    {
        primeira_geracao = true;

        if (i == 1)
        {
            for(int i = 0; i < 10; i++)
            {
                double end1 = clock();

                setup();

                cout << "Numero de solucoes encontradas: " << contSolucao << endl;

                double elapsed = ((double)(end1 - start)) / CLOCKS_PER_SEC;
                cout << "-->"
                     << " Setup time (in seconds): " << elapsed << "s" << endl;

            }
            return 0;
        }
        ga.initialize();
        //saveInitialPop = false;

        while (!ga.done())
        {
            cout << "\rIteracao " << i << " (de " << outIter << ") e Geracao " << ga.generation() + 2;
            ++ga;
        }
        //Muda a(s) Semente(s)
        changeSeed(ga.statistics());

        // Calcula Caminho Critico
        findCriticalPath(S[0]);

        //Faz busca local na semente atual
        //S[0] = localSearch(S[0], bestCurrent, bestCurrent);
        localSearch();

        cout << "\nScore = " << bestCurrent << endl;
    }

    cout << "\n\nComplete!!! Please check the following files in your 'Results' FOLDER: \n\"BestSequenceIdentified.txt\"\n\"Convergence.txt\"\n\n"
         << endl;

    double end2 = clock();
    double elapsed = ((double)(end2 - start)) / CLOCKS_PER_SEC;
    //    elapsed *= -1;
    horas = (elapsed / horas_seg);
    minutos = (elapsed - (horas_seg * horas)) / 60;
    segundos = (elapsed - (horas_seg * horas) - (minutos * 60));
    cout << "-->"
         << " Elapsed time (in seconds): " << elapsed << "s" << endl;
    printf("--> Elapsed time (human readable): %dh:%dm:%ds", horas, minutos, segundos);

    // PROPORÇAO DE SOLUCOES FACTIVEIS
    float prop = (float)numberOfFactibles / numberOfEvaluations;

    ofstream fileOut;
    fileOut.open(ArqBsi.c_str(), ios::app);
    fileOut << "Tempo de processamento = " << elapsed << "s [ " << horas << ":" << minutos << ":" << segundos << " ]" << endl;
    fileOut << "Problema testado = " << PROBLEMA << " (" << JOB << " jobs x " << MACHINE << " m�quinas)" << endl;
    fileOut << "\nPar�metros utilizados:" << endl;
    fileOut << "\tMakespan para indiv�duos n�o-fact�veis: " << ITER << endl;
    fileOut << "\tN�mero de la�os externos: " << outIter << endl;
    fileOut << "\tN�mero de gera��es: " << ga.nGenerations() + 1 << endl;
    fileOut << "\tN�mero m�nimo de gera��es para converg�ncia: " << ga.nConvergence() << endl;
    fileOut << "\tTamanho da popula��o: " << ga.populationSize() << endl;
    fileOut << "\tPercentual de substitui��o da popula��o: " << ga.pReplacement() * 100 << "%" << endl;
    fileOut << "\tProbabilidade de cruzamento: " << ga.pCrossover() * 100 << "%" << endl;
    fileOut << "\tProbabilidade de muta��o: " << ga.pMutation() * 100 << "%" << endl;
    fileOut << "\tN�mero de melhores indiv�duos (para busca local): " << ga.nBestGenomes() << endl;
    fileOut << "\tPropor��o de solu��es factiveis: " << prop * 100 << endl;
    fileOut << "\tN�mero de solu��es factiveis: " << numberOfFactibles << endl;
    fileOut << "\tN�mero total de avalia��es: " << numberOfEvaluations << endl;
    fileOut.close();

    return 0;
}

/******************************************************************************
* 	Function: changeSeed													  *
*	Short Description: Changes the seed of the next runs based on the best	  *
*					   seed from the previous one.							  *
*******************************************************************************/
void changeSeed(const GAStatistics &g)
{
    cout << endl
         << "************* Muda a semente *********** " << endl;

    // Estruturas temporarias para sementes
    vector<int> SS;

    // cout << "TopList " << n+1 << endl;
    // Pega o n-esimo individuo da lista dos melhores
    GA2DBinaryStringGenome &genome = (GA2DBinaryStringGenome &)g.bestIndividual();
    //bestCurrent = ITER;
    vector<int> score(nBestIndividuals);
    //Permuta a semente em S com base no genoma 'genome'
    // Copia S para uma estrutura temporaria que sofrer� permuta��o
    SS = S[0];

    for (int i = 0; i < genome.width(); i++)
    {
        for (int j = 0; j < genome.height(); j++)
        {
            if (genome.gene(i, j) == 1)
            {
                int idx = ((i + 1) % JOB);
                int val = SS[(j * JOB) + i];
                SS[(j * JOB) + i] = SS[(j * JOB) + idx];
                SS[(j * JOB) + idx] = val;
            }
        }
    }

    // Atualiza a semente original S
    S[0] = SS;
    return;
}

/******************************************************************************
* 	Function: gera_sequencia												  *
*	Short Description: Generates the initial seed from the routes given by    *
*					   the problem, based on FIFO dispatching rule.		  *
*******************************************************************************/

vector<int> gera_sequencia(int rota[])
{
    vector<int> s_temp(MACHINE * JOB);

    /// /////////////////////////////////////////////////////////////////////////
    cout << "\t GERA SEQUENCIA BASEADA EM 3 FASES " << endl;
    /// ************** USA AG PARA DETERMINAR O ATRASO OTIMO **************** ///
    /// /////////////////////////////////////////////////////////////////////////
    s_temp = gera_seed_ga();

    return s_temp;
}

/// ////////////////////////////////////////////////
/// Aplica GA para gerar semente otimizando delay
/// ////////////////////////////////////////////////
void init(GAGenome &g)
{
    GARealGenome &genome = (GARealGenome &)g;

    /// ================= IMPLEMENTACAO =============== ///
    for (int i = 0; i < JOB; i++)
    {
        // VERIFICA A MAQUINA DO STEP 0
        int mq = R[i * MACHINE] - 1;
        // COLOCA O JOB NA FILA DA MAQ
        Queue_rseed[mq].push_back(i);
        // TEMPO
        processTime_rseed[mq].push_back(T[i * MACHINE]);
        // NEXT STATION
        nextStation_rseed[i] = R[i * MACHINE + 1] - 1;
    }
    // ATUALIZA IDLE TIME, Lista eventos futuros e STATUS COM O ATRASO
    for (int mq = 0; mq < MACHINE; mq++)
    {
        float delay = genome.gene(mq * JOB + seized_rseed[mq]);
        idleTime_rseed[mq] = delay;

        futureEvents_rseed.push_back(idleTime_rseed[mq]);
        // ATUALIZA status da maquina mq
        if (delay != 0.0)
            status_rseed[mq] = 2;
        else if (delay == 0.0)
            status_rseed[mq] = 0;

        fitness_rseed = delay;
    }
}

/*****************************************************************************
* 	Function: return_seed                                                    *
*	Short Description: Retorna uma semente utilizando o algoritmo das 3 fases*
******************************************************************************/
void return_seed(vector<int> &seed, GAGenome &g)
{
   //Inicializa
    init(g);

    // faz enquanto n�o escalonar todos os jobs
    while (k_seed < MACHINE * JOB)
    {
       phaseA();
        // Em t=0 n�o h� atividade em andamento que possa se encerrar na fase B
        if (t_rseed > 0)
           phaseB(g);
    backtracking:
        phaseC(g);
    }

    /// Termina o processo de um sentido da arvores
    contTotal++;
    cout << "\rNumero de solucoes encontradas: " << ++contTotal;

    /*for(int i = 0; i < MACHINE; i++)
        for (int j = 0; j < JOB; j++)
            cout << S_rseed[i][j] << endl;*/


    /// Passa para temporária
    for (int i = 0; i < MACHINE; i++)
        for (int j = 0; j < JOB; j++)
            seedTemp[i * JOB + j] = S_rseed[i][j];

    SS_rseed[0] = seedTemp;

    //cout << "Fitness antes de refinar: " << score_temp << endl;

    refinaIndividuos();
    //cout << "\n";
    //cout << "\rMelhor MK encontrado: " << score_temp << "Solucoes testadas: " << ++contTotal;
   //cout << "Fitness apos de refinar: " << score_temp << endl;

    if(score_temp == MAKESPAN){
        cout << "Chegou no makespan" << endl;
        getchar();
    }
    //cout << "\rNumero de soluções encontras: " << ++contSolucao;

    if (menorEncontrado > score_temp)
    {
        menorEncontrado = score_temp;
    }

    //Condição para parar o processamento, visto que pode-se gerar muitas possibilidades
    if((indiceVet >= 0 && score_temp > MAKESPAN) && contTotal < 2000)
    {
        //cout << "\n Score: " << score_temp << endl;
        somaTotal += score_temp;
        retornarProcesso();
        rollback = true;
        goto backtracking;
    }

    cout << "\nMedia das solucoes encontradas: " << somaTotal / contTotal << endl;
    cout << "\nMenor valor encontrado: " << menorEncontrado << endl;

    for (int j = 0; j < JOB * MACHINE; j++)
        cout << SSFinal_rseed[0][j] << endl;

    //Depois de tudo estiver feito, passa a semente para a principal
    for (int j = 0; j < JOB * MACHINE; j++)
        seed[j] = SSFinal_rseed[0][j];
}

/*****************************************************************************
* 	Function: phaseA                                                        *
*	Short Description: Retorna uma semente utilizando o algoritmo das 3 fases*
******************************************************************************/
void phaseA()
{
    futureEvents_rseed.sort();
    t_rseed = *futureEvents_rseed.begin();

    // Remove primeiro elemento da lista de eventos futuros
    std::list<float>::iterator iter = futureEvents_rseed.erase(futureEvents_rseed.begin());
}

/*****************************************************************************
* 	Function: phaseB                                                         *
*	Short Description: Retorna uma semente utilizando o algoritmo das 3 fases*
******************************************************************************/
void phaseB(GAGenome &g)
{
    GARealGenome &genome = (GARealGenome &)g;

    // Encontra atividade
    std::vector<float>::iterator it = find(idleTime_rseed.begin(), idleTime_rseed.end(), t_rseed);

    if (it != idleTime_rseed.end())
    {
        // Identifica m�quina a ser liberada (do atendimento ou da espera)
        int mq = it - idleTime_rseed.begin();
        // Guarda status antigo
        previousStatus_rseed[mq] = status_rseed[mq];

        // Verifica se a m�quina estava em atendimento
        if (previousStatus_rseed[mq] == 1)
        {
            // Identifica job
            int job = S_rseed[mq].back() - 1; //S � 1-index
            // atualiza step do job
            step_rseed[job] = step_rseed[job] + 1;
            //
            // VERIFICA DELAY E COLOCA M�QUINA EM ESPERA
            delay = genome.gene(mq * JOB + seized_rseed[mq]);
            if (delay != 0)
            {
                status_rseed[mq] = 2;
                idleTime_rseed[mq] = t_rseed + delay;
                // Lista eventos futuros
                futureEvents_rseed.push_back(idleTime_rseed[mq]);
            }
            else
            {
                status_rseed[mq] = 0;
                idleTime_rseed[mq] = t_rseed - 1;
            }

            /// Verifica se job tem um pr�ximo passo
            if (step_rseed[job] < MACHINE)
            {
                //Pega m�quina atual
                int currentStation = nextStation_rseed[job];
                // Coloca job na pr�xima fila j� definida
                Queue_rseed[currentStation].push_back(job);
                // Tempo de atendimento do job da fila
                processTime_rseed[currentStation].push_back(T[job * MACHINE + step_rseed[job]]);
            }
            // Atualiza pr�xima m�quina do job (se houver)
            if (step_rseed[job] < MACHINE - 1)
                nextStation_rseed[job] = R[job * MACHINE + step_rseed[job] + 1] - 1;
        }
        else if (previousStatus_rseed[mq] == 2)
        {
            // Libera a m�quina
            status_rseed[mq] = 0;
            idleTime_rseed[mq] = t_rseed - 1;
        }
    } // if(it != idleTime.end())
}

/*****************************************************************************
* 	Function: phaseC                                                        *
*	Short Description: inicia as atividades por maquina                     *
******************************************************************************/
void phaseC(GAGenome &g)
{
    GARealGenome &genome = (GARealGenome &)g;

    if(rollback == true)
        goto returnTree;

    numberInQueue_rseed.clear();
    for (int j = 0; j < MACHINE; j++)
        numberInQueue_rseed.push_back(Queue_rseed[j].size());

    // Percorre as filas das maquinas
    for (mq = 0; mq < MACHINE; mq++)
    {
        // Verifica se a maquina est� livre
        if (status_rseed[mq] == 0)
        {
            // N�mero de elementos na fila
            nf = numberInQueue_rseed[mq];

            // H� fila?
            if (nf != 0)
            {
                /// 'Zera' o tamanho da fila visitada (visita apenas uma vez)
                numberInQueue_rseed[mq] = 0;

                delay = genome.gene(mq * JOB + seized_rseed[mq]);
                /// Verifica se m�quina estava em espera. Nesse caso usa LIFO
                if (delay != 0)
                {
                    /// Atende na base LIFO
                    time_rseed = processTime_rseed[mq].back();
                    theJob_rseed = Queue_rseed[mq].back();
                    /// Remove job da fila
                    Queue_rseed[mq].pop_back();
                    /// Remove tempo
                    processTime_rseed[mq].pop_back();
                    runProcess(mq);
                }
                else if (delay == 0)
                {
                returnTree:
                    //Pega o tamanho da fila de prioridade;
                    int sizeInQueue;

                    //Verifica se está retornando
                    if(rollback)
                    {
                        //Realiza a ação de retorno
                        backtracking(numberOfNodes_rseed.size());
                        runProcess(mq);
                    }
                    else
                    {
                        sizeInQueue = processTime_rseed[mq].size();

                        if(sizeInQueue == 1)
                        {
                            atendeJobFIFO();
                            runProcess(mq);
                        }
                        else if(sizeInQueue > 1)
                        {
                            //Guarda todos os jobs na fila
                            numberOfNodes_rseed.clear();
                            for(int i = 0; i < sizeInQueue; i++)
                                numberOfNodes_rseed.push_back(i);

                            //Pega o próxima (primeira) e apaga.
                            int indiceJob = numberOfNodes_rseed.at(0);
                            numberOfNodes_rseed.erase(numberOfNodes_rseed.begin());
                            armazenarProcesso();

                            atendeJob(indiceJob);
                            runProcess(mq);
                        }
                    }
                }
            }
        }
    }
}

/*************************************************************************
* 	Function: backtracking												  *
*	Short Description: Função faz o processo de retorno
************************************************************************/
void backtracking(int sizeInQueue)
{
    int indiceJob = numberOfNodes_rseed.at(0);
    numberOfNodes_rseed.erase(numberOfNodes_rseed.begin());

    if(sizeInQueue > 1) armazenarProcesso();

    atendeJob(indiceJob);

    rollback = false;
}

/*****************************************************************************
* 	Function: atendeJobFIFO                                                    *
*	Short Description:*
******************************************************************************/
//Atende na base FIFO
void atendeJobFIFO()
{
    /// Atende na base FIFO
    time_rseed = processTime_rseed[mq][0];
    theJob_rseed = Queue_rseed[mq][0];
    /// Remove job da fila
    Queue_rseed[mq].erase(Queue_rseed[mq].begin());
    /// Remove tempo
    processTime_rseed[mq].erase(processTime_rseed[mq].begin());
}

/*****************************************************************************
* 	Function: atendeJob                                                      *
*	Short Description: Atende um determinado Job na fila passando um indice  *
* @param indexJob parametro do indice                                        *
******************************************************************************/
void atendeJob(int indexJob)
{
    /// Atende na base FIFO
    time_rseed = processTime_rseed[mq][indexJob];
    theJob_rseed = Queue_rseed[mq][indexJob];
    /// Remove job da fila
    Queue_rseed[mq].erase(Queue_rseed[mq].begin() + indexJob) ;
    /// Remove tempo
    processTime_rseed[mq].erase(processTime_rseed[mq].begin() + indexJob);
}

/*****************************************************************************
* 	Function: runProcess                                                    *
*	Short Description: Após selecionar qual Job será atendido, executa o processo
******************************************************************************/
void runProcess(int mqn)
{
    // Inicia operacao do 'theJob' em mq
    k_seed++;
    // Muda status
    status_rseed[mqn] = 1;
    // Atualiza idleTime
    idleTime_rseed[mqn] = t_rseed + time_rseed;
    // Atualiza makespan
    if (idleTime_rseed[mqn] > fitness_rseed)
        fitness_rseed = idleTime_rseed[mqn];

    seized_rseed[mqn] = seized_rseed[mqn] + 1;
    // Lista eventos futuros
    futureEvents_rseed.push_back(idleTime_rseed[mqn]);
    // Coloca job na semente (1-index)
    S_rseed[mqn].push_back(theJob_rseed + 1);
}

/*****************************************************************************
* 	Function: armazenarProcesso                                              *
*	Short Description:
******************************************************************************/
//Armazena todo o processo feito até encontrar mais de uma tarefa ótima
void armazenarProcesso()
{
    indiceVet++;

    storeStack.push_back(storeProcess());

    storeStack[indiceVet].nf    = nf;
    storeStack[indiceVet].mq    = mq;
    storeStack[indiceVet].t_class   = t_rseed;
    storeStack[indiceVet].S_class   = S_rseed;
    storeStack[indiceVet].k_class   = k_seed;
    storeStack[indiceVet].step_class    = step_rseed;
    storeStack[indiceVet].Queue_class   = Queue_rseed;
    storeStack[indiceVet].status_class  = status_rseed;
    storeStack[indiceVet].seized_class  = seized_rseed;
    storeStack[indiceVet].fitness_class     = fitness_rseed;
    storeStack[indiceVet].idleTime_class    = idleTime_rseed;
    storeStack[indiceVet].nextStation_class     = nextStation_rseed;
    storeStack[indiceVet].processTime_class     = processTime_rseed;
    storeStack[indiceVet].futureEvents_class    = futureEvents_rseed;
    storeStack[indiceVet].numberInQueue_class   = numberInQueue_rseed;
    storeStack[indiceVet].numberOfNodes_class   = numberOfNodes_rseed;
    storeStack[indiceVet].previousStatus_class  = previousStatus_rseed;
    storeStack[indiceVet].numberOfEventsPerMachine_class    = numberOfEventsPerMachine_rseed;
}

/*****************************************************************************
* 	Function: retornarProcesso                                                    *
*	Short Description:
******************************************************************************/
//Retorna todo o processo do ultimo processo armazenado
void retornarProcesso()
{
    numberOfEventsPerMachine_rseed  = storeStack[indiceVet].numberOfEventsPerMachine_class;
    storeStack[indiceVet].numberOfEventsPerMachine_class.clear();
    previousStatus_rseed    = storeStack[indiceVet].previousStatus_class;
    storeStack[indiceVet].previousStatus_class.clear();
    numberOfNodes_rseed     = storeStack[indiceVet].numberOfNodes_class;
    storeStack[indiceVet].numberOfNodes_class.clear();
    numberInQueue_rseed     = storeStack[indiceVet].numberInQueue_class;
    storeStack[indiceVet].numberInQueue_class.clear();
    futureEvents_rseed  = storeStack[indiceVet].futureEvents_class;
    storeStack[indiceVet].futureEvents_class.clear();
    processTime_rseed   = storeStack[indiceVet].processTime_class;
    storeStack[indiceVet].processTime_class.clear();
    nextStation_rseed   = storeStack[indiceVet].nextStation_class;
    storeStack[indiceVet].nextStation_class.clear();
    idleTime_rseed  = storeStack[indiceVet].idleTime_class;
    storeStack[indiceVet].idleTime_class.clear();
    fitness_rseed   = storeStack[indiceVet].fitness_class;
    seized_rseed    = storeStack[indiceVet].seized_class;
    storeStack[indiceVet].seized_class.clear();
    status_rseed    = storeStack[indiceVet].status_class;
    storeStack[indiceVet].status_class.clear();
    Queue_rseed     = storeStack[indiceVet].Queue_class;
    storeStack[indiceVet].Queue_class.clear();
    step_rseed  = storeStack[indiceVet].step_class;
    storeStack[indiceVet].step_class.clear();
    t_rseed     = storeStack[indiceVet].t_class;
    k_seed     = storeStack[indiceVet].k_class;
    S_rseed     = storeStack[indiceVet].S_class;
    storeStack[indiceVet].S_class.clear();
    mq  = storeStack[indiceVet].mq;
    nf  = storeStack[indiceVet].nf;

    indiceVet--;
}

//////////////////////////////////////////////////////////
// otimiza o atraso com o AG
//////////////////////////////////////////////////////////
vector<int> gera_seed_ga()
{
    //float MIN_VALUE = 0.0;
    //float MAX_VALUE = 135.0;
    //float INC = 10.0;

    GARealAlleleSet alleles(MIN_VALUE, MAX_VALUE, INC);

    GARealGenome genome(MACHINE * JOB, alleles, ObjectiveSeed);
    // CRUZAMENTO
    genome.crossover(GARealUniformCrossover);   //LA3, LA4
    //genome.crossover(GARealEvenOddCrossover);
    //genome.crossover(GARealOnePointCrossover);
    //genome.crossover(GARealTwoPointCrossover); //LA2

    // MUTACAO
    genome.mutator(GARealSwapMutator);

    // GASimpleGA ga(genome);
    GASteadyStateGA ga(genome);
    ga.pReplacement(0.9);
    ga.minimize();               // by default we want to minimize the objective
    ga.populationSize(100);     // how many individuals in the population
    ga.nGenerations(100);        // number of generations to evolve
    ga.pMutation(0.05);          // likelihood of mutating new offspring
    ga.pCrossover(0.85);          // likelihood of crossing over parents
    ga.scoreFilename("bog.dat"); // name of file for scores
    ga.scoreFrequency(50);       // keep the scores of every generation
    ga.flushFrequency(50);       // specify how often to write the score to disk
    ga.selectScores(GAStatistics::Minimum);

    ga.initialize();

    while (!ga.done())
    {
        ga.step();
    }

    // Semente
    vector<int> s_temp(MACHINE * JOB);
    // Atrasos
    genome = ga.statistics().bestIndividual();
    // Semente
    return_seed(s_temp, genome);

    // IMPRIME SOLU�AO ENCONTRADA
    bestCurrent = genome.score();
    cout << "Makespan: " << genome.score() << endl;
    cout << "Atrasos e semente:" << endl;
    for (int i = 0; i < MACHINE * JOB; i++)
        cout << genome[i] << " \t" << s_temp[i] << endl;

    return s_temp;
}

////////////////////////////////////////////////////////////////////
// Faz setup da semente na primeira gera��o de todo o processo
////////////////////////////////////////////////////////////////////
void setup()
{
    //Gera um vector com uma semente
    S[0] = gera_sequencia(R);
    // Calcula Caminho Critico
    findCriticalPath(S[0]);

    /////////////////////////////////////////////////////////////////////////////
    ///   I N T E R P R E T A C A O   O R I G I N A L   D A   S E M E N T E   ///
    /////////////////////////////////////////////////////////////////////////////
    //avalia a semente e atualiza valor de aptid�o em bestCurrent
    int P[MACHINE * JOB];

    int maq = -1;
    for (int i = 0; i < MACHINE * JOB; i++)
    {
        int job = S[0][i] - 1; // Semente S[0] � 1-index; job 0-index
        int pos = (i % JOB) + 1;
        if (pos == 1)
            maq++;

        // Valores em P s�o 1-index
        P[maq * JOB + job] = pos;
    }
    // Calcula aptidao da semente
    int score = factivel(P, R, 1);
    bestCurrent = score;

    cout << endl
         << "*** Score: " << bestCurrent << endl;
    //getchar();

    // Aplica busca local
    //contLocalSearchCall = currentCriticalPath.size();
    //S[0] = localSearch(S[0], score, bestCurrent);
    localSearch();

    //Atualiza Caminho Critico
    setCriticalPath(S[0], bestCurrent);

    ///////////////////////////////////////////
    ///       COMPUTA      FACTIVEIS        ///
    ///////////////////////////////////////////
    numberOfEvaluations++;
    if (score < ITER)
        numberOfFactibles++;
    ///////////////////////////////////////////
    firstBest = score;
    // Inicialmente, todos os  valores s�o iguais
    for (int num = 0; num < nBestIndividuals; num++)
        scoreTopList[num] = score;
}

/*****************************************************************************
* 	Function: refinaInidividuos
*	Short Description: Faz um refinamento na solução encontrada, utilizando busca local
******************************************************************************/
void refinaIndividuos()
{
    /// Calcula Caminho Critico
    findCriticalPath(SS_rseed[0]);

    /// //////////////////////////////////////////////////////////////////////////
    ///   I N T E R P R E T A C A O   O R I G I N A L   D A   S E M E N T E   ///
    /// //////////////////////////////////////////////////////////////////////////
    ///avalia a semente e atualiza valor de aptid�o em bestCurrent
    int P[MACHINE * JOB];

    int maq = -1;
    for (int i = 0; i < MACHINE * JOB; i++)
    {
        int job = SS_rseed[0][i] - 1; // Semente S[0] � 1-index; job 0-index
        int pos = (i % JOB) + 1;
        if (pos == 1)
            maq++;

        // Valores em P s�o 1-index
        P[maq * JOB + job] = pos;
    }
    // Calcula aptidao da semente
    int score = factivel(P, R, 1);
    bestCurrent = score;

    // Aplica busca local no método de três fases
    localSearchIn3phases();

    //Atualiza Caminho Critico
    setCriticalPath(SS_rseed[0], bestCurrent);

    //cout << endl
      //   << "*** Score em refina: " << bestCurrent << endl;

    if(score_temp > bestCurrent)
    {
        score_temp = bestCurrent;
        SSFinal_rseed = SS_rseed;
    }
}

//*****************************************************************************
// calculate the makespan based on the critical path method
int projectDuration(vector<int> Seed)
{
    int *p;
    vector<int> startVec;
    vector<int> endVec;
    vector<int> timeVec;

    // Cria arestas do grafo
    for (int j = 0; j < JOB; j++)
    {
        startVec.push_back(0);
        timeVec.push_back(0);

        for (int m = 1; m <= MACHINE; m++)
        {
            endVec.push_back(j * MACHINE + m);         //1-index
            startVec.push_back(j * MACHINE + m);       //1-index
            timeVec.push_back(T[j * MACHINE + m - 1]); //0-index
        }
        endVec.push_back(JOB * MACHINE + 1);
    }
    //
    for (int m = 0; m < MACHINE; m++)
    {
        for (int j = 0; j < JOB; j++)
        {
            int job = Seed[m * JOB + j] - 1;
            // Vai na 'linha' de R referente ao job
            int start = job * MACHINE;
            int end = job * MACHINE + MACHINE - 1;
            p = std::find(R + start, R + end, m + 1);

            //Posicao
            int pos = p - (R + start);
            int no = job * MACHINE + pos + 1;
            //
            // Coloca na lista de arestas
            if (j == 0)
            {
                //seta apenas como n� inicial
                startVec.push_back(no);       //1-index
                timeVec.push_back(T[no - 1]); //0-index
            }
            if (j == JOB - 1)
            {
                //seta apenas como n� final
                endVec.push_back(no); //1-index
            }
            if (j != 0 && j != JOB - 1)
            {
                //seta como inicial e final
                startVec.push_back(no);       //1-index
                endVec.push_back(no);         //1-index
                timeVec.push_back(T[no - 1]); //0-index
            }
        }
    }
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    int numberVertex = MACHINE * JOB + 2; // number  of vertex
    int numberActivities = endVec.size(); // number of edges
    int j;

    /*vertices are of the form (v, u) */
    //indegree of each vertex (that is, count the number of edges entering them)
    vector<int> indegree(numberVertex, 0);

    std::vector<int> v(startVec);
    std::vector<int> u(endVec);
    std::vector<int> d(timeVec);

    int project_duration = 0; //project duration

    /*# Compute the indegree for each vertex v from the graph:
    for each neighbor u of v: indegree[u] += 1*/
    for (j = 0; j < numberActivities; j++)
    {
        indegree[u[j]]++;
    }

    queue<int> Q; //queue Q = empty queue
    int distance[numberVertex];
    memset(distance, 0, sizeof(int) * numberVertex); //distance = array filled with zeroes

    //vetor de conjuntos de nos
    vector< set<int> > vecSetCriticalPath(numberVertex);

    //for each vertex v:
    //if indegree[v] = 0:
    //insert v on Q
    for (j = 0; j < numberVertex; j++)
    {
        if (indegree[j] == 0)
            Q.push(v[j]);
    }

    int first;

    //printf ("first in the queue=%d\n", Q.front());

    /*for each neighbor u of v:
    d istance[u] = max(distance[u], distance[v] + time(v, u))
    indegree[u] -= 1
    if indegree[u] = 0:
    insert u on Q
    */
    while (!Q.empty())
    { //while Q is not empty:

        first = Q.front(); //v = get front element from Q
        Q.pop();           //delete de first from queue

        //Para todo vizinho de first faca
        for (int i = 0; i < numberActivities; i++) //**********MELHORAR ESSA ESTRUTURA **********
        {
            if (v[i] == first) //
            {
                //////////////////////////////////////////////////////////////
                //
                distance[u[i]] = std::max(distance[u[i]], distance[v[i]] + d[i]);
                indegree[u[i]] -= 1;
                ////////////////////////////////////////////////////////////////

                if (indegree[u[i]] == 0)
                {
                    Q.push(u[i]);
                }
            }
        }
    } //fecha while

    /*Now, select the vertex x with the largest distance.
    This is the minimum total project_duration.*/
    project_duration = *std::max_element(distance, distance + numberVertex);
    return project_duration;
}

/******************************************************************************
* Calcula o caminho cr�tico localmente
******************************************************************************/
vector<int> calcCriticalPath(vector<int> s)
{
    int *p;
    vector<int> startVec;
    vector<int> endVec;
    vector<int> timeVec;
    // Cria arestas do grafo
    for (int j = 0; j < JOB; j++)
    {
        startVec.push_back(0);
        timeVec.push_back(0);

        for (int m = 1; m <= MACHINE; m++)
        {
            //
            endVec.push_back(j * MACHINE + m);         //1-index
            startVec.push_back(j * MACHINE + m);       //1-index
            timeVec.push_back(T[j * MACHINE + m - 1]); //0-index
        }
        endVec.push_back(JOB * MACHINE + 1);
    }
    //
    for (int m = 0; m < MACHINE; m++)
    {
        for (int j = 0; j < JOB; j++)
        {
            int job = s[m * JOB + j] - 1;
            // cout << "Job: " << job+1 << endl;
            // Vai na 'linha' de R referente ao job
            int start = job * MACHINE;
            int end = job * MACHINE + MACHINE - 1;
            //cout << "start = " << start << " end = " << end << endl;
            p = std::find(R + start, R + end, m + 1);

            //cout << "Maq: " << m+1 << endl;
            //Posicao
            int pos = p - (R + start);
            //cout << "Operacao: " << pos+1 << endl;
            int no = job * MACHINE + pos + 1;
            //cout << endl << "****No***** " << no << endl;
            //
            // Coloca na lista de arestas
            if (j == 0)
            {
                //seta apenas como n� inicial
                startVec.push_back(no);       //1-index
                timeVec.push_back(T[no - 1]); //0-index
            }
            if (j == JOB - 1)
            {
                //seta apenas como n� final
                endVec.push_back(no); //1-index
            }
            if (j != 0 && j != JOB - 1)
            {
                //seta como inicial e final
                startVec.push_back(no);       //1-index
                endVec.push_back(no);         //1-index
                timeVec.push_back(T[no - 1]); //0-index
            }
        }
    }
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    int numberVertex = MACHINE * JOB + 2; // number  of vertex
    int numberActivities = endVec.size(); // number of edges
    int j;

    /*vertices are of the form (v, u) */
    //indegree of each vertex (that is, count the number of edges entering them)
    vector<int> indegree(numberVertex, 0);

    std::vector<int> v(startVec);
    std::vector<int> u(endVec);
    std::vector<int> d(timeVec);

    //    for (j=0;j<numberActivities; j++)
    //        cout << "***** " << v[j] << "  " << u[j] << "  " << d[j] << endl;
    //    getchar();

    int project_duration = 0; //project duration

    /*# Compute the indegree for each vertex v from the graph:
    for each neighbor u of v: indegree[u] += 1*/
    for (j = 0; j < numberActivities; j++)
    {
        indegree[u[j]]++;
    }

    //    for (j=0;j<numberVertex; j++)
    //        printf ("indegree %d= %d\n",j,indegree[j] );
    //getchar();

    queue<int> Q; //queue Q = empty queue
    int distance[numberVertex];
    memset(distance, 0, sizeof(int) * numberVertex); //distance = array filled with zeroes

    //vetor de conjuntos de nos
    vector< set<int> > vecSetCriticalPath(numberVertex);

    //for each vertex v:
    //if indegree[v] = 0:
    //insert v on Q
    for (j = 0; j < numberVertex; j++)
    {
        if (indegree[j] == 0)
            Q.push(v[j]);
    }

    int first;

    //printf ("first in the queue=%d\n", Q.front());

    /*for each neighbor u of v:
    d istance[u] = max(distance[u], distance[v] + time(v, u))
    indegree[u] -= 1
    if indegree[u] = 0:
    insert u on Q
    */
    while (!Q.empty())
    { //while Q is not empty:

        first = Q.front(); //v = get front element from Q
        Q.pop();           //delete de first from queue

        //Para todo vizinho de first faca
        for (int i = 0; i < numberActivities; i++) //**********MELHORAR ESSA ESTRUTURA **********
        {
            if (v[i] == first) //
            {
                //////////////////////////////////////////////////////////////
                //
                distance[u[i]] = std::max(distance[u[i]], distance[v[i]] + d[i]);
                indegree[u[i]] -= 1;
                ////////////////////////////////////////////////////////////////

                if (indegree[u[i]] == 0)
                {
                    Q.push(u[i]);
                }
            }
        }
        ////////////////////////////////////////////
        //
        //for(int ii=0; ii<numberVertex; ii++)
        //cout << *std::max_element(distance,distance+numberVertex) << endl;
        //getchar();
        ////////////////////////////////////////////
    } //fecha while

    /*Now, select the vertex x with the largest distance.
    This is the minimum total project_duration.*/
    project_duration = *std::max_element(distance, distance + numberVertex);
    //printf ("Total Project Duration %d\n", project_duration);
    //getchar();
    ////////////////////////////////////////////////////////////////////////////
    // Pega posicao do elemento m�ximo = primeiro na 'lista critica'
    int criticalNode = std::find(distance, distance + numberVertex, project_duration) - distance;
    vector<int> criticalPath;
    //criticalPath.push_back(criticalNode);
    ////////////////////////////////////////////////////////////////////////////
    while (criticalNode != 0)
    {
        vector<int> neighborhood;
        vector<int> neighDistances;
        //
        // varre recursivamente os vizinhos e pega o m�ximo
        for (int i = 0; i < numberActivities; i++) //**********MELHORAR ESSA ESTRUTURA **********
        {
            if (u[i] == criticalNode)
            {
                neighborhood.push_back(v[i]);
                neighDistances.push_back(distance[v[i]]);
                //cout << v[i] << " " << distance[v[i]] << endl;
                //getchar();
            }
        }
        int max = *std::max_element(neighDistances.begin(), neighDistances.end());
        int pp = std::find(neighDistances.begin(), neighDistances.end(), max) - neighDistances.begin();
        if (max == 0)
            pp = neighDistances.size() - 1;
        criticalNode = neighborhood[pp];
        criticalPath.push_back(criticalNode);
        //cout << "criticalNode = " << criticalNode << "Max: " << max << endl;
    }
    ////////////////////////////////////////////////////////////////////////////
    /*    cout << "criticalPath: " << endl;
    for (int j=0; j<criticalPath.size(); j++)
        cout << criticalPath[j] << endl;

    cout << '\n';
*/
    return criticalPath;
}

/******************************************************************************
* void findCriticalPath
*******************************************************************************/
void findCriticalPath(vector<int> s)
{
    int *p;
    vector<int> startVec;
    vector<int> endVec;
    vector<int> timeVec;

    // Cria arestas do grafo
    for (int j = 0; j < JOB; j++)
    {
        startVec.push_back(0);
        timeVec.push_back(0);

        for (int m = 1; m <= MACHINE; m++)
        {
            //
            endVec.push_back(j * MACHINE + m);         //1-index
            startVec.push_back(j * MACHINE + m);       //1-index
            timeVec.push_back(T[j * MACHINE + m - 1]); //0-index
        }
        endVec.push_back(JOB * MACHINE + 1);
    }
    //
    for (int m = 0; m < MACHINE; m++)
    {
        for (int j = 0; j < JOB; j++)
        {
            int job = s[m * JOB + j] - 1;
            // cout << "Job: " << job+1 << endl;
            // Vai na 'linha' de R referente ao job
            int start = job * MACHINE;
            int end = job * MACHINE + MACHINE - 1;
            //cout << "start = " << start << " end = " << end << endl;
            p = std::find(R + start, R + end, m + 1);

            //cout << "Maq: " << m+1 << endl;
            //Posicao
            int pos = p - (R + start);
            //cout << "Operacao: " << pos+1 << endl;
            int no = job * MACHINE + pos + 1;
            //cout << endl << "****No***** " << no << endl;
            //
            // Coloca na lista de arestas
            if (j == 0)
            {
                //seta apenas como n� inicial
                startVec.push_back(no);       //1-index
                timeVec.push_back(T[no - 1]); //0-index
            }
            if (j == JOB - 1)
            {
                //seta apenas como n� final
                endVec.push_back(no); //1-index
            }
            if (j != 0 && j != JOB - 1)
            {
                //seta como inicial e final
                startVec.push_back(no);       //1-index
                endVec.push_back(no);         //1-index
                timeVec.push_back(T[no - 1]); //0-index
            }
        }
    }
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    int numberVertex = MACHINE * JOB + 2; // number  of vertex
    int numberActivities = endVec.size(); // number of edges
    int j;

    /*vertices are of the form (v, u) */
    //indegree of each vertex (that is, count the number of edges entering them)
    vector<int> indegree(numberVertex, 0);

    std::vector<int> v(startVec);
    std::vector<int> u(endVec);
    std::vector<int> d(timeVec);

    //    for (j=0;j<numberActivities; j++)
    //        cout << "***** " << v[j] << "  " << u[j] << "  " << d[j] << endl;
    //    getchar();

    int project_duration = 0; //project duration

    /*# Compute the indegree for each vertex v from the graph:
    for each neighbor u of v: indegree[u] += 1*/
    for (j = 0; j < numberActivities; j++)
        indegree[u[j]]++;

    queue<int> Q; //queue Q = empty queue
    int distance[numberVertex];
    memset(distance, 0, sizeof(int) * numberVertex); //distance = array filled with zeroes

    //vetor de conjuntos de nos
    vector< set<int> > vecSetCriticalPath(numberVertex);

    for (j = 0; j < numberVertex; j++)
        if (indegree[j] == 0)
            Q.push(v[j]);

    int first;

    //printf ("first in the queue=%d\n", Q.front());

    /*for each neighbor u of v:
    d istance[u] = max(distance[u], distance[v] + time(v, u))
    indegree[u] -= 1
    if indegree[u] = 0:
    insert u on Q
    */
    while (!Q.empty())
    { //while Q is not empty:
        first = Q.front(); //v = get front element from Q
        Q.pop();           //delete de first from queue

        //Para todo vizinho de first faca
        for (int i = 0; i < numberActivities; i++) //**********MELHORAR ESSA ESTRUTURA **********
        {
            if (v[i] == first) //
            {
                distance[u[i]] = std::max(distance[u[i]], distance[v[i]] + d[i]);
                indegree[u[i]] -= 1;

                if (indegree[u[i]] == 0)
                    Q.push(u[i]);
            }
        }
        //for(int ii=0; ii<numberVertex; ii++)
        //cout << *std::max_element(distance,distance+numberVertex) << endl;
        //getchar();
        ////////////////////////////////////////////
    } //fecha while

    /*Now, select the vertex x with the largest distance.
    This is the minimum total project_duration.*/
    project_duration = *std::max_element(distance, distance + numberVertex);
    //printf ("Total Project Duration %d\n", project_duration);
    //getchar();
    ////////////////////////////////////////////////////////////////////////////
    // Pega posicao do elemento m�ximo = primeiro na 'lista critica'
    int criticalNode = std::find(distance, distance + numberVertex, project_duration) - distance;
    vector<int> criticalPath;
    //criticalPath.push_back(criticalNode);
    ////////////////////////////////////////////////////////////////////////////
    while (criticalNode != 0)
    {
        vector<int> neighborhood;
        vector<int> neighDistances;

        // varre recursivamente os vizinhos e pega o m�ximo
        for (int i = 0; i < numberActivities; i++) //**********MELHORAR ESSA ESTRUTURA **********
        {
            if (u[i] == criticalNode)
            {
                neighborhood.push_back(v[i]);
                neighDistances.push_back(distance[v[i]]);
                //cout << v[i] << " " << distance[v[i]] << endl;
                //getchar();
            }
        }
        int max = *std::max_element(neighDistances.begin(), neighDistances.end());
        int pp = std::find(neighDistances.begin(), neighDistances.end(), max) - neighDistances.begin();

        if (max == 0)
            pp = neighDistances.size() - 1;

        criticalNode = neighborhood[pp];
        criticalPath.push_back(criticalNode);
    }

    currentCriticalPath = criticalPath;

    return;
}

/*********************************************************************
* Seta semente
**********************************************************************/
void setSeed(vector< vector<int> > s)
{
    //cria estrutura da matriz s
    vector<int> ss(MACHINE * JOB);
    S[0] = ss;
    for (int i = 0; i < MACHINE; i++)
        for (int j = 0; j < JOB; j++)
        {
            S[0][i * JOB + j] = s[i][j];
        }
}

/*********************************************************************
* Seta o caminho critico
*********************************************************************/
void setCriticalPath(vector<int> s, int theBestScore)
{
    currentCriticalPath = s;
    bestCurrent = theBestScore;
}

//////////////////////////////////////////////////////////////////////
// Funcao objetivo do problema de otimiza��o continuo para encontrar
// o sequenciamento com atraso otimo
//////////////////////////////////////////////////////////////////////
float ObjectiveSeed(GAGenome &g)
{
    GARealGenome &genome = (GARealGenome &)g;
    /// /////////////////////////////////////////////////////////////////////////
    //cout << "************** USA AG PARA DETERMINAR O ATRASO OTIMO ******" << endl;
    /// /////////////////////////////////////////////////////////////////////////
    // vetor para a semente
    vector< vector<int> > S(MACHINE);
    // filas nas m�quinas
    vector< vector<int> > Queue(MACHINE);
    // tempos de processamento
    vector< vector<int> > processTime(MACHINE);
    // eventos de fim de atividade (idle time)
    vector<float> idleTime(MACHINE, 0);
    // contador de alocacao das maquinas
    vector<int> seized(MACHINE, 0);
    // vetor de pr�ximas m�quinas
    vector<int> nextStation(JOB);
    // vetor de eventos programados para m�quinas
    vector<int> numberOfEventsPerMachine(MACHINE, 0);
    // vetor de passos do job
    vector<int> step(JOB, 0);
    // STATUS DAS M�QUINAS
    vector<int> status(MACHINE, 0);
    // STATUS ANTERIOR DAS M�QUINAS
    vector<int> previousStatus(MACHINE, 0);
    // Rel�gio de simula��o
    float t = 0;
    // Lista de eventos futuros
    list<float> futureEvents(1, 0.0);
    // Contador de opera��es escalonadas
    int k = 0;
    //
    float fitness = 0;
    //
    int theJob;
    int time;
    /// ================= IMPLEMENTACAO =============== ///
    // INICIALIZA��O
    for (int i = 0; i < JOB; i++)
    {
        // VERIFICA A MAQUINA DO STEP 0
        int mq = R[i * MACHINE] - 1;
        // COLOCA O JOB NA FILA DA MAQ
        Queue[mq].push_back(i);
        // TEMPO
        processTime[mq].push_back(T[i * MACHINE]);
        // NEXT STATION
        nextStation[i] = R[i * MACHINE + 1] - 1;
    }

    //string ArqBsi = "Debug_objective_seed."+(string)PROBLEMA+".txt";
    //ofstream fout;
    //fout.open(ArqBsi.c_str(), ios::trunc);

    // ATUALIZA IDLE TIME, Lista eventos futuros e STATUS COM O ATRASO
    for (int mq = 0; mq < MACHINE; mq++)
    {
        float delay = genome.gene(mq * JOB + seized[mq]);
        idleTime[mq] = delay;
        //
        futureEvents.push_back(idleTime[mq]);
        // ATUALIZA status da maquina mq
        if (delay != 0.0)
            status[mq] = 2;
        else if (delay == 0.0)
            status[mq] = 0;
        //
        fitness = delay;
        //
        //fout << "t " << futureEvents.back() << " status " << status[mq] << endl;
    }
    //
    // faz enquanto n�o escalonar todos os jobs
    while (k < MACHINE * JOB)
    {
        /// ///////////////////
        /// // FASE A      ////
        /// ///////////////////
        ///
        //fout << "Fase A" << endl;
        ///
        futureEvents.sort();
        t = *futureEvents.begin();
        //
        // Remove primeiro elemento da lista de eventos futuros
        std::list<float>::iterator iter = futureEvents.erase(futureEvents.begin());

        //fout << "Lista de eventos futuros:" << endl;
        //for (std::list<float>::iterator it = futureEvents.begin(); it != futureEvents.end(); it++)
        //fout << *it << ' ';

        //fout << endl;

        //fout << "****************" << endl;
        //fout << "*** TNOW = " << t << endl;
        //fout << "****************" << endl;
        //
        // Em t=0 n�o h� atividade em andamento que possa se encerrar na fase B
        if (t > 0)
        {
            /// ///////////////////
            /// // FASE B      ////
            /// ///////////////////
            ///
            //fout << "Fase B" << endl;
            ///
            // Encontra atividade
            std::vector<float>::iterator it = find(idleTime.begin(), idleTime.end(), t);
            //
            if (it != idleTime.end())
            {
                // Identifica m�quina a ser liberada (do atendimento ou da espera)
                int mq = it - idleTime.begin();
                // Guarda status anterior
                previousStatus[mq] = status[mq];

                //fout << "Mq: " << mq << " status: " << status[mq] << " Fila: " << Queue[mq].size() << " idle: " << idleTime[mq] << endl;
                //
                ////////////////////////////////////////////////
                // Verifica se a m�quina estava em atendimento
                if (previousStatus[mq] == 1)
                {
                    // Identifica job
                    int job = S[mq].back() - 1; //S � 1-index
                    // atualiza step do job
                    step[job] = step[job] + 1;
                    //
                    // VERIFICA DELAY E COLOCA M�QUINA EM ESPERA
                    float delay = genome.gene(mq * JOB + seized[mq]);
                    if (delay != 0)
                    {
                        status[mq] = 2;
                        idleTime[mq] = t + delay;
                        // Lista eventos futuros
                        futureEvents.push_back(idleTime[mq]);
                    }
                    else
                    {
                        status[mq] = 0;
                        // seta como t-1 para que essa m�quina n�o seja pega novamente
                        idleTime[mq] = t - 1;
                    }
                    //
                    /// Verifica se job tem um pr�ximo passo
                    if (step[job] < MACHINE)
                    {
                        //Pega m�quina atual
                        int currentStation = nextStation[job];
                        // Coloca job na pr�xima fila j� definida
                        Queue[currentStation].push_back(job);
                        // Tempo de atendimento do job da fila
                        processTime[currentStation].push_back(T[job * MACHINE + step[job]]);
                        //
                    }
                    // Atualiza pr�xima m�quina do job (se houver)
                    if (step[job] < MACHINE - 1)
                    {
                        nextStation[job] = R[job * MACHINE + step[job] + 1] - 1;
                    }
                }
                else if (previousStatus[mq] == 2)
                {
                    // Libera a m�quina
                    status[mq] = 0;
                    idleTime[mq] = t - 1;
                }
            } // if(it != idleTime.end())
        }     //if(t>0)
        /// /////////////////////////////////////////////////////
        /// // FASE C: inicia as atividades por maquina      ////
        /// /////////////////////////////////////////////////////
        ///
        //fout << "Fase C" << endl;
        ///
        vector<int> numberInQueue;
        for (int j = 0; j < MACHINE; j++)
        {
            numberInQueue.push_back(Queue[j].size());
        }
        // Percorre as filas das maquinas
        for (int m = 0; m < MACHINE; m++)
        {
            // M�quina
            int mq = m;

            //fout << "Mq: " << mq << " status: " << status[mq] << " Fila: " << Queue[mq].size() << " idle: " << idleTime[mq] << endl;

            // Verifica se a maquina est� livre
            if (status[mq] == 0)
            {
                // N�mero de elementos na fila
                int nf = numberInQueue[mq];

                // H� fila?
                if (nf != 0)
                {
                    /// 'Zera' o tamanho da fila visitada (visita apenas uma vez)
                    numberInQueue[mq] = 0;
                    ///
                    float delay = genome.gene(mq * JOB + seized[mq]);
                    /// Verifica se m�quina estava em espera. Nesse caso usa LIFO
                    if (delay != 0)
                    {
                        /// Atende na base LIFO
                        time = processTime[mq].back();
                        theJob = Queue[mq].back();
                        /// Remove job da fila
                        Queue[mq].pop_back();
                        /// Remove tempo
                        processTime[mq].pop_back();
                    }
                    else if (delay == 0)
                    {
                        /// Atende na base FIFO
                        time = processTime[mq][0];
                        theJob = Queue[mq][0];
                        /// Remove job da fila
                        Queue[mq].erase(Queue[mq].begin());
                        /// Remove tempo
                        processTime[mq].erase(processTime[mq].begin());
                    }
                    // cout << "Inicia operacao do job " << theJob << " em mq: " << mq << ". Tempo: " << time << endl;
                    k++;
                    // cout << "Iteracao: " << k << endl;
                    //
                    /// /////////////////////////
                    //fout << ">>>> k = " << k << endl;
                    /// /////////////////////////
                    // Muda status
                    status[mq] = 1;
                    // Atualiza idleTime
                    //fout << ">>>>" << t << " * " << idleTime[mq] << " * " << t + time << endl;
                    idleTime[mq] = t + time;
                    // Atualiza makespan
                    if (idleTime[mq] > fitness)
                        fitness = idleTime[mq];
                    //
                    seized[mq] = seized[mq] + 1;
                    //
                    ////////////////////////////////////
                    //
                    // Lista eventos futuros
                    futureEvents.push_back(idleTime[mq]);
                    // Coloca job na semente (1-index)
                    S[mq].push_back(theJob + 1);
                    //
                    //fout << "Programa Job " << theJob << " na Maquina " << mq << endl;
                    //getchar();
                }
            }
        } //if m==0
    }     // FECHA WHILE
    //
    //cout << fitness << endl;
    ////////////////////////////////
    return fitness;
}

/******************************************************************************
* 	Function: Objective														  *
*	Short Description: Function called automatically by the GA, where the     *
*					   chromosome are evaluated (and can be manipulated as	  *
*					   well.
*******************************************************************************/
float Objective(GAGenome &g)
{
    GA2DBinaryStringGenome &genome = (GA2DBinaryStringGenome &)g;
    vector<int> score(nBestIndividuals);
    int P[MACHINE * JOB];

    //////////////////////////////////////////////////////////////
    // Na primeira gera��o apenas adiciona a semente na popula��o
    // A aptid�o j� foi calculada no setup ou na changeSeed
    //////////////////////////////////////////////////////////////
    if (primeira_geracao == true)
    {
        // Coloca individuo nulo na populacao
        genome.unset(0, 0, genome.width(), genome.height());

        primeira_geracao = false;

        // Retorna o valor de aptidao da melhor semente
        return bestCurrent;
    }
    else
    {
        /////////////////////////////////////////////
        /// TESTA O GENOMA
        /////////////////////////////////////////////
        //Permuta todas as sementes em S com base no genoma 'genome'
        vector<int> SS;

        //bestCurrent = ITER;
        for (int num = 0; num < nBestIndividuals; num++)
        {
            // Copia S para uma estrutura temporaria SS que sofre MUDAN�A
            SS = S[num];

            for (int i = 0; i < genome.width(); i++)
            {
                for (int j = 0; j < genome.height(); j++)
                {
                    if (genome.gene(i, j) == 1)
                    {
                        int idx = ((i + 1) % JOB);
                        int val = SS[(j * JOB) + i];
                        SS[(j * JOB) + i] = SS[(j * JOB) + idx];
                        SS[(j * JOB) + idx] = val;
                    }
                }
            }

            /////////////////////////////////////////////////////////////////////////////
            ///   I N T E R P R E T A C A O   O R I G I N A L   D A   S E M E N T E   ///
            /////////////////////////////////////////////////////////////////////////////
            int block = -1;
            for (int i = 0; i < MACHINE * JOB; i++)
            {
                int p = (i % JOB) + 1;
                if (p == 1)
                    block++;
                P[(block * JOB) + SS[i] - 1] = p;
            }

            // Apenas o menor score interessa.
            score[num] = factivel(P, R, 1);

            ///////////////////////////////////////////
            ///       COMPUTA      FACTIVEIS        ///
            ///////////////////////////////////////////
            numberOfEvaluations++;
            if (score[0] < ITER)
                numberOfFactibles++;
            ///////////////////////////////////////////

            sort(score.begin(), score.end());
            if (score[0] < bestCurrent)
                bestCurrent = score[0];
            ///////////////////////////////////////////////////
            // Imprime se melhorou
            ///////////////////////////////////////////////////
            if (score[0] <= score_ant)
            {
                ofstream fileOut;
                fileOut.open(ArqBsi.c_str(), ios::out);
                for (int seq = 0; seq < JOB * MACHINE; seq++)
                    fileOut << P[seq] << endl;

                fileOut << "\nMakespan da seq��ncia prioridades acima = " << bestCurrent << " (UT)" << endl;
                fileOut.close();
                score_ant = score[0];
                return (float)score[0];
            }
        } //fecha for
    }     //fecha else

    return (float)score[0];
}

/******************************************************************************
* 	Function: Local Search (multiple seeds)                                    *
*	Short Description:
*******************************************************************************/
void localSearch(const GAStatistics &g)
{
    // Estruturas temporarias para sementes
    vector< vector<int> > Stemp(nBestIndividuals);
    vector<int> SS;
    vector<int> bestSeed;
    int val;
    int P[MACHINE * JOB];
    int bgn, fnl, pos, job, opr, maq, bgn2, fnl2, pos2, job2, opr2, maq2;

    //Faz uma busla local em cada um dos melhores indiv�duos
    for (int n = 0; n < newBestIndividuals; n++)
    {
        if (n % 5 == 0)
        {
            //cout << "\nTopList " << n+1 << endl;
            // Pega o n-esimo individuo da lista dos melhores
            GA2DBinaryStringGenome &genome = (GA2DBinaryStringGenome &)g.bestIndividual(n);
            //bestCurrent = ITER;
            vector<int> score(newBestIndividuals);
            //cout << genome << endl;

            // Copia S para uma estrutura temporaria que sofrer� permuta��o
            SS = S[0];
            // Permuta os valores de SS baseado no cromossomo gerado pelo AG
            for (int i = 0; i < genome.width(); i++)
            {
                for (int j = 0; j < genome.height(); j++)
                {
                    if (genome.gene(i, j) == 1)
                    {
                        int idx = ((i + 1) % JOB);
                        val = SS[(j * JOB) + i];
                        SS[(j * JOB) + i] = SS[(j * JOB) + idx];
                        SS[(j * JOB) + idx] = val;
                    }
                }
            }
            //
            findCriticalPath(SS);
            if (n == 0)
                bestSeed = SS;
            //
            //        // Calcula score do individuo
            //        int block = -1;
            //        for(int i=0; i<MACHINE*JOB; i++)
            //        {
            //            int p = (i%JOB)+1;
            //            if(p == 1) block++;
            //            P[(block*JOB)+SS[i]-1] = p;
            //        }
            //        // Apenas o menor score interessa.
            //        float escore = factivel(P,R,1);
            //        cout << "Escore do individuo: " << escore << endl;
            //        getchar();

            ////////////////////////////////////////////////////////////////////////
            // INICIA BUSCA LOCAL NA SEMENTE PERMUTADA PELO INDIVIDUO
            // BUSCA LOCAL COM BASE NO CAMINHO CR�TICO DA SEMENTE PERMUTADA
            ////////////////////////////////////////////////////////////////////////
            for (unsigned int i = 0; i < currentCriticalPath.size() - 2; i++)
            {
                ////////////////////////////////////////////////////////////////////////////
                //identifica a i-�sima operacao do caminho critico
                job = (int)(currentCriticalPath[i] - 1) / MACHINE;         //0-index
                opr = (int)((currentCriticalPath[i] - 1) - job * MACHINE); //0-index
                // busca na matriz R quem � a maquina
                maq = R[job * MACHINE + opr]; //1-index
                // encontra posicao na semente
                bgn = (maq - 1) * JOB;
                fnl = bgn + (JOB - 1);
                pos = std::find(SS.begin() + bgn, SS.begin() + fnl, (job + 1)) - SS.begin();

                //cout << endl;
                //cout << currentCriticalPath[i] << "  " << job+1 << "  " << opr+1 << "  " << maq << endl;
                //cout << currentCriticalPath[i] << "  " << bgn << "  " << fnl << "  " << pos << endl;
                //getchar();
                ////////////////////////////////////////////////////////////////
                // identifica todas as outras operacoes cr�ticas da mesma maq
                for (unsigned int j = i + 1; j < currentCriticalPath.size() - 1; j++)
                {
                    //
                    job2 = (int)(currentCriticalPath[j] - 1) / MACHINE;          //0-index
                    opr2 = (int)((currentCriticalPath[j] - 1) - job2 * MACHINE); //0-index
                    // busca na matriz R quem � a maquina
                    maq2 = R[job2 * MACHINE + opr2]; //1-index
                    // encontra posicao na semente
                    bgn2 = (maq2 - 1) * JOB;
                    fnl2 = bgn2 + (JOB - 1);
                    pos2 = std::find(SS.begin() + bgn2, SS.begin() + fnl2, (job2 + 1)) - SS.begin();
                    //
                    if (maq2 == maq)
                    { //permuta
                        //
                        vector<int> newSS = SS;
                        std::vector<int> contOcurr(JOB, 0);
                        std::vector<int> contador(JOB, 0);

                        val = newSS[pos];
                        newSS[pos] = newSS[pos2];
                        newSS[pos2] = val;
                        //
                        //////////////////////////////////////////////////////////
                        // Transforma S em prioridades - interpretacao original //
                        //////////////////////////////////////////////////////////
                        int block = -1;
                        for (int i = 0; i < MACHINE * JOB; i++)
                        {
                            int p = (i % JOB) + 1;
                            if (p == 1)
                                block++;
                            P[(block * JOB) + newSS[i] - 1] = p;
                        }

                        /////////////////////////////////////////////////////////////
                        // Apenas o menor score interessa
                        /////////////////////////////////////////////////////////////
                        score[0] = factivel(P, R, 1);

                        ///////////////////////////////////////////
                        ///       COMPUTA      FACTIVEIS        ///
                        ///////////////////////////////////////////
                        numberOfEvaluations++;
                        if (score[0] < ITER)
                            numberOfFactibles++;
                        ///////////////////////////////////////////
                        //
                        int fitness = (int)score[0];
                        ////////gera_sequencia////////////////////////////////////////////
                        //cout << endl << projectDuration(SS) << " " << fitness << endl;
                        ////////////////////////////////////////////////////
                        if (fitness < bestCurrent)
                        {
                            bestCurrent = fitness;
                            //mudo a semente
                            bestSeed = newSS;
                            cout << "\nAtualiza solucao em Local Search: " << bestCurrent << endl;
                            //getchar();
                        }
                        //cout << endl;
                        //cout << currentCriticalPath[j] << "  " << job2+1 << "  " << opr2+1 << "  " << maq2 << endl;
                        //cout << currentCriticalPath[j] << "  " << bgn2 << "  " << fnl2 << "  " << pos2 << endl;
                        //getchar();

                    } //Fecha if maq == maq2
                }
            }
        }
    }

    findCriticalPath(S[0]);
    S[0] = bestSeed;
    return;
}

/******************************************************************************
* 	Function: Local Search	(single seed)									  *
*	Short Description:
*******************************************************************************/
void localSearch()
{
    findCriticalPath(S[0]);

    vector<int> score(currentCriticalPath.size());
    int val = 0;
    int maq, maq2, pos, pos2;
    vector<int> bestInitialSeed;

    int P[MACHINE * JOB];

    /////////////////////////////////////////////
    /// Permuta semente e calcula aptidao
    /////////////////////////////////////////////
    //Permuta todas as sementes em S
    vector<int> SS;
    //bestCurrent = ITER;
    bool melhora = true;
    while (melhora)
    {
        // Copia S para uma estrutura temporaria que sofrer� permuta��o
        SS = S[0];
        melhora = false;
        //
        //cout << "*" << currentCriticalPath.size() << endl;
        //getchar();
        for (unsigned int i = 0; i < currentCriticalPath.size() - 1; i++)
        {
            //identifica a i-�sima operacao do caminho critico
            int job = (int)(currentCriticalPath[i] - 1) / MACHINE;         //0-index
            int opr = (int)((currentCriticalPath[i] - 1) - job * MACHINE); //0-index
            // busca na matriz R quem � a maquina
            maq = R[job * MACHINE + opr]; //1-index
            // encontra posicao na semente
            int bgn = (maq - 1) * JOB;
            int fnl = bgn + (JOB - 1);
            pos = std::find(SS.begin() + bgn, SS.begin() + fnl, (job + 1)) - SS.begin();

            //cout << endl;
            //cout << currentCriticalPath[i] << "  " << job+1 << "  " << opr+1 << "  " << maq << endl;
            //cout << currentCriticalPath[i] << "  " << bgn << "  " << fnl << "  " << pos << endl;
            //getchar();
            ////////////////////////////////////////////////////////////////
            // identifica todas as outras operacoes cr�ticas da mesma maq
            for (unsigned int j = i + 1; j < currentCriticalPath.size() - 1; j++)
            {
                //
                int job2 = (int)(currentCriticalPath[j] - 1) / MACHINE;          //0-index
                int opr2 = (int)((currentCriticalPath[j] - 1) - job2 * MACHINE); //0-index
                // busca na matriz R quem � a maquina
                maq2 = R[job2 * MACHINE + opr2]; //1-index
                // encontra posicao na semente
                int bgn2 = (maq2 - 1) * JOB;
                int fnl2 = bgn2 + (JOB - 1);
                pos2 = std::find(SS.begin() + bgn2, SS.begin() + fnl2, (job2 + 1)) - SS.begin();

                if (maq == maq2)
                { //permuta
                    //
                    vector<int> newSS = SS;

                    val = newSS[pos];
                    newSS[pos] = newSS[pos2];
                    newSS[pos2] = val;
                    //
                    /////////////////////////////////////////////////////////////////////////////
                    ///   I N T E R P R E T A C A O   O R I G I N A L   D A   S E M E N T E   ///
                    /////////////////////////////////////////////////////////////////////////////
                    // Transforma S em prioridades
                    int block = -1;
                    for (int i = 0; i < MACHINE * JOB; i++)
                    {
                        int p = (i % JOB) + 1;
                        if (p == 1)
                            block++;
                        P[(block * JOB) + newSS[i] - 1] = p;
                    }
                    /////////////////////////////////////////////////////////////////////////////
                    // Apenas o menor score interessa.
                    score[0] = factivel(P, R, 1);
                    //
                    int fitness = score[0];
                    //cout << endl << scoreTopList[0] << " " << fitness << endl;
                    if (fitness < bestCurrent)
                    {
                        bestCurrent = fitness;
                        // Guardo solucao melhor para trocar a semente
                        bestInitialSeed = newSS;
                        //S[0] = newSS;
                        //findCriticalPath(S[0]);
                        melhora = true;
                        cout << "Melhora em Local Search: " << fitness << endl;
                        //getchar();
                    }
                }
            }
        }
        if (melhora == true)
        {
            cout << "****Atualiza Semente**** " << endl;
            S[0] = bestInitialSeed;
            findCriticalPath(S[0]);
        }
    } //fecha while
    return;
}

/******************************************************************************
* METODO DE BUSCA LOCAL BASEADO NO CAMINHO CRITICO. GENERALIZA O METODO VOID,
* POIS BUSCA NA VIZINHAN�A DE TODAS AS SOLU��ES MELHORES QUE A ORIGINAL
******************************************************************************/
void localSearchIn3phases()
{
    //cout << "****  EM LOCAL SEARCH *****" << endl;
    findCriticalPath(SS_rseed[0]);

    vector<int> score(currentCriticalPath.size());
    int val = 0;
    int maq, maq2, pos, pos2;
    vector<int> bestInitialSeed;

    int P[MACHINE * JOB];

    /////////////////////////////////////////////
    /// Permuta semente e calcula aptidao
    /////////////////////////////////////////////
    //Permuta todas as sementes em S
    vector<int> SS;
    //bestCurrent = ITER;
    bool melhora = true;
    while (melhora)
    {
        // Copia S para uma estrutura temporaria que sofrer� permuta��o
        SS = SS_rseed[0];
        melhora = false;

        for (unsigned int i = 0; i < currentCriticalPath.size() - 1; i++)
        {
            //identifica a i-�sima operacao do caminho critico
            int job = (int)(currentCriticalPath[i] - 1) / MACHINE;         //0-index
            int opr = (int)((currentCriticalPath[i] - 1) - job * MACHINE); //0-index
            // busca na matriz R quem � a maquina
            maq = R[job * MACHINE + opr]; //1-index
            // encontra posicao na semente
            int bgn = (maq - 1) * JOB;
            int fnl = bgn + (JOB - 1);
            pos = std::find(SS.begin() + bgn, SS.begin() + fnl, (job + 1)) - SS.begin();

            ////////////////////////////////////////////////////////////////
            // identifica todas as outras operacoes cr�ticas da mesma maq
            for (unsigned int j = i + 1; j < currentCriticalPath.size() - 1; j++)
            {
                int job2 = (int)(currentCriticalPath[j] - 1) / MACHINE;          //0-index
                int opr2 = (int)((currentCriticalPath[j] - 1) - job2 * MACHINE); //0-index
                // busca na matriz R quem � a maquina
                maq2 = R[job2 * MACHINE + opr2]; //1-index
                // encontra posicao na semente
                int bgn2 = (maq2 - 1) * JOB;
                int fnl2 = bgn2 + (JOB - 1);
                pos2 = std::find(SS.begin() + bgn2, SS.begin() + fnl2, (job2 + 1)) - SS.begin();

                if (maq == maq2)
                { //permuta
                    //
                    vector<int> newSS = SS;

                    val = newSS[pos];
                    newSS[pos] = newSS[pos2];
                    newSS[pos2] = val;

                    /////////////////////////////////////////////////////////////////////////////
                    ///   I N T E R P R E T A C A O   O R I G I N A L   D A   S E M E N T E   ///
                    /////////////////////////////////////////////////////////////////////////////
                    // Transforma S em prioridades
                    int block = -1;
                    for (int i = 0; i < MACHINE * JOB; i++)
                    {
                        int p = (i % JOB) + 1;
                        if (p == 1)
                            block++;
                        P[(block * JOB) + newSS[i] - 1] = p;
                    }
                    /////////////////////////////////////////////////////////////////////////////
                    // Apenas o menor score interessa.
                    score[0] = factivel(P, R, 1);

                    int fitness = score[0];
                    //cout << endl << scoreTopList[0] << " " << fitness << endl;
                    if (fitness < bestCurrent)
                    {
                        bestCurrent = fitness;
                        // Guardo solucao melhor para trocar a semente
                        bestInitialSeed = newSS;
                        //S[0] = newSS;
                        //findCriticalPath(S[0]);
                        melhora = true;
                        //cout << "Melhora em Local Search: " << fitness << endl;
                        //getchar();
                    }
                }
            }
        }
        if (melhora == true)
        {
            SS_rseed[0] = bestInitialSeed;
            findCriticalPath(SS_rseed[0]);
        }
    } //fecha while
    return;
}
