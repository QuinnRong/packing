/* 
   Packing of hard spheres via molecular dynamics
   Developed by Monica Skoge, 2006, Princeton University
   Contact: Aleksandar Donev (adonev@math.princeton.edu) with questions
   This code may be used, modified and distributed freely.
   Please cite:
   
   "Packing Hyperspheres in High-Dimensional Euclidean Spaces"
   	M. Skoge, A. Donev, F. H. Stillinger and S. Torquato, 2006
	
   if you use these codes.	
*/

//-----------------------------------------------------------------------------
// Box maker
//---------------------------------------------------------------------------

#ifndef  BOX_H
#define  BOX_H


#include <vector>
#include <math.h>

#include "vector.h"
#include "grid_field.h"
#include "event.h"
#include "sphere.h"
#include "heap.h"


#define PI     3.141592653589793238462643
#define SIZE   1.0            // size of box
#define VOLUMESPHERE pow(PI,((double)(DIM))/2.)/exp(lgamma(1+((double)(DIM))/2.)) // volume prefactor for sphere
#define DBL_EPSILON  2.2204460492503131e-016 // smallest # such that 1.0+DBL_EPSILON!=1.0
#define M 1.0


//---------------------------------------------------------------------------
// Class neighbor
//---------------------------------------------------------------------------
class neighbor
{
public:
    neighbor(int i_i);

    virtual void Operation(int j, vector<DIM, int>& pboffset) = 0;

    int i;
};


//---------------------------------------------------------------------------
// Class box
//---------------------------------------------------------------------------
class box
{
public:
    // constructor and destructor
    box(int N_i, double r_i, double growthrate_i, double maxpf_i);
    ~box();

    // Creating configurations
    int Optimalngrids(double maxpf);
    void CreateSpheres(double temp);
    void CreateSphere(int Ncurrent);   
    double Velocity(double temp);
    void VelocityGiver(double temp);  
    void SetInitialEvents();
    void RecreateSpheres(const char* filename, double temp);
    void ReadPositions(const char* filename);
    void AssignCells();

    void Process(int n);
    /**
     * repeat ProcessEvent() n times;
     */
    void ProcessEvent();
    /**
     * i = h.extractmax(), e = s[i].nextevent;
     * if e.j == INF, which means check:
     * s[i].nextevent = FindNextEvent(int i), downheap(1);
     * if 0 <= e.j < N, which means collision between i and j:
     * Collision(e) then give i and j checks;
     * else, which means transfer:
     * Transfer(e) then give i check;
     */
    event FindNextEvent(int i);
    /**
     * t = FindNextTransfer(i);
     * c = FindNextCollision(int i);
     * if c.time < t.time, CollisionChecker(c), return c;
     * else, return t;
     */
    void CollisionChecker(event c);
    /**
     * if c.j's next event is collision with other sphere, set it as check since collision between i and j happens before that;
     * set c.j's next event as c(c.time, j, i, c.v * (-1)), which is symmetric with c.i;
     */
    event FindNextTransfer(int i);
    /**
     * calculate which wall the sphere i will hit first
     * return:
     * an transfer event recording transfer time (ttime + gtime) and wall index
     */
    event FindNextCollision(int i);
    /**
     * collision cc(i, this);
     * vl(all -1s), vr(all 1s);
     * ForAllNeighbors(i, vl, vr, cc);
     * cc.cpartner == i  --> no collisions --> return an event with INF time and check;
     * else --> closest collision time cc.time with cc.cpartner, and stores its pboffset;
     */
    void ForAllNeighbors(int i, vector<DIM, int> vl, vector<DIM,int> vr, neighbor& operation);
    /**
     * for all neighbor cells of i ranging from vl(all -1s) to vr(all 1s), calculate pboffset;
     * for all spheres in the cell, operation.Operation(j, pboffset),
     * which actually is PredictCollision(i, j, pboffset, ctime, cpartner, cpartnerpboffset);
     * ctime initialized as dbINF and cpartner initialized as i;
     */
    void PredictCollision(int i, int j, vector<DIM, int> pboffset, double& ctime, int& cpartner, 
                            vector<DIM, int>& cpartnerpboffset);
    /**
     * if j != i, 
     * ctimej = CalculateCollision(i, j, pboffset.Double()) + gtime;
     * if (ctimej < ctime) && (ctimej < s[j].nextevent.time)
     * ctime = ctimej, cpartner = j; cpartnerpboffset = pboffset; 
     */
    double CalculateCollision(int i, int j, vector<DIM> pboffset);
    /**
     * input:
     * i, j: sphere ID
     * pboffset: offset of j because of periodic boundary in all directions
     * -1 means i & j are neighbors across left boundary;
     * 1 means i & j are neighbors across right boundary;
     * 0 means i & j are neighbors without crossing boundary.
     * return:
     * time of collision between i & j
     * INF if it will not happen.
     */
    double QuadraticFormula(double a, double b, double c);
    /**
     * solve ax^2 + 2bx + c = 0, where c > 0;
     * return:
     * x if a < 0 || b < 0;
     * INF otherwise.
     */
    void Collision(event e);
    /**
     * gtime = e.time;
     * update position and lutime;
     * update velocity according to Lubachevsky–Stillinger algorithm;
     * record momentum change;
     */
    void Transfer(event e);
    /**
     * gtime = e.time;
     * update position and lutime;
     * calculate new cell beased on e.j, and pdateCell(e.i, new cell);
     */
    void UpdateCell(int i, vector<DIM,int>& celli);
    /**
     * delete i from cell array at cell
     * add i to cell array at celli
     */

    void Synchronize(bool rescale);
    /**
     * updates position and lutime of all spheres to gtime;
     * r += gtime*growthrate;
     * gtime = 0;
     * can change growth rate and recale velocity;
     */

    // Debugging
    void TrackPositions();
    void OutputEvents();
    void OutputCells();
    void GetInfo();

    // Statistics
    double Energy();
    double PackingFraction();  
    void PrintStatistics();
    void RunTime();
    void WriteConfiguration(const char* wconfigfile);

    //variables
    const int N;                   // number of spheres
    int ngrids;                    // number of cells in one direction
    double maxpf;
    double growthrate;             // growthrate of the spheres
    double r;                      // radius, defined at gtime = 0
    double gtime;                  // this is global clock
    double rtime;                  // reset time, total time = rtime + gtime

    // statistics
    double pressure;               // pressure
    double xmomentum;              // exchanged momentum
    double pf;                     // packing fraction
    double energy;                 // kinetic energy
    double energychange;
    int ncollisions;               // number of collisions
    int ntransfers;                // number of transfers
    int nchecks;                   // number of checks
    int neventstot;                // total number of events 
    int ncycles;                   // counts # cycles for output

    time_t start, error, end;      // run time of program

    // arrays
    sphere *s;                      // array of spheres
    grid_field<DIM, int> cells; // array that keeps track of spheres in each cell
    int *binlist;                   // linked-list for cells array
    heap h;                         // event heap
    vector<DIM> *x;                 // positions of spheres.used for graphics
};


//---------------------------------------------------------------------------
// Predicts collisions, inherits neighbor operation
//---------------------------------------------------------------------------
class collision : public neighbor 
{
public:
    collision(int i_i, box *b);

    virtual void Operation(int j, vector<DIM, int>& pboffset);

    box *b; 
    double ctime;
    int cpartner;
    vector<DIM,int> cpartnerpboffset;
};

#endif 
