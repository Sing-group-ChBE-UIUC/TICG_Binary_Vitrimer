#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h> 
#include <errno.h>
#include <sys/stat.h> 
#include <stdbool.h>
#include <limits.h> 
#ifndef PARAM_NOTE
#define PARAM_NOTE ""
#endif
#define HAS_NOTE (sizeof(PARAM_NOTE) > 1)
/* 
------------- CHANGE LOG --------------------------------------------------------------------------- 
In Version 5.0, the following changes have been made:
- preallocate memory for MSD calculation to avoid dynamic allocation during BD loop. 
- Refactor force calculation to be pairwise, enforcing Newton's third law explicitly. 
- Dump more info to coord file for offline analysis. 
- Cap multiorigin MSD calculation to MAX_MSD_STARTS and MSD_CHECK_EVERY_BD to avoid excessive memory use. 
------------- HISTORY ------------------------------------------------------------------------------
In Version 4.9, the following changes have been made: 
- Add multiorigin and choose a subset for MSD calculation. 
- Add system drift correction in MSD calculation.
- Cap delta U in Metropolis test to MAX_DU to avoid exp overflow. 
- Refactor mean and sem calculation into a helper function. 
- New `wrap` and `mic` functions for safer PBC handling. 

In Version 4.8, the following changes have been made: 
- Calculating rolling average and standard error of degree of separation during the production run, and output to log file every DUMP_LOG_STRIDE steps.
- Adjust sampling density based on acceptance ratio of swap move. 

In Version 4.7, the following changes have been made:
- New feature: Calculate bond autocorrelation function (for dynamic bond)
  Each new start is dumped out, mean and sem are later calculated in post-processing python script. 

In Version 4.6, the following changes have been made:
- Add non-reactive ratio to control the fraction of chain ends that do not participate in bonding. 
  Take bound ratio and non-reactive ratio as input parameters. 
  The rest of the chain ends are "free" and "reactive" that can participate in swap moves. (BER) 

In Version 4.5, the following changes have been made:
- Do not skip swap move attempts if one is accepted. 
  Instead, loop over all neighbors for swap move within one MC attempt. 
  To avoid messing up with rate calculation, the acceptance ratio is calculated based on if the junction-end bond changes bonding partner at the end of the neighbor loop. 
  Note that the originally bound end is not listed as a candidate, so if no swap is accepted, the bond remains unchanged. 
- Rename counters to show if it's cumulative, window-based, or per-junction. 

In Version 4.4, the following changes have been made:
- Modify the function assignJuncFE() to be more random. 
- Add an equilibration loop before the production run.
- Loop over all neighbors for swap move within one MC attempt, instead of picking one randomly. 
  However, once one swap is accepted, the rest are skipped. 

In Version 4.3, the following changes have been made:
- Output fraction of all A and all B junctions to log file to monitor composition evolution. 

In Version 4.1, the following changes have been made:
- Refactor the parameter part of code to improve readability and maintainability. 

In Version 4.0, the following changes have been made: 
- Define Kuhn mass, Kuhn length, and reference temperature as the three base reference units. 
    \tau_0 = sqrt{m b^2 / k_B T} --- reference time unit
- Bug: Directly adding barrer to energy difference in Metropolis test kills the detailed balance. 
  Fix: Use an activation gate to thin the Markov chain. 
- Fix bug: type of junc.img should be int or long, not double. 

In Version 3.0, the following changes have been made: 
- Bug: in_contact array is stale inside the MC loop, causing incorrect swap rate and acceptance ratio calculation. 
    Fix: use the fresh retrun from `attempt_swap()` as the contact status; no need to check `in_contact` array. 
    Affected: `k_cond` and `acc_ratio_swap_cond` calculations.
- Define a tunable energy barrier for junction swap move to represent the activation energy of the bond exchange reaction. 
- Define spring constant for monomer-monomer bond and monomer-junction bond separately, such that the acceptance ratio of swap move given available neighbors is approximately 0.5. 
- Measure mean r^2 for monomer-junction bond and log it to file. 
- Add file output junction and free end coordinates. 
- Put file output strides into param file as well. 

In Version 2.0, the following changes have been made:
- Perform stride-convergence test
    - k_swap [1/BD time]
        - definition: 
            # of total accecpted swap / (total simulation time * number of junction) 
        - meaning: 
            The overall swap rate per junction = the macroscopic exchange rate 
    - k_cond [1/BD time]
        - definition: 
            # of accepted swap that occured while junction is in contact with free ends / total time duration of contact episodes 
        - meaning:
            Intrinsic swap rate while contact exists. 
            As MC stride decreases, k_cond should converge to a constant value. 
    - f_miss [#]
        Define \Delta t_MC being `SWAP_MOVE_STRIDE` * `DT`
        - definition: 
            (# of contact windows that has duration < \Delta t_MC) / # of contact windows 
        - meaning: 
            The fraction of contact episodes that are so short that they are completely unobservable by a stride-based algorithm. 
            A large f_miss indicates that the stride is too coarse. 

In Version 1.0, the following changes have been made: 
- Modified from ticg_1d_BD_BB_v9.c, now add volumeless junctions to the system. 
- Remove the reseeding of RNG in swap move to prevent correlation. 
*/ 

// ---------- constants ---------------------------------------------------------------------------- 
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NDIV (1+IMM1/NTAB)
#define STR_LEN 256
#define kB 1.380649e-23 // Boltzmann constant [Joule/Kelvin] 
#define CtoK 273.15 // zero Celsius temperature [Kelvin] 
#define N_AVO 6.02214076e23 // Avogadro's number [#] 
#define MAX_DU 700.0 // Threshold for Metropolis test to avoid exp overflow
#define MAX_CAND 128 // Maximum number of candidates for swap move 
#define MAX_ARM 10 
#define MAX_STARTS 1024 // Maximum number of start points to track for bond correlation function
#define MAX_MSD_STARTS 16
#define N_EQBR 3000 // Number of BD steps to run with Ea=0 and chiN=0 for equilibration 
#define SYSTEM "BV"
#define DIM 1 // dimensionality (1, 2, or 3) 
#define VERSION "5.0" 
// --- compile-time parameters ------------------- 
#define N_TYPE 2 
#define N_BEAD 10 
#define N_INT 15 
#define MG 64 
#define NGRID ((int) MG)
// --- constants --------------------------------- 
const double            ALPHA_NB = 0.25; 
const double            BETA_NOISE = 0.25; 
const double            ALPHA_BOND = 0.25; 
const double            fA = 0.5; 
const double            fB = 0.5;
// default values -------------------------------- 
static double           DT = 1e-3; // will be updated by recommend_dt() 
static int              REPLICA = 0; // default replica index for multiple runs 
// ----------------------------------------------------------------------------- 
const int DUMP_PHI_STRIDE        = (int) 5e3;
const int DUMP_COORD_STRIDE      = (int) 5e5;
const int DUMP_LOG_STRIDE        = (int) 5e3;
const int DUMP_AVG_SAMPLES       = (int) 5e4;
const int DUMP_DISTRIB_STRIDE    = (int) 5e3;
const int DUMP_WAIT_STRIDE       = (int) 1e4; 
const int DUMP_SITE_COORD_STRIDE = (int) 1e5; 
const long long nBDmax           = (long long) 1e8; 
static int NUM_BONDS             = 0; // will be updated after filling parameters  

// bond lifetime calculation 
static int CHECK_EVERY_BD        = 1; // how oftern to check the bond correlation after a start time 
static int MAX_LAG_BD            = (int) 1e6; // the largest time to cover for a given start time of bond correlation 
static long long MSD_MAX_LAG_BD  = (long long) 1e9; // the largest time to cover for a given start time of MSD calculation
static int NUM_CHECKS            = 0; 
static int MIN_ORG_GAP           = (int) 5e3; // minimum gap between two start times for bond correlation 
static int ORG_ACC_CNT_THR       = 10; // minimum number of accepted swaps for a start time to be considered 
static double ORG_AR_THR         = 0.1; // minimum acceptance ratio for a start time to be considered 
static int next_start_slot       = 0; 
static int last_start_nBD        = (int) -1e9; 
static int accepts_at_last_start = 0;
static int acf_max_i_seen = 0; // highest lag index that has any samples
// MSD calculation 
static int MSD_CHECK_EVERY_BD    = 1000;
static int MSD_ORIGIN_GAP        = 500; 
static int MSD_NUM_CHECKS        = 0; 
static int next_msd_start_slot   = 0;
static int last_msd_start_nBD    = (int) -1e9;
static int msd_max_i_seen      = 0; 
// ---------- data structures ---------------------------------------------------------------------- 
enum { 
    F_PARAM = 0, 
    F_PHI, 
    F_COORD, 
    F_PHI_AVG, 
    F_LOG, 
    F_DISTRIB,
    F_WAIT_TIME, 
    F_DYNSITE, 
    F_DEGSEP_AVG, 
    F_BOND_ACF, 
    F_MSD, 
    NUM_STR
};
enum { // for command-line input 
    chi0N_idx = 1,
    mc_stride_idx,
    Ea_idx,
    arm_idx,
    bound_ratio_idx,
    inert_ratio_idx,
    replica_idx,
    NUM_CLI
}; 
typedef struct {
    double r[N_BEAD]; // bead positions along x (1‑D)
    double f[N_BEAD]; // force along x (1‑D)
    long   img[N_BEAD]; // cumulative number of box crossings, recorded for unwrapping positions to compute MSD 
    int    type[N_BEAD]; // 0 = A, 1 = B 
} chain;
typedef struct {
    double r; 
    double f;
    long   img;  
    int    type; 
} junction; 

typedef struct {
    const char *prefix;
    const char *ext;
} FileSpec;

typedef struct {
    int active; 
    int start_nBD; 
    int next_check_index; 
    int *bond_at_start; // copy of BonJ at the start time 
    unsigned int *change_at_start; // copy of arm_change_count at the start time 
    long long id; 
} StartPoint; 

typedef struct {
    int active; 
    int start_nBD; 
    int next_check_index;

    // segmental 
    int n_beads; 
    double *rBead0_O;

    // chain COM 
    double *comChain0_O; 

    // system COM 
    double comSys0_O; 

    long long id; 
} MSDOrigin; 

typedef struct {
    // --- model knobs ---------------------------
    int    sqrtNbar; 
    double kappaN; 
    double chi0N; 
    int    mc_stride;
    double bound_ratio; 
    double inert_ratio; 
    double inert_fracA; 
    int    inert_random; // 1: randomly pick ENDs to be inert; 0: pick inert CHAINS with `inert_farcA` of A type 
    int    arm; 

    // --- reference units ----------------------- 
    double b_ref, m_ref, Tref_K, kBTref, tau_ref, gamma_ref;

    // --- geometry/size (reduced, in [b]) ------- 
    double B2_tilde; 
    double Re_b; 
    double box_len_b; 
    int    N_chain_tot; 
    int    N_J; 
    int    N_NR; // non-reactive chain ends 
    int    N_FE; // free, reactive chain ends 

    // density cloud (quasi-1D) ------------------ 
    double rho0_per_Re;         // [#/Re]
    double deltaL_b;            // cloud size in [b]
    double Cw_cloud;            // 1D normalization
    double DPHI;                // density increment per bead for depositing on grid 
    double R_swap_b;            // reaction radius
    double L_cell_sw_b, L_cell_nb_b;
    int    N_cell_sw, N_cell_nb;

    // --- thermodynamics ------------------------ 
    double Tphys_K;             // current T
    double T_tilde;             // T/Tref
    double kBT_tilde;           // == T_tilde
    double gamma_tilde;         // friction in gamma_ref units 

    // springs
    double ks_mon_tilde;
    double ks_junc_tilde;

    // barrier
    double Ea_kJmol; 
    double Ea_tilde;            // Ea/(kBTref)

    // misc
    double fA, fB;              // compositions 
} SimParams;

// Global objects -------------------------------------------------------------- 
static SimParams P; 
// Global pointers ------------------------------------------------------------- 
static chain *poly; 
static junction *junc; 
// bookkeeping arrays 
static int *BonJ; 
static int *FEnd; 
static int *NREnd; 
// linked-cell list 
static int *cellHead_SW; 
static int *cellNext_SW;
static int *cellHead_NB; 
static int *cellNext_NB;
// metrics for BD-MC cadence diagnostics 
static long long     *swaps_J;            // per-junction accepted swaps
static char          *in_contact;         // bool per junction (0/1)
static double        *contact_start;      // time when current contact window opened
static double        *contact_time;       // accumulated contact time per junction
static long long     *windows_total;      // number of contact windows per junction
static long long     *windows_missed;     // windows with duration < Δt_MC
static double        *wait_sum;           // sum of waits per junction
static double        *wait_sumsq;         // sum of waits^2 per junction
static unsigned long *wait_n;             // number of waits per junction
static double        *last_swap_time;     // last swap time per junction (−1 if none yet)
// Bond correlation 
static unsigned int *arm_change_count;    // count whenever the bonding arm of a junction changes partner (used to block double counting if swap and reform happen in the same check interval) 
static StartPoint starts[MAX_STARTS]; 
static double *acf_sum   = NULL;          // sum of C over starts, per lag i
static double *acf_sumsq = NULL;          // sum of C^2 over starts, per lag i (for SEM)
static unsigned long *acf_count = NULL;   // number of starts contributing to lag i
// MSD 
static MSDOrigin msd_starts[MAX_MSD_STARTS];
static double *msd_seg_sum   = NULL; 
static double *msd_seg_sumsq = NULL;
static double *msd_com_sum   = NULL;
static double *msd_com_sumsq = NULL;
static unsigned long *msd_seg_cnt = NULL;
static unsigned long *msd_com_cnt = NULL;
static double Rcom0_sysRef = 0.0; 
// MSD pools
static double *msd_r0_pool = NULL; // size = MAX_MSD_STARTS * P.N_chain_tot * N_BEAD 
static double *msd_com0_pool = NULL; // size = MAX_MSD_STARTS * P.N_chain_tot 
// ----------------------------------------------------------------------------- 
// indices and weights for particle in mesh interpolation (only for outputting density fields) 
typedef struct {
    int    idx[2]; 
    double w  [2]; 
} gridPointInterpolator;

// ---------- helper functions --------------------------------------------------------------------- 
void fill_params(SimParams *P, 
                double chi0_from_cli, int mc_stride_from_cli, double Ea_tilde_from_cli, int arm_from_cli, 
                double bound_ratio_from_cli, double inert_ratio_from_cli) { 
    
    /*
    Extra notes 

    ***** sqrtNbar ***** 
    \sqrt{\bar{N}} is the invariant degree of polymerization, defined as 
            \sqrt{\bar{N}} = \rho_0 R_e^3 / N
    where \rho_0 is the bead number density, R_e is the end-to-end distance of a chain, and N is the degree of polymerization. 

    \rho_0 = nN / V; V = L^3 
    \sqrt{\bar{N}} = \rho_0 R_e^3 / N
                = (nN / V) (R_e^3 / N)
                = nN / (L^3 / R_e^3) / N 
                = n / (L^3 / R_e^3)
                = n / BOX_LENGTH^3
    where BOX_LENGTH is the length of the simulation box in units of R_e. 

    ***** n_int *****
    n_int = c_w N \sqrt{\bar{N}} * (a_int / R_e)^3
    where c_w is a cloud shape dependent constant, N is the degree of polymerisation, \sqrt{\bar{N}} is the invariant degree of polymerisation, a_int is the density cloud size, and R_e is the end-to-end distance of a chain. 

    ***** deltaL ***** 
    substitute \sqrt{\bar{N}} = \rho_0 R_e^3 / N in n_int: 
    n_int = c_w \rho_0 a_int^3 = \rho_0 a_int^3 (let c_w = 1) 

    therefore, \Delta L = a_int = (n_int / \rho_0)^{1/3}
    or more generally, \Delta L = a_int = (n_int / \rho_0)^{1/d}, where d is the dimensionality of the system. 

    ***** Cw_cloud ***** 
    Note that for this quasi-1D system, rho_0 is in units of [R_e^-1], and \Delta L is in units of [R_e]. 
    For 3D system, rho_0 would be in units of [R_e^-3] and C_w_CLOUD should multiply by \Delta L^3.
    */

    // knobs ------------------------------------- 
    P->sqrtNbar     = 32;
    P->kappaN       = 50.0;
    P->chi0N        = chi0_from_cli;
    P->mc_stride    = mc_stride_from_cli;
    P->bound_ratio  = bound_ratio_from_cli;
    P->inert_ratio  = inert_ratio_from_cli;
    P->inert_fracA  = 0.5; // fraction of A type within inert chains 
    P->inert_random = 0; 
    P->arm          = arm_from_cli;
    P->fA = 0.5; 
    P->fB = 0.5;

    // reference units (SI units) ----------------
    P->b_ref       = 7.5e-10; // PEO statistical segment length [m]; estimated from reported persistent length of 3.7 [A] 
    double M0_gmol = 44.05; // PEO monomer weight [g/mol] 
    double monomer_len = 2.8e-10; // estimated monomer length [m] 
    double m_gmol = M0_gmol * (P->b_ref / monomer_len); // Kuhn mass of PEO [g/mol] 

    P->m_ref       = m_gmol * 1e-3 / N_AVO; // reference mass [kg] 
    P->Tref_K      = 70.0 + CtoK; // reference temperature [K]
    P->kBTref      = kB * P->Tref_K; // reference thermal energy [J]
    P->tau_ref     = sqrt(P->m_ref * P->b_ref * P->b_ref / P->kBTref); // derived reference time [s] 
    P->gamma_ref   = P->kBTref * P->tau_ref / (P->b_ref * P->b_ref); // derived reference friction coefficient [kg/s] 

    // geometry (reduced, in [b]) ---------------- 
    P->B2_tilde    = 1.0;
    P->Re_b        = sqrt((double)N_BEAD * P->B2_tilde);
    P->box_len_b   = 15.0 * P->Re_b;
    P->N_chain_tot = (int)(1.0 * 1.0 * P->box_len_b / P->Re_b * P->sqrtNbar); // total number of chains in the system; assuming 1 [Re] in both y and z directions
    P->N_J         = (int)(P->bound_ratio * (double)P->N_chain_tot * 2.0 / (double)P->arm); 
    P->N_NR        = (int)(P->inert_ratio * (double)P->N_chain_tot * 2.0); if (P->N_NR % 2 != 0) P->N_NR += 1;
    P->N_FE        = (int)(2 * P->N_chain_tot - P->N_J * P->arm - P->N_NR); 

    P->rho0_per_Re = (double)(P->N_chain_tot * N_BEAD) / (1.0 * 1.0 * P->box_len_b / P->Re_b); // y and z dimensions are assumed to be 1.0 [Re]; thus rho_0 [#/R_e^3]; Unit in this quasi-1D system is [#/R_e] 
    P->deltaL_b    = ((double)N_INT / P->rho0_per_Re) * P->Re_b; // density cloud size [b] ( = a_int = \Delta_L in paper) NOTE: inverse rho_0 firstly gives [Re]; it is then converted to [b] by multiplying by RE
    P->Cw_cloud    = 1.0 / (P->rho0_per_Re * (P->deltaL_b / P->Re_b)); // normalization factor for quasi-1D density cloud volume, \rho_0 [#/R_e] \times \Delta L [R_e] 
    P->DPHI        = (double)MG / (double)(P->N_chain_tot * N_BEAD); // density increment per bead 

    P->R_swap_b    = 0.5 * P->Re_b; // swap move reaction region radius [b]
    P->L_cell_sw_b = P->R_swap_b; 
    P->N_cell_sw   = (int)ceil(P->box_len_b / P->L_cell_sw_b);

    P->L_cell_nb_b = 1.0 * P->deltaL_b;
    P->N_cell_nb   = (int)ceil(P->box_len_b / P->L_cell_nb_b);

    // thermodynamics ---------------------------- 
    P->Tphys_K     = P->Tref_K; // simulation temperature [K] 
    P->T_tilde     = P->Tphys_K / P->Tref_K; // reduced simulation temperature [Tref] 
    P->kBT_tilde   = P->T_tilde; // reduced thermal energy [kBTref]
    P->gamma_tilde = 1.0; // reduced friction [gamma_ref] 

    // springs
    P->ks_mon_tilde  = P->kBT_tilde * (double)DIM / P->B2_tilde; // bond spring constant for monomer-monomer bond [kBTref/b_ref^2] 
    P->ks_junc_tilde = 1.0 * P->ks_mon_tilde; // monomer-junction bond is set stronger for tuning the swap acceptance ratio

    // barrier
    P->Ea_tilde = Ea_tilde_from_cli;
    P->Ea_kJmol = P->Ea_tilde * N_AVO * P->kBTref / 1000.0; // back-calculate Ea in [kJ/mol] for logging purpose
}

// wrap coordinates back into the box 
static inline void wrap(double *x) {
    double L = P.box_len_b;
    double y = fmod(*x, L);
    if (y < 0.0) y += L;
    *x = y;
}

// minimum‑image convention for periodic boundary conditions 
static inline double mic(double d) {
    return d - P.box_len_b * nearbyint(d / P.box_len_b);
}

static inline void globalToChain(int global_bead_idx, int *chain_idx, int *monomer_idx) {
    *chain_idx = global_bead_idx / N_BEAD; // chain index
    *monomer_idx = global_bead_idx - (*chain_idx) * N_BEAD; // monomer index within the chain
}

static inline int cellIndex_1D(double x, int Ncell, double Lcell) { 
    int cx = (int) floor(x / Lcell);
    cx = (cx % Ncell + Ncell) % Ncell; // wrap around
    return cx;
}

static inline int end_flag(int m) {
    return (m == 0) ? 1 : ((m == N_BEAD - 1) ? 2 : 0);
}

// Build local bead->junction maps on demand 
static void build_bead_to_junc_maps(int *bead_to_J, short *bead_to_arm) {
    int nTot = P.N_chain_tot * N_BEAD;
    for (int g = 0; g < nTot; ++g) { bead_to_J[g] = -1; bead_to_arm[g] = -1; }
    for (int J = 0; J < P.N_J; ++J) {
        for (int a = 0; a < P.arm; ++a) {
            int g = BonJ[P.arm * J + a];
            bead_to_J[g]  = J;
            bead_to_arm[g]= (short)a;
        }
    }
}

// reuse COM computation 
static inline void snapshot_COMs(double *RcomNow_sys, double *com_chain_now) {
    int g = 0;
    double Rsum = 0.0;
    for (int c = 0; c < P.N_chain_tot; ++c) {
        double sum = 0.0;
        for (int m = 0; m < N_BEAD; ++m, ++g) {
            sum += poly[c].r[m] + poly[c].img[m] * P.box_len_b; // unwrapped
        }
        com_chain_now[c] = sum / (double)N_BEAD;
        Rsum += com_chain_now[c];
    }
    *RcomNow_sys = Rsum / (double)P.N_chain_tot;
}

// Given accumulators S = sum(x), S2 = sum(x^2), N = #samples, report mean and SEM 
static inline void mean_sem_from_sums(double S, double S2, double N,
                                      double *mean_out, double *sem_out) {
    if (N <= 0.0) { *mean_out = NAN; *sem_out = NAN; return; }
    double mu  = S / N;
    double var = (N > 1.0) ? (S2 - (S* S)/N) / (N - 1.0) : 0.0;
    if (var < 0.0) var = 0.0;                 // numeric guard
    double sem = (N > 1.0) ? sqrt(var / N) : NAN;
    *mean_out = mu;
    *sem_out  = sem;
}

// build gridPointInterpolator for point x in a grid of size MG with grid spacing BOX_LENGTH/MG 
static gridPointInterpolator mkInterp(double x) {
    double h = P.box_len_b / MG;
    double gx = x / h;  
    int ix0 = (int) floor(gx);
    // fractions 
    double fx = gx - ix0; 
    // wrap 
    ix0 = (ix0 % MG + MG) % MG;
    int ix1 = (ix0 + 1) % MG;

    gridPointInterpolator s;
    // 1D bilinear indices and weights
    s.idx[0] = ix0; s.w[0] = (1.0 - fx);
    s.idx[1] = ix1; s.w[1] = fx;
    // Ensure indices are non-negative
    for (int i = 0; i < 2; ++i) {
        if (s.idx[i] < 0 || s.idx[i] >= NGRID) {
            fprintf(stderr, "Invalid grid index: %d at i=%d\n", s.idx[i], i);
            exit(EXIT_FAILURE);
        }
    }
    return s;
}

static void depositPhi(double *Phi, const gridPointInterpolator *s, double dphi, int sign) {
    // Deposit or remove dphi onto eight grid points (stencil) in the gridPointInterpolator
    for (int i = 0; i < 2; ++i) {
        Phi[s->idx[i]] += sign * dphi * s->w[i];
    }
}

// squared distance under PBC
static double rsq(double x1, double x2, double *dx_out) {
    double dx = mic(x1 - x2); 
    if (dx_out) *dx_out = dx; 
    return dx * dx; 
}

// bonded force on bead i of chain c
static double bonded_force_one_bead(int c, int i, double ks_mon_tilde, double ks_junc_tilde) {
    // force = negative gradient of energy
    const chain *ch = &poly[c]; 
    double f = 0.0; 
    if (i > 0) {
        double dx; rsq(ch->r[i], ch->r[i-1], &dx);
        f += -ks_mon_tilde * dx; // force on bead i from previous bead
    }
    if (i < N_BEAD-1) {
        double dx; rsq(ch->r[i], ch->r[i+1], &dx);
        f += -ks_mon_tilde * dx; // force on bead i from next bead
    }
    if (i == 0 || i == N_BEAD-1) {
        int g = c * N_BEAD + i; // global bead index 
        for (int j = 0; j < P.arm * P.N_J; ++j) {
            if (BonJ[j] == g) { // this bead is bonded to junction j/ARM
                int J = j / P.arm; // junction index
                double dx; rsq(ch->r[i], junc[J].r, &dx);
                f += -ks_junc_tilde * dx; // force on bead i from jth arm of junction J 
                break;
            }
        }
    }
/*
Note that the bonded potential (Gaussian-spring model) has a hidden assumption that the harmonic bond has a rest length of 0. 
(U_bond minimum at r_ij = 0) 
This means when springs get compressed, there is no push-back force. Repulsion comes from other interactions, such as non-bonded forces. 
Only when the position of bead `i` is stretched away from its neighbor, x_i > x_{i+1}, then `i+1` "pulls" `i` "left", back towards `i-1`. 
*/
    return f;
}

static double bonded_force_one_junc(int J, double ks_junc) {
    double f = 0.0; 
    for (int j = P.arm * J; j < P.arm * J + P.arm; ++j) {
        int b = BonJ[j]; // global bead index of junction-bound chain end bead 
        int c, m; globalToChain(b, &c, &m); 
        double dx; rsq(junc[J].r, poly[c].r[m], &dx);
        f += -ks_junc * dx; // force on junction J from all ARM bound beads 
    }
    return f;
}

// overlap of two cubic clouds along one dimension 
static inline double overlap_1d(double d) {
    d = fabs(d);
    if (d >= P.deltaL_b) return 0.0;
    double overlap_len = P.deltaL_b - d; // linear overlap length
    return overlap_len * P.Cw_cloud * P.Cw_cloud; // normalize 
}

static inline double dmin(double a, double b) { return (a < b) ? a : b; }
static inline double dmax(double a, double b) { return (a > b) ? a : b; }

static double recommend_dt(double ks_junc_tilde, double gamma_tilde, double kBT_tilde, 
                           double alpha_nb, double beta_noise, double alpha_bond,
                           double *out_dt_nb, double *out_dt_bond, double *out_dt_noise) {
    /*
    One BD step moves a bead by \Delta x_{det} = \frac{\Delta t}{\gamma}|\mathbf{F}| 
    We want this to be a small fraction `ALPHA` of the smallest feature of the potential, which is DEL`TA_L.
    Therefore, we set:
        \frac{\Delta t}{\gamma}|\mathbf{F}| \leq \alpha \Delta L
    =>  \Delta t \leq \alpha\frac{\Delta L \gamma}{|\mathbf{F}|}
    where |\mathbf{F}| is estimated using the target interaction number N_INT.     
    */                        
    double pref = ((double)P.sqrtNbar / P.Re_b); // prefactor for non-bonded interaction      
    double C2   = P.Cw_cloud * P.Cw_cloud;

    double base = P.kappaN + P.chi0N * (2.0 * fA * fB); // mean unlike contribution
    double Feff_nb = kBT_tilde * pref * base * C2 * N_INT; // estimate of nonbonded force magnitude
    const double dt_nb = alpha_nb * P.deltaL_b * gamma_tilde / dmax(Feff_nb, EPS); 

    // Bond relaxation time tau = gamma / k_spring
    // To be more conservative, we use the stiffer spring constant ks_junc
    const double tau_bond = gamma_tilde / dmax(ks_junc_tilde, EPS);
    const double dt_bond  = alpha_bond * tau_bond;
    
    // Noise drift = sqrt(2 kT DT / GAMMA) <= beta_noise * DELTA_L 
    const double dt_noise = (beta_noise*beta_noise) * gamma_tilde * P.deltaL_b * P.deltaL_b / (2.0 * kBT_tilde);

    if (out_dt_nb)    *out_dt_nb    = dt_nb;
    if (out_dt_bond)  *out_dt_bond  = dt_bond;
    if (out_dt_noise) *out_dt_noise = dt_noise;

    const double dt_rec = dmin(dt_nb, dmin(dt_bond, dt_noise)); // take the most conservative 

    return dt_rec;
}

// pair potential energy 
static double pair_energy(double xi, double xj,
                          int type_i, int type_j) {

    double dx; 
    rsq(xi, xj, &dx); 

    // volume overlap for cubic clouds 
    double v_overlap = overlap_1d(dx); 
    
    // convert overlap to local‐density contribution 
    double e = 0.0;
    // incompatibility (chi term) 
    if (type_i != type_j) {
        e += P.chi0N * v_overlap; 
    } // else no penalty for same type 

    // incompressibility (kappa term) 
    e += P.kappaN * v_overlap; 
    return e;
}

// // non-bonded force on bead i of chain c 
// static double nonbonded_force_one_bead(int c_idx, int i_idx, double kbT_tilde)
// {
//     double xi = poly[c_idx].r[i_idx];
//     int  ti =  poly[c_idx].type[i_idx];

//     const double pref = (double)P.sqrtNbar / (double)P.Re_b; // prefactor for non-bonded interaction 

//     double fi = 0.0;

//     int cx0 = cellIndex_1D(xi, P.N_cell_nb, P.L_cell_nb_b);

//     for (int dcx = -1; dcx <= 1; ++dcx) { // neighboring cells 
//         int cx = (cx0 + dcx + P.N_cell_nb) % P.N_cell_nb; // wrap around

//         // traverse the linked-cell list in cell cx
//         for (int g = cellHead_NB[cx]; g != -1; g = cellNext_NB[g]) { 
//             int c, m; globalToChain(g, &c, &m); // get chain and monomer indices
//             if (c == c_idx && m == i_idx) continue; // skip self 

//             double xj = poly[c].r[m];
//             int    tj = poly[c].type[m];

//             double dx; rsq(xi, xj, &dx);
//             double base = ((ti != tj) ? P.chi0N : 0.0) + P.kappaN; 
//             double coeff = kbT_tilde * pref * base * P.Cw_cloud * P.Cw_cloud; // 2: sign, C^2 

//             if (fabs(dx) < P.deltaL_b) { // within interaction range 
//                 fi += coeff * (dx > 0 ? 1.0 : -1.0); // force direction 
//             }
//         }
//     }
//     return fi;
// }

// // wrapper
// static void compute_conservative_forces(double ks_mon_tilde, double ks_junc_tilde, double kbT_tilde) {
//     for (int c = 0; c < P.N_chain_tot; ++c) {
//         for (int i = 0; i < N_BEAD; ++i) {
//             double fb = bonded_force_one_bead(c, i, ks_mon_tilde, ks_junc_tilde);
//             double fn = nonbonded_force_one_bead(c, i, kbT_tilde);
//             poly[c].f[i] = fb + fn;
//         }
//     } 

//     for (int J = 0; J < P.N_J; ++J) {
//         junc[J].f = bonded_force_one_junc(J, ks_junc_tilde); // already considered sum over three arms 
//     }
// }

static void compute_conservative_forces_pairwise(double ks_mon_tilde, double ks_junc_tilde, double kbT_tilde) {
    // zero forces
    for (int c = 0; c < P.N_chain_tot; ++c)
        for (int i = 0; i < N_BEAD; ++i)
            poly[c].f[i] = 0.0;
    for (int J = 0; J < P.N_J; ++J)
        junc[J].f = 0.0;

    // bonded polymer-polymer & polymer-junction 
    for (int c = 0; c < P.N_chain_tot; ++c)
        for (int i = 0; i < N_BEAD; ++i)
            poly[c].f[i] += bonded_force_one_bead(c, i, ks_mon_tilde, ks_junc_tilde);

    for (int J = 0; J < P.N_J; ++J)
        junc[J].f += bonded_force_one_junc(J, ks_junc_tilde);

    // symmetric nonbonded over cell pairs (1D)
    const double pref   = (double)P.sqrtNbar / (double)P.Re_b;
    const double C2     = P.Cw_cloud * P.Cw_cloud;
    const double coeff0 = kbT_tilde * pref * C2;

    // Iterate each cell and its right neighboring cell 
    for (int cx = 0; cx < P.N_cell_nb; ++cx) {
        for (int g1 = cellHead_NB[cx]; g1 != -1; g1 = cellNext_NB[g1]) { // traverse the linked-cell list to get all beads in cell cx 
            // for each bead `g1` in cell `cx` 
            int c1, m1; globalToChain(g1, &c1, &m1);
            double x1 = poly[c1].r[m1];
            int    t1 = poly[c1].type[m1];

            // find another bead `g2` in the same cell and calculate pairwise non-bonded force 
            for (int g2 = cellNext_NB[g1]; g2 != -1; g2 = cellNext_NB[g2]) {
                // cellNext_NB[g1] means starting g2 after g1, to avoid double counting 
                int c2, m2; globalToChain(g2, &c2, &m2);
                double x2 = poly[c2].r[m2];
                int    t2 = poly[c2].type[m2];

                double dx; rsq(x1, x2, &dx);
                double adx = fabs(dx);
                if (adx >= P.deltaL_b) continue;

                double base  = ((t1 != t2) ? P.chi0N : 0.0) + P.kappaN;
                double fi    = coeff0 * base * (dx >= 0.0 ? 1.0 : -1.0);

                poly[c1].f[m1] += fi;
                poly[c2].f[m2] -= fi;
            }

            // also find beads `g2` in the right neighboring cell 
            int cx2 = (cx + 1) % P.N_cell_nb;
            for (int g2 = cellHead_NB[cx2]; g2 != -1; g2 = cellNext_NB[g2]) {
                int c2, m2; globalToChain(g2, &c2, &m2);
                double x2 = poly[c2].r[m2];
                int    t2 = poly[c2].type[m2];

                double dx; rsq(x1, x2, &dx);
                double adx = fabs(dx);
                if (adx >= P.deltaL_b) continue;

                double base  = ((t1 != t2) ? P.chi0N : 0.0) + P.kappaN;
                double fi    = coeff0 * base * (dx >= 0.0 ? 1.0 : -1.0);

                poly[c1].f[m1] += fi;
                poly[c2].f[m2] -= fi;
            }
        }
    }
}

// initialize random number generator
long initRan()
{
    unsigned long a = clock();
    unsigned long b = time(NULL);
    unsigned long c = getpid();
    a = a - b;  a = a - c;  a = a ^ (c >> 13);
    b = b - c;  b = b - a;  b = b ^ (a << 8);
    c = c - a;  c = c - b;  c = c ^ (b >> 13);
    a = a - b;  a = a - c;  a = a ^ (c >> 12);
    b = b - c;  b = b - a;  b = b ^ (a << 16);
    c = c - a;  c = c - b;  c = c ^ (b >> 5);
    a = a - b;  a = a - c;  a = a ^ (c >> 3);
    b = b - c;  b = b - a;  b = b ^ (a << 10);
    c = c - a;  c = c - b;  c = c ^ (b >> 15);
    return c % 1000000000;
}

// generate uniformly distributed random number between 0 and 1  
float ran1(long *idum)
{
    int j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;
    
    if (*idum <= 0)
    {
        if (-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);
        idum2 = (*idum);
        for (j = NTAB + 7; j >= 0; --j)
        {
            k = (*idum) / IQ1;
            *idum = IA1 * (*idum - k * IQ1) - k * IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ1;
    *idum = IA1 * (*idum - k * IQ1) - k * IR1;
    if (*idum < 0) *idum += IM1;
    k = idum2 / IQ2;
    idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
    if (idum2 < 0) idum2 += IM2;
    j = iy / NDIV;
    iy = iv[j] - idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp = AM * iy) > RNMX) return RNMX;
    else return temp;
}

// generate standard normally distributed random number, N(0,1) 
float gasdev(long *idum)
{
    static int iset = 0;
    static float gset;
    float fac, rsq, v1, v2;
    
    if (*idum < 0) iset = 0;
    if (iset == 0)
    {
        do
        {
            v1 = 2.0 * ran1(idum) - 1.0;
            v2 = 2.0 * ran1(idum) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        }
        while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    }
    else
    {
        iset = 0;
        return gset;
    }
}

// function for creating a directory if it does not exist 
void create_directory(const char* path) {
    if (mkdir(path, 0777) == -1) {
        if (errno != EEXIST) {
            printf("Error creating directory: %s\n", path);
            exit(1);
        }
    }
    printf("Directory created or already exists: %s\n", path);
}

void assignJuncFE(long *seed, int inert_random, double inert_fracA) {
/*
inert_random: 1 = pick inert ENDS uniformly at random
              0 = pick inert CHAINS (paired ends); use inert_fracA for A/B split
*/
    const int N_end   = 2 * P.N_chain_tot;
    const int N_bound = P.arm * P.N_J;
    const int N_nr    = P.N_NR;
    const int N_fe    = P.N_FE;

    if (N_bound + N_fe + N_nr != N_end) {
        fprintf(stderr, "assignJuncFE: counts mismatch (%d + %d + %d != %d)\n",
                N_bound, N_fe, N_nr, N_end);
        exit(EXIT_FAILURE);
    }

    BonJ  = (int*) calloc(N_bound, sizeof(int));
    FEnd  = (int*) calloc(N_fe,   sizeof(int));
    NREnd = (int*) calloc(N_nr,   sizeof(int));
    if (!BonJ || !FEnd || !NREnd) { perror("calloc"); exit(EXIT_FAILURE); }

    // Build all ends (bead ids)
    int *ends = (int*) malloc(N_end * sizeof(int));
    if (!ends) { perror("malloc ends"); exit(EXIT_FAILURE); }
    for (int c = 0; c < P.N_chain_tot; ++c) {
        ends[2*c]     = c * N_BEAD;
        ends[2*c + 1] = c * N_BEAD + (N_BEAD - 1);
    }

    // Boolean by bead id: 1 = inert
    unsigned char *is_inert = (unsigned char*) calloc(P.N_chain_tot * N_BEAD, 1);
    if (!is_inert) { perror("calloc is_inert"); exit(EXIT_FAILURE); }

    if (inert_random) {
        // Randomly assign inert ENDS
        // shuffle ends 
        for (int i = N_end - 1; i > 0; --i) {
            int j = (int)(ran1(seed) * (i + 1));
            int t = ends[i]; ends[i] = ends[j]; ends[j] = t;
        }
        // assign first N_nr ends as inert 
        for (int k = 0; k < N_nr; ++k) { 
            NREnd[k] = ends[k]; 
            is_inert[ends[k]] = 1;
        }
    } else {
        // Paired inert CHAINS with A/B mix
        const int n_pairs = N_nr / 2; 

        int *chains = (int*) malloc(P.N_chain_tot * sizeof(int)); // indices of all chains 
        if (!chains) { perror("malloc chains"); exit(EXIT_FAILURE); }
        for (int c = 0; c < P.N_chain_tot; ++c) chains[c] = c;

        // shuffle chains 
        for (int i = P.N_chain_tot - 1; i > 0; --i) {
            int j = (int)(ran1(seed) * (i + 1));
            int t = chains[i]; chains[i] = chains[j]; chains[j] = t;
        }

        // Target A/B counts, then pick first chains meeting the split
        int nA_inert = (int) round(inert_fracA * n_pairs);
        if (nA_inert < 0) nA_inert = 0; if (nA_inert > n_pairs) nA_inert = n_pairs; // safeguard 
        int nB_inert = n_pairs - nA_inert;

        // First pass: satisfy A/B request
        int put = 0;
        for (int i = 0; i < P.N_chain_tot && put < N_nr; ++i) {
            int c = chains[i], t = poly[c].type[0];
            if ((t == 0 && nA_inert) || (t == 1 && nB_inert)) {
                int s = c*N_BEAD, e = s + (N_BEAD-1);
                NREnd[put++] = s; NREnd[put++] = e;
                is_inert[s] = is_inert[e] = 1;
                if (t == 0) --nA_inert; else --nB_inert;
            }
        }
        // Second pass: fill any shortfall with whatever remains
        for (int i = 0; i < P.N_chain_tot && put < N_nr; ++i) {
            int c = chains[i];
            int s = c*N_BEAD, e = s + (N_BEAD-1);
            if (!is_inert[s]) { NREnd[put++] = s; NREnd[put++] = e; is_inert[s] = is_inert[e] = 1; }
        }
        if (put != N_nr) { fprintf(stderr, "assignJuncFE: Not enough chains to place inert pairs (%d < %d)\n", put, N_nr); exit(EXIT_FAILURE); }
    }

    // Build reactive list (everything not inert), shuffle, then split to BonJ/FEnd
    int *reactiveEnd = (int*) malloc((N_end - N_nr) * sizeof(int));
    if (!reactiveEnd) { perror("malloc react"); exit(EXIT_FAILURE); }
    int rsz = 0; // reactive size
    for (int c = 0; c < P.N_chain_tot; ++c) {
        int s = c * N_BEAD, e = s + (N_BEAD - 1);
        // fill in reactive ends in array
        if (!is_inert[s]) reactiveEnd[rsz++] = s;
        if (!is_inert[e]) reactiveEnd[rsz++] = e;
    }
    if (rsz != N_bound + N_fe) {
        fprintf(stderr, "assignJuncFE: reactive count mismatch (%d != %d + %d)\n",
                rsz, N_bound, N_fe);
        exit(EXIT_FAILURE);
    }
    // shuffle reactive ends 
    for (int i = rsz - 1; i > 0; --i) {
        int j = (int)(ran1(seed) * (i + 1));
        int t = reactiveEnd[i]; reactiveEnd[i] = reactiveEnd[j]; reactiveEnd[j] = t;
    }
    // assign to junction-bound ends and reactive free ends 
    int idx = 0;
    for (int k = 0; k < N_bound; ++k, ++idx) BonJ[k] = reactiveEnd[idx];
    for (int k = 0; k < N_fe;    ++k, ++idx) FEnd[k] = reactiveEnd[idx];


    free(is_inert);
    free(reactiveEnd);
}


static void updateJType(int J) {
    // update junction type based on the number of bound beads
    int num_B = 0;
    for (int i = 0; i < P.arm; ++i) {
        int bead_idx = BonJ[J * P.arm + i];
        int c, m; globalToChain(bead_idx, &c, &m); // get chain and monomer indices
        if (poly[c].type[m] == 1) num_B++;
    }
    junc[J].type = num_B; // update junction type based on bound beads
}

// function to initialize the system 
static void init_sys(void) {
    poly = calloc(P.N_chain_tot, sizeof(chain)); 

    int c = 0; // counter for total number of chains  
    int N_CHAIN = (int) (P.N_chain_tot / N_TYPE); // number of chains per type 
    for (int type = 0; type < N_TYPE; ++type) {
        // initialize random number generator for each type
        long seed = -labs(initRan());
        for (int k = 0; k < N_CHAIN; ++k, ++c) { // total chain counter c increments here 
            poly[c].type[0] = type; // set type for first bead 
            // place first bead at a random position in the box 
            poly[c].r[0] = ran1(&seed) * P.box_len_b;
            poly[c].img[0] = 0; // initialize image counter for first bead
            // place the subsequent beads at a distance DELTA_L from the previous one 
            for (int i = 1; i < N_BEAD; ++i) {
                poly[c].type[i] = type; // set type for subsequent beads 
                poly[c].r[i] = poly[c].r[i-1] + P.deltaL_b * gasdev(&seed);
                poly[c].img[i] = 0; // initialize counter
                // wrap coordinates back into the box
                wrap(&poly[c].r[i]);
            }
        }
    }

    // initialize junctions
    junc = calloc(P.N_J, sizeof(junction)); 

    long seed = -labs(initRan());
    for (int j = 0; j < P.N_J; ++j) {
        junc[j].type = 0; // will be updated later based on bound beads 
        junc[j].img = 0; 
        junc[j].r = ran1(&seed) * P.box_len_b;
    }

    // initialize linked-cell list
    cellHead_SW = (int*) calloc(P.N_cell_sw, sizeof(int));
    cellNext_SW = (int*) calloc(P.N_FE, sizeof(int));
    cellHead_NB = (int*) calloc(P.N_cell_nb, sizeof(int));
    cellNext_NB = (int*) calloc(P.N_chain_tot * N_BEAD, sizeof(int)); 
    for (int i = 0; i < P.N_cell_sw; ++i) cellHead_SW[i] = -1; 
    for (int i = 0; i < P.N_FE; ++i) cellNext_SW[i] = -1; 
    for (int i = 0; i < P.N_cell_nb; ++i) cellHead_NB[i] = -1;
    for (int i = 0; i < P.N_chain_tot * N_BEAD; ++i) cellNext_NB[i] = -1;

    // initialize bookkeeping arrays
    assignJuncFE(&seed, P.inert_random, P.inert_fracA); 
    for (int J = 0; J < P.N_J; ++J) updateJType(J); // initialize junction types 
    arm_change_count = (unsigned int*)calloc(P.arm * P.N_J, sizeof(unsigned int));
}

void init_arr(int *JuncDist, int *CONNECT, int *CIJ, int nUpperTri) {
    for (int i = 0; i < nUpperTri; ++i) {
        CIJ[i] = 0; 
    }
    for (int i = 0; i < P.N_J * (P.arm + 1); ++i) {
        JuncDist[i] = 0; 
    }
    for (int i = 0; i < 2 * P.N_chain_tot; ++i) {
        CONNECT[i] = -1; // default -1 means free end
    }
} 

static void init_metrics(void) {
    swaps_J        = calloc(P.N_J, sizeof(long long));
    in_contact     = calloc(P.N_J, sizeof(char));
    contact_start  = calloc(P.N_J, sizeof(double));
    contact_time   = calloc(P.N_J, sizeof(double));
    windows_total  = calloc(P.N_J, sizeof(long long));
    windows_missed = calloc(P.N_J, sizeof(long long));
    last_swap_time = calloc(P.N_J, sizeof(double));
    wait_sum       = calloc(P.N_J, sizeof(double));
    wait_sumsq     = calloc(P.N_J, sizeof(double));
    wait_n         = calloc(P.N_J, sizeof(unsigned long));

    for (int J = 0; J < P.N_J; ++J) {
        in_contact[J] = 0;
        contact_start[J] = 0.0;
        contact_time[J] = 0.0;
        windows_total[J] = 0;
        windows_missed[J] = 0;
        last_swap_time[J] = -1.0;
    }
}


void computePHI(double *PHIA, double *PHIB) {
    for (int c = 0; c < P.N_chain_tot; ++c) {
        for (int i = 0; i < N_BEAD; ++i) {
            gridPointInterpolator s = mkInterp(poly[c].r[i]); 
            if (poly[c].type[i] == 0) { 
                depositPhi(PHIA, &s, P.DPHI, 1); 
            } else { 
                depositPhi(PHIB, &s, P.DPHI, 1); 
            } 
        }
    }
}

static inline int has_contact_count(int J) {
    int cx0 = cellIndex_1D(junc[J].r, P.N_cell_sw, P.L_cell_sw_b);
    int nCand = 0;
    for (int dcx = -1; dcx <= 1; ++dcx) {
        int cell = (cx0 + dcx + P.N_cell_sw) % P.N_cell_sw;
        for (int fe = cellHead_SW[cell]; fe != -1; fe = cellNext_SW[fe]) {
            int b = FEnd[fe];
            int c, m; globalToChain(b, &c, &m);
            double dx; double r2 = rsq(junc[J].r, poly[c].r[m], &dx);
            if (r2 < P.R_swap_b * P.R_swap_b) {
                if (++nCand >= 1) return nCand; // early-exit; just need a boolean
            }
        }
    }
    return nCand;
}


static int cand_swap(int J, int *cand) {
    // center cell of the junction J 
    int cx0 = cellIndex_1D(junc[J].r, P.N_cell_sw, P.L_cell_sw_b); 
    int nCand = 0; 

    // search the 3 surrounding cells 
    for (int dcx = -1; dcx <= 1; ++dcx) {
        int cellIdx = (cx0 + dcx + P.N_cell_sw) % P.N_cell_sw; // wrap around
        
        // traverse the linked-cell list in cell c 
        for (int fe = cellHead_SW[cellIdx]; fe != -1; fe = cellNext_SW[fe]) {
            int b = FEnd[fe]; // global bead index 
            int c, m; globalToChain(b, &c, &m); // get chain and monomer indices
            double dx; 
            double r2 = rsq(junc[J].r, poly[c].r[m], &dx); // distance from junction to polymer free end bead 
            if (r2 < P.R_swap_b * P.R_swap_b) { // within reaction radius
                cand[nCand++] = fe; // store free-end slot index 
                if (nCand == MAX_CAND) return nCand; // reached maximum candidates
            }
        }
    }
    return nCand; 
}

// Fisher–Yates shuffle for int array
static inline void shuffle_int(int *a, int n, long *seed){
    for (int i=n-1; i>0; --i){
        int j = (int)(ran1(seed) * (i+1));
        int tmp=a[i]; a[i]=a[j]; a[j]=tmp;
    }
}

static int attempt_swap_with_target(int J, int randArm, int fe_slot, long *seed,
                                    double Ea_tilde, double ks_junc_tilde, double T_tilde)
{
    int jBE_idx  = J * P.arm + randArm;

    int b_fe     = FEnd[fe_slot];
    int b_onJ    = BonJ[jBE_idx];

    int c_onJ,m_onJ,c_fe,m_fe;
    globalToChain(b_onJ,&c_onJ,&m_onJ);
    globalToChain(b_fe, &c_fe, &m_fe);

    double dx, r2_old = rsq(junc[J].r, poly[c_onJ].r[m_onJ], &dx);
    double     r2_new = rsq(junc[J].r, poly[c_fe].r[m_fe],   &dx);

    double dU_over_kBT = 0.5 * ks_junc_tilde * (r2_new - r2_old) / T_tilde;

    if (dU_over_kBT > MAX_DU) dU_over_kBT = MAX_DU; // cap to avoid overflow 
    const double p_act = exp(-Ea_tilde / T_tilde);
    double p_acc = (dU_over_kBT <= 0.0) ? p_act : p_act * exp(-dU_over_kBT);
    if (p_acc > 1.0) p_acc = 1.0; if (p_acc < 0.0) p_acc = 0.0;
    if (ran1(seed) >= p_acc) return 0; // reject

    // commit
    BonJ[jBE_idx] = b_fe;
    FEnd[fe_slot] = b_onJ;
    // mark that this arm changed partner
    arm_change_count[jBE_idx] += 1u;
    return 1;
}

static void rebuild_list_SW() {
    for (int i = 0; i < P.N_cell_sw; ++i) cellHead_SW[i] = -1; 
    for (int i = 0; i < P.N_FE; ++i) cellNext_SW[i] = -1; 

    for (int i = 0; i < P.N_FE; ++i) {
        int b_fe = FEnd[i]; 
        int c, m; 
        globalToChain(b_fe, &c, &m); // get chain index and monomer index

        int cellIdx = cellIndex_1D(poly[c].r[m], P.N_cell_sw, P.L_cell_sw_b);

        cellNext_SW[i] = cellHead_SW[cellIdx]; 
        cellHead_SW[cellIdx] = i; 
    }
}

static void rebuild_list_NB() {
    for (int i = 0; i < P.N_cell_nb; ++i) cellHead_NB[i] = -1; 
    for (int i = 0; i < P.N_chain_tot * N_BEAD; ++i) cellNext_NB[i] = -1; 

    for (int g = 0; g < P.N_chain_tot * N_BEAD; ++g) {
        int c, m; globalToChain(g, &c, &m); // get chain index and monomer index
        int cellIdx = cellIndex_1D(poly[c].r[m], P.N_cell_nb, P.L_cell_nb_b);
        cellNext_NB[g] = cellHead_NB[cellIdx];
        cellHead_NB[cellIdx] = g;
    }
}

void compute_MSD_undrifted(double *MSD_seg, double *MSD_com, double *r0, double *Rcom0_chain, double Rcom0_sysRef, double *RcomNow_chain) {
    *MSD_seg = *MSD_com = 0.0;

    // current chain COMs and system COM (unwrapped)
    double RcomNow_sys; 
    snapshot_COMs(&RcomNow_sys, RcomNow_chain); 
    double dR = RcomNow_sys - Rcom0_sysRef;

    // segmental
    double msd_seg_acc = 0.0;
    int bead = 0;
    for (int c = 0; c < P.N_chain_tot; ++c) {
        for (int i = 0; i < N_BEAD; ++i, ++bead) {
            double pos_unw = poly[c].r[i] + poly[c].img[i] * P.box_len_b;
            double dx = (pos_unw - r0[bead]) - dR;
            msd_seg_acc += dx * dx;
        }
    }
    *MSD_seg = msd_seg_acc / (double)(P.N_chain_tot * N_BEAD);

    // chain COM (average over chains)
    double msd_com_acc = 0.0;
    for (int c = 0; c < P.N_chain_tot; ++c) {
        double dxc = (RcomNow_chain[c] - Rcom0_chain[c]) - dR;
        msd_com_acc += dxc * dxc;
    }
    *MSD_com = msd_com_acc / (double)P.N_chain_tot;
}

static void start_msd_origin(int nBD) {
    MSDOrigin *O = &msd_starts[next_msd_start_slot];

    const int nBeads = P.N_chain_tot * N_BEAD;
    const int nChains = P.N_chain_tot; 
    O->n_beads = nBeads; 
    O->rBead0_O    = &msd_r0_pool[next_msd_start_slot * nBeads];
    O->comChain0_O = &msd_com0_pool[next_msd_start_slot * nChains];

    int b = 0; 
    double com_sys_acc = 0.0;
    for (int c = 0; c < P.N_chain_tot; ++c) {
        double sum_chain = 0.0;
        for (int i = 0; i < N_BEAD; ++i, ++b) {
            double pos_unw = poly[c].r[i] + poly[c].img[i] * P.box_len_b;
            O->rBead0_O[b] = pos_unw;
            sum_chain += pos_unw;
        }
        O->comChain0_O[c] = sum_chain / (double)N_BEAD;
        com_sys_acc += O->comChain0_O[c];
    }
    O->comSys0_O = com_sys_acc / (double)P.N_chain_tot;

    O->active = 1;
    O->start_nBD = nBD;
    O->next_check_index = 0;
    static long long msd_id_counter = 0;
    O->id = msd_id_counter++;

    last_msd_start_nBD = nBD;
    next_msd_start_slot = (next_msd_start_slot + 1) % MAX_MSD_STARTS;
}


// Average squared bond length over all 1–2 neighbors
static double mean_bond_r2(void) {
    double sum_r2 = 0.0;
    long long n_bonds = 0;

    for (int c = 0; c < P.N_chain_tot; ++c) {
        for (int i = 0; i < N_BEAD - 1; ++i) {
            double dx; 
            double r2 = rsq(poly[c].r[i+1], poly[c].r[i], &dx); // PBC-aware
            sum_r2 += r2;
            ++n_bonds;
        }
    }
    return (n_bonds > 0) ? (sum_r2 / (double)n_bonds) : 0.0;
}

static double mean_jb_r2(void) {
    double sum_r2 = 0.0; 
    long long n_bonds = 0;
    for (int J = 0; J < P.N_J; ++J)
        for (int a = 0; a < P.arm; ++a) {
            int g = BonJ[P.arm*J + a];
            int c, m; globalToChain(g, &c, &m);
            double dx; double r2 = rsq(junc[J].r, poly[c].r[m], &dx);
            sum_r2 += r2;
            ++n_bonds; 
        }
    return (n_bonds > 0) ? (sum_r2 / (double)n_bonds) : 0.0;
}

// ------------------------------------------------------------------------------------------------- 

// ---------- main function ------------------------------------------------------------------------
int main(int argc, const char *argv[]) {    
    // command line input ------------------------------------------------------ 
    double CHI_0_N_cli      = (argc > chi0N_idx) ? atof(argv[chi0N_idx]) : 0.0;
    int    STRIDE_cli       = (argc > mc_stride_idx) ? atoi(argv[mc_stride_idx]) : 10;
    double Ea_tilde_cli     = (argc > Ea_idx) ? atof(argv[Ea_idx]) : 0.0; 
    int    arm_cli          = (argc > arm_idx) ? atoi(argv[arm_idx]) : 3;
    double bound_ratio_cli  = (argc > bound_ratio_idx) ? atof(argv[bound_ratio_idx]) : 0.95;
    double inert_ratio_cli  = (argc > inert_ratio_idx) ? atof(argv[inert_ratio_idx]) : 0.0;
    REPLICA                 = (argc > replica_idx) ? atoi(argv[replica_idx]) : 0;

    if (argc < NUM_CLI) {
        fprintf(stderr, "Using defaults CHI_0_N=%.1f, STRIDE=%d, Ea_tilde=%.1f, arm=%d, bound_ratio=%.2f, inert_ratio=%.2f, replica=%d\n", 
                CHI_0_N_cli, STRIDE_cli, Ea_tilde_cli, arm_cli, bound_ratio_cli, inert_ratio_cli, REPLICA);
    }
    // parameter setup ---------------------------------------------------------
    fill_params(&P, CHI_0_N_cli, STRIDE_cli, Ea_tilde_cli, arm_cli, bound_ratio_cli, inert_ratio_cli);
    const double chi0N_target = P.chi0N; 
    const double Ea_tilde_target = P.Ea_tilde; 
    // guard 
    if (P.N_cell_nb <= 0 || P.L_cell_nb_b <= 0.0) { 
        fprintf(stderr, "Bad NB grid: N_cell_nb=%d L_cell_nb_b=%.6g\n", P.N_cell_nb, P.L_cell_nb_b);
        exit(EXIT_FAILURE);
    }
    if (P.N_cell_sw <= 0 || P.L_cell_sw_b <= 0.0) { 
        fprintf(stderr, "Bad SW grid: N_cell_sw=%d L_cell_sw_b=%.6g\n", P.N_cell_sw, P.L_cell_sw_b);
        exit(EXIT_FAILURE);
    }

    // step size recommender 
    double dt_nb=0.0, dt_bond=0.0, dt_noise=0.0; 
    DT = recommend_dt(P.ks_junc_tilde, P.gamma_tilde, P.kBT_tilde, ALPHA_NB, BETA_NOISE, ALPHA_BOND, &dt_nb, &dt_bond, &dt_noise);
    const double sigma = sqrt(2.0 * P.kBT_tilde * DT / P.gamma_tilde); 
    // DT should be in unit of [tau_ref]
    int DUMP_DEGSEP_AVG_STRIDE = (DUMP_AVG_SAMPLES + DUMP_LOG_STRIDE - 1) / DUMP_LOG_STRIDE;

    // bond acf setup ----------------------------
    NUM_BONDS = P.arm * P.N_J; 

    // estimate ACF spacing from exponential gate 
    const double ACF_PTS_PER_DECADE = 128.0; // number of ACF points per decade of time
    double p_act = exp(-P.Ea_tilde / P.T_tilde); 
    double tMC = P.mc_stride * DT; 
    double k_est = (p_act > 0.0 && tMC > 0.0) ? (p_act / tMC) : 0.0;     
    double dt_acf_target = (k_est > 0.0) ? (1.0 / k_est) / ACF_PTS_PER_DECADE : 10.0 * tMC; // target dt_acf ~ \tau / ACF_PTS_PER_DECADE, where \tau = 1/k_est 
    if (!isfinite(dt_acf_target) || dt_acf_target <= 0.0) dt_acf_target = 10.0 * tMC;
    
    int stride_from_gate = (int)ceil(dt_acf_target / DT);
    if (stride_from_gate < 1) stride_from_gate = 1;
    
    CHECK_EVERY_BD = stride_from_gate;
    int mult = (CHECK_EVERY_BD + P.mc_stride - 1) / P.mc_stride;   // ceil-div
    CHECK_EVERY_BD = mult * P.mc_stride;                           // >= mc_stride
    // cover at most D decades of tau
    const double ACF_DECADES = 4.0;                 // tweakable
    double tau_est = (k_est > 0.0) ? 1.0 / k_est : 1e99;
    double tmax_tau = tau_est * pow(10.0, ACF_DECADES);

    int num_checks_tau = (int)ceil(tmax_tau / (CHECK_EVERY_BD * DT));
    if (num_checks_tau < 1) num_checks_tau = 1;

    NUM_CHECKS = (int)fmin(NUM_CHECKS, num_checks_tau);

    if (MIN_ORG_GAP < P.mc_stride) MIN_ORG_GAP = P.mc_stride; 
    // guard for ACF sizing 
    NUM_CHECKS = MAX_LAG_BD / CHECK_EVERY_BD;
    if (NUM_CHECKS < 1) NUM_CHECKS = 1;

    const int ACF_NUM_CHECKS_CAP = 500000; 
    if (NUM_CHECKS > ACF_NUM_CHECKS_CAP) {
        fprintf(stderr, "[ACF] Capping NUM_CHECKS from %d to %d.\n",
                NUM_CHECKS, ACF_NUM_CHECKS_CAP);
        NUM_CHECKS = ACF_NUM_CHECKS_CAP;
    }
    acf_sum   = (double*)calloc(NUM_CHECKS, sizeof(double));
    acf_sumsq = (double*)calloc(NUM_CHECKS, sizeof(double));
    acf_count = (unsigned long*)calloc(NUM_CHECKS, sizeof(unsigned long));
    if (!acf_sum || !acf_sumsq || !acf_count) { perror("calloc ACF"); exit(EXIT_FAILURE); }

    // MSD setup ---------------------------------
    
    // Max lag we can possibly see in this run (respect equilibration & user max)
    long long max_observable_lag = (nBDmax > N_EQBR) ? (nBDmax - N_EQBR) : nBDmax;
    if (max_observable_lag > MSD_MAX_LAG_BD) max_observable_lag = MSD_MAX_LAG_BD;

    MSD_NUM_CHECKS = (int)((max_observable_lag + MSD_CHECK_EVERY_BD - 1) / MSD_CHECK_EVERY_BD);
    if (MSD_NUM_CHECKS < 1) MSD_NUM_CHECKS = 1;

    // Hard cap so allocations stay sane 
    const int MSD_NUM_CHECKS_CAP = 2000000;
    if (MSD_NUM_CHECKS > MSD_NUM_CHECKS_CAP) {
        fprintf(stderr, "[MSD] Capping NUM_CHECKS from %d to %d to limit memory.\n",
                MSD_NUM_CHECKS, MSD_NUM_CHECKS_CAP);
        MSD_NUM_CHECKS = MSD_NUM_CHECKS_CAP;
    }
    // origin spacing
    long long span = 1LL * MSD_NUM_CHECKS * MSD_CHECK_EVERY_BD;

    long long gap_ll = (span + MAX_MSD_STARTS - 1) / MAX_MSD_STARTS; // ceil
    if (gap_ll < MSD_CHECK_EVERY_BD) gap_ll = MSD_CHECK_EVERY_BD;
    if (gap_ll > INT_MAX) gap_ll = INT_MAX; // cap before cast
    MSD_ORIGIN_GAP = (int)gap_ll;

    msd_seg_sum = (double*)calloc(MSD_NUM_CHECKS, sizeof(double));
    msd_seg_sumsq = (double*)calloc(MSD_NUM_CHECKS, sizeof(double));
    msd_seg_cnt = (unsigned long*)calloc(MSD_NUM_CHECKS, sizeof(unsigned long));
    msd_com_sum = (double*)calloc(MSD_NUM_CHECKS, sizeof(double));
    msd_com_sumsq = (double*)calloc(MSD_NUM_CHECKS, sizeof(double));
    msd_com_cnt = (unsigned long*)calloc(MSD_NUM_CHECKS, sizeof(unsigned long));
    msd_r0_pool = (double*)calloc(MAX_MSD_STARTS * P.N_chain_tot * N_BEAD, sizeof(double)); 
    msd_com0_pool = (double*)calloc(MAX_MSD_STARTS * P.N_chain_tot, sizeof(double));

    if (!msd_seg_sum || !msd_seg_sumsq || !msd_seg_cnt ||
        !msd_com_sum || !msd_com_sumsq || !msd_com_cnt) {
        perror("calloc MSD accumulators"); exit(EXIT_FAILURE);
    }
    if (!msd_r0_pool || !msd_com0_pool) { perror("calloc MSD pools"); exit(EXIT_FAILURE); }

    // allocate memories -------------------------------------------------------
    // density fields 
    double *PHIA = calloc(NGRID, sizeof(double));
    double *PHIB = calloc(NGRID, sizeof(double));
    double *PHIA_sum = calloc(NGRID, sizeof(double)); 
    double *PHIB_sum = calloc(NGRID, sizeof(double));
    // MSD reference positions 
    double *r0 = calloc(P.N_chain_tot * N_BEAD, sizeof(double)); 
    double *Rcom0_chain = calloc(P.N_chain_tot, sizeof(double)); 
    double *RcomNow_chain = calloc(P.N_chain_tot, sizeof(double));
    // bookkeeping for distributions 
    int *JuncDist = calloc(P.N_J * (P.arm + 1), sizeof(int)); 
    int *CONNECT = calloc(2 * P.N_chain_tot, sizeof(int)); 
    int nUpperTri = ((P.arm + 1) * (P.arm + 2)) / 2; // number of upper triangular elements in the CIJ matrix 
    int *CIJ = calloc(nUpperTri, sizeof(int));     
    for (int i = 0; i < MAX_STARTS; ++i) {
        starts[i].active = 0;
        starts[i].start_nBD = 0;
        starts[i].next_check_index = 0;
        starts[i].bond_at_start = NULL;
        starts[i].change_at_start = NULL;
    }
    next_start_slot = 0;
    last_start_nBD = -1000000000;
    accepts_at_last_start = 0;    
    // file IO ----------------------------------------------------------------- 
    char outputPath[STR_LEN];
    snprintf(outputPath, sizeof(outputPath), "./data_BD_vtrmr_v%s/", VERSION);     

    create_directory(outputPath);

    static const FileSpec specs[NUM_STR] = {
        { "param_",       ".txt"  }, // F_PARAM
        { "Phi_",         ".txt"  }, // F_PHI
        { "Coord_",       ".xyz"  }, // F_COORD
        { "Phi_avg_",     ".txt"  }, // F_PHI_AVG
        { "log_",         ".csv"  }, // F_LOG
        { "distr_",       ".txt"  }, // F_DISTRIB
        { "waitTime_",    ".csv"  }, // F_WAIT_TIME
        { "dynSite_",     ".txt"  }, // F_DYNSITE
        { "degSep_avg_",  ".txt"  }, // F_DEGSEP_AVG
        { "bondCorr_",    ".csv"  }, // F_BOND_ACF
        { "msd_",         ".csv"  }  // F_MSD
    };

    FILE *files[NUM_STR];
    char *strs[NUM_STR];
    for (int i = 0; i < NUM_STR; i++) {
        const char *pfx = specs[i].prefix;
        const char *ext = specs[i].ext;
        int len = snprintf(NULL, 0, "%s%sX0N%.1f_sqN%d_Eatilde%.2f_arm%d_bound%.2f_inert%.2f_rep%d_v%s%s", outputPath, pfx, P.chi0N, P.sqrtNbar, P.Ea_tilde, P.arm, P.bound_ratio, P.inert_ratio, REPLICA, VERSION, ext);
        strs[i] = malloc((len + 1)); 
        snprintf(strs[i], len + 1, "%s%sX0N%.1f_sqN%d_Eatilde%.2f_arm%d_bound%.2f_inert%.2f_rep%d_v%s%s", outputPath, pfx, P.chi0N, P.sqrtNbar, P.Ea_tilde, P.arm, P.bound_ratio, P.inert_ratio, REPLICA, VERSION, ext);
    }

    // Headers ----------------------------------------------------------------- 
    files[F_LOG] = fopen(strs[F_LOG], "w");
    if (files[F_LOG]) {
        fprintf(files[F_LOG], "nBD,nBD_prod,tBD,tBD_prod,in_equil,MeanSqRee,MeanSqRee/(N-1),Ub,Unb,Utot,F_tot,R_cm,R_cm_unw,MSD_seg,MSD_com,mean_bond2,mean_Jbond2,nn_avg,weighted_nn_avg,list_eff,degSep,neighbor_avail,acc_ratio_swap,acc_ratio_swap_raw,k_swap_bulk,k_swap_perJ,k_intrinsic,f_miss,allAjunc_frac,allBjunc_frac,alljunc_frac\n");
        fclose(files[F_LOG]);
    }
    files[F_WAIT_TIME] = fopen(strs[F_WAIT_TIME], "w");
    if (files[F_WAIT_TIME]) {
        fprintf(files[F_WAIT_TIME], "nBD_prod,tBD_prod,events,nJ,mean_wait,sem_wait,k_hat\n"); 
        fclose(files[F_WAIT_TIME]);
    }  
    // ------------------------------------------------------------------------- 
    files[F_PARAM] = fopen(strs[F_PARAM], "w");
    // prefix "param_"
    if (files[F_PARAM]) {
        fprintf(files[F_PARAM], "SYSTEM: %s\n", SYSTEM);
        fprintf(files[F_PARAM], "DIM: %d\n", DIM);
        fprintf(files[F_PARAM], "VERSION: %s\n", VERSION);
        fprintf(files[F_PARAM], "N_TYPE: %d\n", N_TYPE);
        fprintf(files[F_PARAM], "SQRT_N_bar: %d\n", P.sqrtNbar);
        fprintf(files[F_PARAM], "N_INT: %d\n", N_INT);
        fprintf(files[F_PARAM], "N_BEAD: %d\n", N_BEAD);
        fprintf(files[F_PARAM], "N_EQBR: %d\n", N_EQBR);
        fprintf(files[F_PARAM], "B2 [ref]: %lf\n", P.B2_tilde);
        fprintf(files[F_PARAM], "Re_b [b]: %lf\n", P.Re_b);
        fprintf(files[F_PARAM], "BOX_LEN [Re]: %lf\n", P.box_len_b / P.Re_b);
        fprintf(files[F_PARAM], "BOX_LEN [b]: %lf\n", P.box_len_b);
        fprintf(files[F_PARAM], "KAPPA_N: %lf\n", P.kappaN);
        fprintf(files[F_PARAM], "CHI_0_N: %lf\n", P.chi0N);                
        fprintf(files[F_PARAM], "N_CHAIN_TOT: %d\n", P.N_chain_tot);
        fprintf(files[F_PARAM], "BOUND_RATIO: %lf\n", P.bound_ratio);
        fprintf(files[F_PARAM], "INERT_RATIO: %lf\n", P.inert_ratio);
        fprintf(files[F_PARAM], "INERT_RANDOM: %d\n", P.inert_random);
        fprintf(files[F_PARAM], "INERT_FRACA: %lf\n", P.inert_fracA); 
        fprintf(files[F_PARAM], "ARM: %d\n", P.arm);
        fprintf(files[F_PARAM], "N_J: %d\n", P.N_J);
        fprintf(files[F_PARAM], "N_NR: %d\n", P.N_NR);
        fprintf(files[F_PARAM], "N_FE: %d\n", P.N_FE);
        fprintf(files[F_PARAM], "rho_0 [#/Re]: %lf\n", P.rho0_per_Re);
        fprintf(files[F_PARAM], "DELTA_L [Re]: %lf\n", P.deltaL_b / P.Re_b);
        fprintf(files[F_PARAM], "DELTA_L [b]: %lf\n", P.deltaL_b);
        fprintf(files[F_PARAM], "R_SWAP [Re]: %lf\n", P.R_swap_b / P.Re_b);
        fprintf(files[F_PARAM], "R_SWAP [b]: %lf\n", P.R_swap_b);
        fprintf(files[F_PARAM], "L_CELL_SW [Re]: %lf\n", P.L_cell_sw_b / P.Re_b);
        fprintf(files[F_PARAM], "L_CELL_SW [b]: %lf\n", P.L_cell_sw_b);
        fprintf(files[F_PARAM], "N_CELL_SW: %d\n", P.N_cell_sw);
        fprintf(files[F_PARAM], "L_CELL_NB [Re]: %lf\n", P.L_cell_nb_b / P.Re_b);
        fprintf(files[F_PARAM], "L_CELL_NB [b]: %lf\n", P.L_cell_nb_b);
        fprintf(files[F_PARAM], "N_CELL_NB: %d\n", P.N_cell_nb);
        fprintf(files[F_PARAM], "MG: %d\n", MG);
        fprintf(files[F_PARAM], "NGRID: %d\n", NGRID);
        fprintf(files[F_PARAM], "Cw_CLOUD: %lf\n", P.Cw_cloud);
        fprintf(files[F_PARAM], "DPHI: %lf\n", P.DPHI);
        fprintf(files[F_PARAM], "DT: %.3e\n", DT);
        fprintf(files[F_PARAM], "dt_nb: %.3e\n", dt_nb);
        fprintf(files[F_PARAM], "dt_bond: %.3e\n", dt_bond);
        fprintf(files[F_PARAM], "dt_noise: %.3e\n", dt_noise);
        fprintf(files[F_PARAM], "ALPHA_NB: %lf\n", ALPHA_NB);
        fprintf(files[F_PARAM], "BETA_NOISE: %lf\n", BETA_NOISE);
        fprintf(files[F_PARAM], "ALPHA_BOND: %lf\n", ALPHA_BOND);
        fprintf(files[F_PARAM], "DUMP_PHI_STRIDE: %d\n", DUMP_PHI_STRIDE);
        fprintf(files[F_PARAM], "DUMP_COORD_STRIDE: %d\n", DUMP_COORD_STRIDE);
        fprintf(files[F_PARAM], "DUMP_LOG_STRIDE: %d\n", DUMP_LOG_STRIDE);
        fprintf(files[F_PARAM], "DUMP_AVG_SAMPLES: %d\n", DUMP_AVG_SAMPLES);
        fprintf(files[F_PARAM], "DUMP_DISTRIB_STRIDE: %d\n", DUMP_DISTRIB_STRIDE);
        fprintf(files[F_PARAM], "DUMP_WAIT_STRIDE: %d\n", DUMP_WAIT_STRIDE);
        fprintf(files[F_PARAM], "DUMP_SITE_COORD_STRIDE: %d\n", DUMP_SITE_COORD_STRIDE);
        fprintf(files[F_PARAM], "MC_STRIDE: %d\n", P.mc_stride);
        fprintf(files[F_PARAM], "CHECK_EVERY_BD: %d\n", CHECK_EVERY_BD);
        fprintf(files[F_PARAM], "MIN_ORG_GAP: %d\n", MIN_ORG_GAP);
        fprintf(files[F_PARAM], "MAX_LAG_BD: %d\n", MAX_LAG_BD);
        fprintf(files[F_PARAM], "NUM_BONDS: %d\n", NUM_BONDS);
        fprintf(files[F_PARAM], "NUM_CHECKS: %d\n", NUM_CHECKS);
        fprintf(files[F_PARAM], "b_ref [m]: %.2e\n", P.b_ref);
        fprintf(files[F_PARAM], "m_ref [kg]: %.2e\n", P.m_ref);
        fprintf(files[F_PARAM], "Tref [K]: %.2f\n", P.Tref_K);
        fprintf(files[F_PARAM], "kBTref [J]: %.2e\n", P.kBTref);
        fprintf(files[F_PARAM], "tau_ref [s]: %.2e\n", P.tau_ref);
        fprintf(files[F_PARAM], "gamma_ref [kg/s]: %.2e\n", P.gamma_ref);
        fprintf(files[F_PARAM], "Tphys [K]: %.2f\n", P.Tphys_K);
        fprintf(files[F_PARAM], "Ea_kJmol: %.2f\n", P.Ea_kJmol);
        fprintf(files[F_PARAM], "T_tilde: %.2f\n", P.T_tilde);
        fprintf(files[F_PARAM], "kBT_tilde: %.2f\n", P.kBT_tilde);
        fprintf(files[F_PARAM], "Ea_tilde: %.2f\n", P.Ea_tilde);
        fprintf(files[F_PARAM], "gamma_tilde: %.2f\n", P.gamma_tilde);
        fprintf(files[F_PARAM], "ks_mon_tilde: %.2f\n", P.ks_mon_tilde);
        fprintf(files[F_PARAM], "ks_junc_tilde: %.2f\n", P.ks_junc_tilde);
        fprintf(files[F_PARAM], "REPLICA: %d\n", REPLICA);
        fprintf(files[F_PARAM], "nBDmax: %lld\n", nBDmax);         
        if (HAS_NOTE) {
            fprintf(files[F_PARAM], "NOTE: %s\n", PARAM_NOTE);
        }        
        fclose(files[F_PARAM]);
    }

    // initialize random number generator
    long rng = -labs(initRan());
    // initialize system
    init_sys(); 
    init_metrics();
    init_arr(JuncDist, CONNECT, CIJ, nUpperTri); 

    int win_prop_swap = 0, win_tried_swap = 0, win_acc_swap = 0;
    int win_tried_swap_raw = 0, win_acc_swap_raw = 0; 
    
    int cum_tried_swap = 0, cum_acc_swap = 0; 
    
    int sample_count = 0, degSep_count = 0;    
    double degSep_sum = 0.0;     
    int start_id = 0;

    for (long long nBD = 0; nBD < nBDmax; ++nBD) {
        // diagnostics ------------------------------------- 
        // if ((nBD % 10) == 0) { 
        //     printf("nBD = %lld\n", nBD); 
        //     fflush(stdout); 
        // }
        // ------------------------------------------------- 
        // equilibration -----------------------------------
        const int in_equil = (nBD < N_EQBR); 
        P.chi0N = in_equil ? 0.0 : chi0N_target;
        P.Ea_tilde = in_equil ? 0.0 : Ea_tilde_target;
        // absolute and production clock 
        const int nBD_prod = in_equil ? 0 : (nBD - N_EQBR);
        const double tBD_abs = (nBD + 1) * DT;
        const double tNow_prod = in_equil? 0.0: (nBD - N_EQBR) * DT; // start-of-step prod time 
        const double tBD_prod = in_equil? 0.0: (nBD_prod + 1) * DT; // end-of-step prod time 
        // -------------------------------------------------
        if (nBD == N_EQBR) { // after equilibration
            // MSD reference positions 
            int idx = 0; 
            double Rcm_sum_unw = 0.0;
            for (int c = 0; c < P.N_chain_tot; ++c) {
                double sum = 0.0; 
                for (int i = 0; i < N_BEAD; ++i, ++idx) {
                    double pos_unw = poly[c].r[i] + poly[c].img[i] * P.box_len_b; // unwrapped position
                    r0[idx] = pos_unw; // reference position of bead i of chain c (global index idx) 
                    sum += pos_unw; 
                }
                Rcom0_chain[c] = sum / N_BEAD; // reference COM of chain c 
                Rcm_sum_unw += Rcom0_chain[c];
            }
            Rcom0_sysRef = Rcm_sum_unw / (double)P.N_chain_tot; // reference system COM
            // reset metrics
            sample_count = 0;
            start_id = 0; 
            memset(PHIA_sum, 0, NGRID * sizeof(double));
            memset(PHIB_sum, 0, NGRID * sizeof(double));
            degSep_sum = 0.0;
            degSep_count = 0;
            cum_acc_swap = 0;
            for (int J = 0; J < P.N_J; ++J) {
                contact_time[J] = 0.0;
                windows_total[J] = 0;
                windows_missed[J] = 0;
                last_swap_time[J] = -1.0;
                if (in_contact[J]) { // if in contact at the end of equilibration
                    in_contact[J] = 0;
                    contact_start[J] = 0.0;
                }
            }
            memset(JuncDist, 0, P.N_J * (P.arm + 1) * sizeof(int));
            memset(CONNECT, -1, 2 * P.N_chain_tot * sizeof(int));
            memset(wait_sum,   0, P.N_J * sizeof(double));
            memset(wait_sumsq, 0, P.N_J * sizeof(double));
            memset(wait_n,     0, P.N_J * sizeof(unsigned long));
            // bond ACF 
            memset(acf_sum,   0, NUM_CHECKS * sizeof(double));
            memset(acf_sumsq, 0, NUM_CHECKS * sizeof(double));
            memset(acf_count, 0, NUM_CHECKS * sizeof(unsigned long));
            acf_max_i_seen = 0;
            // MSD 
            for (int i = 0; i < MAX_MSD_STARTS; ++i) {
                msd_starts[i].active = 0;
                msd_starts[i].start_nBD = 0;
                msd_starts[i].next_check_index = 0;
            }
            next_msd_start_slot = 0;
            last_msd_start_nBD  = -1000000000;
            msd_max_i_seen = 0;
            memset(msd_seg_sum,   0, MSD_NUM_CHECKS*sizeof(double));
            memset(msd_seg_sumsq, 0, MSD_NUM_CHECKS*sizeof(double));
            memset(msd_seg_cnt,   0, MSD_NUM_CHECKS*sizeof(unsigned long));
            memset(msd_com_sum,   0, MSD_NUM_CHECKS*sizeof(double));
            memset(msd_com_sumsq, 0, MSD_NUM_CHECKS*sizeof(double));
            memset(msd_com_cnt,   0, MSD_NUM_CHECKS*sizeof(unsigned long));            
        }

        rebuild_list_NB();
        // compute conservative forces using old positions 
        // compute_conservative_forces(P.ks_mon_tilde, P.ks_junc_tilde, P.kBT_tilde);
        compute_conservative_forces_pairwise(P.ks_mon_tilde, P.ks_junc_tilde, P.kBT_tilde);
        
        // evolve beads positions ---------------------------------------------- 
        // polymer beads 
        for (int c = 0; c < P.N_chain_tot; ++c) {
            for (int i = 0; i < N_BEAD; ++i) {
                // update positions 
                poly[c].r[i] += (DT / P.gamma_tilde) * poly[c].f[i] + sigma * gasdev(&rng); // update position
                double pos_new = poly[c].r[i];

                if (pos_new < 0.0) { // crossed the left wall 
                    poly[c].img[i] -= 1; 
                    pos_new += P.box_len_b; // wrap around 
                } else if (pos_new >= P.box_len_b) { // crossed the right wall 
                    poly[c].img[i] += 1; 
                    pos_new -= P.box_len_b; // wrap around 
                }

                poly[c].r[i] = pos_new; // store the wrapped position 
            }
        } 
        // junction beads 
        for (int J = 0; J < P.N_J; ++J) {
            // update junction positions
            junc[J].r += (DT / P.gamma_tilde) * junc[J].f + sigma * gasdev(&rng); // update position
            double pos_new = junc[J].r;

            if (pos_new < 0.0) { // crossed the left wall 
                junc[J].img -= 1; 
                pos_new += P.box_len_b; // wrap around 
            } else if (pos_new >= P.box_len_b) { // crossed the right wall 
                junc[J].img += 1; 
                pos_new -= P.box_len_b; // wrap around 
            }
            junc[J].r = pos_new; // store the wrapped position
        }    

    // diagnostics on contacts ------------------------------------------------- 
        rebuild_list_SW();
        double timeNow = tNow_prod; 
        for (int J = 0; J < P.N_J; ++J) {
            int nCand = has_contact_count(J);
            char cnow = (nCand > 0); 

            if (cnow && !in_contact[J]) { // window starts 
                in_contact[J] = 1; 
                contact_start[J] = timeNow; 
                windows_total[J]++;
            } else if (!cnow && in_contact[J]) { // window ends
                double contact_dur = timeNow - contact_start[J];
                contact_time[J] += contact_dur;
                if (contact_dur < (double) P.mc_stride * DT) {
                    windows_missed[J]++;
                }
                in_contact[J] = 0;
            }
        }        
    // Swap MC move ------------------------------------------------------------
        if (nBD % P.mc_stride == 0) {
            // list already rebuilt above 
            for (int J = 0; J < P.N_J; ++J) {
                int cand[MAX_CAND];
                int m = cand_swap(J, cand); // get candidates 
                ++win_prop_swap; 
                if (m > 0) { // at least one candidate 
                    ++win_tried_swap; 
                    ++cum_tried_swap; 

                    int old_state[P.arm]; for (int a = 0; a < P.arm; ++a) {old_state[a] = BonJ[J * P.arm + a];} // snapshot before the sweep 

                    if (m > 1) shuffle_int(cand, m, &rng); // randomize candidate order 

                    int randArm  = (int)(ran1(&rng) * P.arm); 

                    for (int k = 0; k < m; ++k) { // loop over candidates 
                        ++win_tried_swap_raw;
                        int status = attempt_swap_with_target(J, randArm, cand[k], &rng, P.Ea_tilde, P.ks_junc_tilde, P.T_tilde);

                        if (status == 1) { // one microscopic accept occurred 
                            ++win_acc_swap_raw; 
                            updateJType(J);
                        }
                    } // end loop over candidates

                    int new_state[P.arm]; for (int a = 0; a < P.arm; ++a) {new_state[a] = BonJ[J * P.arm + a];} // snapshot after the sweep 
                    
                    // Determine if there was a net change in junction connectivity 
                    int net_change = 0; for (int a = 0; a < P.arm; ++a) {if (new_state[a] != old_state[a]) { net_change = 1; break; }}

                    if (net_change) {
                        // Count ONE net accept for this junction this tick
                        ++win_acc_swap;            
                        ++cum_acc_swap; 
                        swaps_J[J]++;

                        // Log wait time ONCE per net event (skip during equilibration)
                        if (!in_equil) {
                            if (last_swap_time[J] >= 0.0) { 
                                double wait = tNow_prod - last_swap_time[J];
                                wait_sum[J]   += wait;
                                wait_sumsq[J] += wait * wait;
                                wait_n[J]++;
                            }
                            last_swap_time[J] = tNow_prod; // update reference for next interval
                        }
                    }
                }
            }
        } // end of swap MC move
        // Decide whether to start a new origin for ACF now (adaptive to acceptance activity) ------ 
        if (!in_equil && NUM_CHECKS > 0 && NUM_BONDS > 0) {
            int waited_enough = (nBD - last_start_nBD) >= MIN_ORG_GAP;

            // recent acceptance ratio in the current window (if none, check cumulative instead)
            double recent_ar = (win_tried_swap > 0) ? ((double)win_acc_swap / (double)win_tried_swap)
                                : ((cum_tried_swap > 0) ? ((double)cum_acc_swap / (double)cum_tried_swap) : 0.0); 


            // accepted swaps since the last origin
            int delta_accepts = cum_acc_swap - accepts_at_last_start;

            if (waited_enough && (delta_accepts >= ORG_ACC_CNT_THR || recent_ar >= ORG_AR_THR)) {
                StartPoint *S = &starts[next_start_slot];
                if (!S->bond_at_start) {
                    S->bond_at_start = (int*)malloc(NUM_BONDS * sizeof(int));
                    if (!S->bond_at_start) { perror("malloc bond_at_start"); exit(EXIT_FAILURE); }
                }
                if (!S->change_at_start) {
                    S->change_at_start = (unsigned int*)malloc(NUM_BONDS * sizeof(unsigned int));
                    if (!S->change_at_start) { perror("malloc change_at_start"); exit(EXIT_FAILURE); }
                }

                // snapshot both bead ids and change counters
                for (int v = 0; v < NUM_BONDS; ++v) {
                    S->bond_at_start[v]   = BonJ[v];
                    S->change_at_start[v] = arm_change_count[v];
                }

                S->active = 1;
                S->start_nBD = nBD;
                S->next_check_index = 0;
                S->id = start_id++;

                last_start_nBD = nBD;
                accepts_at_last_start = cum_acc_swap;

                next_start_slot = (next_start_slot + 1) % MAX_STARTS; // reuse oldest when full
            }
        } // end of new origin decision
        // Evaluate correlation points that are due at this step --------------- 
        if (!in_equil && NUM_CHECKS > 0 && NUM_BONDS > 0) {
            for (int si = 0; si < MAX_STARTS; ++si) {
                StartPoint *S = &starts[si];
                if (!S->active) continue;

                while (S->next_check_index < NUM_CHECKS) {
                    int i = S->next_check_index; 
                    int target_nBD = S->start_nBD + i * CHECK_EVERY_BD;
                    if (nBD < target_nBD) break; 
                    if (nBD > target_nBD) { S->next_check_index++; continue; } 

                    // strict survivor: same bead AND no change since the start
                    int survivors = 0;
                    for (int v = 0; v < NUM_BONDS; ++v) {
                        int same_pair = (BonJ[v] == S->bond_at_start[v]);
                        int unchanged = (arm_change_count[v] == S->change_at_start[v]);
                        if (same_pair && unchanged) survivors++;
                    }

                    // write one row for this start and this lag
                    double C_this = (double)survivors / (double)NUM_BONDS;     // fraction surviving
                    acf_sum[i]   += C_this;
                    acf_sumsq[i] += C_this * C_this;
                    acf_count[i] += 1;
                    if (i + 1 > acf_max_i_seen) acf_max_i_seen = i + 1;

                    S->next_check_index++;
                    if (S->next_check_index >= NUM_CHECKS) {
                        S->active = 0; // finished with this start
                    }
                }
            }
        }
        // Evaluate MSD points that are due at this step ----------------------- 
        if (!in_equil && MSD_NUM_CHECKS > 0) {
            if ((nBD - last_msd_start_nBD) >= MSD_ORIGIN_GAP) {
                MSDOrigin *slot = &msd_starts[next_msd_start_slot];
                if (slot->active && slot->next_check_index < MSD_NUM_CHECKS) {
                    fprintf(stderr, "[MSD] WARN: origin overwrite risk at slot %d (nBD=%lld)\n",
                            next_msd_start_slot, nBD);
                    // skip starting a new one this step
                } else {
                    start_msd_origin(nBD);
                }
            }

            double RcomNow_sys = 0.0;
            snapshot_COMs(&RcomNow_sys, RcomNow_chain); 

            for (int si = 0; si < MAX_MSD_STARTS; ++si) {
                MSDOrigin *O = &msd_starts[si];
                if (!O->active) continue;

                while (O->next_check_index < MSD_NUM_CHECKS) {
                    int i = O->next_check_index;
                    int target_nBD = O->start_nBD + i * MSD_CHECK_EVERY_BD;
                    if (nBD < target_nBD) break;
                    if (nBD > target_nBD) { O->next_check_index++; continue; }

                    const double dR = RcomNow_sys - O->comSys0_O;

                    // segmental MSD over sampled beads (drift-removed)
                    double msd_seg_this = 0.0;
                    int g_all = 0;
                    for (int c = 0; c < P.N_chain_tot; ++c) {
                        for (int i = 0; i < N_BEAD; ++i, ++g_all) {
                            double pos_unw = poly[c].r[i] + poly[c].img[i] * P.box_len_b;
                            double dx = (pos_unw - O->rBead0_O[g_all]) - dR;
                            msd_seg_this += dx * dx;
                        }
                    }
                    msd_seg_this /= (double)(P.N_chain_tot * N_BEAD);

                    // chain COM MSD (average over chains), drift-removed
                    double msd_com_this = 0.0;
                    for (int c = 0; c < P.N_chain_tot; ++c) {
                        double dxc = (RcomNow_chain[c] - O->comChain0_O[c]) - dR;
                        msd_com_this += dxc * dxc;
                    }
                    msd_com_this /= (double)P.N_chain_tot;

                    // accumulate
                    msd_seg_sum[i]   += msd_seg_this;
                    msd_seg_sumsq[i] += msd_seg_this * msd_seg_this;
                    msd_seg_cnt[i]   += 1UL;

                    msd_com_sum[i]   += msd_com_this;
                    msd_com_sumsq[i] += msd_com_this * msd_com_this;
                    msd_com_cnt[i]   += 1UL;

                    if (i + 1 > msd_max_i_seen) msd_max_i_seen = i + 1;

                    O->next_check_index++;
                    if (O->next_check_index >= MSD_NUM_CHECKS) { O->active = 0; }
                }
            }
        }

// Bookkeeping junction types -------------------------------------------------- 
        // cumulative count across simulation 
        for (int i = 0; i < P.N_J; ++i) {
            JuncDist[(P.arm + 1) * i + junc[i].type]++; // count junction types
        }
        // undapte CONNECT array 
        for (int i = 0; i < 2 * P.N_chain_tot; ++i) CONNECT[i] = -1; // reset first
        for (int j = 0; j < P.arm*P.N_J; ++j) {
            int b = BonJ[j]; // global bead index 
            int c, m; globalToChain(b, &c, &m); // get chain and monomer indices
            int J = j / P.arm; 
            if (m == 0) { // first bead of the chain
                CONNECT[2*c] = junc[J].type; 
            }
            if (m == N_BEAD - 1) { // last bead of the chain
                CONNECT[2*c + 1] = junc[J].type; 
            }
        }
        // update CIJ matrix 
        for (int c = 0; c < P.N_chain_tot; ++c) {
            if (poly[c].type[0] == 0) { // check the first bead type 
                int end0 = CONNECT[2*c]; 
                int end1 = CONNECT[2*c + 1]; 
                if (end0 == 0 && end1 == 0) CIJ[0]++;                               // n11
                if ((end0 == 0 && end1 == 1) || (end0 == 1 && end1 == 0)) CIJ[1]++; // n12
                if ((end0 == 0 && end1 == 2) || (end0 == 2 && end1 == 0)) CIJ[2]++; // n13
                if ((end0 == 0 && end1 == 3) || (end0 == 3 && end1 == 0)) CIJ[3]++; // n14 
                if (end0 == 1 && end1 == 1) CIJ[4]++;                               // n22
                if ((end0 == 1 && end1 == 2) || (end0 == 2 && end1 == 1)) CIJ[5]++; // n23
                if ((end0 == 1 && end1 == 3) || (end0 == 3 && end1 == 1)) CIJ[6]++; // n24
                if (end0 == 2 && end1 == 2) CIJ[7]++;                               // n33
                if ((end0 == 2 && end1 == 3) || (end0 == 3 && end1 == 2)) CIJ[8]++; // n34
                if (end0 == 3 && end1 == 3) CIJ[9]++;                               // n44 
            }
        }

// Output ---------------------------------------------------------------------- 
        if (nBD % 1 == 0) {
            // Compute PhiA and PhiB
            memset(PHIA, 0, NGRID * sizeof(double));
            memset(PHIB, 0, NGRID * sizeof(double));
            computePHI(PHIA, PHIB); 

            // Accumulate time-averaged density fields and degree of separation
            for (int i = 0; i < NGRID; ++i) {
                PHIA_sum[i] += PHIA[i];
                PHIB_sum[i] += PHIB[i];
            }
            sample_count++;
            // Write to file 
            if (sample_count == DUMP_AVG_SAMPLES) {
                double inv_count = 1.0 / (double) sample_count; 
                int n_start = nBD - DUMP_AVG_SAMPLES + 1;
                int n_end = nBD;
                int n_mid = (n_start + n_end) / 2;    

                // prefix: "Phi_avg_" 
                files[F_PHI_AVG] = fopen(strs[F_PHI_AVG], "a");
                if (files[F_PHI_AVG]) {
                    fprintf(files[F_PHI_AVG], "%d %d %d\n", n_start, n_mid, n_end); 
                    for (int i = 0; i < NGRID; ++i) {
                        fprintf(files[F_PHI_AVG], "%lf %lf %lf\n", 
                                PHIA_sum[i] * inv_count, 
                                PHIB_sum[i] * inv_count, 
                                (PHIA_sum[i] + PHIB_sum[i]) * inv_count);
                    }
                    fclose(files[F_PHI_AVG]);
                }    
                // reset accumulators
                sample_count = 0;
                memset(PHIA_sum, 0, NGRID * sizeof(double));
                memset(PHIB_sum, 0, NGRID * sizeof(double)); 
            }

            if (nBD % DUMP_PHI_STRIDE == 0){ 
                files[F_PHI] = fopen(strs[F_PHI], "w");
                // prefix: "Phi_" 
                if (files[F_PHI]) {
                    fprintf(files[F_PHI], "%lld\n", nBD);
                    for (int i = 0; i < MG; ++i) {
                        fprintf(files[F_PHI], "%lf %lf %lf\n", 
                                PHIA[i], PHIB[i], PHIA[i] + PHIB[i]);
                    }
                    fclose(files[F_PHI]);
                }
            }
        }

        // Coordinates 
        double box_x = P.box_len_b, box_y = 1.0, box_z = 1.0; // y and z dimensions are set to 1.0 [b] for 1D simulation

        // OPTIONAL: dump only ends to keep the file light. Comment this out to dump all beads.
        // #define DUMP_ONLY_ENDS 1

        if (nBD % DUMP_COORD_STRIDE == 0) {
            // prefix: "coord_"
            
            // Build local inert bead map
            int nTot = P.N_chain_tot * N_BEAD;
            unsigned char *is_inert = (unsigned char*)calloc(nTot, 1);
            if (!is_inert) { perror("calloc is_inert"); exit(EXIT_FAILURE); }
            for (int k = 0; k < P.N_NR; ++k) {
                int b = NREnd[k];
                if (b >= 0 && b < nTot) is_inert[b] = 1;
            }

            // Build local bead->junction maps
            int   *bead_to_J   = (int*)  malloc(nTot * sizeof(int));
            short *bead_to_arm = (short*)malloc(nTot * sizeof(short));
            if (!bead_to_J || !bead_to_arm) { perror("malloc maps"); exit(EXIT_FAILURE); }
            build_bead_to_junc_maps(bead_to_J, bead_to_arm);

            // Decide how many atoms to write
        #ifdef DUMP_ONLY_ENDS
            int nToWrite = 2 * P.N_chain_tot;           // only chain ends
        #else
            int nToWrite = nTot;                         // all beads
        #endif

            files[F_COORD] = fopen(strs[F_COORD], "a");
            if (files[F_COORD]) {
                fprintf(files[F_COORD], "ITEM: TIMESTEP\n%lld\n", nBD);
                fprintf(files[F_COORD], "ITEM: NUMBER OF ATOMS\n%d\n", nToWrite);
                fprintf(files[F_COORD], "ITEM: BOX BOUNDS pp pp pp\n");
                fprintf(files[F_COORD], "0.0 %lf\n", box_x);
                fprintf(files[F_COORD], "0.0 %lf\n", box_y);
                fprintf(files[F_COORD], "0.0 %lf\n", box_z);

                // id type x y z chain bead is_end is_inert is_bound junc_id junc_arm img x_unw
                fprintf(files[F_COORD],
                    "ITEM: ATOMS id type x y z chain bead is_end is_inert is_bound junc_id junc_arm img x_unw\n");

                int atom_id = 1;

        #ifdef DUMP_ONLY_ENDS
                for (int c = 0; c < P.N_chain_tot; ++c) {
                    for (int e = 0; e < 2; ++e, ++atom_id) {
                        int m = (e == 0) ? 0 : (N_BEAD - 1);
                        int g = c * N_BEAD + m;
                        int type   = poly[c].type[m];
                        long img   = poly[c].img[m];
                        double x   = poly[c].r[m];
                        double xunw= x + img * P.box_len_b;
                        int isEnd  = end_flag(m);
                        int inert  = is_inert[g] ? 1 : 0;
                        int Jid    = bead_to_J[g];
                        int Jarm   = bead_to_arm[g];
                        int bound  = (Jid >= 0) ? 1 : 0;

                        fprintf(files[F_COORD],
                            "%d %d %.9g 0.0 0.0 %d %d %d %d %d %d %d %ld %.9g\n",
                            atom_id, type, x, c, m, isEnd, inert, bound, Jid, Jarm, img, xunw);
                    }
                }
        #else
                for (int c = 0; c < P.N_chain_tot; ++c) {
                    for (int m = 0; m < N_BEAD; ++m, ++atom_id) {
                        int g = c * N_BEAD + m;
                        int type   = poly[c].type[m];
                        long img   = poly[c].img[m];
                        double x   = poly[c].r[m];
                        double xunw= x + img * P.box_len_b;
                        int isEnd  = end_flag(m);
                        int inert  = is_inert[g] ? 1 : 0;
                        int Jid    = bead_to_J[g];
                        int Jarm   = bead_to_arm[g];
                        int bound  = (Jid >= 0) ? 1 : 0;

                        fprintf(files[F_COORD],
                            "%d %d %.9g 0.0 0.0 %d %d %d %d %d %d %d %ld %.9g\n",
                            atom_id, type, x, c, m, isEnd, inert, bound, Jid, Jarm, img, xunw);
                    }
                }
        #endif

                fclose(files[F_COORD]);
            }

            free(is_inert);
            free(bead_to_J);
            free(bead_to_arm);
        } // end of coordinate dump 

        if (nBD % DUMP_LOG_STRIDE == 0) {
            // Total contact time (include ongoing windows’ partial time)
            double Tcontact_sum = 0.0;
            long long windows_total_sum = 0, windows_missed_sum = 0;
            for (int J = 0; J < P.N_J; ++J) {
                double T = contact_time[J];
                if (in_contact[J]) { // if junction J has a window still open 
                    T += (tNow_prod - contact_start[J]); // add the open window's partial duration 
                }
                Tcontact_sum += T;
                windows_total_sum  += windows_total[J];
                windows_missed_sum += windows_missed[J];
            }

            double k_swap_bulk = (tBD_prod > 0.0 && P.N_J > 0) ? (double)cum_acc_swap / tBD_prod : 0.0; // how many accepted swaps per unit time in the bulk (across all junctions) 
            double k_swap_perJ = (P.N_J > 0) ? k_swap_bulk / (double) P.N_J : 0.0; 
            double k_intrinsic = (Tcontact_sum > 0.0) ? (double)cum_acc_swap / Tcontact_sum : 0.0; // how many accepted swaps per junction per time given that the junction is in contact with a candidate 
            double f_miss = (windows_total_sum > 0) ? (double)windows_missed_sum / (double)windows_total_sum : 0.0; 

            int Ntot = P.N_chain_tot * N_BEAD; 
            double neighbor_avail = win_prop_swap ? (double) win_tried_swap / win_prop_swap  : 0.0;
            double acc_ratio_swap = win_tried_swap ? (double) win_acc_swap / win_tried_swap : 0.0;                        
            double acc_ratio_swap_raw = win_tried_swap_raw ? (double) win_acc_swap_raw / win_tried_swap_raw : 0.0;
            double MeanSqRee = 0.0, F_tot = 0.0, R_cm = 0.0, R_cm_unw = 0.0; 
            for (int c = 0; c < P.N_chain_tot; ++c) {
                for (int b = 0; b < N_BEAD; ++b) {
                    F_tot += poly[c].f[b]; 
                    R_cm += poly[c].r[b];
                    R_cm_unw += poly[c].r[b] + poly[c].img[b] * P.box_len_b; 
                }
                double dx = 0.0; 
                MeanSqRee += rsq(poly[c].r[0], poly[c].r[N_BEAD-1], &dx);
            }
            for (int J = 0; J < P.N_J; ++J) {
                F_tot += junc[J].f; 
            }

            MeanSqRee /= (double)P.N_chain_tot; 
            R_cm /= (double)Ntot;
            R_cm_unw /= (double)Ntot;

            rebuild_list_NB(); 
            double Ub = 0.0, Unb = 0.0, Utot = 0.0; 
            double pref = (double) P.sqrtNbar / (double) P.Re_b; 

            long long cand_pairs = 0, near_pairs = 0; 
            double weighted_pairs = 0.0; 

            double W_unlike = 0.0, W_all = 0.0; 

            for (int g1 = 0; g1 < Ntot; ++g1) {
                int c1, m1; globalToChain(g1, &c1, &m1); 
                double xi = poly[c1].r[m1];
                int ti = poly[c1].type[m1];

                // bonded neighbor intra-chain energy 
                if (m1 < N_BEAD - 1) {
                    double dx; double r2 = rsq(xi, poly[c1].r[m1 + 1], &dx);
                    Ub += 0.5 * P.ks_mon_tilde * r2; 
                }

                // non-bonded energy
                int cx0 = cellIndex_1D(xi, P.N_cell_nb, P.L_cell_nb_b);
                for (int dcx = -1; dcx <= 1; ++dcx) {
                    int cx = (cx0 + dcx + P.N_cell_nb) % P.N_cell_nb; // wrap around
                    
                    // traverse the linked-cell list in cell cx
                    for (int g2 = cellHead_NB[cx]; g2 != -1; g2 = cellNext_NB[g2]) {
                        if (g2 <= g1) continue; // avoid double counting and self-interaction 

                        cand_pairs++; // count candidate pairs

                        int c2, m2; globalToChain(g2, &c2, &m2);     
                        double xj = poly[c2].r[m2];
                        int    tj = poly[c2].type[m2];

                        double dx; rsq(xi, xj, &dx);
                        double w = overlap_1d(dx); 
                        W_all += w; 
                        if (ti != tj) {W_unlike += w;}

                        double adx = fabs(dx); 
                        if (adx < P.deltaL_b) {
                            near_pairs++; // count near pairs
                            double weight = (P.deltaL_b - adx) / P.deltaL_b;
                            weighted_pairs += weight;
                        }

                        Unb += P.kBT_tilde * pref * pair_energy(xi, xj, ti, tj);
                    }
                }
            }
            Utot = Ub + Unb; // total energy

            double nn_avg = (Ntot > 0) ? 2.0 * (double) near_pairs / (double) Ntot : 0.0; 
            double weighted_nn_avg = (Ntot > 0) ? 2.0 * weighted_pairs / (double) Ntot : 0.0;
            double list_eff = (cand_pairs > 0) ? (double) near_pairs / (double) cand_pairs : 0.0; 
            double degMix = (W_all > 0.0) ? (W_unlike / (2.0 * P.fA * P.fB * W_all)) : 0.0; 
            double degSep = 1.0 - degMix; 

            if (degSep >= 0.0) {
                degSep_sum += degSep; 
                degSep_count++; 
            }

            if (degSep_count == DUMP_DEGSEP_AVG_STRIDE) {
                int n_end = nBD; 
                int span = degSep_count * DUMP_LOG_STRIDE;
                int n_start = n_end - span + 1;
                int n_mid = (n_start + n_end) / 2;

                double degSep_avg = degSep_sum / (double) degSep_count;
                files[F_DEGSEP_AVG] = fopen(strs[F_DEGSEP_AVG], "a");
                // prefix: "degSep_avg_" 
                if (files[F_DEGSEP_AVG]) {
                    fprintf(files[F_DEGSEP_AVG], "%d %d %d %lf\n", n_start, n_mid, n_end, degSep_avg);
                    fclose(files[F_DEGSEP_AVG]);
                }
                degSep_sum = 0.0;
                degSep_count = 0;
            }


            int allAjunc = 0, allBJunc = 0;
            for (int J = 0; J < P.N_J; ++J) {
                if (junc[J].type == 0) allAjunc += 1; 
                if (junc[J].type == P.arm) allBJunc += 1; 
            }
            double allAjunc_frac = (P.N_J > 0) ? (double) allAjunc / (double) P.N_J : 0.0;
            double allBJunc_frac = (P.N_J > 0) ? (double) allBJunc / (double) P.N_J : 0.0;
            double alljunc_frac = allAjunc_frac + allBJunc_frac; 

            // prefix: "log_" 
            files[F_LOG] = fopen(strs[F_LOG], "a");
            if (files[F_LOG]) {
                double MSD_seg = 0.0, MSD_com = 0.0;
                if (!in_equil) { // reference frame has been set
                    compute_MSD_undrifted(&MSD_seg, &MSD_com, r0, Rcom0_chain, Rcom0_sysRef, RcomNow_chain); 
                } else {
                    MSD_seg = MSD_com = 0.0; 
                }

                double mean_bond2 = mean_bond_r2(); 
                double mean_Jbond2 = mean_jb_r2(); 
                // Remember to modify header in `Headers` section ----------------------------------
                fprintf(files[F_LOG], 
                    "%lld,"  "%d,"       "%lf,"     "%lf,"      "%d,"        "%lf,"        "%lf,"                            "%lf,"    "%lf,"    "%lf,"    "%lf,"     "%lf,"     "%lf,"     "%lf,"      "%lf,"     "%lf,"        "%lf,"        "%lf,"      "%lf,"             "%lf,"      "%lf,"     "%lf,"            "%e,"            "%e,"                 "%e,"          "%e,"          "%e,"          "%lf,"      "%lf,"           "%lf,"           "%lf\n",                                                  
                    nBD,   nBD_prod,   tBD_abs,   tBD_prod,   in_equil,    MeanSqRee,    MeanSqRee/(double)(N_BEAD-1.0),   Ub,       Unb,      Utot,     F_tot,     R_cm,      R_cm_unw,   MSD_seg,   MSD_com,   mean_bond2,   mean_Jbond2,   nn_avg,    weighted_nn_avg,    list_eff,   degSep,   neighbor_avail,   acc_ratio_swap,  acc_ratio_swap_raw,   k_swap_bulk,   k_swap_perJ,   k_intrinsic,   f_miss,     allAjunc_frac,   allBJunc_frac,   alljunc_frac);                                                                                            
                fclose(files[F_LOG]);
                // --------------------------------------------------------------------------------- 
            }
            // zero for the next window 
            win_prop_swap = win_tried_swap = win_acc_swap = 0; 
            win_tried_swap_raw = 0, win_acc_swap_raw = 0;
        } // end of nBD % DUMP_LOG_STRIDE == 0 

        if (!in_equil && nBD % DUMP_LOG_STRIDE == 0 && acf_max_i_seen > 0) {
            files[F_BOND_ACF] = fopen(strs[F_BOND_ACF], "w"); 
            if (files[F_BOND_ACF]) {
                fprintf(files[F_BOND_ACF], "nBD_prod,tBD_prod,i,tlag,mean,sem,N\n");

                for (int i = 0; i < acf_max_i_seen; ++i) {
                    double N   = (double)acf_count[i];
                    double mu, sem; 
                    mean_sem_from_sums(acf_sum[i], acf_sumsq[i], N, &mu, &sem); 
                    double tl  = (double)(i * CHECK_EVERY_BD) * DT;

                    fprintf(files[F_BOND_ACF], "%d,%.9e,%d,%.9e,%.9e,%.9e,%lu\n",
                            nBD_prod, tBD_prod, i, tl, mu, sem, (unsigned long)acf_count[i]);
                }
                fclose(files[F_BOND_ACF]); 
            }
        }

        if (!in_equil && nBD % DUMP_LOG_STRIDE == 0 && msd_max_i_seen > 0) {
            files[F_MSD] = fopen(strs[F_MSD], "w");
            if (files[F_MSD]) {
                fprintf(files[F_MSD], "nBD_prod,tBD_prod,i,tlag,msd_seg_mean,msd_seg_sem,N_seg,msd_com_mean,msd_com_sem,N_com\n");
                for (int i = 0; i < msd_max_i_seen; ++i) {
                    double Nseg = (double)msd_seg_cnt[i];
                    double seg_mu, seg_sem;
                    mean_sem_from_sums(msd_seg_sum[i], msd_seg_sumsq[i], Nseg, &seg_mu, &seg_sem);

                    double Ncom = (double)msd_com_cnt[i];
                    double com_mu, com_sem;
                    mean_sem_from_sums(msd_com_sum[i], msd_com_sumsq[i], Ncom, &com_mu, &com_sem);

                    double tl = (double)(i * MSD_CHECK_EVERY_BD) * DT;
                    fprintf(files[F_MSD], "%d,%.9e,%d,%.9e,%.9e,%.9e,%lu,%.9e,%.9e,%lu\n",
                            nBD_prod, tBD_prod, i, tl,
                            seg_mu, seg_sem, (unsigned long)msd_seg_cnt[i],
                            com_mu, com_sem, (unsigned long)msd_com_cnt[i]);
                }
                fclose(files[F_MSD]);
            }
        }

        if (nBD % DUMP_DISTRIB_STRIDE == 0) {
            // prefix: "distr_"
            files[F_DISTRIB] = fopen(strs[F_DISTRIB], "w");
            if (files[F_DISTRIB]) {
                fprintf(files[F_DISTRIB], "%lld\n", nBD);
                fprintf(files[F_DISTRIB], "Junction type distribution:\n"); 
                for (int i = 0; i < (P.arm + 1); ++i) {
                    for (int j = 0; j < P.N_J; ++j) {
                        fprintf(files[F_DISTRIB], "%d ", JuncDist[i + j * (P.arm + 1)]);
                    }
                    fprintf(files[F_DISTRIB], "\n");
                }
                fprintf(files[F_DISTRIB], "CIJ matrix:\n");
                for (int i = 0; i < nUpperTri; ++i) {
                    fprintf(files[F_DISTRIB], "%d ", CIJ[i]);
                }
                fclose(files[F_DISTRIB]);
            }
        }

        if (nBD % DUMP_SITE_COORD_STRIDE == 0) {
            // prefix: "dynSite_" 
            files[F_DYNSITE] = fopen(strs[F_DYNSITE], "w");
            if (files[F_DYNSITE]) {
                fprintf(files[F_DYNSITE], "# nBD\n"); 
                fprintf(files[F_DYNSITE], "%lld\n", nBD);
                // Free end snapshots 
                fprintf(files[F_DYNSITE], "# Free ends (index, type, x):\n");
                for (int fe = 0; fe < P.N_FE; ++fe) {
                    int b = FEnd[fe]; // global bead index of the free end 
                    int c, m; globalToChain(b, &c, &m); // get chain and monomer indices
                    fprintf(files[F_DYNSITE], "%d %d %lf\n", 
                            fe, poly[c].type[m], poly[c].r[m]);
                }
                fprintf(files[F_DYNSITE], "\n");
                // Junction snapshots 
                fprintf(files[F_DYNSITE], "# Junctions (index, type, x):\n");
                for (int J = 0; J < P.N_J; ++J) {
                    fprintf(files[F_DYNSITE], "%d %d %lf\n", 
                            J, junc[J].type, junc[J].r);
                } 
                fprintf(files[F_DYNSITE], "\n");
                // inert end snapshots 
                fprintf(files[F_DYNSITE], "\n# Non-reactive ends (index, type, x):\n");
                for (int nr = 0; nr < P.N_NR; ++nr) {
                    int b = NREnd[nr]; // global bead index of the non-reactive end 
                    int c, m; globalToChain(b, &c, &m); // get chain and monomer indices
                    fprintf(files[F_DYNSITE], "%d %d %lf\n", 
                            nr, poly[c].type[m], poly[c].r[m]);
                }
                fclose(files[F_DYNSITE]);
            }
        }

        if (!in_equil) {
            if (nBD % DUMP_WAIT_STRIDE == 0) {
                // prefix: "waitTime_"
                const double t_anchor = (nBD_prod + 1) * DT; // end-of-step prod time 
                // overall event-weighted mean
                unsigned long long Nswap_tot = 0ULL; 
                double S = 0.0, S2 = 0.0; 
                for (int J = 0; J < P.N_J; ++J) { 
                    S += wait_sum[J]; 
                    S2 += wait_sumsq[J];
                    Nswap_tot += wait_n[J];
                }
                double mean_all, sem_wait; 
                mean_sem_from_sums(S, S2, (double)Nswap_tot, &mean_all, &sem_wait);
                double khat_all = (mean_all > 0.0 ? 1.0 / mean_all : NAN);
                if (Nswap_tot > 1ULL) {
                    double N   = (double)Nswap_tot;
                    double var = (S2 - (S*S)/N) / (N - 1.0);   // unbiased sample variance
                    if (var < 0.0) var = 0.0;                  // numeric guard
                    sem_wait = sqrt(var / N);
                }                

                files[F_WAIT_TIME] = fopen(strs[F_WAIT_TIME], "a");
                if (files[F_WAIT_TIME]) {
                    fprintf(files[F_WAIT_TIME], 
                            "%d,"        "%.9e,"    "%llu,"                          "%d,"    "%.9e,"     "%.9e,"     "%.9e\n",
                            nBD_prod,    tBD_prod,  (unsigned long long)Nswap_tot,   P.N_J,   mean_all,   sem_wait,   khat_all);
                    fclose(files[F_WAIT_TIME]);
                }
                for (int J = 0; J < P.N_J; ++J) {
                    last_swap_time[J] = t_anchor; 
                }
                // reset for the next window 
                memset(wait_sum, 0, P.N_J*sizeof(double));
                memset(wait_sumsq, 0, P.N_J*sizeof(double));
                memset(wait_n, 0, P.N_J*sizeof(unsigned long));
            }
        }

        // Reset CIJ every 10000 timesteps
        if (nBD % 10000 == 0) {
            for (int i = 0; i < nUpperTri; ++i) {
                CIJ[i] = 0;
            }
        }         
         
    } // end of time loop 

    // cleanup ----------------------------------------------------------------- 
    for (int i = 0; i < NUM_STR; ++i) {
        free(strs[i]);
    }
    free(poly); 
    free(junc); 
    free(PHIA);
    free(PHIB);
    free(PHIA_sum);
    free(PHIB_sum);
    free(BonJ); 
    free(FEnd); 
    free(NREnd);
    free(r0);
    free(Rcom0_chain);
    free(cellHead_SW); 
    free(cellNext_SW);
    free(cellHead_NB); 
    free(cellNext_NB);
    free(JuncDist);
    free(CONNECT); 
    free(CIJ); 
    free(swaps_J);
    free(in_contact);
    free(contact_start);
    free(contact_time);
    free(windows_total);
    free(windows_missed);
    free(last_swap_time);
    free(wait_sum);
    free(wait_sumsq);
    free(wait_n);
    for (int i = 0; i < MAX_STARTS; ++i) {
        free(starts[i].bond_at_start);
        free(starts[i].change_at_start);
    }
    free(arm_change_count);
    free(acf_sum);
    free(acf_sumsq);
    free(acf_count);
    free(msd_seg_sum);
    free(msd_seg_sumsq);
    free(msd_seg_cnt);
    free(msd_com_sum);
    free(msd_com_sumsq);
    free(msd_com_cnt);
    free(RcomNow_chain); 
    free(msd_r0_pool);
    free(msd_com0_pool);

    return 0;
}