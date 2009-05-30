REAL S_matrix[NX+5][NY+4] ;
REAL Stmp_matrix[NX+5][NY+4] ;
REAL e_matrix[NX+4][NY+4] ;
REAL vx_matrix[NX+4][NY+4] ;
REAL vy_matrix[NX+4][NY+4] ;
REAL pot_matrix[NX+4][NY+4] ;
REAL fl_matrix[NX+4][NY+4] ;
REAL Sneq_matrix[NX+4][NY+4] ;
REAL vyneq_matrix[NX+4][NY+4] ;
REAL eneq_matrix[NX+4][NY+4] ;
REAL N_matrix[NX+4][NY+4] ;
REAL work1_matrix[NX+4][NY+4] ;
REAL work2_matrix[NX+4][NY+4] ;
REAL work3_matrix[NX+4][NY+4] ;
REAL work4_matrix[NX+4][NY+4] ;
REAL work5_matrix[NX+4][NY+4] ;
REAL work6_matrix[NX+4][NY+4] ;

REAL (*S)[NY+4] ;
REAL (*Stmp)[NY+4] ;
REAL (*e)[NY+4] ;
REAL (*vx)[NY+4] ;
REAL (*vy)[NY+4] ;
REAL (*pot)[NY+4] ;
REAL (*fl)[NY+4] ;
REAL (* Sneq)[NY+4] ;
REAL (* eneq)[NY+4] ;
REAL (* vyneq)[NY+4] ;
REAL (*N)[NY+4] ;
REAL (*work1)[NY+4] ;
REAL (*work2)[NY+4] ;
REAL (*work3)[NY+4] ;
REAL (*work4)[NY+4] ;
REAL (*work5)[NY+4] ;
REAL (*work6)[NY+4] ;

double cour ;
double cs ;
double gam,g2,deltag ;
double dx,dy ;
double dt ;
double t,tf ;
double nu_vnr ;
double nu_l ;
double G ;
double cool_fac ;

int i_crit_dt ;
int j_crit_dt ;
int i_crit_tc ;
int j_crit_tc ;
int dt_crit ;

double DTd ;
double DTi ;
double DTl ;
double res ;
double nu_sh ;
double mu ;

/* shearing sheet terms */
double W ;
double omega ;
double q ;
double Lx ;
double Ly ;
double K ;
double S0,cs0,u0,alphag,alphah ;
double kxeq,kyeq,kx0,ky0,m0,n0,k0 ;
double kmax,ksqmax ;
double fx,fy ;
double Qeq ;
double pamp ;
int is_restart ;
double h0,phi,eps ;
double N2min,Xmin ;
double dampinit,damp ;
double a,ar ;
double vxr,vxi,vyr,vyi,sr,si,ur,ui ;

int dump_cnt ;
int im_cnt ;

