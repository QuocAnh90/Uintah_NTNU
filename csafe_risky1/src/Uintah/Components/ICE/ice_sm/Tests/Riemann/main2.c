 #include "Riemann.c"

#define switchDebug_main_custom 1 
/* 
======================================================================*/
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <stdlib.h>
#include "nrutil+.h"
#include "functionDeclare.h"
#include "parameters.h"
#include "switches.h"
#include "macros.h"
#include "cpgplot.h"            /*must have this for plotting to work   */

#include <ieeefp.h>            /* needed by Steve Parker's malloc Library*/
/* ---------------------------------------------------------------------
GENERAL INFORMATION
 Function:  main--Main program
 Filename:  main.c 
 Purpose:    This is the main program for the Uintah ICE cfd code. 

History: 
Version   Programmer         Date       Description                      
     -------   ----------         ----       -----------                 
        1.0     Todd Harman       02/22/99                               
                                                                    
    Programming Conventions
        i, j, k         Loop indices for the x, y, z directions respectively
        f               is a loop index for face-centered values.
        m               Loop index for the different materials

                                 ________ 
                                /  1    /|
                               /_______/ |
                              |       | ______(3)
                       (4)____| I,J,K |  |     
                              |       | /      
                              |_______|/
                                  |               (6) = back face
                                 (2)              (5) = front face

 STEPS:
    - Set some eviromnental variables required for PGPLOT
    - Initialize some variables that are mainly used in testing
    - MEMORY SECTION: Allocate the memory needed for all of the arrays
      For all of the face-centered arrays set equate the common face addresses
      [i][j][k][RIGHT][m] = [i-1][j][k][LEFT][m]
    - PROBLEM INITIALIZATION SECTION: Read in the input file, test the inputs,
      set the boundary condtions, generate the grid
    - MAIN LOOP
        to be filled in
_____________________________________________________________________*/ 
main()
{
    int i, j, k, m, 
        xLoLimit,                       /* x array lower limits             */
        yLoLimit,                       /* y array lower limits             */
        zLoLimit,
        xHiLimit,
        yHiLimit,
        zHiLimit,
        printSwitch,
        should_I_write_output,          /* flag for dumping output          */
        fileNum,                        /* tecplot file number              */
        stat,                           /* status of putenv and getenv      */
        **BC_inputs,                    /* BC_types[wall][m] that contains  */
                                        /* the users boundary condition     */
                                        /* selection for each wall          */
        ***BC_types,                    /* each variable can have a Neuman, */
                                        /* or Dirichlet type BC             */
                                        /* BC_types[wall][variable][m]=type */
        ***BC_float_or_fixed,           /* BC_float_or_fixed[wall][variable][m]*/
                                        /* Variable on boundary is either   */
                                        /* fixed or it floats during the    */
                                        /* compuation                       */
        nMaterials;                     /* Number of materials              */

                                        
/* ______________________________   
*  Geometry                        
* ______________________________   */     
     double  delX,                      /* Cell width                       */
             delY,                      /* Cell Width in the y dir          */
             delZ,                      /* Cell width in the z dir          */
             delt,                      /* time step                        */
             CFL,                       /* Courant-Friedrichs and Lewy      */
             t_final,                   /* final problem time               */
            *t_output_vars,             /* array holding output timing info */
                                        /* t_output_vars[1] = t_initial     */
                                        /* t_output_vars[2] = t final       */
                                        /* t_output_vars[3] = delta t       */
            *delt_limits,               /* delt_limits[1]   = delt_minimum  */
                                        /* delt_limits[2]   = delt_maximum  */
                                        /* delt_limits[3]   = delt_initial  */
             t,                         /* current time                     */
             ***x_CC,                   /* x-coordinate of cell center      */
             ***y_CC,                   /* y-coordinate of cell center      */
             ***z_CC,                   /* z-coordinate of cell center      */
             ***Vol_CC,                 /* vol of the cell at the cellcenter*/
                                        /* (x, y, z)                        */
            /*------to be treated as pointers---*/
             *****x_FC,                 /* x-coordinate of face center      */
                                        /* x_FC(i,j,k,face)                 */
                                        /* cell i,j,k                       */
             *****y_FC,                 /* y-coordinate of face center      */
                                        /* y_FC(i,j,k,face)                 */
                                        /* of cell i,j,k                    */
             *****z_FC;                 /* z-coordinate of face center      */
                                        /* z_FC(i,j,k,face)                 */
                                        /* of cell i,j,k                    */
            /*----------------------------------*/ 

/* ______________________________   
*  Cell-centered and Face centered                      
* ______________________________   */ 
    double                              /* (x,y,z,material                  */
            ****uvel_CC,                /* u-cell-centered velocity         */
            ****vvel_CC,                /*  v-cell-centered velocity        */
            ****wvel_CC,                /* w cell-centered velocity         */
            ****delPress_CC,            /* cell-centered change in pressure */                                                                                       
            ****press_CC,               /* Cell-centered pressure           */
            ****Temp_CC,                /* Cell-centered Temperature        */
            ****rho_CC,                 /* Cell-centered density            */
            ****viscosity_CC,           /* Cell-centered Viscosity          */ 
            ****thermalCond_CC,         /* Cell-centered thermal conductivity*/
            ****cv_CC,                  /* Cell-centered specific heat      */ 
            ****mass_CC,                /* total mass, cell-centered        */
            ****xmom_CC,                /* x-dir momentum cell-centered     */
            ****ymom_CC,                /* y-dir momentum cell-centered     */
            ****zmom_CC,                /* z-dir momentum cell-centered     */
            ****int_eng_CC,             /* Internal energy cell-centered    */
            ****total_eng_CC,           /* Total energy cell-centered       */  
            ****div_velFC_CC,           /* Divergence of the face centered  */
                                        /* velocity that lives at CC        */    
            ****scalar1_CC,             /* Cell-centered scalars            */   
            ****scalar2_CC,             /* (x, y, z, material)              */
            ****scalar3_CC,
            /*------to be treated as pointers---*/
                                        /*______(x,y,z,face, material)______*/
            ******uvel_FC,              /* u-face-centered velocity         */
            ******vvel_FC,              /* *v-face-centered velocity        */
            ******wvel_FC,              /* w face-centered velocity         */
            ******press_FC,             /* face-centered pressure           */
            ******tau_X_FC,             /* *x-stress component at each face */
            ******tau_Y_FC,             /* *y-stress component at each face */
            ******tau_Z_FC,             /* *z-stress component at each face */
            /*----------------------------------*/                                              
           *grav,                       /* gravity (dir)                    */
                                        /* x-dir = 1, y-dir = 2, z-dir = 3  */
            *gamma,
            ****speedSound;             /* speed of sound (x,y,z, material) */    
             
/* ______________________________   
*  Lagrangian Variables            
* ______________________________   */ 
    double                              /*_________(x,y,z,material)_________*/
            ****rho_L_CC,               /* Lagrangian cell-centered density */
            ****mass_L_CC,              /* Lagrangian cell-centered mass    */
            ****Temp_L_CC,              /* Lagrangian cell-centered Temperature */
            ****press_L_CC,             /* Lagrangian cell-centered pressure*/            
            ****xmom_L_CC,              /* Lagrangian cell-centered momentum*/
            ****ymom_L_CC,              /* Lagrangian cell-centered momentum*/
            ****zmom_L_CC,              /* Lagrangian cell-centered momentum*/
            ****int_eng_L_CC,           /* Lagrangian cc internal energy    */
            ****Vol_L_CC,               /* Lagrangian cell-centered volume  */
/*__________________________________
* source terms
*___________________________________*/
            ****mass_source,            /* Mass source term (x,y,z, material */
            ****xmom_source,            /* momentum source terms            */
                                        /* (x, y, z, material)              */
            ****ymom_source,
            ****zmom_source,
            ****int_eng_source,         /* internal energy source           */
/*__________________________________
*   MISC Variables
*___________________________________*/            
            ***BC_Values,                /* BC values BC_values[wall][variable][m]*/  
            *R;                         /* gas constant R[material]          */ 

    char    output_file_basename[30],   /* Tecplot filename description     */
            output_file_desc[50];       /* Title used in tecplot stuff      */
/*__________________________________
*   Variables specific to this the Riemann
*   test problem
*___________________________________*/
/**/    double       
/**/            u0,
/**/            u1,     u4,
/**/            rho1,   rho4,
/**/            p1,     p4,
/**/            a1,     a4,
/**/            *u_Rieman,
/**/            *a_Rieman,
/**/            *p_Rieman,
/**/            *rho_Rieman,
/**/            *T_Rieman,
/**/            delQ;
/**/    int     xLo, yLo, zLo,
/**/            xHi, yHi, zHi,
/**/            qLoLimit, qHiLimit,
/**/            Q_MAX_LIM;
/**/            
/**/    int     counter;

/*__________________________________
* Variables specific to the blast
* wave problem
* Note the units of E are
*   M L^2 T^-2  Spherical case
*   M L   T^-2  Cylindrical case
*   M     T^-2  Plane case. 
*   This test problem is valid if pinitial = 0.0 
*
* From the graph on pg 231
*   alpha = 1.175   Spherical case
*         = 1.375  Cylindrical case
*         = 1.5    Plane case. 
*___________________________________*/                             
    double  residual,
/**/        E,
/**/        E0,
/**/        alpha,                      /* from graph on pg231             */
/**/        rho0,
/**/        speedSound0,
/**/        ***vvel_exact,              /* exact velocity ydir              */
/**/        ***uvel_exact,              /* exact velocity xdir              */
/**/        ***rho_exact,               /* exact density                    */
/**/        ***press_exact,             /* exact pressure                   */
/**/        ***Temp_exact,              /* exact Temperature K              */
/**/        blastWave_radius;
/**/double  uvel_initial,
            vvel_initial,
            rho_initial,
            press_initial,
            Temp_initial;
            
/**/int     ignition_pt_xdir,           /* x location of where E was added  */
/**/        ignition_pt_ydir,           /* y location of where E was added  */
/**/        BW_radius_cell;             /* cell location where the blast wave*/
/**/                                    /* is in.  Used mainly in plotting  */
/*__________________________________
*   Plotting variables
*___________________________________*/
#if (switchDebug_main_custom == 1|| switchDebug_main_custom == 2 || switchDebug_main_input == 1)
    #include "plot_declare_vars.h"   
#endif
    stat = putenv("PGPLOT_DIR=" PGPLOT_DIR);
    stat = putenv("PGPLOT_I_AM_HERE=0");              
                                        /* tell the plotting routine that  */
                                        /* you're at the top of main       */      

    stat = putenv("PGPLOT_PLOTTING_ON_OFF=1");
    stat = putenv("PGPLOT_OPEN_NEW_WINDOWS=1");  

/*______________________________________________________________________
*   Initialize variables
*_______________________________________________________________________*/ 
/**/    printSwitch = 1;    
/**/    t           = 0.0;  
/**/    m           = 1;
/**/    fileNum     = 1;
/**/    
/**/    u4          = 0.0;                      /* left chamber                     */
/**/    p4          = 100.0*1000.0;             /* Pa                               */
/**/    rho4        = 1.0;
/**/    a4          = sqrt(1.4 * p4 /rho4);
/**/    
/**/    u1          = 0.0;                      /* Right chamber                    */
/**/    p1          = 10.0 * 1000.0;            /* Pa                               */    
/**/    rho1        = 0.125;
/**/    a1          = sqrt(1.4 * p1 /rho1);  
/**/    Q_MAX_LIM   = 105;
       
/*______________________________________________________________________ 
*    M  E  M  O  R  Y     S  E  C  T  I  O  N 
*   - Allocate memory for the arrays                                          
*_______________________________________________________________________*/
/**/    u_Rieman    = dvector_nr(0, Q_MAX_LIM);
/**/    a_Rieman    = dvector_nr(0, Q_MAX_LIM);
/**/    p_Rieman    = dvector_nr(0, Q_MAX_LIM);
/**/    rho_Rieman  = dvector_nr(0, Q_MAX_LIM);
/**/    T_Rieman    = dvector_nr(0, Q_MAX_LIM);
#include "allocate_memory.i"
 /*__________________________________
* Needed by Steve Parkers Malloc library
*___________________________________*/    
    /* fpsetmask(FP_X_UFL|FP_X_OFL|FP_X_DZ|FP_X_INV); */
    /* audit(); */                      /* Steve Parkers memory tool        */
/*______________________________________________________________________
*
*  P  R  O  B  L  E  M     I  N  I  T  I  A  L  I  Z  A  T  I  O  N  
*  - read input file
*   - test the input variables
*   - Equate the address of the face centered variables
*   - Generate a grid
*   - zero all of the face-centered arrays
*   
*                  
* -----------------------------------------------------------------------  */
                                        
       readInputFile(   &xLoLimit,      &yLoLimit,      &zLoLimit,     
                        &xHiLimit,      &yHiLimit,      &zHiLimit,
                        &delX,          &delY,          &delZ,
                        uvel_CC,        vvel_CC,        wvel_CC, 
                        Temp_CC,        press_CC,       rho_CC,
                        scalar1_CC,     scalar2_CC,     scalar3_CC,
                        viscosity_CC,   thermalCond_CC, cv_CC,
                        R,              gamma,
                        &t_final,       t_output_vars,  delt_limits,
                        output_file_basename,           output_file_desc,       
                        grav,           speedSound,
                        BC_inputs,      BC_Values,      &CFL,
                        &nMaterials);      
    
    testInputFile(      xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delX,           delY,           delZ,
                        Temp_CC,        press_CC,       rho_CC,
                        viscosity_CC,   thermalCond_CC, cv_CC,
                        speedSound,      
                        t_final,        t_output_vars,  delt_limits,
                        BC_inputs,      printSwitch,    CFL,
                        nMaterials); 
                   
    definition_of_different_physical_boundary_conditions(              
                        BC_inputs,      BC_types,       BC_float_or_fixed,
                        BC_Values,      nMaterials  );  
                        
/*__________________________________
* Now make sure that the face centered
* values know about each other.
* for example 
* [i][j][k][RIGHT][m] = [i-1][j][k][LEFT][m]
*___________________________________*/  

    equate_ptr_addresses_adjacent_cell_faces(              
                        x_FC,           y_FC,           z_FC,
                        uvel_FC,        vvel_FC,        wvel_FC,
                        press_FC,
                        tau_X_FC,       tau_Y_FC,       tau_Z_FC,
                        nMaterials);   

    /*__________________________________
    * Generate a grid
    *___________________________________*/ 
    generateGrid(       xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delX,           delY,           delZ,
                        x_CC,           y_CC,           z_CC,   Vol_CC,  
                        x_FC,           y_FC,           z_FC );
    /*__________________________________
    *   zero the face-centered arrays
    *___________________________________*/
    zero_arrays_6d(
                        xLoLimit,       yLoLimit,       zLoLimit,             
                        xHiLimit,       yHiLimit,       zHiLimit,
                        1,              N_CELL_FACES,
                        1,              nMaterials,     
                        7,             
                        uvel_FC,        vvel_FC,        wvel_FC,
                        press_FC,
                        tau_X_FC,       tau_Y_FC,       tau_Z_FC);                         
    stat = putenv("PGPLOT_PLOTTING_ON_OFF=1");
                            
    /*  audit();  */                    /* Steve Parkers memory tool    */   
/*______________________________________________________________________
*  TESTING: HARDWIRE SOME OF THE INPUTS
*       HARDWIRE FOR NOW
*   Comment this out 
*_______________________________________________________________________*/

/*`==========TESTING==========*/ 
/**/    if(xHiLimit > xLoLimit && yHiLimit == yLoLimit)
/**/    {
/**/        qLoLimit = xLoLimit;
/**/        qHiLimit = xHiLimit;
/**/        delQ     = delX;
/**/    }
/**/    if(yHiLimit > yLoLimit && xHiLimit == xLoLimit)
/**/    {
/**/        qLoLimit = yLoLimit;
/**/        qHiLimit = yHiLimit;
/**/        delQ     = delY;
/**/    }
 /*==========TESTING==========`*/
 #if switchOveride_Initial_Conditions                               
  #include "overide_initial_conds.i"
#endif 
/*______________________________________________________________________
*   Plot the inputs (MUST HARDWIRE WHAT YOU WANT TO VIEW)
*   To keep the code clean I moved the code to another file
*_______________________________________________________________________*/
#if switchDebug_main_input
    #define switchInclude_main_1 1
    #include "debugcode.i"
    #undef switchInclude_main_1
#endif     
/*__________________________________
*   For the first time through
*   set some variables
*___________________________________*/
    delt    = delt_limits[3];              
    t       = delt;
    fprintf(stderr,"\nInitial time %f, timestep is %f\n",t,delt);
    
    
    
/*______________________________________________________________________
*   M  A  I  N     A  D  V  A  N  C  E     L  O  O  P 
*_______________________________________________________________________*/                       
    while( t <= t_final)
    {
         should_I_write_output = Is_it_time_to_write_output( t, t_output_vars  );
        /* fprintf(stderr, "should _ I write_output %i\n",should_I_write_output);      */
        

    /*__________________________________
    * update the physical boundary conditions
    * and initialize some arrays
    *___________________________________*/                        
    update_CC_FC_physical_boundary_conditions( 
                        xLoLimit,       yLoLimit,       zLoLimit,             
                        xHiLimit,       yHiLimit,       zHiLimit,             
                        delX,           delY,           delZ,
                        BC_types,       BC_float_or_fixed,
                        BC_Values, 
                        nMaterials,     3,                 
                        uvel_CC,        UVEL,           uvel_FC,
                        vvel_CC,        VVEL,           vvel_FC,
                        wvel_CC,        WVEL,           wvel_FC);
                        
    update_CC_physical_boundary_conditions( 
                        xLoLimit,       yLoLimit,       zLoLimit,             
                        xHiLimit,       yHiLimit,       zHiLimit,             
                        delX,           delY,           delZ,
                        BC_types,       BC_float_or_fixed,
                        BC_Values, 
                        nMaterials,     3,                 
                        Temp_CC,TEMP,   rho_CC,DENSITY, press_CC,PRESS);
                        
    zero_arrays_4d(
                        xLoLimit,       yLoLimit,       zLoLimit,             
                        xHiLimit,       yHiLimit,       zHiLimit,
                        1,              nMaterials,     8,             
                        mass_source,    delPress_CC,    int_eng_source,  
                        xmom_source,    ymom_source,    zmom_source,
                        Vol_L_CC,       mass_CC);
/*`==========TESTING==========*/ 
/*______________________________________________________________________
*   Specific to solving the Riemann problem
*_______________________________________________________________________*/
/**/        for ( i = 1; i <= Q_MAX_LIM; i++)
/**/        {
/**/           u_Rieman[i]     = 0.0;
/**/           a_Rieman[i]     = 0.0;
/**/           p_Rieman[i]     = 0.0;
/**/           rho_Rieman[i]   = 0.0; 
/**/           T_Rieman[i]     = 0.0;          
/**/        }
/**/
/**/        Solve_Riemann_problem(  
/**/                            qLoLimit,       qHiLimit,       
/**/                            delQ,           t,              gamma[m],
/**/                            p1,             rho1,           u1, a1,
/**/                            p4,             rho4,           u4, a4,
/**/                            u_Rieman,       a_Rieman,       p_Rieman,
/**/                            rho_Rieman,     T_Rieman,       R[m]);                      
 /*==========TESTING==========`*/
 /*__________________________________
 *  STEP 1
 *  Use the equation of state to get
 *  P at the cell center
 *___________________________________*/
#if switch_step1_OnOff
        equation_of_state(
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        R,
                        press_CC,       rho_CC,         Temp_CC,
                        cv_CC,          nMaterials   );
                        
        speed_of_sound(
                        xLoLimit,       yLoLimit,       zLoLimit,       
                        xHiLimit,       yHiLimit,       zHiLimit,       
                        gamma,          R,              Temp_CC,     
                        speedSound,     nMaterials   );
    #endif

    /*__________________________________
    *    S  T  E  P     2 
    *   Use Euler's equation thingy to solve
    *   for the n+1 Lagrangian press (CC)
    *   and the n+1 face centered fluxing
    *   velocity
    *___________________________________*/ 
     /*__________________________________
    *   Take (*)vel_CC and interpolate it to the 
    *   face-center.  Advection operator needs
    *   uvel_FC and so does the pressure solver
    *___________________________________*/ 
        stat = putenv("PGPLOT_PLOTTING_ON_OFF=1"); 
        compute_face_centered_velocities( 
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delX,           delY,           delZ,
                        delt,           
                        BC_types,       BC_float_or_fixed,
                        BC_Values,
                        rho_CC,         grav,           press_CC,
                        uvel_CC,        vvel_CC,        wvel_CC,
                        uvel_FC,        vvel_FC,        wvel_FC,
                        nMaterials ); 
                        
        divergence_of_face_centered_velocity(  
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delX,           delY,           delZ,
                        uvel_FC,        vvel_FC,        wvel_FC,
                        div_velFC_CC,   nMaterials); 
        stat = putenv("PGPLOT_PLOTTING_ON_OFF=1");

 
#if switch_step2_OnOff                        

    explicit_delPress
             (  
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delX,           delY,           delZ,
                        div_velFC_CC,
                        delPress_CC,    press_CC,
                        rho_CC,         delt,           speedSound,
                        nMaterials );
                
    update_CC_physical_boundary_conditions( 
                        xLoLimit,       yLoLimit,       zLoLimit,             
                        xHiLimit,       yHiLimit,       zHiLimit,             
                        delX,           delY,           delZ,
                        BC_types,       BC_float_or_fixed,
                        BC_Values, 
                        nMaterials,     1,                 
                        delPress_CC,    DELPRESS);           
#endif 
 
    
   
    /* ______________________________   
    *    S  T  E  P     3    
    *   Compute the face-centered pressure
    *   using the "continuity of acceleration"
    *   principle                     
    * ______________________________   */
    #if switch_step3_OnOff                                  
        press_face(         
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delX,           delY,           delZ,
                        BC_types,       BC_float_or_fixed, BC_Values,
                        press_CC,       press_FC,       rho_CC, 
                        nMaterials );
    #endif
        /*__________________________________
        *
        *___________________________________*/
    update_CC_FC_physical_boundary_conditions( 
                        xLoLimit,       yLoLimit,       zLoLimit,             
                        xHiLimit,       yHiLimit,       zHiLimit,             
                        delX,           delY,           delZ,
                        BC_types,       BC_float_or_fixed,
                        BC_Values, 
                        nMaterials,     1,                 
                        press_CC,       PRESS,          press_FC);


/* ______________________________  
*    S  T  E  P     4                           
*   Compute ssources of mass, momentum and energy
*   For momentum, there are sources
*   due to mass conversion, gravity
*   pressure, divergence of the stress
*   and momentum exchange
* ______________________________   */
#if (switch_step4_OnOff == 1 && switch_Compute_burgers_eq == 0) 
    accumulate_momentum_source_sinks(
                        xLoLimit,       yLoLimit,       zLoLimit,                  
                        xHiLimit,       yHiLimit,       zHiLimit,                  
                        delt,                      
                        delX,           delY,           delZ,                      
                        grav,                  
                        mass_CC,        rho_CC,         press_FC,            
                        Temp_CC,        cv_CC,
                        uvel_CC,        vvel_CC,        wvel_CC,
                        tau_X_FC,       tau_Y_FC,       tau_Z_FC,               
                        viscosity_CC,              
                        xmom_source,    ymom_source,    zmom_source,           
                        nMaterials   ); 

 
   accumulate_energy_source_sinks(
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delt,            
                        delX,           delY,           delZ,    
                        grav,           mass_CC,        rho_CC,          
                        press_CC,       delPress_CC,    Temp_CC,         
                        cv_CC,          speedSound,     
                        uvel_CC,        vvel_CC,        wvel_CC,
                        div_velFC_CC,         
                        int_eng_source,  
                        nMaterials   );

    #endif


    /*__________________________________
    *    S  T  E  P     5                        
    *   Compute Lagrangian values for the volume 
    *   mass, momentum and energy.
    *   Lagrangian values are the sum of the time n
    *   values and the sources computed in 4
    *___________________________________*/
    #if switch_step5_OnOff 
    lagrangian_vol(     xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delX,           delY,           delZ,
                        delt,           
                        Vol_L_CC,       Vol_CC,
                        uvel_FC,        vvel_FC,        wvel_FC,
                        nMaterials);
                        
    calc_flux_or_primitive_vars(    -1,           
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        rho_CC,         Vol_CC,         
                        uvel_CC,        vvel_CC,        wvel_CC,        
                        xmom_CC,        ymom_CC,        zmom_CC,
                        cv_CC,          int_eng_CC,     Temp_CC,
                        nMaterials );                       
                        
    lagrangian_values(  
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        Vol_L_CC,       Vol_CC,         rho_CC,
                        rho_L_CC,
                        xmom_CC,        ymom_CC,        zmom_CC,
                        uvel_CC,        vvel_CC,        wvel_CC,
                        xmom_L_CC,      ymom_L_CC,      zmom_L_CC,
                        mass_L_CC,      mass_source,    
                        xmom_source,    ymom_source,    zmom_source,
                        int_eng_CC,     int_eng_L_CC,   int_eng_source,
                        nMaterials);
    #endif  
                                     
    /*_________________________________   
    *    S  T  E  P     6                            
    *   Compute the advection of mass,
    *   momentum and energy.  These
    *   quantities are advected using the face
    *   centered velocities velocities from 2
    *                  
    *    S  T  E  P     7 
    *   Compute the time advanced values for
    *   mass, momentum and energy.  "Time advanced"
    *   means the sum of the "Lagrangian" values,
    *   found in 5 and the advection contribution
    *   from 6                      
    *______________________________ */  
    #if (switch_step7_OnOff== 1 || switch_step6_OnOff == 1)
     advect_and_advance_in_time(   
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        delX,           delY,           delZ,
                        Vol_CC,         rho_CC,
                        xmom_CC,        ymom_CC,        zmom_CC,
                        Vol_L_CC,       rho_L_CC,       mass_L_CC,
                        xmom_L_CC,      ymom_L_CC,      zmom_L_CC,
                        int_eng_CC,     int_eng_L_CC,
                        uvel_FC,        vvel_FC,        wvel_FC,
                        delt,           nMaterials);

         
    /*__________________________________
    *   Backout the velocities from the 
    *   the momentum
    *___________________________________*/                        
    calc_flux_or_primitive_vars(    1,           
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        rho_CC,         Vol_CC,         
                        uvel_CC,        vvel_CC,        wvel_CC,        
                        xmom_CC,        ymom_CC,        zmom_CC,
                        cv_CC,          int_eng_CC,     Temp_CC,
                        nMaterials ); 
#endif

    /*__________________________________
    *    T  E  C  P  L  O  T  
    *___________________________________*/     
     
    #if tecplot
    if ( should_I_write_output == YES)
    {                     
        tecplot_CC(         
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        x_CC,           y_CC,           z_CC,
                        uvel_CC,        vvel_CC,        wvel_CC,
                        press_CC,       Temp_CC,        rho_CC,
                        scalar1_CC,     scalar2_CC,     scalar3_CC,
                        fileNum,        output_file_basename,       output_file_desc,
                        nMaterials);

        tecplot_FC(         
                        xLoLimit,       yLoLimit,       zLoLimit,
                        xHiLimit,       yHiLimit,       zHiLimit,
                        x_FC,           y_FC,           z_FC,
                        uvel_FC,        vvel_FC,        wvel_FC,
                        fileNum,        output_file_basename,       output_file_desc,
                        nMaterials );
                            
        fileNum ++;
    } 
    #endif 


    /*__________________________________
    *  P  L  O  T  T  I  N  G     S  E  C  T  I  O  N 
    *___________________________________*/
    #define switchDebug_main_custom 1
    #if switchDebug_main_custom
    if ( should_I_write_output == YES)
    {
         #define switchInclude_main_custom 1
         #include "debugcode.i"
         #undef switchInclude_main_custom 
    }
    #endif
         /*__________________________________
         *  Clean up the plotting windows 
         *___________________________________*/
         putenv("PGPLOT_I_AM_HERE=1");              
                                         /* tell the plotting routine that   */
                                         /* you're at the bottom of main     */
         putenv("PGPLOT_OPEN_NEW_WINDOWS=1"); 
         
         
    /*__________________________________
    *    A  D  V  A  N  C  E     I  N     T  I  M  E 
    *___________________________________*/
                        
        find_delta_time_based_on_CC_vel(
                        xLoLimit,        yLoLimit,      zLoLimit,
                        xHiLimit,        yHiLimit,      zHiLimit,
                        &delt,           delt_limits,
                        delX,            delY,          delZ,
                        uvel_CC,         vvel_CC,       wvel_CC,
                        speedSound,      CFL,           nMaterials );
           
        t = t + delt;
        fprintf(stderr,"\nTime is %f, timestep is %f\n",t,delt);
 
 }
/* -----------------------------------------------------------------------  
*   F  R  E  E     T  H  E     M  E  M  O  R  Y                                                     
* -----------------------------------------------------------------------  */
    fprintf(stderr,"Now deallocating memory");
    #include "free_memory.i"
    free_dvector_nr( u_Rieman,     0, Q_MAX_LIM);
    free_dvector_nr( a_Rieman,     0, Q_MAX_LIM);
    free_dvector_nr( p_Rieman,     0, Q_MAX_LIM);
    free_dvector_nr( rho_Rieman,   0, Q_MAX_LIM); 
    free_dvector_nr( T_Rieman,     0, Q_MAX_LIM);    
    
    #if switch_explicit_implicit 
        PetscFinalize();
    #endif    
/*__________________________________
*   Quite fullwarn compiler remarks
*___________________________________*/
    i = i;      j = j;      k = k;    
    residual    = residual;
    xLo = xLo;  xHi = xHi;
    yLo = yLo;  yHi = yHi;
    zLo = zLo;  zHi = zHi;
    u0  = u0;
    QUITE_FULLWARN(stat);                       
    QUITE_FULLWARN(fileNum); 

    return(1);
/*STOP_DOC*/
}
