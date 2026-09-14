//*******************************************************************
//          Основные структуры

class plant_t
{
  
  double   noise[][]     =      {
                                 {1.0e-04, 0, 0, 0},
                                 {0, 1.0e-04, 0, 0},
                                 {0, 0, 1.0e-04, 0},
                                 {0, 0, 0, 1.0e-04}
                                }; 
                                
       
  double h      = 0.05,
         ut,
         up,
         
         K11,
         K12,
         K21,
         K22,
         K31,
         K32,
         K41,
         K42,
         
         //dtheta1:  [rad/s] angular velocity of inner pendulum
         //dtheta2:  [rad/s] angular velocity of outer pendulum
         //theta1:   [rad]   angle of inner pendulum
         //theta2:   [rad]   angle of outer pendulum
         
         x1_i_1,
         x1_i,
         x2_i_1,
         x2_i,
         x3_i_1,
         x3_i,
         x4_i_1,
         x4_i,
         
         m1 = 0.5, // % [kg]     mass of 1st link
         m2 = 0.5, //  % [kg]     mass of 2nd link
         b1 = 0.0, //  % [Ns/m]  coefficient of friction (1st joint)
         b2 = 0.0, //  % [Ns/m]  coefficient of friction (2nd joint)
         l1 = 0.5, //  % [m]      length of 1st pendulum
         l2 = 0.5, //  % [m]      length of 2nd pendulum
         g  = 9.82, // % [m/s^2]  acceleration of gravity
         I1 = m1 * l1 * l1 / 12.0, //  % moment of inertia around pendulum midpoint (inner link)
         I2 = m2 * l2 * l2 / 12.0, //  % moment of inertia around pendulum midpoint (outer link)
         
         a11 = l1 * l1 * (0.25 * m1 + m2) + I1,
         a12 = 0.5 * m2 * l1 * l2,
         a22 = l2 * l2 * 0.25 * m2 + I2,
         c1 = g * l1 * (0.5 * m1 + m2),
         c11 = - 0.5 * m2 * l1 * l2,
         c2 = 0.5 * m2 * l2;
         
  double  A[][]          =      {
                                  {a11, 0},
                                  {0, a22}
                                },
                                
          B[]          =        {0, 0},
          x[]          =        {0, 0};  

  void  calcKoeff()
  {
   
         
  }

  void  calcShema(double u)
  {
        
         //A = [l1^2*(0.25*m1+m2) + I1,      0.5*m2*l1*l2*cos(z(3)-z(4));
         //      0.5*m2*l1*l2*cos(z(3)-z(4)), l2^2*0.25*m2 + I2          ];
        
         //b = [g*l1*sin(z(3))*(0.5*m1+m2) - 0.5*m2*l1*l2*z(2)^2*sin(z(3)-z(4))...
         //                                               + f(t) - b1*z(1);
         //0.5*m2*l2*( l1*z(1)^2*sin(z(3)-z(4)) + g*sin(z(4)) )    - b2*z(2)];
       
        x1_i = x1_i_1;
        x2_i = x2_i_1;
        x3_i = x3_i_1;
        x4_i = x4_i_1;
        
        A[0][1] = a12 * cos_db(x3_i - x4_i) ;
        A[1][0] = A[0][1];
        
        B[0] = c1 * sin_db(x3_i) + c11 * x2_i * x2_i * sin_db(x3_i - x4_i) + u - b1 * x1_i;
        B[1] = c2 * ( l1 * x1_i * x1_i * sin_db(x3_i - x4_i) + g * sin_db(x4_i) ) - b2 * x2_i;
        
        x = operInvSlash(A, B);
        
        K11 = x[0];
        K21 = x[1];  
        
        K31 = x1_i;
        K41 = x2_i;
        
        A[0][1] = a12 * cos_db((x3_i + h * K31) - (x4_i + h * K41)) ;
        A[1][0] = A[0][1];
        
        B[0] = c1 * sin_db(x3_i + h * K31) + c11 * (x2_i + h * K21) * (x2_i + h * K21) * 
                  sin_db((x3_i + h * K31) - (x4_i + h * K41)) + u - b1 * (x1_i + h * K11);
        B[1] = c2 * ( l1 * (x1_i + h * K11) * (x1_i + h * K11) * 
                sin_db((x3_i + h * K31) - (x4_i + h * K41)) + g * sin_db(x4_i + h * K41) ) 
                - b2 * (x2_i + h * K21);
        
        x = operInvSlash(A, B);
        
        K12 = x[0];
        K22 = x[1];  
        
        K32 = x1_i + h * K11;
        K42 = x2_i + h * K21;
                
        x1_i_1 = x1_i + 0.5 * h * (K11 + K12);
        x2_i_1 = x2_i + 0.5 * h * (K21 + K22);
        x3_i_1 = x3_i + 0.5 * h * (K31 + K32);
        x4_i_1 = x4_i + 0.5 * h * (K41 + K42);
        
  }
  
  double[]  dynamics(double stateL[], double uL)
  {
   double ret[] = new double[4];
   double stateW[] = stateL;//vectorCoef(stateL, 1.0);
   x1_i_1 = stateW[0];
   x2_i_1 = stateW[1];
   x3_i_1 = stateW[2];
   x4_i_1 = stateW[3];
   h = 0.025;
   calcShema(uL);
   calcShema(uL);
   calcShema(uL);
   
   ret[0] = x1_i_1;
   ret[1] = x2_i_1;
   ret[2] = x3_i_1;
   ret[3] = x4_i_1;
   
   return ret;
  }
 
  plant_t(double hn) 
  {      
    h = hn;
    x1_i_1 = 0.0;
    x2_i_1 = 0.0;
    x3_i_1 = 3.0;
    x4_i_1 = 0.0;
  }
  
}

class dynmodel_t 
{
  
  double  inputs[][],
          targets[][],
          hyp[][];
  
  double  iK[][][],
          _K[][][],
          beta[][];
      
  int     n,
          mIn,
          mTg,
          dyni[]  = {0, 1, 4, 5, 6, 7};//{0, 1, 4, 5, 6, 7};
          
  dynmodel_t(int ni, int mi, int mt) 
  {  
    
    inputs   = new double[ni][mi];
    targets  = new double[ni][mt];
    n = ni;
    mIn = mi;
    mTg = mt;
    
  }
  
  void prt() 
  {
    
    println();
    println(" Inputs " + n + " x " + mIn);
    printMat(inputs);
    println();
    println(" Targets " + n + " x " + mTg);
    printMat(targets);
    
  }
  
}

class rollout_r 
{
  
  double[][]    xLock;
  double[][]    yLock;
  double[][]    latentLock;
  
  double        L[];
  
  int          i,
               j,
               Hlock;
  
  rollout_r(int count) 
  { 
    
    xLock       = new double[count + 1][9];//new double[count + 1][5];
    yLock       = new double[count][4];//new double[count][2];
    latentLock  = new double[count + 1][5];//new double[count + 1][3];
    Hlock       = count;
    L           = new double[count];
    
  }
  
  void printVal() 
  {
    
    println("---- x data ----");
    for(i = 0; i < Hlock ; i++)
    {
     
      for(j = 0; j < 9 ; j++)
        print(xLock[i][j] + "  ");
      println();
    }
    
    println("---- y data ----");
    for(i = 0; i < Hlock ; i++)
    {
     
      for(j = 0; j < 4 ; j++)
        print(yLock[i][j] + "  ");
      println();
    }
    println("---- latent data ----");
    for(i = 0; i < Hlock ; i++)
    {
     
      for(j = 0; j < 3 ; j++)
        print(latentLock[i][j] + "  ");
      println();
    }
     println("---- L data ----");
    for(i = 0; i < Hlock ; i++)
    {
        print(L[i] + "  ");
     
    }
    println();
  }
  
}

class tp_p
{
  // GP
  double  inputs[][],
          targets[],
          hyp[];
  // Lin       
  double  w[],
          b;
  tp_p()
  {
  };
};

class nargout_valueSh_t 
{
  
  double J;
      double[] L;      // ← этого поля нет? 
  nargout_valueSh_t() 
  {   
  
  }
  
}

class nargout_value_t 
{
  
  double J;
  double[] L; 
  tp_p   dJdp;
        
  nargout_value_t() 
  {   
  
  }
  
}

class policy_t 
{
  
  double  maxU     = 3.5;
          
  int     angle[]  = {2, 3},
          e        = 1;   
  tp_p    p;
  
  int     poli[]   = {0, 1, 4, 5, 6, 7};
  
  policy_t() 
  {  
 
  }
  
  
}

class nargout_SE_A_t // Выходные параметры для функции квадратичной экспоненциальной ковариации, 
                     // функции суммы ковариации и независимой ковариационной функциям. Начало
{
  
  int     n;
  double   A[][];
    
  nargout_SE_A_t(int n) 
  {   
    
    A = new double[n][n];
      
  }
  
} // Выходные параметры для функции квадратичной экспоненциальной ковариации,
  // функции суммы ковариации и независимой ковариационной функциям. Конец

//  Выходные параметры для функции, которая возвращает минус логарифм правдоподобия и его частные производные по гиперпараметрам;
class nargout_gpr_t 
{
  
  int     n;
  double  out1;
  double  out2[];
    
  nargout_gpr_t(int n_1) 
  {   
    
    n = n_1;
    out2 = new double[n];
      
  }
  
  void prt()
  {
    
    println(" gpr result ");
    println(" out1 = " + out1);
    print(" out2 = ");
    for(int i = 0; i < n ; i++)
      print(out2[i] + "  ");
      
    println();
    
  }
  
}

// Параметры для штрафования экстремальных гиперпараметров
class curb_t 
{
  
  int     m,
          i,
          D,
          n;
  double  snr,
          ls;
  double  _std[];
    
  curb_t(double sn, double l, double inp[][]) 
  {   
    D = inp[0].length;
    n = inp.length;
    m = D;
    _std = std(inp);
    snr = sn;
    ls = l;
    
  }
  
  void prt()
  {
    
    println(" curb result ");
    println(" snr =" + snr);
    println(" ls = " + ls);
    println(" std = ");
    
    for(i = 0; i < m ; i++) 
        print(_std[i] + "  ");
      
   println();
    
  }
  
}

// Выходные параметры для функции прогнозирования gp0. Настройка
class nargout_gp0_t 
{
  
  double M[],
         V[][],
         S[][];
    
  nargout_gp0_t() 
  {   
  
  }
  
}

// Выходные параметры для функции прогнозирования gp0. Прогноз с производными
class nargout_gp0d_t 
{
  
  double M[],
         V[][],
         S[][],
         dMdm[][],
         dSdm[][],
         dVdm[][],
         dMds[][],
         dSds[][],
         dVds[][];
    
  nargout_gp0d_t() 
  {   
  
  }
  
}

class nargout_gp2_t 
{
  
  double M[],
         V[][],
         S[][],
         dMdm[][],
         dMds[][],
         dSdm[][],
         dSds[][],
         dVdm[][],
         dVds[][],
         dMdi[][],
         dMdt[][],
         dMdX[][],
         dSdi[][],
         dSdt[][],
         dSdX[][],
         dVdi[][],
         dVdt[][],
         dVdX[][];
    
  nargout_gp2_t() 
  {   
  
  }
  
}

class cost_t 
{
  
  double  target[] = {0, 0, 0, 0},
          gamma    = 1.0f,
          p[]      = {0.5f, 0.5f},
          width_c  = 0.25,
          expl     = 0;
  
  double  z[],
          W[][];
          
  int     angle[]  = {2, 3};
  cost_t() 
  {  
    
  }
  
  void fcn() 
  {
    
        
  }
  
}

class nargout_loss_cp_t 
{
  
  double  L,
          dLdm[],  // 1xsize(m)
          dLds[],  // 1xsize(m^2)
          S2;
  
  nargout_loss_cp_t() 
  {  
    
  }
  
}

class nargout_lossSat_t 
{
  
  double L, 
         dLdm[], 
         dLds[][], 
         S, 
         dSdm[], 
         dSds[][];
  
  nargout_lossSat_t() 
  {  
    
  
  }
  
}

class gTrig_r 
{
  
  double  M[];
  double  V[][];
  double  C[][];

  gTrig_r() 
  {   
  }
    
}

class gTrig_Full_r 
{
  
  double  M[];
  double  V[][];
  double  C[][];
  
  double  dMdm[][],
          dMdv[][],
          dVdm[][],
          dVdv[][],
          dCdv[][],
          dCdm[][];
          
  gTrig_Full_r() 
  {   
  }
    
}

//=================================================================================  
 // Структуры, содержащие наборы входных и выходных параметров функций оптимизации
 
  class p_linesearch_r // Для работы функции линейного поиска. Начало
  {
  
    double  x       = 0,
            f       = 0,
            s       = 0;
    int     D       = 0;
            
    char    fs,
            ss;
            
    double  df[];
    
    p_linesearch_r(double x0, double f0, double df0[], double s0) 
    {   
      
      x = x0;
      f = f0;
      D = df0.length;
      df = new double [D];
      for(int i = 0; i < D; i++)
        df[i] = df0[i];
      s = s0;
      
    }
    
  }  // Для работы функции линейного поиска. Конец
  
   class nargout_linesearch2_r // Набор выходных параметров функции линейного поиска. Начало
  {
  
    double  a       = 0,
            fx      = 0;
    int     i       = 0;
    double  x[],
            df[];
    
    nargout_linesearch2_r() 
    {   
      
    }
    
  }    // Набор выходных параметров функции линейного поиска. Конец
    
  class p_wp_t // Для работы функции wp проверка Wolfe-Powell условий. Начало
  {
  double a,
         b,
         c,
         sig,
         rho;
  p_wp_t()
  {
    
  }
  }; // Для работы функции wp проверка Wolfe-Powell условий. Конец
          
  class pt_struct_t // Набор параметров для работы функции поиска минимума. Начало
  {
    
    String  S = "               #",
            S2 = "function evaluation #";
    int     lgth        = 150,
            MFEPLS      = 30,
            verbosity   = 1,
            fh          = 1;
    double  MSR         = 100,
            SIG         = 0.5,
            H[][];  
    
    pt_struct_t() 
    {   
      
    }
    
  }    // Набор параметров для работы функции поиска минимума. Конец
  
  class nargout_BFGS2_t // Набор выходных параметров функции поиска минимума функции методом BFGS. Начало
  {
  
    double         x[],
                   fx[]; 
    int            i;
    pt_struct_t    p;
  
    nargout_BFGS2_t() 
    {   
     
    }
    
  }    // Набор выходных параметров функции поиска минимума функции методом BFGS. Начало
  
  class nargout_f_t 
  {
     
    double fx,
           dfx[];
   
    nargout_f_t() 
    {   
   
    }
   
  }
  
// class nargout_congp_t 
// {
  
//  //function [M, S, C]
//  double M[],
//         S[][],
//         C[][];
    
//  nargout_congp_t() 
//  {   
  
//  }
  
//}

//class nargout_congpd_t 
//{
  
//  //function [M, S, C, dMdm, dSdm, dCdm, dMds, dSds, dCds, dMdp, dSdp, dCdp]
//  double M[],
//         S[][],
//         C[][],
//         dMdm[][],
//         dSdm[][],
//         dCdm[][],
//         dMds[][],
//         dSds[][],
//         dCds[][],
//         dMdp[],
//         dSdp[],
//         dCdp[][];
    
//  nargout_congpd_t() 
//  {   
  
//  }
  
//}

class nargout_conpols_t 
 {
  
  //function [M, S, C]
  double M[],
         S[][],
         C[][];
    
  nargout_conpols_t() 
  {   
  
  }
  
}

class nargout_conpolsd_t 
{
  
  //function [M, S, C, dMdm, dSdm, dCdm, dMds, dSds, dCds, dMdp, dSdp, dCdp]
  double M[],
         S[][],
         C[][],
         dMdm[][],
         dSdm[][],
         dCdm[][],
         dMds[][],
         dSds[][],
         dCds[][],
         dMdp[],
         dSdp[],
         dCdp[][];
    
  nargout_conpolsd_t() 
  {   
  
  }
  
}

class nargout_gSat_t 
{
  
  //function [M, S, C, dMdm, dSdm, dCdm, dMdv, dSdv, dCdv] = gSat(m, v, i, e)
  
  double M,
         S,
         C[],
         dMdm[],
         dSdm[],
         dCdm[][],
         dMdv[],
         dSdv[],
         dCdv[][];
    
  nargout_gSat_t() 
  {   
  
  }
  
}

class nargout_gSin_t 
{
  
  //function [M, V, C, dMdm, dVdm, dCdm, dMdv, dVdv, dCdv] = gSin(m, v, i, e)
  
  double M[],
         V[][],
         C[][],
         dMdm[][],
         dVdm[][],
         dCdm[][],
         dMdv[][],
         dVdv[][],
         dCdv[][];
    
  nargout_gSin_t() 
  {   
  
  }
  
}

class nargout_propagate_t 
{
  
  double Mnext[],
        Snext[][];
    
  nargout_propagate_t() 
  {   
  
  }
  
}

class nargout_propagated_t 
{
  
  double Mnext[],
         Snext[][],
         dMdm[][],
         dSdm[][],
         dMds[][],
         dSds[][],
         dMdp[][],
         dSdp[][];
        
  nargout_propagated_t() 
  {   
  
  }
  
}

class nargout_fillInFl_t 
{
  
  double S[][],
         Mdm[][],
         Mds[][],
         Sdm[][],
         Sds[][],
         Mdp[][],
         Sdp[][];
        
  nargout_fillInFl_t() 
  {   
  
  }
  
}

class nargin_fillInFl_t 
{
  
  double S[][],
         C[][],
         mdm[][],
         sdm[][],
         Cdm[][],
         mds[][],
         sds[][],
         Cds[][],
         Mdm[][],
         Sdm[][],
         Mds[][],
         Sds[][],
         Mdp[][],
         Sdp[][],
         Cdp[][];
  
  int    i[],
         j[],
         k[],
         D;
        
  nargin_fillInFl_t() 
  {   
  
  }
  
}

class nargout_fillInSh_t 
{
  
  double S[][],
         Mdm[][],
         Mds[][],
         Sdm[][],
         Sds[][];
        
  nargout_fillInSh_t() 
  {   
  
  }
  
}

class nargin_fillInSh_t 
{
  
  double S[][],
         C[][],
         mdm[][],
         sdm[][],
         Cdm[][],
         mds[][],
         sds[][],
         Cds[][],
         Mdm[][],
         Sdm[][],
         Mds[][],
         Sds[][];
  
  int    i[],
         j[],
         k[],
         D;
        
  nargin_fillInSh_t() 
  {   
  
  }
  
}

class nargout_fillInMd_t 
{
  
  double S[][],
         Mdm[][],
         Mds[][],
         Sdm[][],
         Sds[][],
         Mdp[][],
         Sdp[][];
        
  nargout_fillInMd_t() 
  {   
  
  }
  
}

class nargin_fillInMd_t 
{
  
  double S[][],
         C[][],
         mdm[][],
         sdm[][],
         Cdm[][],
         mds[][],
         sds[][],
         Cds[][],
         Mdm[][],
         Sdm[][],
         Mds[][],
         Sds[][],
         Mdp[][],
         Sdp[][];
  
  int    i[],
         j[],
         k[],
         D;
        
  nargin_fillInMd_t() 
  {   
  
  }
  
}

class nargout_minimize_gp_t // Набор выходных параметров функции минимизации. Начало
  {
  
    double hyp[];
    double v[];
                   
    nargout_minimize_gp_t() 
    {   
                 
    }
        
  }  // Набор выходных параметров функции минимизации. Конец
  
class nargout_minimize_pol_t // Набор выходных параметров функции минимизации. Начало
  {
  
    tp_p pol;
    double fX[];
                   
    nargout_minimize_pol_t() 
    {   
          
    }
        
  }  // Набор выходных параметров функции минимизации. Конец
  
 class sparse_t
 {
   int m,
       n,
       elements,
       i[],
       j[];
       
   double val[];
   
   sparse_t(int n1, int m1, int elms)
   {
     m = m1;
     n = n1;
     elements = elms;
     i = new int [elements];
     j = new int [elements];
     val = new double [elements];
   }
   
   sparse_t(double inp[][])
   {
     m = inp[0].length;
     n = inp.length;
     int k = 0;
     i = new int [1];
     j = new int [1];
     for(int i1 = 0; i1 < n; i1++)
       for(int j1 = 0 ; j1 < m; j1++)
         if(inp[i1][j1] != 0)
         {
           i[k] = i1;
           j[k] = j1;
           i = append(i, 1);
           j = append(j, 1);
           k++;
         }
       elements = k;
       val = new double [elements];

       for(int i1 = 0; i1 < elements; i1++)
       {
           val[i1] = inp[i[i1]][j[i1]];
       }     
   }
   
   void print()
   {
     println(" Прореженная матрица размером " + n + " на " + m + " элементов");
     for(int i1 = 0; i1 < elements; i1++)
     {
       println(" (" + i[i1] + ", " + j[i1] + ")        " + val[i1]);
     }
     
   }
   
 }
  
