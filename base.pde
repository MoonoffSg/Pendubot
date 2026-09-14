
//%  1  dtheta1        angular velocity of inner pendulum
//%  2  dtheta2        angular velocity of outer pendulum
//%  3  theta1         angle inner pendulum
//%  4  theta2         angle outer pendulum
//%  5  sin(theta1)    complex representation ...
//%  6  cos(theta1)    ... of angle of inner pendulum
//%  7  sin(theta2)    complex representation ...
//%  8  cos(theta2)    ... of angle of outer pendulum
//%  9  u              torque applied to the inner joint

rollout_r rollout(double start[], int H, plant_t plantL, cost_t cost)
{

  rollout_r    ret            = new rollout_r(H);
  gTrig_r      getgTrig       = new gTrig_r();
  
  double[][]   xLock          = new double[1][8]; // dtheta1 dtheta2 theta1 theta12 sin(theta1) cos(theta1) sin(theta2) cos(theta2)
  double[][]   latentLock     = new double[1][5]; // dtheta1 dtheta2 theta1 theta12 u
  double[]     uLock          = new double[1];    // u
  
  double[]     randn          = new double[4];
  double[][]   chol_noise     = new double[4][4];
  double       zeros[][]      = {
                                  {0, 0, 0, 0},
                                  {0, 0, 0, 0},
                                  {0, 0, 0, 0},
                                  {0, 0, 0, 0}
                                 };
  int          indices[]      = {2, 3};
  double       e              = 1;
  
  double[]     stateLock      = vectorCoef(start, 1.0);// В исходном коде State(simi) = start в массив записываем начальное состояние
  double[]     nextLock       = new double[4];
  double[]     L              = new double[H];
  int          i              = 0,
               j              = 0;
   
    //  Вычисляем randn(size(simi))*chol(plant.noise)
    randn[0]  = (double)randomGaussian() ;
    randn[1]  = (double)randomGaussian() ;
    randn[2]  = (double)randomGaussian() ;
    randn[3]  = (double)randomGaussian() ;
    
    Chol_S(plant.noise, 4); 
    chol_noise = getSch();
    randn = matrixMultiplyC_T2(randn, chol_noise);
    
  // Вычисляем x(1,simi) = start' + randn(size(simi))*chol(plant.noise);
    xLock[0][0]  = start[0] + randn[0];
    xLock[0][1]  = start[1] + randn[1];
    xLock[0][2]  = start[2] + randn[2];
    xLock[0][3]  = start[3] + randn[3];
    
    for(i = 0; i < H ; i++)
    {
     
      // Вычисляем тригонометрию угла для текущего состояния из массива x и дополняем массив x
      getgTrig = gTrig(xLock[i], zeros, indices, e);
      xLock[i][4] = getgTrig.M[0];
      xLock[i][5] = getgTrig.M[1];
      xLock[i][6] = getgTrig.M[2];
      xLock[i][7] = getgTrig.M[3];
      
      // Apply policy ... or random actions
      
      uLock[i] = policy.maxU * (double)random(-1.0, 1.0);
             
      uLock = append_db(uLock, 0); 
          
      // Заполняем не зашумленными данными массив latent и добавляем управляющий сигнал
      latentLock[i][0] = stateLock[0];
      latentLock[i][1] = stateLock[1];
      latentLock[i][2] = stateLock[2];
      latentLock[i][3] = stateLock[3];
      latentLock[i][4] = uLock[i];
      latentLock = append_db(latentLock); 
      // Вычисляем следующее состояние системы при переходе из текущего не зашумленного
      nextLock = simulate(stateLock, uLock[i], plantL);
            
      stateLock[0] = nextLock[0];
      stateLock[1] = nextLock[1];
      stateLock[2] = nextLock[2];
      stateLock[3] = nextLock[3];
      
      // Вычисляем randn(size(simi))*chol(plant.noise)
      randn[0]  = (double)randomGaussian() ;
      randn[1]  = (double)randomGaussian() ;
      randn[2]  = (double)randomGaussian() ;
      randn[3]  = (double)randomGaussian() ;
    
      randn = matrixMultiplyC_T2(randn, chol_noise);
      // Вычисляем x(i+1,simi) = state(simi) + randn(size(simi))*chol(plant.noise);
      xLock = append_db(xLock); 
      xLock[i + 1][0]  = stateLock[0] + randn[0];
      xLock[i + 1][1]  = stateLock[1] + randn[1];
      xLock[i + 1][2]  = stateLock[2] + randn[2];
      xLock[i + 1][3]  = stateLock[3] + randn[3];
      
      // Compute Cost 
      double zerosS[][] = new double[stateLock.length][stateLock.length];
      nargout_loss_cp_t r = loss_pendubot(cost, stateLock, zerosS);
      
      L[i] = r.L;
      
    }
    
    ret.yLock = new double[H][4];
    for(int k = 1; k < (H + 1) ; k++)
    {
      ret.yLock[k - 1][0] = xLock[k][0];
      ret.yLock[k - 1][1] = xLock[k][1];
      ret.yLock[k - 1][2] = xLock[k][2];
      ret.yLock[k - 1][3] = xLock[k][3];
    }
    
    ret.xLock = new double[H][9];
    for(int k = 0; k < H ; k++)
    {
      ret.xLock[k][0] = xLock[k][0];
      ret.xLock[k][1] = xLock[k][1];
      ret.xLock[k][2] = xLock[k][2];
      ret.xLock[k][3] = xLock[k][3];
      ret.xLock[k][4] = xLock[k][4];
      ret.xLock[k][5] = xLock[k][5];
      ret.xLock[k][6] = xLock[k][6];
      ret.xLock[k][7] = xLock[k][7];
      ret.xLock[k][8] = uLock[k];
    }
    
    latentLock[H][0] = stateLock[0];
    latentLock[H][1] = stateLock[1];
    latentLock[H][2] = stateLock[2];
    latentLock[H][3] = stateLock[3];
   
    ret.latentLock = latentLock;//coeffProdMat(1.0, latentLock);
    ret.L = L;//vectorCoef(L, 1.0);
    
  return ret;
  
}

rollout_r rollout(double start[], policy_t _policy, int H, plant_t plantL, cost_t cost)
{
  
  rollout_r    ret            = new rollout_r(H);
  gTrig_r      getgTrig       = new gTrig_r();
  
  double[][]   xLock          = new double[1][8]; // dtheta1 dtheta2 theta1 theta12 sin(theta1) cos(theta1) sin(theta2) cos(theta2)
  double[][]   latentLock     = new double[1][5]; // dtheta1 dtheta2 theta1 theta12 u
  double[]     uLock          = new double[1];    // u
  
  double[]     randn          = new double[4];
  double[][]   chol_noise     = new double[4][4];
  double       zeros[][]      = {
                                  {0, 0, 0, 0},
                                  {0, 0, 0, 0},
                                  {0, 0, 0, 0},
                                  {0, 0, 0, 0}
                                 };
  int          indices[]      = {2, 3};
  double       e              =  1;
  
  double[]     stateLock      = vectorCoef(start, 1.0);// В исходном коде State(simi) = start в массив записываем начальное состояние
  double[]     nextLock       = new double[4];
  double[]     L              = new double[H];
  int          i              = 0,
               j              = 0;
   
    //  Вычисляем randn(size(simi))*chol(plant.noise)
    randn[0]  = (double)randomGaussian() ;
    randn[1]  = (double)randomGaussian() ;
    randn[2]  = (double)randomGaussian() ;
    randn[3]  = (double)randomGaussian() ;
    
    Chol_S(plant.noise, 4); 
    chol_noise = getSch();
    randn = matrixMultiplyC_T2(randn, chol_noise);
    
  // Вычисляем x(1,simi) = start' + randn(size(simi))*chol(plant.noise);
    xLock[0][0]  = start[0] + randn[0];
    xLock[0][1]  = start[1] + randn[1];
    xLock[0][2]  = start[2] + randn[2];
    xLock[0][3]  = start[3] + randn[3];
          
    for(i = 0; i < H ; i++)
    {
     
      // Вычисляем тригонометрию угла для текущего состояния из массива x и дополняем массив x
      getgTrig = gTrig(xLock[i], zeros, indices, e);
      xLock[i][4] = getgTrig.M[0];
      xLock[i][5] = getgTrig.M[1];
      xLock[i][6] = getgTrig.M[2];
      xLock[i][7] = getgTrig.M[3];
    
      // Apply policy ... or random actions
      //u(i,:) = policy.fcn(policy,s(poli),zeros(length(poli)));
      double zer[][] = new double[_policy.poli.length][_policy.poli.length];
      nargout_conpols_t con = conCat(_policy, getVec(xLock[i], _policy.poli), zer);
      uLock[i] = con.M[0];
      uLock = append_db(uLock, 0); 
      
     // Заполняем не зашумленными данными массив latent и добавляем управляющий сигнал
      latentLock[i][0] = stateLock[0];
      latentLock[i][1] = stateLock[1];
      latentLock[i][2] = stateLock[2];
      latentLock[i][3] = stateLock[3];
      latentLock[i][4] = uLock[i];
      latentLock = append_db(latentLock); 
      
      // Вычисляем следующее состояние системы при переходе из текущего не зашумленного
      nextLock = simulate(stateLock, uLock[i], plantL);
            
      stateLock[0] = nextLock[0];
      stateLock[1] = nextLock[1];
      stateLock[2] = nextLock[2];
      stateLock[3] = nextLock[3];
      
      // Вычисляем randn(size(simi))*chol(plant.noise)
      randn[0]  = (double)randomGaussian() ;
      randn[1]  = (double)randomGaussian() ;
      randn[2]  = (double)randomGaussian() ;
      randn[3]  = (double)randomGaussian() ;
      
      randn = matrixMultiplyC_T2(randn, chol_noise);
      // Вычисляем x(i+1,simi) = state(simi) + randn(size(simi))*chol(plant.noise);
      xLock = append_db(xLock); 
      xLock[i + 1][0]  = stateLock[0] + randn[0];
      xLock[i + 1][1]  = stateLock[1] + randn[1];
      xLock[i + 1][2]  = stateLock[2] + randn[2];
      xLock[i + 1][3]  = stateLock[3] + randn[3];
      
      // Compute Cost 
      double zerosS[][] = new double[stateLock.length][stateLock.length];
      nargout_loss_cp_t r = loss_pendubot(cost, stateLock, zerosS);
      
      L[i] = r.L;
      
    }
    
    ret.yLock = new double[H][4];
    for(int k = 1; k < (H + 1) ; k++)
    {
      ret.yLock[k - 1][0] = xLock[k][0];
      ret.yLock[k - 1][1] = xLock[k][1];
      ret.yLock[k - 1][2] = xLock[k][2];
      ret.yLock[k - 1][3] = xLock[k][3];
    }
    
    ret.xLock = new double[H][9];
    for(int k = 0; k < H ; k++)
    {
    ret.xLock[k][0] = xLock[k][0];
      ret.xLock[k][1] = xLock[k][1];
      ret.xLock[k][2] = xLock[k][2];
      ret.xLock[k][3] = xLock[k][3];
      ret.xLock[k][4] = xLock[k][4];
      ret.xLock[k][5] = xLock[k][5];
      ret.xLock[k][6] = xLock[k][6];
      ret.xLock[k][7] = xLock[k][7];
      ret.xLock[k][8] = uLock[k];
    }
    
    latentLock[H][0] = stateLock[0];
    latentLock[H][1] = stateLock[1];
    latentLock[H][2] = stateLock[2];
    latentLock[H][3] = stateLock[3];
   
   
    ret.latentLock = latentLock;//coeffProdMat(1.0, latentLock);
    ret.L = L;//vectorCoef(L, 1.0);
    
  return ret;
  
}

double[] simulate(double stateL[], double uL, plant_t plantL)
{
  
  double[]   ret = new double[4];
  
  ret = plantL.dynamics(stateL, uL);
  
  return ret;
  
}

void trainDynModel(double xL[][], double yL[][])
{
  
 // В input не записваем велиину угла, а только косинус и синус3
 // угловая скорость / косинус угла / синус угла / управляющий сигнал
 // определяется массивом dyni
 // также используется в функции propagate()

 dynmodel = new dynmodel_t(xL.length, xL[0].length - 2, yL[0].length);
 
 for(int i = 0; i < xL.length; i++)
 {
   
   dynmodel.inputs[i][0] = xL[i][0];
   dynmodel.inputs[i][1] = xL[i][1];
   dynmodel.inputs[i][2] = xL[i][4];
   dynmodel.inputs[i][3] = xL[i][5];
   dynmodel.inputs[i][4] = xL[i][6];
   dynmodel.inputs[i][5] = xL[i][7];
   dynmodel.inputs[i][6] = xL[i][8];
   
 }

 for(int i = 0; i < yL.length; i++)
 {
   
   dynmodel.targets[i][0] = yL[i][0] - xL[i][0];
   dynmodel.targets[i][1] = yL[i][1] - xL[i][1];
   dynmodel.targets[i][2] = yL[i][2] - xL[i][2];
   dynmodel.targets[i][3] = yL[i][3] - xL[i][3];
  
 }
 
  dynmodel = train(dynmodel, trainOpt);

  printMat(dynmodel.hyp);
  
  print("Learned noise std: ");
  for(int i1 = 0; i1 < dynmodel.hyp[0].length ; i1++)
  {
    print("  " + (float)exp_db(dynmodel.hyp[dynmodel.hyp.length - 1][i1]));
  }
  
  println();
  print("SNRs             : ");
  for(int i1 = 0; i1 < dynmodel.hyp[0].length ; i1++)
  {
    print("  " + (float)exp_db(dynmodel.hyp[dynmodel.hyp.length - 2][i1] - dynmodel.hyp[dynmodel.hyp.length - 1][i1]));
  }
  println();
 
}

dynmodel_t train(dynmodel_t gpmodel, int iter)
{
  
  int    D     = gpmodel.inputs[0].length,
         E     = gpmodel.targets[0].length;
  curb_t cr    = new curb_t(1000, 100, gpmodel.inputs);
  
  gpmodel.hyp  = new double[D + 2][E];
  nlml         = new double[E];      
  double lh[][] = new double[D + 2][E];
  double lv[]   = new double[D + 2];
  double lvD[]  = new double[D];
  double lvE[]  = new double[E];
  double gp_trg[]   = new double[gpmodel.n];
  // При первом запуске иниализируем матрицв гиперпараметров
  
  lvD = std(gpmodel.inputs);
  
  for(int i = 0; i < D ; i++)
  for(int j = 0; j < E ; j++)
    lh[i][j] = log_db(lvD[i]);
  
  lvE = std(gpmodel.targets);
  for(int i = 0; i < E ; i++)
  {
    lh[D][i]     = log_db(lvE[i]);
    lh[D + 1][i] = log_db(lvE[i] / 10.0);
  }  
  
  // Для повторного не добавил пока
  
  println("Train hyper-parameters of full GP ...\n");
  
  for(int i = 0; i < E; i++)
  {
    println("GP " + (i + 1) + "/" + E);
    
    for(int j = 0; j < gpmodel.n ; j++)
      gp_trg[j] = gpmodel.targets[j][i];
 
    for(int j = 0; j < (D + 2); j++)
      lv[j] = lh[j][i];

    nargout_minimize_gp_t mOUT = minimize_gp(lv, gpmodel.inputs, gp_trg, cr);
       
    for(int j = 0; j < (D + 2) ; j++)
      gpmodel.hyp[j][i] = mOUT.hyp[j];
      
  }
  return gpmodel;
  
}

//function [J, dJdp] = value(p, m0, S0, dynmodel, policy, plant, cost, H)


nargout_valueSh_t valueSh(tp_p p1, double m0[], double S0[][], dynmodel_t _gpmodel, policy_t _policy, cost_t cost, int H)
{
 println("valueSh: m0 length = " + m0.length);
  println("valueSh: S0 size = " + S0.length + " x " + S0[0].length);
  
  double m[] = m0.clone();
  double S[][] = S0;
  double L[] = new double[H];
  
  for(int t = 0; t < H; t++)
  {
    nargout_propagate_t pd = propagate(m, S, _gpmodel, _policy);
    
    if(t == 0) {
      println("После propagate t=0:");
      println("  Mnext length = " + pd.Mnext.length);
      println("  Mnext = "); printVec(pd.Mnext);
      println("  Snext size = " + pd.Snext.length);
    }
    
    m = pd.Mnext;
    S = pd.Snext;
    
    nargout_loss_cp_t loss = loss_pendubot(cost, m, S);
    
    if(t == 0) {
      println("loss t=0: L = " + loss.L);
    }
    
    L[t] = pow_db(cost.gamma, t) * loss.L;
  }
  
  nargout_valueSh_t ret = new nargout_valueSh_t();
  ret.J = SUM(L);
  ret.L = L;
  return ret;
 //nargout_valueSh_t ret = new nargout_valueSh_t();
 
 ////policy.p = p;             overwrite policy.p with new parameters from minimize
 ////p = unwrap(policy.p) 
 
 // _policy.p = let_p(p1);
 // double[] p = unwrap(_policy.p);
 
 ////dp = 0*p
 ////m = m0; S = S0; L = zeros(1,H)
 
 //double m[] = m0; //vectorCoef(m0, 1.0);
 //double S[][] = S0;//coeffProdMat(1.0, S0);
 //double L[] = new double[H];
 
 ////[m, S, dMdm] = plant.prop(m, S, plant, dynmodel, policy) 
 //nargout_propagate_t  pd2 = new nargout_propagate_t();
 
 //nargout_loss_cp_t loss = new nargout_loss_cp_t();
 //////************************************
 
 ////for t = 1:H 
 //for(int t = 0; t < H ; t++)
 //{
   
 //  //[m, S] = plant.prop(m, S, plant, dynmodel, policy)      get next state
 //  pd2 = propagate(m, S, _gpmodel, _policy);
 //  m =pd2.Mnext; //vectorCoef(pd2.Mnext, 1.0);
 //  S =pd2.Snext; //coeffProdMat(1.0, pd2.Snext);
  
 //  // L(t) = cost.gamma^t.*cost.fcn(cost, m, S);     % expected discounted cost
 //  loss = loss_pendubot(cost, m, S);
 //  L[t] = pow_db(cost.gamma, t) * loss.L;
 //}

 // printMat(S);
 // ret.J = SUM(L);
 // ret.L = L;

 //return ret;
 
}

nargout_value_t valueSh2(tp_p p1, double m0[], double S0[][], dynmodel_t _gpmodel, policy_t _policy, cost_t cost, int H)
{
 
 nargout_value_t ret = new nargout_value_t();
 
 //policy.p = p;             overwrite policy.p with new parameters from minimize
 //p = unwrap(policy.p) 
 
  _policy.p = let_p(p1);
  double[] p = unwrap(_policy.p);
 
 //dp = 0*p
 //m = m0; S = S0; L = zeros(1,H)
 
 double m[] = m0;//vectorCoef(m0, 1.0);
 double S[][] = S0;//coeffProdMat(1.0, S0);
 double L[] = new double[H];
 
 //[m, S, dMdm] = plant.prop(m, S, plant, dynmodel, policy) 
 nargout_propagated_t  pd = new nargout_propagated_t();
 
 nargout_loss_cp_t loss = new nargout_loss_cp_t();
 ////************************************

 //for t = 1:H 
 for(int t = 0; t < H ; t++)
 {
   
   //[m, S] = plant.prop(m, S, plant, dynmodel, policy)      get next state
   pd = propagated(m, S, _gpmodel, _policy);
   m = pd.Mnext;//vectorCoef(pd.Mnext, 1.0);
   S = pd.Snext; //coeffProdMat(1.0, pd.Snext);
  
   // L(t) = cost.gamma^t.*cost.fcn(cost, m, S);     % expected discounted cost
   loss = loss_pendubot(cost, m, S);
   L[t] = pow_db(cost.gamma, t) * loss.L;
 }

  ret.J = SUM(L);
  ret.L = L;
 return ret;
 
}

nargout_value_t value(tp_p p1, double m0[], double S0[][], dynmodel_t _gpmodel, policy_t _policy, cost_t cost, int H)
{
 
 nargout_value_t ret = new nargout_value_t();
 
 //policy.p = p;             overwrite policy.p with new parameters from minimize
 //p = unwrap(policy.p) 
 
  _policy.p = let_p(p1);
  double[] p = unwrap(_policy.p);
  
 //dp = 0*p
 //m = m0; S = S0; L = zeros(1,H)
 
 double dp[] = new double[p.length];
 double m[] = m0;//vectorCoef(m0, 1.0);
 double S[][] = S0;//coeffProdMat(1.0, S0);
 double L[] = new double[H];
 
 double dmdmO[][], 
        dSdmO[][], 
        dmdSO[][], 
        dSdSO[][], 
        dmdp[][], 
        dSdp[][];
 
 // dmOdp = zeros([size(m0,1), length(p)]);
 // dSOdp = zeros([size(m0,1)*size(m0,1), length(p)]);
 double dmOdp[][] = new double[m0.length][p.length];
 double dSOdp[][] = new double[m0.length * m0.length][p.length];
 
 nargout_propagated_t pd = new nargout_propagated_t();
 nargout_loss_cp_t loss  = new nargout_loss_cp_t();
 
 for(int t = 0; t < H ; t++)
 {
  // [m, S, dmdmO, dSdmO, dmdSO, dSdSO, dmdp, dSdp] = ...
  //    plant.prop(m, S, plant, dynmodel, policy); % get next state
  
   pd = propagated(m, S, _gpmodel, _policy);
   m = pd.Mnext;//vectorCoef(pd.Mnext, 1.0);
   S = pd.Snext;//coeffProdMat(1.0, pd.Snext);
   
   //println("Итерация №", t);
   //printVec(m);
   //printMat(S);
   
   dmdmO =pd.dMdm ;//coeffProdMat(1.0, pd.dMdm);
   dSdmO = pd.dSdm;//coeffProdMat(1.0, pd.dSdm);
   dmdSO = pd.dMds;//coeffProdMat(1.0, pd.dMds);
   dSdSO =pd.dSds; //coeffProdMat(1.0, pd.dSds);
   dmdp =pd.dMdp; //coeffProdMat(1.0, pd.dMdp);
   dSdp =pd.dSdp ;//coeffProdMat(1.0, pd.dSdp);
   
   //dmdp = dmdmO*dmOdp + dmdSO*dSOdp + dmdp;
   //dSdp = dSdmO*dmOdp + dSdSO*dSOdp + dSdp;
   dmdp = matrixADD(dmdp, matrixMultiply(dmdmO, dmOdp));
   dmdp = matrixADD(dmdp, matrixMultiply(dmdSO, dSOdp));
      
   dSdp = matrixADD(dSdp, matrixMultiply(dSdmO, dmOdp));
   dSdp = matrixADD(dSdp, matrixMultiply(dSdSO, dSOdp));
    
   // [L(t), dLdm, dLdS] = cost.fcn(cost, m, S);              % predictive cost
   // L(t) = cost.gamma^t*L(t)
   loss = loss_pendubot(cost, m, S);
   double koef = pow_db(cost.gamma, t);
   L[t] =  koef * loss.L;
  
   double dLdm[] = loss.dLdm;//vectorCoef(loss.dLdm, 1.0);
   double dLdS[] = loss.dLds;//vectorCoef(loss.dLds, 1.0);

   // dp = dp + cost.gamma^t*( dLdm(:)'*dmdp + dLdS(:)'*dSdp )'
   
   dp = vectorADD(dp, 
       vectorCoef(vectorADD(matrixMultiplyC_T2(dLdm, dmdp), matrixMultiplyC_T2(dLdS, dSdp)), koef));
  
   // dmOdp = dmdp; dSOdp = dSdp
   dmOdp =dmdp; //coeffProdMat(1.0, dmdp);
   dSOdp = dSdp;//coeffProdMat(1.0, dSdp);
   
 }
 
 //J = sum(L)
 ret.J = SUM(L);

 //dJdp = rewrap(policy.p, dp)
 ret.dJdp = rewrap(_policy.p, dp);

 return ret; 
 
}
