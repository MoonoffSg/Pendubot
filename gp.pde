 //<>// //<>// //<>// //<>//
//=================================================================================  
// Структуры, содержащие наборы входных и выходных параметров функций для описания GP
 
double  K_SE[][];
int     K_SE_m   = 0,
        K_SE_n   = 0;

double  nlml[];


//=====================================================================================================================
// Функция квадратичной экспоненциальной ковариации с Автоматическим определением релевантности (ARD) меры расстояния. 
// Функция ковариации параметризуется как:
// k(x^p,x^q) = sf2 * exp(-(x^p - x^q)'*inv(P)*(x^p - x^q)/2)
// где матрица P диагональна с параметрами ARD ell_1^2,...,ell_D^2, 
// где D — размерность входного пространства, а sf2 — дисперсия сигнала. 
// Гиперпараметры:
// loghyper = [ log(ell_1)
//              log(ell_2)
//               .
//              log(ell_D)
//              log(sqrt(sf2)) ]
              
nargout_SE_A_t covSEardA(double loghypers[], double x[][])
{
  
  int            n   = x.length,
                 D   = x[0].length;
  nargout_SE_A_t ret = new nargout_SE_A_t(n);
  double [][]    ell = new double[D][D];
  double [][]    tem = new double[D][n];
  double         sf2 = 0;
  
  K_SE_n = n;
  K_SE_m = n;
  K_SE = new double[K_SE_n][K_SE_m];
  
  for(int i = 0; i < D; i++)
    ell[i][i] = 1.0 / exp_db(loghypers[i]);  
    
  sf2 = exp_db(2.0 * loghypers[D]);   
  
  tem = matrixMultiply(ell, transMat(x)); // D * n
  
  ret.A = sq_dist(D, tem, n, tem, n);
  for(int i = 0; i < n ; i++)
    for(int j = 0; j < n ; j++)
    {
      ret.A[i][j] = sf2 * exp_db(-ret.A[i][j] * 0.5);
      K_SE[i][j] = ret.A[i][j];
    }
     
  return ret;
  
}

nargout_SE_A_t covSEardA(double loghyper[], double x[][], int z)
{
  
  int            n   = x.length,
                 D   = x[0].length;
  nargout_SE_A_t ret = new nargout_SE_A_t(n);
  double [][]    ell = new double[D][D];
  double [][]    tem = new double[D][n];
  double         sf2 = 0;
  
  for(int i = 0; i < D; i++)
    ell[i][i] = 1.0 / exp_db(loghyper[i]);  
    
  sf2 = exp_db(2.0 * loghyper[D]);   
  
  tem = matrixMultiply(ell, transMat(x)); // D * n
  
  ret.A = sq_dist(D, tem, n, tem, n);
  
  boolean any_size = (K_SE_n == n)&&(K_SE_m == n);
  
  if(!any_size)
  {
    
    K_SE_n = n;
    K_SE_m = n;
    
    K_SE = new double[K_SE_n][K_SE_m];
    
    for(int i = 0; i < n ; i++)
      for(int j = 0; j < n ; j++)
        K_SE[i][j] = sf2 * exp_db(-ret.A[i][j] * 0.5);
        
  }
  
  if(z < D)
  {
    //A = K.*sq_dist(x(:,z)'/ell(z));
    double [] vectr = new double[n];
    
    for(int i = 0; i < n ; i++)
      vectr[i] = x[i][z] / exp_db(loghyper[z]);
      
    ret.A = sq_dist(vectr, n);
    for(int i = 0; i < n ; i++)
      for(int j = 0; j < n ; j++)
       ret.A[i][j] = K_SE[i][j] * ret.A[i][j];
    
  }
  else
  {
    
    for(int i = 0; i < n ; i++)
      for(int j = 0; j < n ; j++)
      {
         ret.A[i][j] = 2.0 * K_SE[i][j];
         K_SE[i][j] = 0;
      }
    
  }
  
  return ret;
  
}

// Формирует ковариационную функцию как сумму других ковариационных функций. 
// Функция ыполняет некоторую обработку и вызывает другие ковариационные функции 
// для выполнения фактической работы.

nargout_SE_A_t covSum(char covfunc, double loghyper[], double x[][])
{
  
  int   D = loghyper.length,
        n = x.length;
  
  double loghyp[]; 
  nargout_SE_A_t ret = new nargout_SE_A_t(n);
  
  loghyp = new double[D - 1];
  for(int i = 0; i < (D - 1) ; i++) loghyp[i] = loghyper[i];
  ret = covSEardA(loghyp, x);
    
  if(covfunc == 'D') // Вызов covSEard + covNoise
  {
    
   ret.A = matrixADD(ret.A, covNoise(loghyper[D - 1], n).A);
    
  }
     
  return ret;
  
}

nargout_SE_A_t covSum(char covfunc, double loghyper[], double x[][], int z)
{
  
  int   D = loghyper.length,
        n = x.length;
        
  double loghyp[]; 
  nargout_SE_A_t ret = new nargout_SE_A_t(n);
  
  loghyp = new double[D - 1];
  for(int i = 0; i < (D - 1); i++) loghyp[i] = loghyper[i];
  ret = covSEardA(loghyp, x, z);
    
  if(covfunc == 'D') // Вызов covSEard + covNoise
  {
    
   ret.A = matrixADD(ret.A, covNoise(loghyper[D - 1], n).A);
    
  }
     
  return ret;
  
}

// Независимая ковариационная функция, т. е. «белый шум», с заданной дисперсией.
// Ковариационная функция задается как:

// k(x^p,x^q) = s2 * \delta(p,q)
//
// где s2 — это шумовая дисперсия, а \delta(p,q) — это дельта-функция Кронекера
// которая равна 1, если p=q, и нулю в противном случае. 
// Гиперпараметр — это
// logtheta = [ log(sqrt(s2)) ]

nargout_SE_A_t covNoise(double logtheta, int n)
{
  
  nargout_SE_A_t ret = new nargout_SE_A_t(n);
  
  double s2  = exp_db(2.0 * logtheta);
  
  for(int i = 0; i < n ; i++)
    for(int j = 0; j < n ; j++)
    {
      
      if(i == j)
        ret.A[i][j] = s2;
      else
        ret.A[i][j] = 0;
      
    }
    
  return ret;
  
}

nargout_SE_A_t covNoise(double logtheta, int n, int z)
{
  
  nargout_SE_A_t ret = new nargout_SE_A_t(n);
  
  double s2  = exp_db(2.0 * logtheta);
  
  // A = 2*s2*eye(size(x,1));
    
  for(int i = 0; i < n ; i++)
    for(int j = 0; j < n ; j++)
    {
      
      if(i == j)
        ret.A[i][j] = 2.0 * s2;
      else
        ret.A[i][j] = 0;
      
    }
    
  return ret;
  
}

// Регрессия гауссовского процесса с именованной функцией ковариации. 
// Функция возвращает минус логарифм правдоподобия и его частные производные по гиперпараметрам; 
// этот режим используется для подгонки гиперпараметров. 

nargout_gpr_t gpr(char covfunc, double loghyper[], double x[][], double y[])
{
  
  int   n                    = x.length;                
  double Klock[][]           = new double [n][n],
         //Llock[][]           = new double [n][n],
         W1[][]              = new double [n][n],
         W[][]               = new double [n][n],
         eye[][]             = new double [n][n],
         alpha[]             = new double [n],
         sumTm[]             = new double [n];
  nargout_SE_A_t  ret        = new nargout_SE_A_t(n);
  nargout_gpr_t   out        = new nargout_gpr_t(loghyper.length);
  
  // Вычисляем матрицу K
  //Klock = covSum(covfunc, loghyper, x).A;
  Klock = covSum(covfunc, loghyper, x).A;
  // Вычисляем матрицу L - разложение Холецкого для матрицы K
  Chol_S(Klock, n);
  alpha = solve_hol(y, n);
  
  // Вычиления значения логарифма функции правдоподобия
  out.out1 = 0;
  for(int i = 0; i < n ; i++)
      out.out1 += log_db(Sch[i + 1][i + 1]);
  out.out1 = 0.5 * matrixMultiplyC_V(y, alpha) + out.out1 + 0.5 * n * log_db(2.0 * PI); 
  // Вычисление градиента по гиперпараметрам лог функции правдоподобия
 
  eye = eye(n);
  
  // Вычисление выражения  W = L'\(L\eye(n))-alpha*alpha'
  // L'\(L\eye(n))
  W1 = solve(eye, n);
  // alpha*alpha'
  W  = Outer_Product(alpha, alpha);
  // Вычитание получившихся матриц L'\(L\eye(n))-alpha*alpha'
  W  = matrixSUB(W1, W);
  
  // Вычисление суммарной ковариации sum(sum(W.*feval(covfunc{:}, logtheta, x, i)))/2
  // Применяется matrixDOT - вычисление поэлементное умножение матриц W.*feval(covfunc{:}, logtheta, x, i)
  // Применяется SUM - вычисление суммы элементов слобцов матрицы и запись их в вектор
  for(int i = 0; i < (loghyper.length - 1) ; i++)
  {
    // feval(covfunc{:}, logtheta, x, i)
    ret = covSum('S', loghyper, x, i);
    
    out.out2[i] = 0;
    sumTm = SUM(matrixDOT(W, ret.A));
    for(int j = 0; j < n ; j++)
      out.out2[i] += sumTm[j];
      
     out.out2[i] /= 2.0;
       
  }
  
  // Теже действия для второй ковариационной функции covNoise
  ret = covNoise(loghyper[loghyper.length - 1], n, 1);
  sumTm = SUM(matrixDOT(W, ret.A));
  
  for(int j = 0; j < n ; j++)
      out.out2[loghyper.length - 1] += sumTm[j];
      
   out.out2[loghyper.length - 1] /= 2.0;
  
  return out;
  
}

// Оболочка для обучения GP (через gpr), штрафующая за большое отношение сигнал/шум и экстремальные масштабы длины, 
// чтобы избежать числовых нестабильностей
// На вход D + 2 гиперпараметров
// Входные аргументы:

// lh log-hyper-parameters [D+2 x E ]
// covfunc функция ковариации, например,
// covfunc = {'covSum', {'covSEard', 'covNoise'}};
// x обучающие входы [n x D ]
// y обучающие цели [n x E ]
// curb (необязательные) параметры для штрафования экстремальных гиперпараметров
//  .ls масштабы длин
//  .snr отношение сигнал/шум (старайтесь, чтобы оно было ниже 500)
//  std дополнительный параметр, необходимый для штрафа за масштаб длины

// Выходные аргументы:

// f штрафное отрицательное логарифмическое предельное правдоподобие
// df производная штрафного отрицательного логарифмического предельного правдоподобия wrt
// GP логарифмические гиперпараметры

nargout_f_t hypCurb2(double loghyper[], double x[][], double y[], curb_t curb)
{
 
  
  double          sum     = 0,
                  dif     = 0;
  int             p       = 30,
                  D       = x[0].length;
  nargout_gpr_t   out     = new nargout_gpr_t(loghyper.length);
  nargout_f_t     ret     = new nargout_f_t();
  
  out = gpr('D', loghyper, x, y);
  
  // Обработка гиперпараметров масштаба-длины
  sum = 0;
  for(int i = 0; i < D ; i++)
  {
    
    double div = loghyper[i] - log_db(curb._std[i]);
    sum += pow_db(div / log_db(curb.ls), p);
    
    dif = (double)p * pow_db(div, p - 1) / pow_db(log_db(curb.ls), p);
    out.out2[i] += dif;
    
  }
  
  // Добавление штрафа и изменение производных
  out.out1 += sum;
  
  // Обработка сигнал / шум. Оставшиеся гиперпараметры
  
  sum = pow_db((loghyper[D] - loghyper[D + 1]) / log_db(curb.snr), p);
  out.out1 += sum;
  
  dif = (double)p * pow_db(loghyper[D] - loghyper[D + 1], p - 1) / pow_db(log_db(curb.snr), p);
  out.out2[D] += dif;
  
  out.out2[D + 1] -= dif;
  
  ret.fx = out.out1;
  ret.dfx = out.out2;//vectorCoef(out.out2, 1.0);
  return ret;
  
}

nargout_gpr_t hypCurb(char covfunc, double loghyper[], double x[][], double y[], curb_t curb)
{
 
  
  double          sum     = 0,
                  dif     = 0;
  int             p       = 30,
                  D       = x[0].length;
  nargout_gpr_t   out     = new nargout_gpr_t(loghyper.length);
  out = gpr(covfunc, loghyper, x, y);
  
  // Обработка гиперпараметров масштаба-длины
  sum = 0;
  for(int i = 0; i < D ; i++)
  {
    
    double div = loghyper[i] - log_db(curb._std[i]);
    sum += pow_db(div / log_db(curb.ls), p);
    
    dif = (double)p * pow_db(div, p - 1) / pow_db(log_db(curb.ls), p);
    out.out2[i] += dif;
    
  }
  
  // Добавление штрафа и изменение производных
  out.out1 += sum;
  
  // Обработка сигнал / шум. Оставшиеся гиперпараметры
  
  sum = pow_db((loghyper[D] - loghyper[D + 1]) / log_db(curb.snr), p);
  out.out1 += sum;
  
  dif = (double)p * pow_db(loghyper[D] - loghyper[D + 1], p - 1) / pow_db(log_db(curb.snr), p);
  out.out2[D] += dif;
  
  out.out2[D + 1] -= dif;
  
  return out;
  
}

// Вычисляет совместные прогнозы для нескольких GP с неопределенными входными данными. 
// Настройка

dynmodel_t gp0(dynmodel_t _gpmodel)
{
   dynmodel_t ret    = new dynmodel_t(0, 0, 0);
   int       n       = _gpmodel.inputs.length,
             D       = _gpmodel.inputs[0].length,
             E       = _gpmodel.targets[0].length,
             i       = 0,
             j       = 0,
             g       = 0;
             
   double    X[][]   = new double[D + 2][E],
             L[][],
             inp[][] = new double[n][D];  
             
   double    iK[][][],
             _K[][][],
             beta[][];
   
   iK   = new double[E][n][n];
   _K   = new double[E][n][n];
   beta = new double[n][E];
   
   ret.iK = new double[E][n][n];
   ret._K = new double[E][n][n];
   ret.beta = new double[n][E];
   
   X = _gpmodel.hyp;//coeffProdMat(1.0, _gpmodel.hyp);
  
   for(i = 0; i < E; i++)
   {
     
     double X_D[] = new double[D];
     for(j = 0 ; j < D ; j++)
        X_D[j] = exp_db(X[j][i]);
     for(j = 0 ; j < n ; j++)
       for(g = 0; g < D ; g++)
        inp[j][g] = _gpmodel.inputs[j][g] / X_D[g];
    
    double MH[][] = maha(inp, inp);
    for(j = 0 ; j < n ; j++)
       for(g = 0; g < n ; g++)
       _K[i][j][g] = exp_db(2.0 * X[D][i] - MH[j][g] * 0.5);
   
     // L = chol(K(:,:,i) + exp(2*X(D+2,i))*eye(n))'
     L = eye(n);
     for(j = 0 ; j < n ; j++)
       for(g = 0; g < n ; g++)
        L[j][g] *= exp_db(2.0 * X[D + 1][i]);
        L = matrixADD(_K[i], L);
        
        L = matrixMultiply(L, eye(n));
        Chol_S(L, n);
        L = getSch();
     
     //beta(:,i) = L'\(L\gpmodel.targets(:,i));
     double tg[] = new double[n];
     for(j = 0 ; j < n ; j++)
       tg[j] = _gpmodel.targets[j][i];
       
     tg = solve_hol(tg, n);
     // tg = operInvSlash(transMat(L), operInvSlash(L, tg));// Альтернативный вариант, результат тот же
     
     //iK(:,:,i) = L'\(L\eye(n))
     
     //iK[i] = solve(eye(n), n);
     iK[i] = operInvSlash(transMat(L), operInvSlash(L, eye(n)));// Альтернативный вариант, результат тот же
     
     for(j = 0 ; j < n ; j++)
     {
       
       beta[j][i] = tg[j];
       
     }
  
   }
  
  for(int i1 = 0; i1 < iK.length ; i1 ++)
    ret.iK[i1] = iK[i1];//coeffProdMat(1.0, iK[i1]);
  for(int i1 = 0; i1 < _K.length ; i1 ++)
    ret._K[i1] = _K[i1];//coeffProdMat(1.0, _K[i1]);
  ret.beta = beta;//coeffProdMat(1.0, beta);

  ret.inputs = _gpmodel.inputs;//coeffProdMat(1.0, _gpmodel.inputs);
  ret.targets =_gpmodel.targets;// coeffProdMat(1.0, _gpmodel.targets);
  ret.hyp =  _gpmodel.hyp;//coeffProdMat(1.0, _gpmodel.hyp);
   
  return ret;
}

// Вычисляет совместные прогнозы для нескольких GP с неопределенными входными данными. 
// Прогнозирование

nargout_gp0_t _gp0(dynmodel_t _gpmodel, double m[], double s[][])
{
  
  int        n       = _gpmodel.inputs.length,
             D       = _gpmodel.inputs[0].length,
             E       = _gpmodel.targets[0].length,
             i       = 0,
             j       = 0,
             g       = 0;
  double     k[][]   = new double[n][E],
             M[]     = new double[E],
             V[][]   = new double[D][E],
             S[][]   = new double[E][E],
             inp[][] = new double[n][D],
             ii[][]  = new double[n][D],
             X[][]   = new double[D + 2][D + 2];
  
  nargout_gp0_t ret = new nargout_gp0_t();
  
  X = _gpmodel.hyp;//coeffProdMat(1.0, _gpmodel.hyp);
  
  //inp = bsxfun(@minus,gpmodel.inputs,m');  
  for(i = 0; i < _gpmodel.inputs.length; i++)
  {
   for(j = 0; j < _gpmodel.inputs[0].length; j++)
     inp[i][j] = _gpmodel.inputs[i][j] - m[j];
  }
  
//   2) compute predicted mean and inv(s) times input-output covariance
  for(i = 0; i < E; i++)
  {
   
     double iL[][] = new double [D][D];
     
     for(j = 0 ; j < D ; j++)
       iL[j][j] = exp_db(-X[j][i]);
     
      
     double iN[][] = new double[n][D];
     iN = matrixMultiply(inp, iL);
    
     double B[][];
     B = matrixMultiply(iL, s);
     B = matrixMultiply(B, iL);
     for(j = 0 ; j < D ; j++)
       B[j][j] += 1.0;
   
    
    double t[][] = solveMatM(iN, B);
    
    //l = exp(-sum(in.*t,2)/2)
    double l[] = EXP_F(vectorCoef(SUM(matrixDOT(iN, t), 2), -0.5));
    
    //lb = l.*beta(:,i);
    double lb[] = vectorDOT(l, Mat_to_Vec(_gpmodel.beta, i));
    //printVecTst(Mat_to_Vec(_gpmodel.beta, i));
    //println(SUM(lb));
    //tiL = t*iL
    double tiL[][] = matrixMultiply(t, iL);
   
    //c = exp(2*X(D+1,i))/sqrt(det(B))
    double c = exp_db(2.0 * X[D][i]) / sqrt_db(detMatrix(B));
     
    //M(i) = sum(lb)*c
    M[i] = SUM(lb) * c;
  
    //V(:,i) = tiL'*lb*c                    inv(s) times input-output covariance    
    double tiLT_lb_c[] = vectorCoef(matrixMultiplyC(transMat(tiL), lb), c);
    V = Vec_to_Mat(V, tiLT_lb_c, i);
  
    //k(:,i) = 2*X(D+1,i)-sum(in.*in,2)/2
    double sum2[] = vectorCoef(SUM(matrixDOT(iN, iN), 2), 0.5);
    for(j = 0 ; j < sum2.length ; j++)
    {
      k[j][i] = 2.0 * X[D][i] - sum2[j];
    }
 
  }
  
//  3) ompute predictive covariance, non-central moments
  for(i = 0; i < E; i++)
  {
    // ii = bsxfun(@rdivide,inp,exp(2*X(1:D,i)'));
    double X_D[] = new double[D];
     for(j = 0 ; j < D ; j++)
        X_D[j] = exp_db(2.0 * X[j][i]);
     for(j = 0 ; j < n ; j++)
       for(g = 0; g < D ; g++)
        ii[j][g] = inp[j][g] / X_D[g];
     
     for(int l = 0; l <= i ; l++) //  Для корректной работы i + 1 !!!
     {
       //R = s*diag(exp(-2*X(1:D,i))+exp(-2*X(1:D,j)))+eye(D); 
       
       double R[][] = new double [D][D];
       for(j = 0 ; j < D ; j++)
        X_D[j] = exp_db(-2.0 * X[j][i]) + exp_db(-2.0 * X[j][l]);
       R = matrixMultiply(s, diag(X_D));
       R = matrixADD(R, eye(D));
    
        // t = 1/sqrt(det(R));
        
        double t_d = 1.0;
        t_d = 1.0 / sqrt_db(detMatrix(R));
     
        // ij = bsxfun(@rdivide,inp,exp(2*X(1:D,j)'));
        double ij[][] = new double[inp.length][inp[0].length];
        for(j = 0 ; j < D ; j++)
          X_D[j] = exp_db(2.0 * X[j][l]);
        for(j = 0 ; j < n ; j++)
         for(g = 0; g < D ; g++)
          ij[j][g] = inp[j][g] / X_D[g];
       
        // L = exp(bsxfun(@plus,k(:,i),k(:,j)')+maha(ii,-ij,R\s/2))
        double L[][] = new double [k.length][k.length];
        double sb[] = new double [k.length],
               st[] = new double [k.length];
        sb = getSlb(k, i);
        st = getSlb(k, l);
        
        L = matrixADD(bsxfun_plus(sb, st), maha(ii, matrixCoef(ij, -1.0), matrixCoef(operInvSlash(R, s), 0.5)));
        L = EXP_F(L);
     
        if(i == l)
        {
          //S(i,i) = t*(beta(:,i)'*L*beta(:,i) - sum(sum(iK(:,:,i).*L)));
          double sbb[] = new double [_gpmodel.beta.length];
          sbb = getSlb(_gpmodel.beta, i);
         
          S[i][i] = t_d * (matrixMultiplyC_V(matrixMultiplyC_T2(sbb, L), sbb) -
                          SUM(SUM(matrixDOT(_gpmodel.iK[i], L), 1))); 
        }
        else
        {
          //S(i,j) = beta(:,i)'*L*beta(:,j)*t; 
          //S(j,i) = S(i,j);
          double sbI[] = new double [_gpmodel.beta.length],
                 sbJ[] = new double [_gpmodel.beta.length];
          sbI = getSlb(_gpmodel.beta, i);
          sbJ = getSlb(_gpmodel.beta, l);
          S[i][l] = t_d * matrixMultiplyC_V(matrixMultiplyC_T2(sbI, L), sbJ);
          S[l][i] = S[i][l];
        }
        
     }
    
    
     S[i][i] += exp_db(2.0 * X[D][i]);
    
  }
  
  // 4) centralize moments
  //S = S - M*M';   
  S = matrixSUB(S, (Outer_Product(M, M)));
    
  ret.M =M; //vectorCoef(M, 1.0);
  ret.V = V;//coeffProdMat(1.0, V);
  ret.S = S;//coeffProdMat(1.0, S);
 
  return ret;
  
}

//function [M, S, V, dMdm, dSdm, dVdm, dMds, dSds, dVds] = gp0d(gpmodel, m, s)

nargout_gp0d_t _gp0d(dynmodel_t _gpmodel, double m[], double s[][])
{
  
  int        n       = _gpmodel.inputs.length,
             D       = _gpmodel.inputs[0].length,
             E       = _gpmodel.targets[0].length;
             
  //k = zeros(n,E); M = zeros(E,1); V = zeros(D,E); S = zeros(E)      

  double     k[][]   = new double[n][E],
             M[]     = new double[E],
             V[][]   = new double[D][E],
             S[][]   = new double[E][E];
             
  //dMds = zeros(E,D,D); dSdm = zeros(E,E,D)
  //dSds = zeros(E,E,D,D); dVds = zeros(D,E,D,D); T = zeros(D)
  double  dMds[][][]   = new double[E][D][D],
          dSdm[][][]   = new double[E][E][D],
          dSds[][][][] = new double[E][E][D][D],
          dVds[][][][] = new double[D][E][D][D], 
          T[][]        = new double[D][D],
          inp[][]      = new double[_gpmodel.inputs.length][_gpmodel.inputs[0].length];
       
          
  nargout_gp0d_t ret = new nargout_gp0d_t();
  
  double input[][]   =  _gpmodel.inputs;//coeffProdMat(1.0, _gpmodel.inputs);  
//  double target[][]  = coeffProdMat(1.0, _gpmodel.targets); 
  double X[][]       =  _gpmodel.hyp;//coeffProdMat(1.0, _gpmodel.hyp);
 
  //inp = bsxfun(@minus,gpmodel.inputs,m');  
  for(int i = 0; i < input.length; i++)
  {
   for(int j = 0; j < input[0].length; j++)
     inp[i][j] = input[i][j] - m[j];
  }

//   2) compute predicted mean and inv(s) times input-output covariance
  
  for(int i = 0; i < E; i++)
  {
     //iL = diag(exp(-X(1:D,i)));  inverse length-scales
     double iL[][] = new double[D][D];
     for(int j = 0 ; j < D ; j++)
       iL[j][j] = exp_db(-X[j][i]);
     
     //in = inp*iL
     double iN[][] = matrixMultiply(inp, iL);
     
     //B = iL*s*iL+eye(D)
     double B[][] = new double[D][D];
     B = matrixMultiply(iL, s); B = matrixMultiply(B, iL);
     B = matrixADD(B, eye(D));
     
     //LiBL = iL/B*iL
     double LiBL[][] = matrixMultiply(solveMat(iL, B), iL);
     
     //t = in/B
     double t[][] = solveMatM(iN, B);
    
    //l = exp(-sum(in.*t,2)/2)
    double l[] = EXP_F(vectorCoef(SUM(matrixDOT(iN, t), 2), -0.5));
       
    //lb = l.*beta(:,i);
    double lb[] = vectorDOT(l, Mat_to_Vec(_gpmodel.beta, i));
  
    //tL = t*iL
    double tL[][] = matrixMultiply(t, iL);
   
    //tlb = bsxfun(@times,tL,lb)
    double tlb[][] = matrixDOT_T(tL, lb);
    
    //c = exp(2*X(D+1,i))/sqrt(det(B))
    double c = exp_db(2.0 * X[D][i]) / sqrt_db(detMatrix(B));
   
    //M(i) = sum(lb)*c; 
    M[i] = SUM(lb) * c;
    
    //V(:,i) = tL'*lb*c;                    inv(s) times input-output covariance
    double tiLT_lb_c[] = vectorCoef(matrixMultiplyC(transMat(tL), lb), c);
    V = Vec_to_Mat(V, tiLT_lb_c, i);
  
    //dMds(i,:,:) = c*tL'*tlb/2 - LiBL*M(i)/2
    dMds[i] = coeffProdMat(c * 0.5, matrixMultiply(transMat(tL), tlb));
    dMds[i] = matrixSUB(dMds[i], coeffProdMat(M[i] * 0.5, LiBL));
 
    //for d = 1:D
    //dVds(d,i,:,:) = c*bsxfun(@times,tL,tL(:,d))'*tlb/2 - LiBL*V(d,i)/2 ...
    //  - (V(:,i)*LiBL(d,:) + LiBL(:,d)*V(:,i)')/2
    //end
    for(int d = 0; d < D; d++)
    {
       double dVds_[][] = matrixMultiply(transMat(matrixCoef(matrixDOT_T(tL, Mat_to_Vec(tL, d)), c * 0.5)), tlb);
       dVds_ = matrixADD(dVds_, matrixCoef(LiBL, -V[d][i] * 0.5));
       dVds_ = matrixADD(dVds_, matrixCoef(Outer_Product(Mat_to_Vec(V, i), Mat_to_VecT(LiBL, d)), -0.5));
       dVds_ = matrixADD(dVds_, matrixCoef(Outer_Product(Mat_to_Vec(LiBL, d), Mat_to_Vec(V, i)), -0.5));
       for(int j11 = 0; j11 < dVds_.length; j11++)
       for(int j12 = 0; j12 < dVds_[0].length; j12++)
         dVds[d][i][j11][j12] = dVds_[j11][j12];    
    }
    
    //k(:,i) = 2*X(D+1,i)-sum(in.*in,2)/2;
    double sum[] = vectorCoef(SUM(matrixDOT(iN, iN), 2), 0.5);
    for(int j1 = 0 ; j1 < k.length ; j1++)
    {
      k[j1][i] = 2.0 * X[D][i] - sum[j1];
    }
   
  }
  
  //dMdm = V'; dVdm = 2*permute(dMds,[2 1 3])            derivatives wrt m
  double dMdm[][] =  transMat(V);//coeffProdMat(1.0, transMat(V));
  double dVdm[][][] = new double[dMds[0].length][dMds.length][dMds[0][0].length];
     for(int j13 = 0; j13 < dMds[0][0].length; j13++)
     for(int j11 = 0; j11 < dMds.length; j11++)
       for(int j12 = 0; j12 < dMds[0].length; j12++)
         dVdm[j12][j11][j13] = 2.0 * dMds[j11][j12][j13]; 
  
  //iell2 = exp(-2*gpmodel.hyp(1:D,:))
  int ind[] = createInd(0, D);
  double iell2[][] = EXP_F(coeffProdMat(-2.0, getMat(_gpmodel.hyp, ind)));
  
  //inpiell2 = bsxfun(@times,inp,permute(iell2,[3,1,2])); % N-by-D-by-E
  double inpiell2[][][] = new double[inp.length][inp[0].length][iell2[0].length];
  
     for(int j13 = 0; j13 < iell2[0].length; j13++)// 4
     for(int j11 = 0; j11 < inp.length; j11++)//40
       for(int j12 = 0; j12 < inp[0].length; j12++)//6
         inpiell2[j11][j12][j13] = inp[j11][j12] * iell2[j12][j13]; 
         
  //3) compute predictive covariance matrix, non-central moments    
  
  double ii[][];
  for(int i = 0; i < E ; i++)
  {
    //ii = inpiell2(:,:,i)
    ii = new double[inpiell2.length][inpiell2[0].length];
    for(int j11 = 0; j11 < inpiell2.length; j11++)//40
       for(int j12 = 0; j12 < inpiell2[0].length; j12++)//6
         ii[j11][j12] = inpiell2[j11][j12][i]; 
        
    for(int j = 0; j <= i ; j++)
    {
      //R = s*diag(iell2(:,i)+iell2(:,j))+eye(D)
      double R[][] = new double [D][D];
      R = matrixMultiply(s, diag(vectorADD(Mat_to_Vec(iell2, i), Mat_to_Vec(iell2, j))));
      R = matrixADD(R, eye(D));
      
      //t = 1/sqrt(det(R));
      double t_d = detMatrix(R);
        t_d = 1.0 / sqrt_db(t_d);
        
      //ij = inpiell2(:,:,j)
      double ij[][] = new double[inpiell2.length][inpiell2[0].length];
      for(int j11 = 0; j11 < inpiell2.length; j11++)
       for(int j12 = 0; j12 < inpiell2[0].length; j12++)
         ij[j11][j12] = inpiell2[j11][j12][j]; 
      
       //L = exp(bsxfun(@plus,k(:,i),k(:,j)')+maha(ii,-ij,R\s/2))
       double L[][] = new double [k.length][k.length];
       double sb[]  = new double [k.length],
              st[] = new double [k.length];
        sb = getSlb(k, i);
        st = getSlb(k, j);
        L = matrixADD(bsxfun_plus(sb, st), maha(ii, matrixCoef(ij, -1.0), matrixCoef(operInvSlash(R, s), 0.5)));
        L = EXP_F(L);
        double r[];
        if(i == j)
        {
          //iKL = iK(:,:,i).*L
          double iKL[][] = matrixDOT(_gpmodel.iK[i], L);
          
          //s1iKL = sum(iKL,1)
          double s1iKL[] = SUM(iKL,1);
          
          //s2iKL = sum(iKL,2)
          double s2iKL[] = SUM(iKL,2);
          
          //S(i,j) = t*(beta(:,i)'*L*beta(:,i) - sum(s1iKL))
          double sbb[] = getSlb(_gpmodel.beta, i);
          S[i][j] = t_d * (matrixMultiplyC_V(matrixMultiplyC_T2(sbb, L), sbb) -
                          SUM(s1iKL)); 
          
          //zi = ii/R
          double zi[][] = solveMatM(ii, R);
      
          //bibLi = L'*beta(:,i).*beta(:,i); cbLi = L'*bsxfun(@times, beta(:,i), zi)
          double bibLi[] = vectorDOT(matrixMultiplyC(transMat(L), sbb), sbb); 
          double cbLi[][] = matrixMultiply(transMat(L), bsxfun_times(sbb, zi));
          
          //r = (bibLi'*zi*2 - (s2iKL' + s1iKL)*zi)*t
          r = vectorCoef(matrixMultiplyC_T2(bibLi, zi), 2.0);
          r = vectorSUB(r, matrixMultiplyC_T2(vectorADD(s2iKL, s1iKL), zi));
          r = vectorCoef(r, t_d);
          
          //float zi_1_d[][];
          //float zi_d[];
          for(int d = 0; d < D ; d++)
          {
          //  T(d,1:d) = 2*(zi(:,1:d)'*(zi(:,d).*bibLi) + ...
          //cbLi(:,1:d)'*(zi(:,d).*beta(:,i)) - zi(:,1:d)'*(zi(:,d).*s2iKL) ...
          //- zi(:,1:d)'*(iKL*zi(:,d)))
            double zi_1_d[][] = getMat_to_Mat(zi, 0, d); 
            double zi_d[][] = getMat_to_Mat(zi, d);
            double T_1_d[][] = matrixMultiply(transMat(zi_1_d), matrixDOT_T(zi_d, bibLi));
            T_1_d = matrixADD(T_1_d, matrixMultiply(transMat(getMat_to_Mat(cbLi, 0, d)), 
                matrixDOT(zi_d, getMat_to_Mat(_gpmodel.beta, i))));
            T_1_d = matrixSUB(T_1_d, matrixMultiply(transMat(zi_1_d), 
                matrixDOT_T(zi_d, s2iKL)));
            T_1_d = matrixSUB(T_1_d, matrixMultiply(transMat(zi_1_d), 
                matrixMultiply(iKL, zi_d)));
            T_1_d = coeffProdMat(2.0, T_1_d);
            for(int j1 = 0; j1 < (d + 1) ; j1++)
            {
              T[d][j1] = T_1_d[j1][0];
              //T(1:d,d) = T(d,1:d)';
              T[j1][d] = T_1_d[j1][0];
            }
          }
        }
        else
        {
          //zi = ii/R; zj = ij/R
          double zi[][] = solveMatM(ii, R);
          double zj[][] = solveMatM(ij, R);
          
          //S(i,j) = beta(:,i)'*L*beta(:,j)*t
          //S(j,i) = S(i,j)
          double sbI[] = getSlb(_gpmodel.beta, i),
                 sbJ[] = getSlb(_gpmodel.beta, j);
          S[i][j] = t_d * matrixMultiplyC_V(matrixMultiplyC_T2(sbI, L), sbJ);
          S[j][i] = S[i][j];
          
          //bibLj = L*beta(:,j).*beta(:,i)
          double bibLj[] = vectorDOT(matrixMultiplyC(L, sbJ), sbI);
          
          //bjbLi = L'*beta(:,i).*beta(:,j)
          double bjbLi[] = vectorDOT(matrixMultiplyC(transMat(L), sbI), sbJ);
          
          //cbLi = L'*bsxfun(@times, beta(:,i), zi)
          double cbLi[][] = matrixMultiply(transMat(L), bsxfun_times(sbI, zi));
          
          //cbLj = L*bsxfun(@times, beta(:,j), zj)
          double cbLj[][] = matrixMultiply(L, bsxfun_times(sbJ, zj));
          
          //r = (bibLj'*zi+bjbLi'*zj)*t
          r = matrixMultiplyC_T2(bibLj, zi);
          r = vectorADD(r, matrixMultiplyC_T2(bjbLi, zj));
          r = vectorCoef(r, t_d);
          
          for(int d = 0; d < D ; d++)
          {
                    
          //  T(d,1:d) = zi(:,1:d)'*(zi(:,d).*bibLj) + ...
          //  cbLi(:,1:d)'*(zj(:,d).*beta(:,j)) + zj(:,1:d)'*(zj(:,d).*bjbLi) + ...
          //  cbLj(:,1:d)'*(zi(:,d).*beta(:,i))
          //  T(1:d,d) = T(d,1:d)'
          
            double zi_1_d[][] = getMat_to_Mat(zi, 0, d); 
            double zi_d[][] = getMat_to_Mat(zi, d);
            double zj_1_d[][] = getMat_to_Mat(zj, 0, d); 
            double zj_d[][] = getMat_to_Mat(zj, d);
            
            double T_1_d[][] = matrixMultiply(transMat(zi_1_d), matrixDOT_T(zi_d, bibLj));
            T_1_d = matrixADD(T_1_d, matrixMultiply(transMat(getMat_to_Mat(cbLi, 0, d)),
            matrixDOT_T(zj_d, Mat_to_Vec(_gpmodel.beta, j))));
            T_1_d = matrixADD(T_1_d, matrixMultiply(transMat(zj_1_d), 
                matrixDOT_T(zj_d, bjbLi)));
            T_1_d = matrixADD(T_1_d, matrixMultiply(transMat(getMat_to_Mat(cbLj, 0, d)), 
                matrixDOT(zi_d, getMat_to_Mat(_gpmodel.beta, i))));
           
            for(int j1 = 0; j1 < (d + 1) ; j1++)
            {
              T[d][j1] = T_1_d[j1][0];
              //T(1:d,d) = T(d,1:d)';
              T[j1][d] = T_1_d[j1][0];
            }
            
          }
          
        }
        
      //dSdm(i,j,:) = r - M(i)*dMdm(j,:)-M(j)*dMdm(i,:)
      double dSdm_i_j[] = vectorSUB(r, vectorCoef(Mat_to_VecT(dMdm, j), M[i]));
      dSdm_i_j = vectorSUB(dSdm_i_j, vectorCoef(Mat_to_VecT(dMdm, i), M[j]));

      //dSdm(j,i,:) = dSdm(i,j,:)
      dSdm[i][j] = vectorCoef(dSdm_i_j, 1.0);
      dSdm[j][i] = vectorCoef(dSdm_i_j, 1.0);
      
      //T = (t*T-S(i,j)*diag(iell2(:,i)+iell2(:,j))/R)/2   
      double T_[][] = coeffProdMat(S[i][j], diag(vectorADD(Mat_to_Vec(iell2, i), Mat_to_Vec(iell2, j))));
      T = coeffProdMat(t_d, T);
      T = matrixSUB(T, solveMatM(T_, R));
      T = coeffProdMat(0.5, T);
      
      //T = T - reshape(M(i)*dMds(j,:,:) + M(j)*dMds(i,:,:),D,D)
      T_ = coeffProdMat(M[i], dMds[j]);
      T_ = matrixADD(T_, coeffProdMat(M[j], dMds[i]));
      T = matrixSUB(T, T_);
      
      //dSds(i,j,:,:) = T 
      //dSds(j,i,:,:) = T
      dSds[i][j] = T;//coeffProdMat(1.0, T);
      dSds[j][i] = T;//coeffProdMat(1.0, T);
     
    }
    
    // S(i,i) = S(i,i) + exp(2*X(D+1,i))
    S[i][i] += exp_db(2.0 * X[D][i]);
    
  }
  
  //4) centralize moments
  //S = S - M*M'
  S = matrixSUB(S, bsxfun_times(M, M));
  
  //5) vectorize derivatives
  //dMds = reshape(dMds,[E D*D]);
  //dSds = reshape(dSds,[E*E D*D]); dSdm = reshape(dSdm,[E*E D]);
  //dVds = reshape(dVds,[D*E D*D]); dVdm = reshape(dVdm,[D*E D]);
  //dMds=reshape(dMds,[E D*D]);

  ret.dMds = reshape(dMds, E, D * D);
  
  //dSdm=reshape(dSdm,[E*E D]);
  ret.dSdm = reshape(dSdm, E * E, D);
  
  //dSds=reshape(dSds,[E*E D*D]);
  ret.dSds = reshape(dSds, E * E, D * D);
  
  //dVdm=reshape(dVdm,[D*E D]);
  ret.dVdm = reshape(dVdm, D * E, D);
  
  //dVds=reshape(dVds,[D*E D*D]);
  ret.dVds = reshape(dVds, D * E, D * D);
 
  ret.dMdm = dMdm;
  ret.M = M;//vectorCoef(M, 1.0);
  ret.V = V;//coeffProdMat(1.0, V);
  ret.S = S;//coeffProdMat(1.0, S);

  return ret;
  
}

dynmodel_t gp2(dynmodel_t _gpmodel)
{
   dynmodel_t ret    = new dynmodel_t(_gpmodel.inputs.length, _gpmodel.inputs[0].length, _gpmodel.targets.length);
   int       n       = _gpmodel.targets.length,
             D       = _gpmodel.inputs[0].length,
             E       = _gpmodel.targets[0].length,
             i       = 0,
             j       = 0,
             g       = 0;
             
   double    X[][]   = new double[_gpmodel.hyp.length][E],
             L[][],
             inp[][] = new double[n][D];  
   
   double iK[][][]   = new double[E][n][n];
   double _K[][][]   = new double[E][n][n];
   double beta[][]   = new double[n][E];
   ret.iK = new double[E][n][n];
   ret._K = new double[E][n][n];
   ret.beta = new double[n][E];
   X = _gpmodel.hyp;//coeffProdMat(1.0, _gpmodel.hyp);

   for(i = 0; i < E; i++)
   {
     
     double X_D[] = new double[D];
     for(j = 0 ; j < D ; j++)
        X_D[j] = exp_db(X[j][i]);
     for(j = 0 ; j < n ; j++)
       for(g = 0; g < D ; g++)
        inp[j][g] = _gpmodel.inputs[j][g] / X_D[g];
   
    
     double MH[][] = maha(inp, inp);
     for(j = 0 ; j < n ; j++)
       for(g = 0; g < n ; g++)
       _K[i][j][g] = exp_db(2.0 * X[D][i] - MH[j][g] * 0.5);
     
     // L = chol(K(:,:,i) + exp(2*X(D+2,i))*eye(n))'
     L = eye(n);
     for(j = 0 ; j < n ; j++)
       for(g = 0; g < n ; g++)
        L[j][g] *= exp_db(2.0 * X[D + 1][i]);
        L = matrixADD(_K[i], L);
        
        //L = matrixMultiply(L, eye(n));
        Chol_S(L, n);
        L = getSch();
       
     //  iK(:,:,i) = L'\(L\eye(n));
     iK[i] = solve(eye(n), n);
     //iK[i] = operInvSlash(transMat(L), operInvSlash(L, eye(n)));// Альтернативный вариант, результат тот же
 
     //beta(:,i) = L'\(L\gpmodel.targets(:,i));
     double tg[] = Mat_to_Vec(_gpmodel.targets, i);
    
     //tg = solve_hol(tg, n);
     tg = operInvSlash(transMat(L), operInvSlash(L, tg));// Альтернативный вариант, результат тот же
           
     for(j = 0 ; j < n ; j++)
     {
       
       beta[j][i] = tg[j];
       
     }

   }
   
  for(int i1 = 0; i1 < iK.length ; i1 ++)
    ret.iK[i1] = iK[i1];//coeffProdMat(1.0, iK[i1]);
  for(int i1 = 0; i1 < _K.length ; i1 ++)
    ret._K[i1] =_K[i1];// coeffProdMat(1.0, _K[i1]);
  ret.beta =beta;// coeffProdMat(1.0, beta);
  
  ret.inputs = _gpmodel.inputs;//coeffProdMat(1.0, _gpmodel.inputs);
  ret.targets =_gpmodel.targets;// coeffProdMat(1.0, _gpmodel.targets);
  ret.hyp =_gpmodel.hyp;// coeffProdMat(1.0, _gpmodel.hyp);
     
  return ret;

}

nargout_gp0_t _gp2(dynmodel_t _gpmodel, double m[], double s[][])
{
  
  int        n       = _gpmodel.inputs.length,
             D       = _gpmodel.inputs[0].length,
             E       = _gpmodel.targets[0].length,
             i       = 0,
             j       = 0,
             g       = 0;
  double     k[][]   = new double[n][E],
             M[]     = new double[E],
             V[][]   = new double[D][E],
             S[][]   = new double[E][E],
             inp[][] = new double[n][D],
             ii[][]  = new double[n][D],
             X[][]   = new double[_gpmodel.hyp.length][E];
  
  nargout_gp0_t ret = new nargout_gp0_t();
  
  X = _gpmodel.hyp;//coeffProdMat(1.0, _gpmodel.hyp);
  
  //inp = bsxfun(@minus,gpmodel.inputs,m');  
  for(i = 0; i < _gpmodel.inputs.length; i++)
  {
   for(j = 0; j < _gpmodel.inputs[0].length; j++)
     inp[i][j] = _gpmodel.inputs[i][j] - m[j];
  }

//   2) compute predicted mean and inv(s) times input-output covariance
  for(i = 0; i < E; i++)
  {
   
     double iL[][] = new double[D][D];
     
     for(j = 0 ; j < D ; j++)
       iL[j][j] = exp_db(-X[j][i]);
     double iN[][] = new double[n][D];
     iN = matrixMultiply(inp, iL);
      
     double B[][] = new double[D][D];
     B = matrixMultiply(iL, s);
     B = matrixMultiply(B, iL);
     for(j = 0 ; j < D ; j++)
       B[j][j] += 1.0;
    
     double t[][] = new double[n][D];
     t = solveMatM(iN, B);
     
    //l = exp(-sum(in.*t,2)/2)
    double l[] = new double[n];
    double sum = 0;
    for(j = 0 ; j < n ; j++)
    {
       sum = 0;
       for(g = 0 ; g < D ; g++)
        sum += iN[j][g] * t[j][g];
        l[j] = exp_db(- sum / 2.0);
    }
    
    //lb = l.*beta(:,i);
    double lb[] = new double[n];
    for(j = 0 ; j < n ; j++)
    {
       lb[j] = l[j] * _gpmodel.beta[j][i];
    }
 
    //tiL = t*iL
    double tiL[][] = matrixMultiply(t, iL); 

    //c = exp(2*X(D+1,i))/sqrt(det(B))
    double detB = detMatrix(B);
    double c = exp_db(2.0 * X[D][i]) / sqrt_db(detB);

    //M(i) = sum(lb)*c; 
    sum = 0;
    for(j = 0 ; j < n ; j++)
      sum += lb[j];
    M[i] = sum * c;

    //V(:,i) = tiL'*lb*c;                    % inv(s) times input-output covariance 
    double tiLT_lb_c[] = vectorCoef(matrixMultiplyC(transMat(tiL), lb), c);
    
    for(j = 0 ; j < tiLT_lb_c.length ; j++)
    {
      V[j][i] = tiLT_lb_c[j];
    }

    //k(:,i) = 2*X(D+1,i)-sum(in.*in,2)/2;
    double sum2[] = vectorCoef(SUM_d(matrixDOT(iN, iN), 2), 0.5);
    
    for(j = 0 ; j < k.length ; j++)
    {
      k[j][i] = 2.0 * X[D][i] - sum2[j];
    }
  }

  // 3) ompute predictive covariance, non-central moments
  for(i = 0; i < E; i++)
  {
    // ii = bsxfun(@rdivide,inp,exp(2*X(1:D,i)'));
    double X_D[] = new double[D];
     for(j = 0 ; j < D ; j++)
        X_D[j] = exp_db(2.0 * X[j][i]);
     for(j = 0 ; j < n ; j++)
       for(g = 0; g < D ; g++)
        ii[j][g] = inp[j][g] / X_D[g];
  
  /////////////////////
     for(int l = 0; l <= i; l++) //  Для корректной работы i + 1 !!!
     {
       //R = s*diag(exp(-2*X(1:D,i))+exp(-2*X(1:D,j)))+eye(D); 
       
       double R[][] = new double [D][D];
       for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = exp_db(-2.0 * X[j1][i]) + exp_db(-2.0 * X[j1][l]);
        R = matrixMultiply(s, diag(X_D));
        R = matrixADD(R, eye(D));
         
        // t = 1/sqrt(det(R));
      
        double t_d = detMatrix(R);
        t_d = 1.0 / sqrt_db(t_d);

        // ij = bsxfun(@rdivide,inp,exp(2*X(1:D,j)'));
        double ij[][] = new double[inp.length][inp[0].length];
        for(j = 0 ; j < D ; j++)
          X_D[j] = exp_db(2.0 * X[j][l]);
        for(int j1 = 0 ; j1 < n ; j1++)
         for(int g1 = 0; g1 < D ; g1++)
          ij[j1][g1] = inp[j1][g1] / X_D[g1];
       
        // L = exp(bsxfun(@plus,k(:,i),k(:,j)')+maha(ii,-ij,R\s/2))
        double L[][] = new double [k.length][k.length];
        double sb[] = new double [k.length],
               st[] = new double [k.length];
        sb = getSlb(k, i);
        st = getSlb(k, l);
        
        L = matrixADD(bsxfun_plus(sb, st), maha(ii, matrixCoef(ij, -1.0), 
           matrixCoef(operInvSlash(R, s), 0.5)));
        L = EXP_F(L);
        
        // S(i,j) = t*beta(:,i)'*L*beta(:,j); S(j,i) = S(i,j);
        
        double sbI[] = new double [_gpmodel.beta.length],
               sbJ[] = new double [_gpmodel.beta.length];
          sbI = getSlb(_gpmodel.beta, i);
          sbI = vectorCoef(sbI, t_d);
          sbJ = getSlb(_gpmodel.beta, l);
          S[i][l] = matrixMultiplyC_V(matrixMultiplyC_T2(sbI, L), sbJ);
          S[l][i] = S[i][l];
         
     }
 
     //S(i,i) = S(i,i) + 1e-6;
     S[i][i] += 1e-6;
      
  }
  
  // 4) centralize moments
  //S = S - M*M';   
  
  S = matrixSUB(S, Outer_Product(M, M));
  
  ret.M = M;// vectorCoef(M, 1.0);
  ret.V = V;//coeffProdMat(1.0, V);
  ret.S = S;//coeffProdMat(1.0, S);
  
  return ret;
  
}

nargout_gp2_t _gp2d(dynmodel_t _gpmodel, double m[], double s[][])
{
  
  int        n       = _gpmodel.inputs.length,
             D       = _gpmodel.inputs[0].length,
             E       = _gpmodel.targets[0].length;
  double     k[][]   = new double[n][E],
             M[]     = new double[E],
             V[][]   = new double[D][E],
             S[][]   = new double[E][E],
             R[][]   = new double[D][D],
             L[][]   = new double[k.length][k.length],
             B[][]   = new double[k.length][k.length],
             iR[][]  = new double[k.length][k.length],
             t[][]   = new double[n][D],
             l[]     = new double[n],
             lb[]    = new double[n],
             inp[][] = new double[n][D],
             inp2[][] = new double[n][D],
             ii[][]  = new double[n][D],
             X[][]   = new double[D + 2][D + 2],
             K2[][]  = new double[n][n],
             dslb[]  = new double[D],
             tliK[][],
             liK[],
             tlb[][],
             detdX[],
             cdX[],
             dldX[][],
             iK2beta[],
             dsi[][],
             sqdi[][],
             sqdiBi[],
             tlbdi2[][],
             dsqdX[][],
             dKdX[][],
             dKdXbeta[],
             dlb[],
             dtdX[][],
             dlbt[],
             dMdm[][],
             tdX[],
             tdXi[],
             tdXj[],
             bLiKi[],
             bLiKj[],
             Q2[][],
             aQ[][],
             bQ[][],
             Q[][],
             RTi[][],
             RTj[][],
             diRi[][],
             diRj[][],
             QdXi[][],
             QdXj[][] = new double[D][D],
             daQi[][],
             dsaQi[],
             dsaQj[],
             dsbQi[],
             dsbQj[],
             dm2i[][],
             dm2j[][],
             dbQi[][],
             dbQj[][],
             daQj[][],
             dm1i[][],
             dm1j[][],
             dmahai[][],
             dmahaj[][],
             LdXi[][],
             LdXj[][];
      //       r[]      = new float[D];
             
  double    c,
            input[][];
//  dMds = zeros(E,D,D); dSdm = zeros(E,E,D); r = zeros(1,D);
//  dSds = zeros(E,E,D,D); dVds = zeros(D,E,D,D); T = zeros(D);
//  tlbdi = zeros(n,D); dMdi = zeros(E,n,D); dMdt = zeros(E,n,E);
//  dVdt = zeros(D,E,n,E); dVdi = zeros(D,E,n,D); dSdt = zeros(E,E,n,E);
//  dSdi = zeros(E,E,n,D); dMdX = zeros(E,D+2,E); dSdX = zeros(E,E,D+2,E);
//  dVdX = zeros(D,E,D+2,E); Z = zeros(n,D);
//  bdX = zeros(n,E,D); kdX = zeros(n,E,D+1);
  double  dMds[][][]   = new double[E][D][D],
          dSdm[][][]   = new double[E][E][D],
          r[]          = new double[D],
          dSds[][][][] = new double[E][E][D][D],
          dVds[][][][] = new double[D][E][D][D], 
          T[][]        = new double[D][D],
          tlbdi[][]    = new double[n][D],
          dMdi[][][]   = new double[E][n][D],
          dMdt[][][]   = new double[E][n][E],
          dVdt[][][][] = new double[D][E][n][E], 
          dVdi[][][][] = new double[D][E][n][D], 
          dSdt[][][][] = new double[E][E][n][E], 
          dSdi[][][][] = new double[E][E][n][D], 
          dMdX[][][]   = new double[E][D + 2][E], 
          dSdX[][][][] = new double[E][E][D + 2][E],
          dVdX[][][][] = new double[D][E][D + 2][E], 
          
          Z[][]        = new double[n][D],
          bdX[][][]    = new double[n][E][D], 
          kdX[][][]    = new double[n][E][D + 1],
          kK[][]       = new double[n][D],
          dVdm[][][]; 
          
  nargout_gp2_t ret = new nargout_gp2_t();
  
  input = _gpmodel.inputs;//coeffProdMat(1.0, _gpmodel.inputs);  
 // target = coeffProdMat(1.0, _gpmodel.targets); 
  X = _gpmodel.hyp;//coeffProdMat(1.0, _gpmodel.hyp);

  //inp = bsxfun(@minus,gpmodel.inputs,m');  
  for(int i = 0; i < input.length; i++)
  {
   for(int j = 0; j < input[0].length; j++)
     inp[i][j] = input[i][j] - m[j];
  }

//   2) compute predicted mean and inv(s) times input-output covariance
  for(int i = 0; i < E; i++)
  {
   
     // K2 = K(:,:,i)+exp(2*X(D+2,i))*eye(n); 
     K2 = _gpmodel._K[i];
     K2 = matrixADD(K2, matrixCoef(eye(n), exp_db(2.0 * X[D + 1][i])));
    
     // inp2 = bsxfun(@rdivide,input,exp(X(1:D,i)'));
     double X_D[] = new double[D];
     for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = exp_db(X[j1][i]);
     for(int j1 = 0 ; j1 < n ; j1++)
       for(int j2 = 0; j2 < D ; j2++)
        inp2[j1][j2] = input[j1][j2] / X_D[j2];
     
     //ii = bsxfun(@rdivide,input,exp(2*X(1:D,i)'))   
     for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = exp_db(2.0 * X[j1][i]);
     for(int j1 = 0 ; j1 < n ; j1++)
       for(int j2 = 0; j2 < D ; j2++)
        ii[j1][j2] = input[j1][j2] / X_D[j2];
       
     //R = s+diag(exp(2*X(1:D,i)))
     R = matrixADD(s, diag(X_D));
       
     //L = diag(exp(-X(1:D,i)))
     for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = exp_db(-X[j1][i]);
     L = diag(X_D);
     
     
     //B = L*s*L+eye(D) 
     B = matrixMultiply(matrixMultiply(L, s), L);
     B = matrixADD(B, eye(D));
         
     //iR = L/B*L;
     iR = matrixMultiply(solveMat(L, B), L);
        
     //t = inp*iR
     t = matrixMultiply(inp, iR);
     
     //l = exp(-sum(t.*inp,2)/2)
     l = EXP_F(vectorCoef(SUM(matrixDOT(t, inp), 2), -0.5));

     //lb = l.*beta(:,i);
     lb = vectorDOT(l, Mat_to_Vec(_gpmodel.beta, i));
     
     //for(int j1 = 0 ; j1 < n ; j1++)
     //{
     //  lb[j1] = l[j1] * _gpmodel.beta[j1][i];
     //}
      
     //tliK = t'*bsxfun(@times,l,iK(:,:,i))
     tliK = matrixMultiply(transMat(t), matrixDOT(l, _gpmodel.iK[i]));
 
     //liK = K2\l
     liK = operInvSlash(K2, l);
     
     //tlb = bsxfun(@times,t,lb);
     tlb = matrixDOT(lb, t);
     
     //c = exp(2*X(D+1,i))/sqrt(det(R))*exp(sum(X(1:D,i)))
     // Диаганальная матрица R
     double detR = detMatrix(R);
     //for(int j1 = 0 ; j1 < R.length ; j1++)
     //  detR *= R[j1][j1];
 
     for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = X[j1][i];
     c = exp_db(2.0 * X[D][i]) / sqrt_db(detR) * exp_db(SUM(X_D));
 
     //detdX = diag(bsxfun(@times,det(R)*iR',2.*exp(2.*X(1:D,i))))
     for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = 2.0 * exp_db(2.0 * X[j1][i]);
     detdX = diag(matrixDOT(X_D, coeffProdMat(detR, transMat(iR))));
      
     //cdX = -0.5*c/det(R).*detdX'+ c.*ones(1,D)
     cdX = vectorCoef(detdX, -0.5 * c / detR);
     for(int j1 = 0 ; j1 < D ; j1++)
       cdX[j1] += c;
      
     //dldX = bsxfun(@times,l,bsxfun(@times,t,2.*exp(2*X(1:D,i)')).*t./2)
     for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = 2.0 * exp_db(2.0 * X[j1][i]);
     dldX = matrixDOT(l, coeffProdMat(0.5, matrixDOT(matrixDOT(t, X_D), t)));
      
     //M(i) = sum(lb)*c; 
     M[i] = SUM(lb) * c;
     
     //iK2beta = K2\beta(:,i)
     double beta_i[] = Mat_to_Vec(_gpmodel.beta, i);
     iK2beta = operInvSlash(K2, beta_i);

     //dMds(i,:,:) = c*t'*tlb/2-iR*M(i)/2
     dMds[i] = coeffProdMat(0.5, matrixMultiply(coeffProdMat(c, transMat(t)), tlb));
     
     dMds[i] = matrixSUB(dMds[i], matrixCoef(iR, M[i] * 0.5));
     
     //dMdX(i,D+2,i) = -c*sum(l.*(2*exp(2*X(D+2,i))*iK2beta))
     dMdX[i][D + 1][i] = -c * SUM(vectorDOT(l, vectorCoef(iK2beta, 2.0 * exp_db(2.0 * X[D + 1][i]))));
    
     //dMdX(i,D+1,i) = -dMdX(i,(i-1)*(D+2)+D+2)
     
     dMdX[i][D][i] = -dMdX[0][i * (D + 1) + D + 1][i];
     
     //dVdX(:,i,D+2,i) = -((l.*(2*exp(2*X(D+2,i))*iK2beta))'*t*c)'
     double dVdX_1D[] = vectorCoef(matrixMultiplyC_T2(
                                  vectorDOT(l, vectorCoef(iK2beta, 2.0 * exp_db(2.0 * X[D + 1][i]))),
                                  t)
                                  , -c);
     for(int j1 = 0 ; j1 < dVdX.length ; j1++)
       dVdX[j1][i][D + 1][i] = dVdX_1D[j1];
        
     //dVdX(:,i,D+1,i) = -dVdX(:,i,D+2,i)
     for(int j1 = 0 ; j1 < dVdX.length ; j1++)
       dVdX[j1][i][D][i] = -dVdX_1D[j1];
       
     //dsi = -bsxfun(@times,inp2,2.*inp2);
     dsi = matrixCoef(matrixDOT(inp2, matrixCoef(inp2, 2.0)), -1.0);
 
     for(int d = 0; d < D; d++)
     {
       //sqdi = K(:,:,i).*bsxfun(@minus,ii(:,d),ii(:,d)')    : 15 x 15
       double ii_d[] = new double [ii.length];
       for(int j1 = 0 ; j1 < ii.length ; j1++)
        ii_d[j1] = ii[j1][d];
       sqdi = matrixDOT(_gpmodel._K[i], bsxfun_minus(ii_d, ii_d));
       
       //sqdiBi = sqdi*beta(:,i)  : 15
       sqdiBi = matrixMultiplyC(sqdi, beta_i);
           
       //tlbdi(:,d) = sqdi*liK.*beta(:,i) + sqdiBi.*liK : 15 2
       double sqdi_d[] = new double [sqdi.length];
       sqdi_d = vectorDOT(matrixMultiplyC(sqdi, liK), beta_i);
       sqdi_d = vectorADD(sqdi_d, vectorDOT(sqdiBi, liK));
       for(int j1 = 0 ; j1 < sqdi_d.length ; j1++)
         tlbdi[j1][d] = sqdi_d[j1];
     
       //tlbdi2 = -tliK*(-bsxfun(@times,sqdi,beta(:,i))'-diag(sqdiBi)) : 2 15
       tlbdi2 = matrixMultiply(matrixCoef(tliK, -1.0), matrixSUB(matrixCoef(transMat(matrixDOT(beta_i, sqdi)), -1.0), diag(sqdiBi)));
              
       //dVdi(:,i,:,d) = c*(iR(:,d)*lb' - bsxfun(@times,t,tlb(:,d))' + tlbdi2)
       double iR_d[] = new double [iR.length];
       for(int j1 = 0; j1 < iR.length; j1++)
         iR_d[j1] = iR[j1][d];
       double tlb_d[] = new double [tlb.length];
       for(int j1 = 0; j1 < tlb.length; j1++)
         tlb_d[j1] = tlb[j1][d];
       double dVdi_i_d[][] = Outer_Product(iR_d, lb);
       dVdi_i_d = matrixSUB(dVdi_i_d, transMat(matrixDOT(tlb_d, t)));
       dVdi_i_d = matrixCoef(matrixADD(dVdi_i_d, tlbdi2), c);
       for(int j11 = 0; j11 < dVdi_i_d.length; j11++)
       for(int j12 = 0; j12 < dVdi_i_d[0].length; j12++)
         dVdi[j11][i][j12][d] = dVdi_i_d[j11][j12];
       
       //dsqdX = bsxfun(@plus,dsi(:,d),dsi(:,d)') + 4.*inp2(:,d)*inp2(:,d)'  : 15 15
       double dsi_d[] = new double [dsi.length];
       for(int j1 = 0; j1 < dsi_d.length; j1++)
         dsi_d[j1] = dsi[j1][d];
       double inp2_d[] = new double [inp2.length];
       for(int j1 = 0; j1 < inp2_d.length; j1++)
         inp2_d[j1] = inp2[j1][d];
       dsqdX = bsxfun_plus(dsi_d, dsi_d);
       dsqdX = matrixADD(dsqdX, matrixCoef(Outer_Product(inp2_d, inp2_d), 4.0));
             
       //dKdX = -K(:,:,i).*dsqdX./2 : 15 15
       dKdX = matrixCoef(matrixDOT(_gpmodel._K[i], dsqdX), -0.5);
              
       //dKdXbeta = dKdX*beta(:,i) : 15
       dKdXbeta = matrixMultiplyC(dKdX, beta_i);
            
       //bdX(:,i,d) = -K2\dKdXbeta : 
       double bdX_i_d[] = new double [bdX.length];
       bdX_i_d = vectorCoef(operInvSlash(K2, dKdXbeta), -1.0);
       for(int j1 = 0; j1 < bdX.length; j1++)
         bdX[j1][i][d] = bdX_i_d[j1];
        
       //dslb(d) = -liK'*dKdXbeta + beta(:,i)'*dldX(:,d) 
       dslb[d] = -matrixMultiplyC_V(liK, dKdXbeta) + matrixMultiplyC_V(beta_i, Mat_to_Vec(dldX, d));
     
       //dlb = dldX(:,d).*beta(:,i) + l.*bdX(:,i,d)
       dlb = vectorDOT(Mat_to_Vec(dldX, d), beta_i);
       dlb = vectorADD(dlb, vectorDOT(l, bdX_i_d));
       
       
       //dtdX = inp*(-bsxfun(@times,iR(:,d),2.*exp(2*X(d,i))*iR(d,:)))
       dtdX = matrixMultiply(inp, matrixCoef(Outer_Product(Mat_to_Vec(iR, d), 
                                   vectorCoef(Mat_to_VecT(iR, d), 2.0 * exp_db(2.0 * X[d][i]))), -1.0));
                                  
       //dlbt = lb'*dtdX + dlb'*t
       dlbt = matrixMultiplyC_T2(lb, dtdX); 
       dlbt = vectorADD(dlbt, matrixMultiplyC_T2(dlb, t));
        
        
       //dVdX(:,i,d,i) = (dlbt'*c + cdX(d)*(lb'*t)');
       double dVdX_[] = vectorADD(vectorCoef(dlbt, c), vectorCoef(matrixMultiplyC_T2(lb, t), cdX[d]));
       for(int j1 = 0; j1 < dVdX.length; j1++)
         dVdX[j1][i][d][i] = dVdX_[j1];
   
     } // d end
     
     //dMdi(i,:,:) = c*(tlbdi - tlb)
     double dMdi_[][] = matrixCoef(matrixSUB(tlbdi, tlb), c);
     for(int j11 = 0; j11 < dMdi_.length; j11++)
      for(int j12 = 0; j12 < dMdi_[0].length; j12++)
         dMdi[i][j11][j12] = dMdi_[j11][j12];
     
     //dMdt(i,:,i) = c*liK'
     double dMdt_[] = vectorCoef(liK, c);
     for(int j1 = 0; j1 < dMdt_.length; j1++)
       dMdt[i][j1][i] = dMdt_[j1];
    
     //dMdX(i,1:D,i) = cdX.*sum(beta(:,i).*l) + c.*dslb
     double dMdX_[] = vectorCoef(cdX, SUM(vectorDOT(beta_i,l)));
     dMdX_ = vectorADD(dMdX_, vectorCoef(dslb, c));
     for(int j1 = 0; j1 < D; j1++)
       dMdX[i][j1][i] = dMdX_[j1];
      
     //v = bsxfun(@rdivide,inp,exp(X(1:D,i)'));
     double v[][] = new double [inp.length][inp[0].length];
     for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = exp_db(X[j1][i]);
     for(int j1 = 0 ; j1 < n ; j1++)
       for(int j2 = 0; j2 < D ; j2++)
        v[j1][j2] = inp[j1][j2] / X_D[j2];
         
     //k(:,i) = 2*X(D+1,i)-sum(v.*v,2)/2
     double kKV[] = vectorCoef(SUM_d((matrixDOT(v, v)), 2), 0.5);
     for(int j1 = 0; j1 < v.length; j1++)
       kK[j1][i] = 2.0 * X[D][i] - kKV[j1];

     //V(:,i) = t'*lb*c
     double vV[] = vectorCoef(matrixMultiplyC(transMat(t), lb), c);
     for(int j1 = 0; j1 < vV.length; j1++)
       V[j1][i] = vV[j1];
     
     for(int d = 0; d < D; d++)
     {
       //dVds(d,i,:,:) = c*bsxfun(@times,t,t(:,d))'*tlb/2 - iR*V(d,i)/2 ...
       //- V(:,i)*iR(d,:)/2 -iR(:,d)*V(:,i)'/2;
       double dVds_[][] = matrixMultiply(transMat(matrixCoef(matrixDOT_T(t, Mat_to_Vec(t, d)), c * 0.5)), tlb);
       dVds_ = matrixADD(dVds_, matrixCoef(iR, -V[d][i] * 0.5));
       dVds_ = matrixADD(dVds_, matrixCoef(Outer_Product(Mat_to_Vec(V, i), Mat_to_VecT(iR, d)), -0.5));
       dVds_ = matrixADD(dVds_, matrixCoef(Outer_Product(Mat_to_Vec(iR, d), Mat_to_Vec(V, i)), -0.5));
       for(int j11 = 0; j11 < dVds_.length; j11++)
       for(int j12 = 0; j12 < dVds_[0].length; j12++)
         dVds[d][i][j11][j12] = dVds_[j11][j12];     
        
       //kdX(:,i,d) = bsxfun(@times,v(:,d),v(:,d))
       double kdX_[] = vectorDOT(Mat_to_Vec(v, d), Mat_to_Vec(v, d));
       for(int j1 = 0; j1 < kdX_.length; j1++)
         kdX[j1][i][d] = kdX_[j1];    
     } // d end
     
     //dVdt(:,i,:,i) = c*tliK
     double dVdt_[][] = matrixCoef(tliK, c);
     for(int j11 = 0; j11 < dVdt_.length; j11++)
       for(int j12 = 0; j12 < dVdt_[0].length; j12++)
         dVdt[j11][i][j12][i] = dVdt_[j11][j12];   
     
     //kdX(:,i,D+1) = 2*ones(1,n)
     for(int j1 = 0; j1 < n; j1++)
         kdX[j1][i][D] = 2.0; 
     }
    
     // dMdm = V'
     dMdm = transMat(V);
     
     //dVdm = 2*permute(dMds,[2 1 3])
     dVdm = new double[dMds[0].length][dMds.length][dMds[0][0].length];
     for(int j13 = 0; j13 < dMds[0][0].length; j13++)
     for(int j11 = 0; j11 < dMds.length; j11++)
       for(int j12 = 0; j12 < dMds[0].length; j12++)
         dVdm[j12][j11][j13] = 2.0 * dMds[j11][j12][j13];

//  3) ompute predictive covariance, non-central moments
  for(int i = 0; i < E; i++)
  {
    // K2 = K(:,:,i)+exp(2*X(D+2,i))*eye(n);
    K2 = _gpmodel._K[i];//coeffProdMat(1.0, _gpmodel._K[i]);
    K2 = matrixADD(K2, matrixCoef(eye(n), exp_db(2.0 * X[D+1][i])));
        
    // ii = bsxfun(@rdivide,inp,exp(2*X(1:D,i)'));
    double X_D[] = new double[D];
     for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = exp_db(2.0 * X[j1][i]);
        
     for(int j1 = 0 ; j1 < n ; j1++)
       for(int j2 = 0; j2 < D ; j2++)
        ii[j1][j2] = inp[j1][j2] / X_D[j2];
   
     for(int j = 0; j <= i; j++) //  Для корректной работы i + 1 !!!
     {
  
      //R = s*diag(exp(-2*X(1:D,i))+exp(-2*X(1:D,j)))+eye(D);    
      for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = exp_db(-2.0 * X[j1][i]) + exp_db(-2.0 * X[j1][j]);
      R = matrixMultiply(s, diag(X_D));
      R = matrixADD(R, eye(D));
    
      // t = 1/sqrt(det(R));
        
       double t_d = detMatrix(R);
        t_d = 1.0 / sqrt_db(t_d);
       
       //if rcond(R) < 1e-15; fprintf('R-matrix in gp2d ill-conditioned'); keyboard; end
       // Проверка на обусловленность матрицы 0 - плохо, 1 - хорошо
       // Пока не буду реализовыать
       
       //iR = R\eye(D)
       iR = operInvSlash(R, eye(D));

       // ij = bsxfun(@rdivide,inp,exp(2*X(1:D,j)'));
       double ij[][] = new double[inp.length][inp[0].length];
       for(int j1 = 0 ; j1 < D ; j1++)
         X_D[j1] = exp_db(2.0 * X[j1][j]);
       for(int j1 = 0 ; j1 < n ; j1++)
        for(int j2 = 0; j2 < D ; j2++)
         ij[j1][j2] = inp[j1][j2] / X_D[j2];
       
       // L = exp(bsxfun(@plus,k(:,i),k(:,j)')+maha(ii,-ij,R\s/2))
       double sb[] = new double [kK.length],
              st[] = new double [kK.length];
       sb = Mat_to_Vec(kK, i);
       st = Mat_to_Vec(kK, j);
       L = matrixADD(bsxfun_plus(sb, st), maha(ii, matrixCoef(ij, -1.0), matrixCoef(operInvSlash(R, s), 0.5)));
       L = EXP_F(L);
            
       //A = beta(:,i)*beta(:,j)'; A = A.*L 
       double A[][] = Outer_Product(Mat_to_Vec(_gpmodel.beta, i), Mat_to_Vec(_gpmodel.beta, j));
       A = matrixDOT(A, L);
       
       //ssA = sum(sum(A))
       double ssA = SUM(SUM(A));
       
       // S(i,j) = t*ssA; S(j,i) = S(i,j);
       S[i][j] = t_d * ssA;
       S[j][i] = S[i][j];
     
       //zzi = ii*(R\s)
       double R_s[][] = operInvSlash(R, s);
       double zzi[][] = matrixMultiply(ii, R_s);
     
       //zzj = ij*(R\s);
       double zzj[][] = matrixMultiply(ij, R_s);
   
       //zi = ii/R; zj = ij/R;
       double zi[][] = solveMatM(ii, R);
       double zj[][] = solveMatM(ij, R);

       //tdX  = -0.5*t*sum(iR'.*bsxfun(@times,s,-2*exp(-2*X(1:D,i)')-2*exp(-2*X(1:D,i)')))
       for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = -4.0 * exp_db(-2.0 * X[j1][i]);
       tdX = vectorCoef(SUM(matrixDOT(transMat(iR), matrixDOT(s, X_D))), -0.5 * t_d);
             
       // tdXi = -0.5*t*sum(iR'.*bsxfun(@times,s,-2*exp(-2*X(1:D,i)')))
       for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = X_D[j1] * 0.5;
       tdXi = vectorCoef(SUM(matrixDOT(transMat(iR), matrixDOT(s, X_D))), -0.5 * t_d);
                
       //tdXj = -0.5*t*sum(iR'.*bsxfun(@times,s,-2*exp(-2*X(1:D,j)')));
       for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = -2.0 * exp_db(-2.0 * X[j1][j]);
       tdXj = vectorCoef(SUM(matrixDOT(transMat(iR), matrixDOT(s, X_D))), -0.5 * t_d);
          
       //bLiKi = iK(:,:,j)*(L'*beta(:,i))
       double beta_i[] = Mat_to_Vec(_gpmodel.beta, i);
       double iK_j[][] = _gpmodel.iK[j];
       bLiKi = matrixMultiplyC(iK_j, matrixMultiplyC(transMat(L), beta_i));
              
       //bLiKj = iK(:,:,i)*(L*beta(:,j))
       double beta_j[] = Mat_to_Vec(_gpmodel.beta, j);
       double iK_i[][] = _gpmodel.iK[i];
       bLiKj = matrixMultiplyC(iK_i, matrixMultiplyC(L, beta_j));
  
       //Q2 = R\s/2
       Q2 = matrixCoef(operInvSlash(R, s), 0.5);
      
       //aQ = ii*Q2
       aQ = matrixMultiply(ii, Q2);
       
       //bQ = ij*Q2
       bQ = matrixMultiply(ij, Q2);
       
       Z = new double[n][D];
       
    
       for(int d = 0; d < D; d++)
       {
        //Z(:,d) = exp(-2*X(d,i))*(A*zzj(:,d) + sum(A,2).*(zzi(:,d) - inp(:,d)))...
        // + exp(-2*X(d,j))*((zzi(:,d))'*A + sum(A,1).*(zzj(:,d) - inp(:,d))')';
        
        double veC[] = vectorDOT(SUM(A, 2),vectorSUB(Mat_to_Vec(zzi, d), Mat_to_Vec(inp, d)));
        veC = vectorADD(matrixMultiplyC(A, Mat_to_Vec(zzj, d)), veC);
        veC = vectorCoef(veC, exp_db(-2.0 * X[d][i]));
        
        double veC2[] = vectorDOT(SUM(A, 1),vectorSUB(Mat_to_Vec(zzj, d), Mat_to_Vec(inp, d)));
        veC2 = vectorADD(matrixMultiplyC_T2(Mat_to_Vec(zzi, d), A), veC2);
        veC2 = vectorCoef(veC2, exp_db(-2.0 * X[d][j]));
        veC = vectorADD(veC, veC2);
        
        for(int j1 = 0 ; j1 < veC.length ; j1++)
          Z[j1][d] = veC[j1];
       
        //Q = bsxfun(@minus,inp(:,d),inp(:,d)');
        Q = bsxfun_minus(Mat_to_Vec(inp, d), Mat_to_Vec(inp, d));
  
        //B = K(:,:,i).*Q;
        B = matrixDOT(_gpmodel._K[i], Q);
      
        //Z(:,d) = Z(:,d)+exp(-2*X(d,i))*(B*beta(:,i).*bLiKj+beta(:,i).*(B*bLiKj));
        double veC3[] = vectorDOT(matrixMultiplyC(B, Mat_to_Vec(_gpmodel.beta, i)), bLiKj);
        double veC4[] = vectorDOT(Mat_to_Vec(_gpmodel.beta, i), matrixMultiplyC(B, bLiKj));
        veC3 = vectorADD(veC3, veC4);
        veC3 = vectorCoef(veC3, exp_db(-2.0 * X[d][i]));
        for(int j1 = 0 ; j1 < veC3.length ; j1++)
          Z[j1][d] += veC3[j1];
  
        //if i~=j; B = K(:,:,j).*Q; end      
        if(i != j)
        {
          B = matrixDOT(_gpmodel._K[j], Q);
        }
       
        //Z(:,d) = Z(:,d)+exp(-2*X(d,j))*(bLiKi.*(B*beta(:,j))+B*bLiKi.*beta(:,j));
        double veC5[] = vectorDOT(bLiKi, matrixMultiplyC(B, Mat_to_Vec(_gpmodel.beta, j)));
        veC5 = vectorADD(veC5, vectorDOT(matrixMultiplyC(B, bLiKi), Mat_to_Vec(_gpmodel.beta, j)));
        veC5 = vectorCoef(veC5, exp_db(-2.0 * X[d][j]));
        for(int j1 = 0 ; j1 < veC5.length ; j1++)
          Z[j1][d] += veC5[j1];
       
        //B = bsxfun(@plus,zi(:,d),zj(:,d)').*A;
        B = bsxfun_plus(Mat_to_Vec(zi, d), Mat_to_Vec(zj, d));
        B = matrixDOT(B, A);
      
        //r(d) = sum(sum(B))*t;
        r[d] = SUM(SUM(B)) * t_d;
            
        //T(d,1:d) = sum(zi(:,1:d)'*B,2) + sum(B*zj(:,1:d))';
//        double zi_1_d[] = new double [zi.length];
//        double zj_1_d[] = new double [zj.length];
//        double zi_1_d2[][] = new double [zi.length][d + 1];
//        double zj_1_d2[][] = new double [zj.length][d + 1];
        
        double T_d_d[] = vectorADD(SUM(matrixMultiply(transMat(getMat(zi, 0, d)) , B), 2), 
                            SUM(matrixMultiply(B, getMat(zj, 0, d)) , 1));
        for(int j0 = 0; j0 < (d + 1); j0++)
        {
          T[d][j0] = T_d_d[j0];
          T[j0][d] = T_d_d[j0];
        }
         
        if(i == j)
        {
          
          //RTi =  bsxfun(@times,s,(-2*exp(-2*X(1:D,i)')-2*exp(-2*X(1:D,j)')));
          for(int j1 = 0 ; j1 < D ; j1++)
            X_D[j1] = -2.0 * exp_db(-2.0 * X[j1][i]) - 2.0 * exp_db(-2.0 * X[j1][j]);
          RTi = matrixDOT(s, X_D);
                   
          //diRi = -R\bsxfun(@times,RTi(:,d),iR(d,:));
          diRi = operInvSlash(matrixCoef(R, -1), bsxfun_times(Mat_to_Vec(RTi, d), Mat_to_VecT(iR, d)));

        }
        else
        {
          
          //RTi = bsxfun(@times,s,-2*exp(-2*X(1:D,i)'));
          for(int j1 = 0 ; j1 < D ; j1++)
            X_D[j1] = -2.0 * exp_db(-2.0 * X[j1][i]);
          RTi = matrixDOT(s, X_D);
          
          //RTj = bsxfun(@times,s,-2*exp(-2*X(1:D,j)'));
          for(int j1 = 0 ; j1 < D ; j1++)
            X_D[j1] = -2.0 * exp_db(-2.0 * X[j1][j]);
          RTj = matrixDOT(s, X_D);
          
          //diRi = -R\bsxfun(@times,RTi(:,d),iR(d,:));
          diRi = operInvSlash(matrixCoef(R, -1), bsxfun_times(Mat_to_Vec(RTi, d), Mat_to_VecT(iR, d)));
          
          //diRj = -R\bsxfun(@times,RTj(:,d),iR(d,:));
          diRj = operInvSlash(matrixCoef(R, -1), bsxfun_times(Mat_to_Vec(RTj, d), Mat_to_VecT(iR, d)));
          
          
          //QdXj = diRj*s/2;
          QdXj = matrixCoef(matrixMultiply(diRj, s), 0.5);
          
        }
        
        //QdXi = diRi*s/2;
        QdXi = matrixCoef(matrixMultiply(diRi, s), 0.5);
  
        if(i == j)
        {
          
         //daQi = ii*QdXi + bsxfun(@times,-2*ii(:,d),Q2(d,:));
         daQi = matrixMultiply(ii, QdXi);
         daQi = matrixADD(daQi, bsxfun_times(vectorCoef(Mat_to_Vec(ii, d), -2.0), Mat_to_VecT(Q2, d)));
     
         //dsaQi = sum(daQi.*ii,2) - 2.*aQ(:,d).*ii(:,d);
         dsaQi = SUM(matrixDOT(daQi, ii), 2);
         dsaQi = vectorADD(dsaQi, vectorDOT(vectorCoef(Mat_to_Vec(aQ, d), - 2), Mat_to_Vec(ii, d)));
         
         dsaQj = dsaQi;//vectorCoef(dsaQi, 1);
          
         //dsbQi = dsaQi; dsbQj = dsbQi;
         dsbQi =dsaQi;// vectorCoef(dsaQi, 1);
         dsbQj =dsbQi;// vectorCoef(dsbQi, 1);
         
         //dm2i = -2*daQi*ii' + 2*(bsxfun(@times,aQ(:,d),ii(:,d)')...
         //      +bsxfun(@times,ii(:,d),aQ(:,d)')); dm2j = dm2i; 
         dm2i = matrixCoef(matrixMultiply(daQi, transMat(ii)), -2.0);    
         dm2i = matrixADD(dm2i,
                               matrixCoef(
                                           matrixADD(bsxfun_times(Mat_to_Vec(aQ, d), Mat_to_Vec(ii, d)),
                                           bsxfun_times(Mat_to_Vec(ii, d), Mat_to_Vec(aQ, d))),
                                           2.0));
         //dm2j = dm2i
         dm2j = matrixCoef(dm2i, 1);
        }
        else
        {
          
          //dbQi = ij*QdXi;  
          dbQi = matrixMultiply(ij, QdXi);
          
          //dbQj = ij*QdXj + bsxfun(@times,-2*ij(:,d),Q2(d,:));
          dbQj = matrixMultiply(ij, QdXj);
          dbQj = matrixADD(dbQj, bsxfun_times(vectorCoef(Mat_to_Vec(ij, d), -2.0), Mat_to_VecT(Q2, d)));
          
          //daQi = ii*QdXi + bsxfun(@times,-2*ii(:,d),Q2(d,:));
          daQi = matrixMultiply(ii, QdXi);
          daQi = matrixADD(daQi, bsxfun_times(vectorCoef(Mat_to_Vec(ii, d), -2.0), Mat_to_VecT(Q2, d)));
          
          //daQj = ii*QdXj;
          daQj = matrixMultiply(ii, QdXj);
          
          //dsaQi = sum(daQi.*ii,2) - 2.*aQ(:,d).*ii(:,d);
          dsaQi = SUM(matrixDOT(daQi, ii), 2);
          dsaQi = vectorADD(dsaQi, vectorDOT(vectorCoef(Mat_to_Vec(aQ, d), - 2), Mat_to_Vec(ii, d)));
          
          //dsaQj = sum(daQj.*ii,2);
          dsaQj = SUM(matrixDOT(daQj, ii), 2);
          
          //dsbQi = sum(dbQi.*ij,2);
          dsbQi = SUM(matrixDOT(dbQi, ij), 2);
          
          //dsbQj = sum(dbQj.*ij,2) - 2.*bQ(:,d).*ij(:,d);
          dsbQj = SUM(matrixDOT(daQi, ij), 2);
          dsbQj = vectorADD(dsbQj, vectorDOT(vectorCoef(Mat_to_Vec(bQ, d), - 2), Mat_to_Vec(ij, d)));
          
          //dm2i = -2*daQi*ij'; 
          dm2i = matrixCoef(matrixMultiply(daQi, transMat(ij)), -2);
          
          //dm2j = -2*ii*(dbQj)';
          dm2j = matrixCoef(matrixMultiply(ii, transMat(dbQj)), -2);
         
        };
        
        //dm1i = bsxfun(@plus,dsaQi,dsbQi');
        dm1i = bsxfun_plus(dsaQi, dsbQi);
               
        //dm1j = bsxfun(@plus,dsaQj,dsbQj');
        dm1j = bsxfun_plus(dsaQj, dsbQj);
        
        //dmahai = dm1i-dm2i;
        dmahai = matrixSUB(dm1i, dm2i);
      
        //dmahaj = dm1j-dm2j;
        dmahaj = matrixSUB(dm1j, dm2j);
        
        if(i == j)
        {
          //LdXi = L.*(dmahai + bsxfun(@plus,kdX(:,i,d),kdX(:,j,d)'));
          double kdX_i_d[] = new double[kdX.length];
          for(int j1 = 0 ; j1 < kdX.length; j1++)
            kdX_i_d[j1] = kdX[j1][i][d];
          LdXi = matrixDOT(L, matrixADD(dmahai, bsxfun_plus(kdX_i_d, kdX_i_d)));
    
          //dSdX(i,i,d,i) = beta(:,i)'*LdXi*beta(:,j);
          dSdX[i][i][d][i] = matrixMultiplyC_V(Mat_to_Vec(_gpmodel.beta, i), matrixMultiplyC(LdXi, Mat_to_Vec(_gpmodel.beta, j)));
        }
        else
        {
          
          //LdXi = L.*(dmahai + bsxfun(@plus,kdX(:,i,d),zeros(n,1)'));
          double zeros_n[] = new double [n];
          double kdX_i_d[] = new double[kdX.length];
          for(int j1 = 0 ; j1 < kdX.length; j1++)
            kdX_i_d[j1] = kdX[j1][i][d];
          LdXi = matrixDOT(L, matrixADD(dmahai, bsxfun_plus(kdX_i_d, zeros_n)));
          
          //LdXj = L.*(dmahaj + bsxfun(@plus,zeros(n,1),kdX(:,j,d)'));
          double kdX_j_d[] = new double[kdX.length];
          for(int j1 = 0 ; j1 < kdX.length; j1++)
            kdX_j_d[j1] = kdX[j1][j][d];
          LdXj = matrixDOT(L, matrixADD(dmahaj, bsxfun_plus(zeros_n, kdX_j_d)));
         
          //dSdX(i,j,d,i) = beta(:,i)'*LdXi*beta(:,j);
          dSdX[i][j][d][i] = matrixMultiplyC_V(Mat_to_Vec(_gpmodel.beta, i), matrixMultiplyC(LdXi, Mat_to_Vec(_gpmodel.beta, j)));
          
          //dSdX(i,j,d,j) = beta(:,i)'*LdXj*beta(:,j);
          dSdX[i][j][d][j] = matrixMultiplyC_V(Mat_to_Vec(_gpmodel.beta, i), matrixMultiplyC(LdXj, Mat_to_Vec(_gpmodel.beta, j)));
         
        }
        
      }// d
      
      // dSdX [E][E][D + 2][E],
    
    
      if(i == j)
      {
        
        //dSdX(i,i,1:D,i) = reshape(dSdX(i,i,1:D,i),D,1) + reshape(bdX(:,i,:),n,D)'*(L+L')*beta(:,i);
        double bdX_n_D[][] = new double [n][D];
        for(int i1 = 0 ; i1 < n ; i1++)
         for(int i2 = 0 ; i2 < D ; i2++)
          bdX_n_D[i1][i2] = bdX[i1][i][i2];
        
        double dSdX_D[] = new double [D];
        for(int i1 = 0 ; i1 < D ; i1++)
          dSdX_D[i1] = dSdX[i][i][i1][i];
        
        double dSdX_D2[] = new double [D];
        dSdX_D2 = vectorADD(dSdX_D, matrixMultiplyC(matrixMultiply(transMat(bdX_n_D), matrixADD(L, transMat(L))), Mat_to_Vec(_gpmodel.beta, i)));
        for(int i1 = 0 ; i1 < D ; i1++)
          dSdX[i][i][i1][i] = dSdX_D2[i1];
        for(int i1 = 0 ; i1 < D ; i1++)
          dSdX_D[i1] = dSdX[i][i][i1][i];
          
        //dSdX(i,i,1:D,i) = reshape(t*dSdX(i,i,1:D,i),D,1)' + tdX*ssA;
        dSdX_D2 = vectorCoef(dSdX_D, t_d);
        dSdX_D2 = vectorADD(dSdX_D2, vectorCoef(tdX, ssA));
        for(int i1 = 0 ; i1 < D ; i1++)
          dSdX[i][i][i1][i] = dSdX_D2[i1];
       
        //dSdX(i,i,D+2,i) = 2*exp(2*X(D+2,i))*t*(-sum(beta(:,i).*bLiKi)-sum(beta(:,i).*bLiKi));
        dSdX[i][i][D + 1][i] = 2.0 * exp_db(2.0 *X[D + 1][i]) * t_d * (- 2.0 * SUM(vectorDOT(Mat_to_Vec(_gpmodel.beta, i), bLiKi)));
        
      }
      else
      {
       
        //dSdX(i,j,1:D,i) = reshape(dSdX(i,j,1:D,i),D,1) + reshape(bdX(:,i,:),n,D)'*(L*beta(:,j));
        double bdX_n_D[][] = new double [n][D];
        for(int i1 = 0 ; i1 < n ; i1++)
         for(int i2 = 0 ; i2 < D ; i2++)
          bdX_n_D[i1][i2] = bdX[i1][i][i2];
        
        double dSdX_D[] = new double [D];
        for(int i1 = 0 ; i1 < D ; i1++)
          dSdX_D[i1] = dSdX[i][j][i1][i];
        
        double dSdX_D2[] = new double [D];
        dSdX_D2 = vectorADD(dSdX_D, matrixMultiplyC(transMat(bdX_n_D), matrixMultiplyC(L, Mat_to_Vec(_gpmodel.beta, j))));
              
        for(int i1 = 0 ; i1 < D ; i1++)
          dSdX[i][j][i1][i] = dSdX_D2[i1];
        
        //dSdX(i,j,1:D,j) = reshape(dSdX(i,j,1:D,j),D,1) + reshape(bdX(:,j,:),n,D)'*(L'*beta(:,i))
        for(int i1 = 0 ; i1 < n ; i1++)
         for(int i2 = 0 ; i2 < D ; i2++)
          bdX_n_D[i1][i2] = bdX[i1][j][i2];
          
        for(int i1 = 0 ; i1 < D ; i1++)
          dSdX_D[i1] = dSdX[i][j][i1][j];
        dSdX_D2 = vectorADD(dSdX_D, matrixMultiplyC(transMat(bdX_n_D), matrixMultiplyC(transMat(L), Mat_to_Vec(_gpmodel.beta, i))));
    
        for(int i1 = 0 ; i1 < D ; i1++)
          dSdX[i][j][i1][j] = dSdX_D2[i1];
          
        //dSdX(i,j,1:D,i) = reshape(t*dSdX(i,j,1:D,i),D,1)' + tdXi*ssA;
        double dSdX_D_i[] = new double [D];
        for(int i1 = 0 ; i1 < D ; i1++)
          dSdX_D_i[i1] = dSdX[i][j][i1][i];
          
        dSdX_D2 = vectorCoef(dSdX_D_i, t_d);
        dSdX_D2 = vectorADD(dSdX_D2, vectorCoef(tdXi, ssA));
        
        for(int i1 = 0 ; i1 < D ; i1++)
          dSdX[i][j][i1][i] = dSdX_D2[i1];
          
        //dSdX(i,j,1:D,j) = reshape(t*dSdX(i,j,1:D,j),D,1)' + tdXj*ssA;
        double dSdX_D_j[] = new double [D];
        for(int i1 = 0 ; i1 < D ; i1++)
          dSdX_D_j[i1] = dSdX[i][j][i1][j];
          
        dSdX_D2 = vectorCoef(dSdX_D_j, t_d);
        dSdX_D2 = vectorADD(dSdX_D2, vectorCoef(tdXj, ssA));
    
        for(int i1 = 0 ; i1 < D ; i1++)
          dSdX[i][j][i1][j] = dSdX_D2[i1];
        
        //dSdX(i,j,D+2,i) = 2*exp(2*X(D+2,i))*t*(-beta(:,i)'*bLiKj);
        dSdX[i][j][D + 1][i] = 2.0 * exp_db(2.0 * X[D + 1][i]) * (-t_d) * matrixMultiplyC_V(Mat_to_Vec(_gpmodel.beta, i), bLiKj);
        
        //dSdX(i,j,D+2,j) = 2*exp(2*X(D+2,j))*t*(-beta(:,j)'*bLiKi);
        dSdX[i][j][D + 1][j] = 2.0 * exp_db(2.0 *X[D + 1][j]) * (-t_d) * matrixMultiplyC_V(Mat_to_Vec(_gpmodel.beta, j), bLiKi);
       
      }
      

      //dSdm(i,j,:) = r - M(i)*dMdm(j,:)-M(j)*dMdm(i,:); dSdm(j,i,:) = dSdm(i,j,:);
      double dSdm_i_j[] = new double[r.length];
      dSdm_i_j = r;//vectorCoef(r, 1);
      dSdm_i_j = vectorSUB(dSdm_i_j, vectorCoef(Mat_to_VecT(dMdm, j), M[i]));
      dSdm_i_j = vectorSUB(dSdm_i_j, vectorCoef(Mat_to_VecT(dMdm, i), M[j]));
      for(int i1 = 0 ; i1 < r.length ; i1++)
      {
       dSdm[i][j][i1] = dSdm_i_j[i1];
       dSdm[j][i][i1] = dSdm_i_j[i1];
      }
    
      //T = (t*T-S(i,j)*diag(exp(-2*X(1:D,i))+exp(-2*X(1:D,j)))/R)/2;
      for(int j1 = 0 ; j1 < D ; j1++)
        X_D[j1] = exp_db(-2.0 * X[j1][i]) + exp_db(-2.0 * X[j1][j]);
     
     T = matrixCoef(T, t_d);
     T = matrixSUB(T, solveMat(matrixCoef(diag(X_D), S[i][j]), R));
     T = matrixCoef(T, 0.5);
           
     //T = T - reshape(M(i)*dMds(j,:,:) + M(j)*dMds(i,:,:),D,D);
     T = matrixSUB(T, matrixADD(matrixCoef(dMds[j], M[i]), matrixCoef(dMds[i], M[j])));
     
     
     //dSds(i,j,:,:) = T; dSds(j,i,:,:) = T;
     for(int i1 = 0 ; i1 < T.length ; i1++)
      for(int i2 = 0 ; i2 < T[0].length ; i2++)
      {
       dSds[i][j][i1][i2] = T[i1][i2];
       dSds[j][i][i1][i2] = T[i1][i2];
      }
      
   
                 if(i == j)
                 {
                   
                   //dSdt(i,i,:,i) = (beta(:,i)'*(L+L'))/(K2)*t - 2*dMdt(i,:,i)*M(i);    
                   double dSdt_i[] = vectorCoef(solveMatV(matrixMultiplyC_T2(beta_i, matrixADD(L, transMat(L))), K2), t_d);
                   double dMdt_i[] = new double [dSdt_i.length];
                   for(int j1 = 0; j1 < dMdt_i.length; j1++)
                     dMdt_i[j1] = dMdt[i][j1][i];
                     
                   dSdt_i = vectorSUB(dSdt_i, vectorCoef(dMdt_i, 2.0 * M[i]));
                   for(int i1 = 0 ; i1 < dSdt_i.length ; i1++)
                    dSdt[i][i][i1][i] = dSdt_i[i1];
               
                   //dSdX(i,j,:,i) = reshape(dSdX(i,j,:,i),1,D+2) - M(i)*dMdX(j,:,j)-M(j)*dMdX(i,:,i);
                   double dSdX_i_j_i[] = new double[D + 2];
                   for(int j1 = 0 ; j1 < (D + 2); j1++)
                    dSdX_i_j_i[j1] = dSdX[i][j][j1][i];
                    
                   double dMdX_j_j[] = new double [D + 2],
                          dMdX_i_i[] = new double [D + 2];
                   for(int j1 = 0 ; j1 < (D + 2); j1++)
                   {
                    dMdX_j_j[j1] = dMdX[j][j1][j];  
                    dMdX_i_i[j1] = dMdX[i][j1][i];
                   }
                   dSdX_i_j_i = vectorSUB(dSdX_i_j_i, vectorCoef(dMdX_j_j, M[i]));
                   dSdX_i_j_i = vectorSUB(dSdX_i_j_i, vectorCoef(dMdX_i_i, M[j]));
                   
                   for(int i1 = 0 ; i1 < dSdX_i_j_i.length ; i1++)
                    dSdX[i][j][i1][i] = dSdX_i_j_i[i1];
                 }
                 else
                 {
                   
                   //dSdt(i,j,:,i) = (beta(:,j)'*L')/(K2)*t - dMdt(i,:,i)*M(j);
                   double dSdt_i_j[] = vectorCoef(solveMatV(matrixMultiplyC_T2(beta_j, transMat(L)), K2), t_d);
                   double dMdt_i[] = new double [dSdt_i_j.length];
                   for(int j1 = 0; j1 < dMdt_i.length; j1++)
                     dMdt_i[j1] = dMdt[i][j1][i];
                     
                   dSdt_i_j = vectorSUB(dSdt_i_j, vectorCoef(dMdt_i, M[j]));
                   
                   for(int i1 = 0 ; i1 < dSdt_i_j.length ; i1++)
                    dSdt[i][j][i1][i] = dSdt_i_j[i1];
                    
                   //dSdt(i,j,:,j) = beta(:,i)'*L/(K(:,:,j)+exp(2*X(D+2,j))*eye(n))*t - dMdt(j,:,j)*M(i);
                   dSdt_i_j = matrixMultiplyC_T2(beta_i, L);
                   dSdt_i_j = vectorCoef(solveMatV(dSdt_i_j, matrixADD(_gpmodel._K[j], 
                               matrixCoef(eye(n), exp_db(2.0 * X[D + 1][j])))), t_d);
                   double dMdt_j[] = new double [dSdt_i_j.length];
                   for(int j1 = 0; j1 < dMdt_j.length; j1++)
                     dMdt_j[j1] = dMdt[j][j1][j];
                   dSdt_i_j = vectorSUB(dSdt_i_j, vectorCoef(dMdt_j, M[i]));
                   for(int i1 = 0 ; i1 < dSdt_i_j.length ; i1++)
                    dSdt[i][j][i1][j] = dSdt_i_j[i1];
                    
                   //dSdt(j,i,:,:) = dSdt(i,j,:,:);
                   for(int i1 = 0 ; i1 < n ; i1++)
                   for(int i2 = 0 ; i2 < E ; i2++)
                   {
                     dSdt[j][i][i1][i2] = dSdt[i][j][i1][i2];
                   }
                   
                   //dSdX(i,j,:,j) = reshape(dSdX(i,j,:,j),1,D+2) - M(i)*dMdX(j,:,j);
                   double dSdX_i_j_j[] = new double[D + 2];
                   for(int j1 = 0 ; j1 < (D + 2); j1++)
                    dSdX_i_j_j[j1] = dSdX[i][j][j1][j];
                  
                   double dMdX_j_j[] = new double [D + 2],
                          dMdX_i_i[] = new double [D + 2];
                   for(int j1 = 0 ; j1 < (D + 2); j1++)
                   {
                    dMdX_j_j[j1] = dMdX[j][j1][j];  
                    dMdX_i_i[j1] = dMdX[i][j1][i];
                   }
                   
                   dSdX_i_j_j = vectorSUB(dSdX_i_j_j, vectorCoef(dMdX_j_j, M[i]));
                 
                   for(int i1 = 0 ; i1 < dSdX_i_j_j.length ; i1++)
                    dSdX[i][j][i1][j] = dSdX_i_j_j[i1];
                    
                   //dSdX(i,j,:,i) = reshape(dSdX(i,j,:,i),1,D+2) - M(j)*dMdX(i,:,i);
                   double dSdX_i_j_i[] = new double[D + 2];
                   for(int j1 = 0 ; j1 < (D + 2); j1++)
                    dSdX_i_j_i[j1] = dSdX[i][j][j1][i];
                   
                   dSdX_i_j_i = vectorSUB(dSdX_i_j_i, vectorCoef(dMdX_i_i, M[j]));
                   
                   for(int i1 = 0 ; i1 < dSdX_i_j_i.length ; i1++)
                    dSdX[i][j][i1][i] = dSdX_i_j_i[i1];
                   
                 }
                    
                   //dSdi(i,j,:,:) = Z*t - reshape(M(i)*dMdi(j,:,:) + dMdi(i,:,:)*M(j),n,D);
                   //dSdi(j,i,:,:) = dSdi(i,j,:,:);
                   for(int i1 = 0 ; i1 < Z.length ; i1++)
                   for(int j1 = 0 ; j1 < Z[0].length ; j1++)
                   {
                     
                     dSdi[i][j][i1][j1] = Z[i1][j1] * t_d - dMdi[j][i1][j1] * M[i] - dMdi[i][i1][j1] * M[j];
                     dSdi[j][i][i1][j1] = dSdi[i][j][i1][j1];
                   
                   }
            
                   //dSdX(j,i,:,:) = dSdX(i,j,:,:);
                   for(int i1 = 0 ; i1 < (D + 2) ; i1++)
                   for(int j1 = 0 ; j1 < E ; j1++)
                   {
                     
                     dSdX[j][i][i1][j1] = dSdX[i][j][i1][j1];
                     
                   }
          
          
    }
    
    //S(i,i) = S(i,i) + 1e-6;
    S[i][i] += 1.0e-6;
      
  }
  
  //dSdX(:,:,D+1,:) = -dSdX(:,:,D+2,:);
  for(int i1 = 0 ; i1 < E ; i1++)
  for(int j1 = 0 ; j1 < E ; j1++)
  for(int k1 = 0 ; k1 < E ; k1++)
  {  
    dSdX[i1][j1][D][k1] = -dSdX[i1][j1][D + 1][k1];
  }
  
  // 4) centralize moments
  //S = S - M*M';   
  
  S = matrixSUB(S, Outer_Product(M, M));
  
  // 5) Vectorize derivatives
  //dMds=reshape(dMds,[E D*D]);
  ret.dMds = reshape(dMds, E, D * D);
  
  //dSdm=reshape(dSdm,[E*E D]);
  ret.dSdm = reshape(dSdm, E * E, D);
  
  //dSds=reshape(dSds,[E*E D*D]);
  ret.dSds = reshape(dSds, E * E, D * D);
  
  //dVdm=reshape(dVdm,[D*E D]);
  ret.dVdm = reshape(dVdm, D * E, D);
  
  //dVds=reshape(dVds,[D*E D*D]);
  ret.dVds = reshape(dVds, D * E, D * D);
  
  //dMdi=reshape(dMdi,E,[]);
  ret.dMdi = reshape(dMdi, E, n * D);
  
  //dMdt=reshape(dMdt,E,[]);
  ret.dMdt = reshape(dMdt, E, n * E);
  
  //dMdX=reshape(dMdX,E,[]);
  ret.dMdX = reshape(dMdX, E, (D + 2) * E);
  
  //dSdi=reshape(dSdi,E*E,[]);
  ret.dSdi = reshape(dSdi, E * E, n * D);
  
  //dSdt=reshape(dSdt,E*E,[]);
  ret.dSdt = reshape(dSdt, E * E, n * E);
  
  //dSdX=reshape(dSdX,E*E,[]);
  ret.dSdX = reshape(dSdX, E * E, (D + 2) * E);
  
  //dVdi=reshape(dVdi,D*E,[]);
  ret.dVdi = reshape(dVdi, D * E, n * D);
  
  //dVdt=reshape(dVdt,D*E,[]);
  ret.dVdt = reshape(dVdt, D * E, n * E);
  
  //dVdX=reshape(dVdX,D*E,[]);
  ret.dVdX = reshape(dVdX, D * E, (D + 2) * E);
  
  ret.dMdm = dMdm;
  ret.M = M;//vectorCoef(M, 1.0);
  ret.V = V;//coeffProdMat(1.0, V);
  ret.S = S;//coeffProdMat(1.0, S);
  
  return ret;
  
}
