
////***********************************************************
////***********************************************************
////***********************************************************
//double sin_db(double x)
//double cos_db(double x)
//double log_db(double x)
//double exp_db(double x)
//double sqrt_db(double x)
//double pow_db(double x, double y)
//double pow_db(double x, int y)
//double abs_db(double x)
//double max_db(double x, double y)
//double min_db(double x, double y)
//double constrain_db(double x, double a, double b)
//double[] append_db(double A[], double newEl)

//void printMat(double matrixA[][])
//void printMat(double matrixA[][], int i)
//void printMat(int matrixA[][])
//void printVecT(double matrixA[])
//void printVec(double matrixA[])
//void printVecTst(double matrixA[])
//void printVec(int matrixA[])
//void printMass3D(double matrixA[][][])
//void printMass4D(double matrixA[][][][])

//double[][] transMat(double matrixA[][])
//int[][] transMat(int matrixA[][])
////******************************************************************
////          Произведение матриц A строки Rows столбцы Cols и B строки Rows столбцы bCol
////
//double[][] matrixMultiply(double matrixA[][], double matrixB[][])
//double[] matrixMultiplyC(double matrixA[][], double matrixB[])
////*************************************************************
////          Произведение транспонировнный вектор * матрица (yT * K)
////          Результат транспонированный вектор (вектор строка)
//double[] matrixMultiplyC_T2(double vectorB[], double matrixA[][])// Эта версия правильная проверить по всему тексту
////*************************************************************
////          Произведение вектор * вектор (x * y)
////          Результат скаляр
//double matrixMultiplyC_V(double vectorA[], double vectorB[])

//double [][] Vec_to_Mat(double A[][], double b[], int i) // Вектор в столец матрицы копирем
//double [] Mat_to_Vec(double A[][], int d)// Вырезаем из матрицы столбец d и в вектор
//double [][][] Mat_to_3Dmass(double A[][][], double mat[][], int i) // Матрицу в i-ю позицию 3-го измерения массива копирем
//double [][] getMat_to_Mat(double A[][], int b, int e)// Вырезаем из матрицы столбецы от b до e и в вектор
//double [][] getMat_to_Mat(double A[][], int d)// Вырезаем из матрицы столбецы от b до e и в вектор
//double [] Mat_to_VecT(double A[][], int d) // Вырезаем строку d и в вектор

//double[][] eye(int n)
//double[] eyeV(int n)

////***************************************************************
////        Внешнее произведение векторов A * B_T
////        Результат матрица
//double[][] Outer_Product(double vectorA[], double vectorB[])

//double[] vectorSUM(double s1, double vectorA[], double s2, double vectorB[])
//double[] vectorADD(double vectorA[], double vectorB[])
//double[][] matrixADD(double matrixA[][], double matrixB[][])
//double[][] matrixSUB(double matrixA[][], double matrixB[][])
//double[] vectorSUB(double vectorA[], double vectorB[])
//double[][] matrixDOT(double matrixA[][], double matrixB[][]) // bsxfun(@times,t,t) t - матрица
//double[][] matrixDOT(double vectorB[], double matrixA[][])// bsxfun(@times, lb, t) t - матрица, lb - вектор
//double[][] matrixDOT(double matrixA[][], double vectorB[])// bsxfun(@times, lb, t) t - матрица, lb - транспонированный вектор
//double[][] matrixDOT_T(double matrixA[][], double vectorB[])// bsxfun(@times, t, lb) t - матрица, lb - транспонированный вектор
//double[] vectorDOT(double vectorA[], double vectorB[])// 
//double[][] matrixCoef(double matrixA[][], double cf)
//double[][] vectorCoef2(double vectorA[], double cf)
//double[] vectorCoef(double vectorA[], double cf)
//int[] vectorCoef(int vectorA[], int cf)
//double[] SUM(double matrixA[][])
//double SUM(double vectorA[])
//double[] SUM(double matrixA[][], int d)
//double[] SUM_d(double matrixA[][], int d)

//double sign(double x)
//float sign(float x)
        
//boolean Chol_S(double matA[][], int n)
//double [][] solve(double A[][], int n)
//double [][] getSch()
//double [] solve_hol(double B[], int n)
//double [][] sparse(double t[][], double a[][], int m, int n)

//double [][] coeffProdMat(double cof, double a[][])
//double [][] coeffProdMat(double cof, double a[])
//double [][][] coeffProd3DMass(double cof, double a[][][])

//double [][] kron(double A[][], double B[][])
//double [][] reshape(double A[][][], int m, int n) // m - новое количество строк, n - новое количество столбцов
//double [][] reshape(double A[][][][], int m, int n) // m - новое количество строк, n - новое количество столбцов

//double [][] Chainrule(double A[][], double B[][], int indB[], double C[][], double D[][], int indD[])
//double [][] Chainrule(double A[][], double B[][], double C[][], double D[][], int indD[])
//double [][] Chainrule(double A[], double B[][], double C[], double D[][])
//double [][] Chainrule(double A[][], double B[][], double C[][], double D[][])

//int Gauss(double matrixA[][], double vectorB[], double x[])
//double[][] solveMat(double matrixB[][], double matrixA[][]) // Оператор MatLab /  для одинаковых квадратных матриц
//double[] solveMatV(double vector[], double matrixA[][]) // Оператор MatLab /  для одинаковых квадратных матриц
//double detMatrix(double[][] matrix)

//double [][] diag(double v[])
//double [] diag(double m[][])

//double [] getStr(double A[][], int k)
//double [] getSlb(double A[][], int k)

//double [][] bsxfun_plus(double A[], double B[])
//double [][] bsxfun_minus(double A[], double B[])
//double [][] bsxfun_times(double A[], double B[])
//double [][] bsxfun_times(double A[], double B[][])

//double[][] solveMatM(double matrixA[][], double matrixB[][])
//double [] operInvSlash(double matrixA[][], double vectorB[])
//double [][] operInvSlash(double matrixA[][], double matrixB[][])

//double[] EXP_F(double A[])
//double[][] EXP_F(double A[][])
//double[][] COS_F(double A[][])
//double[][] SIN_F(double A[][])

//double[][] ForDig(double A[][])

//int [] createInd(int bg, int ed)
//double [] getVec(double M[], int ind[])
//double [][] getMat(double Ss[][], int indR[], int indC[])
//double [][] getMat(double S[][], int indR[])
//double [][] getMat(double S[][], int indB, int indE) // zi(:,1:d)
//double[] setVec(double RCv[], int ind[], double TRm[])
//double[][] setMat(double RCv[][], int indR[], int indC[], double TRm[][])
//double [][] setMat(double RCv[][], int indR[], double TRm[][])
//double [][][] set3DMass(double RCv[][][], int indR[], int indC[], int ind, double TRm[][])

//double[][] addTo(double RCv[][], int ind_i[], int ind_j[], int ind_k[], double TRm[][])
//double[][] addToSimpl(double RCv[][], int ind_i[], int ind_k[], double TRm[][])


//double [][] symmetrizeA(int D, double A[][]) // X = reshape(1:D3*D3,[D3 D3]); XT = X'; Sds = (Sds + Sds(XT(:),:))/2
//double [][] symmetrizeB(int D, double A[][]) // X = reshape(1:D0*D0,[D0 D0]); XT = X'; Sds = (Sds + Sds(:,XT(:)))/2

////  Новая версия оператора \
//double[][] matrix_minor(double A[][], int d)
//double [] vmadd(double a[], double b[], double s, int n) /* c = a + b * s */
//double [][] vmul(double v[], int n) /* m = I - v v^T */
//double vnorm(double A[], int n) /* ||x|| */
//double[] vdiv(double A[], double d, int n) /* y = x / d */
//double[] mcol(double A[][], int c) /* take c-th column of m, put in v */
//double [][] matrix_mul(double x[][], double y[][])
//double [][] matrix_transpose(double A[][])
//double [][] QR_Solver(double A[][], double B[][])
////**********************************************
//// Вычисление решения Ax = B. Если A квадратная верхне триугольная матрица, B вектор
//double[] triMatSolve(double matrixA[][], double vectorB[])

double sin_db(double x)
{
  if(expr_sin == null)
  {
   exprStr = "sin(x)";
   expr_sin = Compile.expression(exprStr, false);
  }
  return(expr_sin.eval(x).answer().toDouble());
}

double cos_db(double x)
{
  if(expr_cos == null)
  {
  exprStr = "cos(x)";
  expr_cos = Compile.expression(exprStr, false);
  }
  return(expr_cos.eval(x).answer().toDouble());
}

double log_db(double x)
{
  if(expr_log == null)
  {
  exprStr = "log(x)";
  expr_log = Compile.expression(exprStr, false);
  }
  return(expr_log.eval(x).answer().toDouble());
}

double exp_db(double x)
{
  if(expr_exp == null)
  {
  exprStr = "exp(x)";
  expr_exp = Compile.expression(exprStr, false);
  }
  return(expr_exp.eval(x).answer().toDouble());
  
}

double sqrt_db(double x)
{ 
  if(expr_sqrt==null)
  {
  exprStr = "sqrt(x)";
  expr_sqrt = Compile.expression(exprStr, false);
  }
  return(expr_sqrt.eval(x).answer().toDouble());
}

double pow_db(double x, double y)
{
  if(expr_pow == null)
  {
  exprStr = "x^y";
  expr_pow = Compile.expression(exprStr, false);
  }
  return(expr_pow.eval(x, y).answer().toDouble());
}

double pow_db(double x, int y)
{
  if(expr_pow == null)
  {
  exprStr = "x^y";
  expr_pow = Compile.expression(exprStr, false);
  }
  return(expr_pow.eval(x, y).answer().toDouble());
}

double abs_db(double x)
{
  if(x >= 0) return x;
  else
             return -x;
}

double max_db(double x, double y)
{
  if(x >= y) return x;
  else return y;
}

double min_db(double x, double y)
{
  if(x <= y) return x;
  else return y;
}

double constrain_db(double x, double a, double b)
{
 if(x <= a) return a;
 if(x >= b) return b;
 else
   return x;
}

double[] append_db(double A[], double newEl)
{
  
  double ret[] = new double [A.length + 1];
  for(int i = 0; i < A.length ; i++)
   ret[i] = A[i];
  ret[A.length] = newEl;
  
  return ret;
  
}

double[][] append_db(double A[][]) // Увеличение матрицы на одну строку
{
  
  double ret[][] = new double [A.length + 1][A[0].length];
  for(int i = 0; i < A.length ; i++)
   ret[i] = A[i];
    
  return ret;
  
}

void printMatFor(double matrixA[][])
{
  int     n    = matrixA.length,
          m    = matrixA[0].length,
          i    = 0,
          j    = 0;  
  println();   
  println("********************");
  
  println("{");
  for(i = 0; i < n; i++)
  {
    print(" {");
    for(j = 0; j < (m - 1); j++)
    { 
     
      print(matrixA[i][j] + ", ");
    
    }
    println(matrixA[i][m - 1] + "},");
    
  }
  println("}");
  println("++++++++++++++++++++");
}

void printMat(double matrixA[][])
{
  int     n    = matrixA.length,
          m    = matrixA[0].length,
          i    = 0,
          j    = 0;  
  println();   
  println("********************");
  println("Матрица размером: строк " + n + ", столбцов " + m);
  println("--------------------");
  for(i = 0; i < n; i++)
  {
    print(" строка =  " + (i + 1));
    for(j = 0; j < m; j++)
    { 
     
      print("    " + matrixA[i][j]);
    
    }
    println(";");
  }
  println("++++++++++++++++++++");
}

void printMat(double matrixA[][], int i)
{
  int     n    = matrixA.length,
          m    = matrixA[0].length,
          j    = 0;  
      
  println("********************");
  println("Матрица размером: строк " + n + ", столбцов " + m);
  
    for(j = 0; j < m; j++)
    { 
     
      println(j + 1 + " " + matrixA[i - 1][j]);
    
    }
    println(";");
  
  println("++++++++++++++++++++");
}

void printMat_f(double matrixA[][])
{
  int     n    = matrixA.length,
          m    = matrixA[0].length,
          i    = 0,
          j    = 0;  
  println();   
  println("********************");
  println("Матрица размером: строк " + n + ", столбцов " + m);
  println("--------------------");
  for(i = 0; i < n; i++)
  {
    print(" строка =  " + (i + 1));
    for(j = 0; j < m; j++)
    { 
     
      print("    " + (float)matrixA[i][j]);
    
    }
    println(";");
  }
  println("++++++++++++++++++++");
}

void printMat(int matrixA[][])
{
  int     n    = matrixA.length,
          m    = matrixA[0].length,
          i    = 0,
          j    = 0;  
      
  println("********************");
  for(i = 0; i < n; i++)
  {
    for(j = 0; j < m; j++)
    { 
     
      print(" " + matrixA[i][j]);
    
    }
    println(";");
  }
  println("++++++++++++++++++++");
}

void printVecT(double matrixA[])
{
  int  i  = 0,
       n  = matrixA.length;
  
  println("********************");
  for(i = 0; i < n; i++)
  {
      print((float)matrixA[i] + "  ");
      
  }
  println("++++++++++++++++++++");
}

void printVecFor(double matrixA[])
{
  int  i  = 0,
       n  = matrixA.length;
  
  println("********************");
  println(" Вектор длиной " + n + " элементов");
  println("{");
  for(i = 0; i < n; i++)
  {
      println(matrixA[i] + " , ");
      
  }
  println("}");
  
  println("++++++++++++++++++++");
}

void printVec(double matrixA[])
{
  int  i  = 0,
       n  = matrixA.length;
  
  println("********************");
  println(" Вектор длиной " + n + " элементов");
  for(i = 0; i < n; i++)
  {
      println(" i = " + i + " : " + matrixA[i] + "  ");
      
  }
  println("++++++++++++++++++++");
}

void printVecTst(double matrixA[])
{
  int  i  = 0,
       n  = matrixA.length;
       
  double B[] = {
                10.8671937687587,
                -17.2055331303966,
                20.8951668019699,
                8.02185468500911,
                -44.4728822921037,
                27.4555433157264,
                -23.8044928435181,
                -16.4640183538914,
                -0.101455457934118,
                -14.2101514920029,
                10.1153826802789,
                -53.5486716782853,
                40.8407556636396,
                40.9979815972716,
                -36.4798909242450,
                32.0156771050565,
                -46.7057222723958,
                62.9161963648390,
                -56.1846505271066,
                14.0249894076769,
                6.54785400113180,
                41.9576835422792,
                9.57134754926921,
                2.41517030955686,
                -28.5361931457791,
                24.8002043015007,
                -9.31923383314560,
                25.2772534106467,
                -2.36278591561255,
                -14.8332588945852,
                20.6122709010438,
                21.9267282372268,
                -12.8975204492073,
                -16.2974424360494,
                -16.7237801402406,
                29.8339347768603,
                -15.3159399173969,
                -16.0192149882764,
                7.17398336569997,
                -16.7961848331781
                };
  println();
  println("============================");
  println("Вектор длиной " + matrixA.length  + " элементов");
  println("============================");
  for(i = 0; i < (n - 1); i++)
  {
      println(B[i] - matrixA[i] + ",");
      
  }
  println(B[n - 1] - matrixA[n - 1] + ";");
  println("*******************");
}

void printVec(int matrixA[])
{
  int  i  = 0,
       n  = matrixA.length;
  
  println("********************");
  for(i = 0; i < n; i++)
  {
      print(matrixA[i] + "  ");
      
  }
  println("++++++++++++++++++++");
}

void printMass3D(double matrixA[][][])
{
  int   sz1 = matrixA.length,                         // 1 е измерение,
        sz2 = matrixA[0].length,                      // 2 е измерение
        sz3 = matrixA[0][0].length,                   // 3 е измерение
        i   = 0,
        j   = 0,
        k   = 0;
      
  println("********* 3D Mass ***********");
  
  for(k = 0; k < sz3 ; k++)
  {
    println("Mass[][][" + k + "]");
    
    for(i = 0; i < sz1 ; i++)
    {
      for(j = 0; j < sz2; j++)
      { 
       
        print(" " + (float)matrixA[i][j][k]);
      
      }
    println(";");
    }
  }
  println("++++++++++++++++++++");
}

void printMass4D(double matrixA[][][][])
{
  int   sz1 = matrixA.length,                         // 1 е измерение,
        sz2 = matrixA[0].length,                      // 2 е измерение
        sz3 = matrixA[0][0].length,                   // 3 е измерение
        sz4 = matrixA[0][0][0].length,
        i   = 0,
        j   = 0,
        k   = 0,
        z   = 0;
      
  println("********* 4D Mass ***********");
  
  for(z = 0; z < sz4 ; z++)
  for(k = 0; k < sz3 ; k++)
  {
    println("Mass[][][" + (k + 1) + "][" + (z + 1) +"]");
    
    for(i = 0; i < sz1 ; i++)
    {
      for(j = 0; j < sz2; j++)
      { 
       
        print(" " + (float)matrixA[i][j][k][z]);
      
      }
    println(";");
    }
  }
  println("++++++++++++++++++++");
}

double[][] transMat(double matrixA[][])
{
  
  int        aRows    = matrixA.length,
             aCols    = matrixA[0].length;
             
  double [][]   temp_mat = new double [aCols][aRows];
 
  int        i, 
             j;
    
    for(i = 0; i < aRows ; i++)
    for(j = 0; j < aCols ; j++)
      temp_mat[j][i] = matrixA[i][j];
      
  return temp_mat;
  
}

sparse_t transMat_SP(sparse_t matrixA)
{
  
  int        aRows    = matrixA.n,
             aCols    = matrixA.m,
             elems    = matrixA.elements;
             
  sparse_t temp_mat = new sparse_t(aCols, aRows, elems);
  
  for(int i = 0; i < elems ; i++)
  {
    temp_mat.i[i] = matrixA.j[i];
    temp_mat.j[i] = matrixA.i[i];
    temp_mat.val[i] = matrixA.val[i];
  }   
  return temp_mat;
   
}

int[][] transMat(int matrixA[][])
{
  
  int        aRows    = matrixA.length,
             aCols    = matrixA[0].length;
             
  int [][]   temp_mat = new int [aCols][aRows];
 
  int        i, 
             j;
    
    for(i = 0; i < aRows ; i++)
    for(j = 0; j < aCols ; j++)
      temp_mat[j][i] = matrixA[i][j];
      
  return temp_mat;
  
}


//******************************************************************
// Произведение матриц A строки Rows столбцы Cols и B строки Rows столбцы bCol
//

double[][] matrixMultiply(double matrixA[][], double matrixB[][])
{
  int          i;
  int          j;
  int          k;
  
  int          aRows = matrixA.length, 
               Cols  = matrixA[0].length, 
               bCols = matrixB[0].length;
  
  double [][] matrixC = new double [aRows][bCols];
  
  for (i = 0; i < aRows; i++)
  {
    
    for (j = 0; j < bCols; j++)
      matrixC[i][j] = 0.0;
      
  }

  for (i = 0; i < aRows; i++)
  {
    for(j = 0; j < bCols; j++)
    {
      for(k = 0;  k < Cols; k++)
      {
       matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
      }
    }
  }
  return matrixC;
}

double[][] matrixMultiply_SP(sparse_t matrixA, double matrixB[][])
{
  
  int          aRows = matrixA.n, 
               bCols = matrixB[0].length;
  
  double [][] matrixC = new double [aRows][bCols];
  for(int j = 0; j < bCols; j++)
    for(int i = 0; i < matrixA.elements; i++)
    {
        matrixC[matrixA.i[i]][j] += matrixB[matrixA.j[i]][j] * matrixA.val[i];
    }
  
  return matrixC;
}

double[][] matrixMultiply_SP(double matrixA[][], sparse_t matrixB)
{
  
  int          aRows = matrixA.length, 
               bCols = matrixB.m;
  
  double [][] matrixC = new double [aRows][bCols];
  for(int j = 0; j < aRows; j++)
    for(int i = 0; i < matrixB.elements; i++)
    {
        matrixC[j][matrixB.j[i]] += matrixA[j][matrixB.i[i]] * matrixB.val[i];
    }
  
  return matrixC;
}

double[] matrixMultiplyC(double matrixA[][], double matrixB[])
{
  int          aRows = matrixA.length,
               aCols = matrixA[0].length;
  int          i;
  int          j;
  int          k;
  double [] matrixC = new double [aRows];
  
  for (i = 0; i < aRows; i++)
  {
    matrixC[i] = 0.0;
  }
 
  for(j = 0; j < aRows; j++)
  {
   for(k = 0;  k < aCols; k++)
   {
    matrixC[j] += matrixB[k] * matrixA[j][k];
   }
  }
  
  return matrixC;
}

double[] matrixMultiplyC_SP(sparse_t matrixA, double vectorB[])
{
  int          aRows = matrixA.n;
  double [] vectorC = new double [aRows];
 
  for(int i = 0; i < matrixA.elements; i++)
  {
      vectorC[matrixA.i[i]] += vectorB[matrixA.j[i]] * matrixA.val[i];
  }
  
  return vectorC;
}

//*************************************************************
//          Произведение транспонировнный вектор * матрица (yT * K)
//          Результат транспонированный вектор (вектор строка)


double[] matrixMultiplyC_T2(double vectorB[], double matrixA[][])// Эта версия правильная проверить по всему тексту
{
  int          aRows = matrixA.length, 
               aCols = matrixA[0].length;
  int          i;
  int          j;
  int          k;
  double []     vectorC = new double [aCols];
  
  for (i = 0; i < aCols; i++)
  {
    vectorC[i] = 0.0;
  }

  for(j = 0; j < aCols; j++)
  {
   for(k = 0;  k < aRows; k++)
   {
    vectorC[j] += vectorB[k] * matrixA[k][j];
   }
  }
  
  return vectorC;
}

//*************************************************************
//          Произведение вектор * вектор (x * y)
//          Результат скаляр

double matrixMultiplyC_V(double vectorA[], double vectorB[])
{
  int          i            = 0;
  double       result       = 0;
  int          aRows        = vectorA.length;   
  for(i = 0; i < aRows; i++)
  {
   result += vectorA[i] * vectorB[i];
  }
  
  return result;
}

double [][][] Mat_to_3Dmass(double A[][][], double mat[][], int ind) 
{
  double ret [][][] = new double [A.length][A[0].length][A[0][0].length];
  
  int  i,
       j,
       k;
  for(i = 0; i < A.length ; i++)
    for(j = 0; j < A[0].length ; j++)
      for(k = 0; k < A[0][0].length ; k++)
        ret[i][j][k] = A[i][j][k];
   
   for(i = 0; i < A.length ; i++)
    for(j = 0; j < A[0].length ; j++)
       ret[i][j][ind] = mat[i][j];
 
   return ret;
  
}

double [][] Vec_to_Mat(double A[][], double b[], int i) // Вектор в столец матрицы копирем
{
 
  double ret [][] = new double [A.length][A[0].length];
  
  ret = coeffProdMat(1.0, A);
  
  for(int j = 0 ; j < b.length ; j++)
  {
    ret[j][i] = b[j];
  }
   return ret;
}

double [] Mat_to_Vec(double A[][], int d)// Вырезаем из матрицы столбец d и в вектор
{
  int n = A.length,
      i = 0;
  double ret [] = new double [n];
  
  for(i = 0; i < n ;i++)
   ret[i] = A[i][d];
  
  return ret;
  
}

double [][] getMat_to_Mat(double A[][], int b, int e)// Вырезаем из матрицы столбецы от b до e и в вектор
{
  int n = A.length,
      i,
      j;
  double ret [][] = new double [n][e - b + 1];
  
  for(i = 0; i < n ; i++)
  for(j = b; j <= e ; j++)
   ret[i][j] = A[i][j];
  
  return ret;
  
}

double [][] getMat_to_Mat(double A[][], int d)// Вырезаем из матрицы столбецы от b до e и в вектор
                                            // d 1 и больше !!!
{
  int n = A.length,
      i,
      j;
  double ret [][] = new double [n][1];
  
  for(i = 0; i < n ; i++)
    ret[i][0] = A[i][d];
  
  return ret;
  
}

double [] Mat_to_VecT(double A[][], int d) // Вырезаем строку d и в вектор
{
  int n = A[0].length,
      i = 0;
  double ret [] = new double [n];
  
  for(i = 0; i < n ;i++)
   ret[i] = A[d][i];
  
  return ret;
  
}

double [][] getMatFromMass3D(double matrixA[][][], int k)
{
  int   sz1 = matrixA.length,                         // 1 е измерение,
        sz2 = matrixA[0].length,                      // 2 е измерение
        sz3 = matrixA[0][0].length,                   // 3 е измерение
        i   = 0,
        j   = 0;
      
 double ret[][] = new double[sz1][sz2];
  
    for(i = 0; i < sz1 ; i++)
    {
      for(j = 0; j < sz2; j++)
      { 
       
        ret[i][j] = matrixA[i][j][k];
      
      }
    
    }

  return ret;
}

double [][][] coeffProd3DMass(double cof, double matrixA[][][])
{
  int   sz1 = matrixA.length,                         // 1 е измерение,
        sz2 = matrixA[0].length,                      // 2 е измерение
        sz3 = matrixA[0][0].length,                   // 3 е измерение
        i   = 0,
        j   = 0,
        k   = 0;
      
 double ret[][][] = new double[sz1][sz2][sz3];
    for(k = 0; k < sz3 ; k++)
    {
      for(i = 0; i < sz1 ; i++)
      {
        for(j = 0; j < sz2; j++)
        { 
         
          ret[i][j][k] = matrixA[i][j][k];
        
        }
      }
    }

  return ret;
}

double[][] eye(int n)
{
 
  int          i;
  int          j;
  
  double [][] matrixC = new double [n][n];
  
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
      matrixC[i][j] = 0.0f;
    matrixC[i][i] = 1.0f;
  }
    
  return matrixC;
  
}

double[] eyeV(int n)
{
 
  int          i;
  int          j;
  
  double [] vector = new double [n];
  
  for (i = 0; i < n; i++)
  {
    vector[i] = 1.0f;
  }
    
  return vector;
  
}
//***************************************************************
//    Внешнее произведение векторов A * B_T
//    Результат матрица

double[][] Outer_Product(double vectorA[], double vectorB[])
{
 
  int          aCols = vectorA.length,
               bCols = vectorB.length;
  int          i;
  int          j;
  
  double [][] matrixC = new double [aCols][bCols];
  
  for (i = 0; i < aCols; i++)
  {
    for(j = 0; j < bCols; j++)
    {
 
       matrixC[i][j] = vectorA[i] * vectorB[j];
      
    }
  }
  
  return matrixC;
  
}


double[] vectorSUM(double s1, double vectorA[], double s2, double vectorB[])
{
  int          i,
               n = vectorA.length;
  
  double [] vectorC = new double [n];
  
  
  for (i = 0; i < n; i++)
  {
   
   vectorC[i] = s1 * vectorA[i] + s2 * vectorB[i];

  }
  return vectorC;
}

double[] vectorADD(double vectorA[], double vectorB[])
{
  int          i,
               n = vectorA.length;
  
  double [] vectorC = new double [n];
  
  
  for (i = 0; i < n; i++)
  {
   
   vectorC[i] = vectorA[i] + vectorB[i];

  }
  return vectorC;
}

double[][] matrixADD(double matrixA[][], double matrixB[][])
{
  int          i;
  int          j;
  int          aRows = matrixA.length,
               aCols = matrixA[0].length;
  
  double [][] matrixC = new double [aRows][aCols];
  
  
  for (i = 0; i < aRows; i++)
  {
    for(j = 0; j < aCols; j++)
    {
     
       matrixC[i][j] = matrixA[i][j] + matrixB[i][j];
     
    }
  }
  return matrixC;
}

double[][] matrixSUB(double matrixA[][], double matrixB[][])
{
  int          i;
  int          j;
  int          aRows = matrixA.length,
               aCols = matrixA[0].length;
  
  double [][] matrixC = new double [aRows][aCols];
  
  
  for (i = 0; i < aRows; i++)
  {
    for(j = 0; j < aCols; j++)
    {
     
       matrixC[i][j] = matrixA[i][j] - matrixB[i][j];
     
    }
  }
  return matrixC;
}

double[] vectorSUB(double vectorA[], double vectorB[])
{
  int          i,
               n = vectorA.length;
  
  double [] vectorC = new double [n];
  
  
  for (i = 0; i < n; i++)
  {
   
   vectorC[i] = vectorA[i] - vectorB[i];

  }
  return vectorC;
}

double[][] matrixDOT(double matrixA[][], double matrixB[][]) // bsxfun(@times,t,t) t - матрица
{
  int          aRows = matrixA.length, 
               aCols = matrixA[0].length;
  int          i;
  int          j;

  double [][] matrixC = new double [aRows][aCols];
  
  
  for (i = 0; i < aRows; i++)
  {
    for(j = 0; j < aCols; j++)
    {
     
       matrixC[i][j] = matrixA[i][j] * matrixB[i][j];
     
    }
  }
  return matrixC;
}

double[][] matrixDOT(double vectorB[], double matrixA[][])// bsxfun(@times, lb, t) t - матрица, lb - вектор
{
  int          aRows = matrixA.length, 
               aCols = matrixA[0].length;
  int          i;
  int          j;

  double [][] matrixC = new double [aRows][aCols];
  
  
  for (i = 0; i < aRows; i++)
  {
    for(j = 0; j < aCols; j++)
    {
     
       matrixC[i][j] = matrixA[i][j] * vectorB[i];
     
    }
  }
  return matrixC;
}

double[][] matrixDOT(double matrixA[][], double vectorB[])// bsxfun(@times, lb, t) t - матрица, lb - транспонированный вектор
{
  int          aRows = matrixA.length, 
               aCols = matrixA[0].length;
  int          i;
  int          j;

  double [][] matrixC = new double [aRows][aCols];
  
  
  for (i = 0; i < aRows; i++)
  {
    for(j = 0; j < aCols; j++)
    {
     
       matrixC[i][j] = matrixA[i][j] * vectorB[j];
     
    }
  }
  return matrixC;
}

double[][] matrixDOT_T(double matrixA[][], double vectorB[])// bsxfun(@times, t, lb) t - матрица, lb - транспонированный вектор
{
  int          aRows = matrixA.length, 
               aCols = matrixA[0].length;
  int          i;
  int          j;

  double [][] matrixC = new double [aRows][aCols];
  
  
  for (i = 0; i < aRows; i++)
  {
    for(j = 0; j < aCols; j++)
    {
     
       matrixC[i][j] = matrixA[i][j] * vectorB[i];
     
    }
  }
  return matrixC;
}

double[] vectorDOT(double vectorA[], double vectorB[])// 
{
  int          aRows = vectorA.length;
  int          i;


  double [] vectorC = new double [aRows];
  
  
  for (i = 0; i < aRows; i++)
  {
   
    vectorC[i] = vectorA[i] * vectorB[i];

  }
  return vectorC;
}

double[][] matrixCoef(double matrixA[][], double cf)
{
  int          aRows = matrixA.length,
               aCols = matrixA[0].length,
               i,
               j;

  double [][] matrixC = new double [aRows][aCols];
  
  
  for (i = 0; i < aRows; i++)
  {
    for(j = 0; j < aCols; j++)
    {
     
       matrixC[i][j] = cf * matrixA[i][j];
     
    }
  }
  return matrixC;
}

double[][] vectorCoef2(double vectorA[], double cf)
{
  int          aRows = vectorA.length,
               i;

  double [][] vectorC = new double [aRows][1];
  
  
  for (i = 0; i < aRows; i++)
  {
    vectorC[i][0] = cf * vectorA[i];
  }
  return vectorC;
}

double[] vectorCoef(double vectorA[], double cf)
{
  int          aRows = vectorA.length,
               i;

  double [] vectorC = new double [aRows];
  
  
  for (i = 0; i < aRows; i++)
  {
    vectorC[i] = cf * vectorA[i];
  }
  return vectorC;
}

int[] vectorCoef(int vectorA[], int cf)
{
  int          aRows = vectorA.length,
               i;

  int [] vectorC = new int [aRows];
  
  
  for (i = 0; i < aRows; i++)
  {
    vectorC[i] = cf * vectorA[i];
  }
  return vectorC;
}

double[] SUM(double matrixA[][])
{
  
  int          aRows = matrixA.length, 
               aCols = matrixA[0].length;
  int          i,
               j;

  double ret[] = new double [aCols];
    
    for(j = 0; j < aCols ; j++)
    {
      
      ret[j] = 0;
      
      for(i = 0; i < aRows; i++)
      {
       
         ret[j] += matrixA[i][j];
       
      }
    }
  
  return ret;
  
}

double SUM(double vectorA[])
{
  
  int          aRows = vectorA.length;
  int          i;

  double ret;
    
      ret = 0;
      
      for(i = 0; i < aRows; i++)
      {
       
         ret += vectorA[i];
       
      }
  
  return ret;
  
}

double[] SUM(double matrixA[][], int d)
{
  
  int          aRows = matrixA.length,   // строки
               aCols = matrixA[0].length,// столбцы
               i,
               j;

  double ret[];
    
    if(d == 1)
    {
      ret = new double [aCols];
      for(j = 0; j < aCols ; j++)
      {
        
        ret[j] = 0;
        
        for(i = 0; i < aRows; i++)
        {
         
           ret[j] += matrixA[i][j];
         
        }
      }
    }
    else
   
    {
     ret = new double [aRows];
      for(j = 0; j < aRows ; j++)
      {
        
        ret[j] = 0;
        
        for(i = 0; i < aCols; i++)
        {
         
           ret[j] += matrixA[j][i];
         
        }
      }
    }
  
  return ret;
}

double[] SUM_d(double matrixA[][], int d)
{
  
  int          aRows = matrixA.length,   // строки
               aCols = matrixA[0].length,// столбцы
               i,
               j;

  double ret[];
    
    if(d == 1)
    {
      ret = new double [aCols];
      for(j = 0; j < aCols ; j++)
      {
        
        ret[j] = 0;
        
        for(i = 0; i < aRows; i++)
        {
         
           ret[j] += matrixA[i][j];
         
        }
      }
    }
    else
   
    {
     ret = new double [aRows];
      for(j = 0; j < aRows ; j++)
      {
        
        ret[j] = 0;
        
        for(i = 0; i < aCols; i++)
        {
         
           ret[j] += matrixA[j][i];
         
        }
      }
    }
  
  return ret;
}
//******************************************************************

double sign(double x)
{
    if (x == 0.0) return 0;
    if (x > 0.0) return 1;
    else return -1;
}

float sign(float x)
{
    if (x == 0.0) return 0;
    if (x > 0.0) return 1;
    else return -1;
}

double  Sch[][],
        Dch[][];
        
boolean Chol_S(double matA[][], int n)
{
  
 double A[][]    = new double[n + 1][n + 1];
 Sch             = new double[n + 1][n + 1];
 Dch             = new double[n + 1][n + 1];
        
 for(int i = 0; i < n; i++)
 for(int j = 0; j < n; j++)
   A[i + 1][j + 1] = matA[i][j];
 
 /*------------------------Инициализация первых элементов матриц S и D------------------------*/
 Sch[1][1] = sqrt_db(abs_db(A[1][1]));
 Dch[1][1] = sign(A[1][1]);
 for(int j = 2; j <= n; ++j) 
 {
   Sch[1][j] = A[1][j] / (Sch[1][1] * Dch[1][1]);
 } 
 /*------------------------Инициализация первых элементов матриц S и D------------------------*/
 
 /*------------------------Расчитывамем все оставшиеся значения S и D------------------------*/
    for(int i = 2; i <= n; ++i) {
        double sum = 0;
        for (int l = 1; l <= i - 1; ++l)
            sum = sum + Sch[l][i] * Sch[l][i] * Dch[l][l];
        Sch[i][i] = sqrt_db(abs_db(A[i][i] - sum));
        Dch[i][i] = sign(A[i][i] - sum);
        for (int j = i + 1; j <= n; ++j) {
            sum = 0;
            for (int l = 1; l <= i - 1; l++) {
                sum = sum + Sch[l][i] * Sch[l][j] * Dch[l][l];
            }
            Sch[i][j] = (A[i][j] - sum) / (Sch[i][i]* Dch[i][i]);
        }
    }
    
/*------------------------Расчитывамем все оставшиеся значения S и D------------------------*/

/*Выводим матрицу после прямого хода, чтобы проверить, что она была приведена
     к ступенчатому виду
  Так как счет идет с ошибкой вычисления, то для наглядности занулим то,
   что находится ниже главной диагонали, а также проверим вектор правой части на наличие
    цифр с ошибкой вычисления                                                       */

    //for (int i = 1; i <= n; i++)
    //    for (int j = 1; j <= n; j++)
    //        if (abs_db(Sch[i][j]) < 0.000001) Sch[i][j] = 0;

    //println("Полученная матрица S:");
    //for (int i = 1; i <= n; i++)
    //{
    // for (int j = 1; j <= n; j++)
    //        print(Sch[j][i] + " ");
    // println();
    //}
   
    //for(int i = 0; i < n ; i++)
    //  for(int j = 0; j < n ; j++)
    //    U[j][i] = (float)S[i + 1][j + 1];
    
  return true;
}

double [][] solve(double A[][], int n)
{
 
  double ret[][] = new double[n][n];
  double y[]     = new double[n];
  
  for(int i = 0; i < n ; i++)
  {
    
    for(int j = 0; j < n; j++)
      y[j] = A[i][j];
      
    y = solve_hol(y, n);
    
    for(int j = 0; j < n; j++)
      ret[i][j] = y[j];
      
  }
  return ret;
  
}

double [][] getSch()
{
  int n = Sch[0].length - 1;
  double ret[][] = new double[n][n];
  
   for(int i = 0; i < n; i++)
      for(int j = 0; j < n ; j++)
        ret[j][i] = Sch[i + 1][j + 1];
    
  return ret;
  
}

double [] solve_hol(double B[], int n)
{
  
   double ret[] = new double[n];
   double y[]   = new double[n + 1],
          x[]   = new double[n + 1],
          f[]   = new double[n + 1];
   
   for(int i = 0; i < n; i++)
    f[i + 1] = B[i];
   
    y[1] = f[1] / Sch[1][1] * Dch[1][1];        //y[1] всегда равен f[1]/s[1,1]
    for(int i = 2; i <= n; ++i) {
        double sum = 0;
        for (int l = 1; l <= (i - 1); ++l)
            sum = sum + Sch[l][i] * y[l]* Dch[l][l];
        y[i] = (f[i] - sum) / (Sch[i][i]* Dch[i][i]);
    }
    
    x[n] = y[n] / Sch[n][n];           //Последний x[n] всегда равен u[n]/s[n,n]
    for(int i = (n - 1); i >= 1; --i) {
        double sum = 0;
        for (int l = (i + 1); l <= n; ++l)
            sum = sum + Sch[i][l] * x[l];
        x[i] = (y[i] - sum) / Sch[i][i];
    }
 
    for(int i = 0; i < n; i++)
    {
      ret[i] = x[i + 1];
    }
    
   return ret;
   
}

double [][] sparse(double t[][], double a[][], int m, int n)
{
  int      ma = a.length,
           na = a[0].length;
  
  for(int i = 0; i < ma ; i++)
  for(int j = 0; j < na ; j++)
    t[m + i][n + j] = a[i][j];
  return t;
}

double [][] coeffProdMat(double cof, double a[][])
{
  int      ma = a.length,
           na = a[0].length;
  double   r[][] = new double[ma][na];
  
  for(int i = 0; i < ma ; i++)
  for(int j = 0; j < na ; j++)
    r[i][j] = cof * a[i][j];
  return r;
}

double [][] coeffProdMat(double cof, double a[])
{
  int      ma = a.length,
           na = 1;
  double   r[][] = new double[na][ma];
  
  for(int i = 0; i < ma ; i++)
      r[0][i] = cof * a[i];
  return r;
}

double [][] kron(double A[], double B[])
{
  int      ma = A.length,
           mb = B.length;
  double   K[][]   = new double[1][ma * mb];
  int      i,
           j,
           l;
   l = 0;       
   for(i = 0; i < ma ; i++)
     {
       for(j = 0; j < ma ; j++)
        K[0][l++] = A[i] * B[j];
     }
       
  return K;
  
}

double [][] kron(double A[][], double B[][])
{
  int      ma = A.length,
           na = A[0].length,
           mb = B.length,
           nb = B[0].length;
  double   K[][]   = new double[ma * mb][na * nb],
           T[][]   = new double[mb][nb];
  int      i,
           j,
           k;
   for(i = 0; i < ma ; i++)
     for(j = 0; j < na ; j++)
     {
       T = coeffProdMat(A[i][j], B);
       K = sparse(K, T, i * mb, j * nb);
     }
       
  return K;
  
}

double [][] reshape(double A[][][], int m, int n) // m - новое количество строк, n - новое количество столбцов
{
  
  double B[][] = new double [m][n];
  
  int   sz1 = A.length,                         // 1 е измерение,
        sz2 = A[0].length,                      // 2 е измерение
        sz3 = A[0][0].length;                   // 3 е измерение
  
  
  double buf[] = new double [sz1 * sz2 * sz3];
  
  int   i,
        j,
        k,
        index;
   
   index = 0;
   
   for(k = 0 ; k < sz3 ; k++)
   for(j = 0 ; j < sz2 ; j++)
   for(i = 0 ; i < sz1 ; i++)
   {
    buf[index] = A[i][j][k];
    index++;
   }
   
   index = 0;
   for(i = 0 ; i < n ; i++)
   for(j = 0 ; j < m ; j++)
   {
    B[j][i] = buf[index];
    index++;
   }
  return B;
}

double [][] reshape(double A[][][][], int m, int n) // m - новое количество строк, n - новое количество столбцов
{
  
  double B[][] = new double [m][n];
  
  int   sz1 = A.length,                         // 1 е измерение,
        sz2 = A[0].length,                      // 2 е измерение
        sz3 = A[0][0].length,                   // 3 е измерение
        sz4 = A[0][0][0].length;
  
  
  double buf[] = new double [sz1 * sz2 * sz3 * sz4];
  
  int   i,
        j,
        k,
        l,
        index;
   
   index = 0;
   
   for(k = 0 ; k < sz4 ; k++)
   for(l = 0 ; l < sz3 ; l++)
   for(j = 0 ; j < sz2 ; j++)
   for(i = 0 ; i < sz1 ; i++)
   {
    buf[index] = A[i][j][l][k];
    index++;
   }
   
   index = 0;
   for(i = 0 ; i < n ; i++)
   for(j = 0 ; j < m ; j++)
   {
    B[j][i] = buf[index];
    index++;
   }
  return B;
}

double [][] Chainrule(double A[][], double B[][], int indB[], double C[][], double D[][], int indD[])
{
  
  int rowA = A.length,
      colA = A[0].length,
      rowB = indB.length,
      colB = B[0].length,
      rowC = C.length,
      colC = C[0].length,
      rowD = indD.length,
      colD = D[0].length;
      
  double BT[][] = new double [rowB][colB];
  double DT[][] = new double [rowD][colD];
  
  for(int i = 0; i < colB ; i++)
  for(int j = 0; j < rowB ; j++)
    BT[j][i] = B[indB[j]][i];
   
  for(int i = 0; i < colD ; i++)
  for(int j = 0; j < rowD ; j++)
    DT[j][i] = D[indD[j]][i];
    
  double R1[][] = new double [rowA][colB];
  double R2[][] = new double [rowA][colB];
  
  R1 = matrixMultiply(A, BT);
  R2 = matrixMultiply(C, DT);

  R1 =  matrixADD(R1, R2);

  return R1;
  
}

double [][] Chainrule(double A[][], double B[][], double C[][], double D[][], int indD[])
{
  
  int rowA = A.length,
      colA = A[0].length,
      colB = B[0].length,
      rowD = indD.length,
      colD = D[0].length;
      
  double DT[][] = new double [rowD][colD];
  
  for(int i = 0; i < colD ; i++)
  for(int j = 0; j < rowD ; j++)
    DT[j][i] = D[indD[j]][i];
    
  double R1[][] = new double [rowA][colB];
  double R2[][] = new double [rowA][colB];
  
  R1 = matrixMultiply(A, B);
  R2 = matrixMultiply(C, DT);

  R1 =  matrixADD(R1, R2);

  return R1;
  
}

double [][] Chainrule(double A[], double B[][], double C[], double D[][])
{
  
  double R1[][],
         R2[][],
         R3[][],
         A1[][] = new double[1][A.length],
         C1[][] = new double[1][C.length];
  
  A1[0] =A; //vectorCoef(A, 1.0);
  C1[0] = C;//vectorCoef(C, 1.0);
  R1 = matrixMultiply(A1, B);
  R2 = matrixMultiply(C1, D);

  R3 =  matrixADD(R1, R2);

  return R3;
  
}

double [][] Chainrule(double A[][], double B[][], double C[][], double D[][])
{
  
  double R1[][],
         R2[][],
         R3[][];
  
  R1 = matrixMultiply(A, B);
  R2 = matrixMultiply(C, D);

  R3 =  matrixADD(R1, R2);

  return R3;
  
}

//*********************************************************
//  Вычисление обратной матрицы с помощью оператора / matlab

int Gauss(double matrixA[][], double vectorB[], double x[])
{
  
  int     aRows = matrixA.length;
  int     i,
          j,
          k;
  double [][] matrixT = new double [aRows][aRows];
  double []   temp    = new double [aRows];
  
  for (i = 0; i < aRows; i++)
  {
    
    for (j = 0; j < aRows; j++)
      matrixT[i][j] = matrixA[i][j];
    
  }
    
  // Прямой ход
  for(k = 1; k < aRows; k++)
  {
  for(i = k; i < aRows; i++)
  {
    
    if(matrixT[k-1][k-1] == 0)
    {
      int i1 = k - 1;
      int j1;
      while((i1 < aRows)&&(matrixT[i1][k-1] == 0)) i1++;
      if(i1 >= aRows)return 0;
      
      for(j1 = 0; j1 < aRows ; j1++)
        temp[j1] = matrixT[i1][j1];
        
      for(j1 = 0; j1 < aRows ; j1++)
        matrixT[i1][j1] = matrixT[k - 1][j1];
        
      for(j1 = 0; j1 < aRows ; j1++)
        matrixT[k - 1][j1] = temp[j1];  
    }
    
    for(j = k; j < aRows; j++)
    matrixT[i][j] = matrixT[i][j] - matrixT[i][k-1] * matrixT[k-1][j] / (matrixT[k-1][k-1]);
    vectorB[i] = vectorB[i] - vectorB[k-1] * matrixT[i][k-1] / (matrixT[k-1][k-1]);
    
    matrixT[i][k-1] = 0;
  }
  }
  // Обратный ход
  
  for(k = (aRows-1); k >= 0; k--)
  {
    x[k] = vectorB[k];
    
    for(i = (aRows-1); i > k; i--)
    x[k] -= matrixT[k][i] * x[i];
    
    x[k] /= matrixT[k][k];
  }
    
  return 1;
}

double[][] solveMat(double matrixB[][], double matrixA[][]) // Оператор MatLab /  для одинаковых квадратных матриц
{
  return (transMat(operInvSlash(transMat(matrixA), transMat(matrixB))));
}

double[] solveMatV(double vector[], double matrixA[][]) // Оператор MatLab /  для одинаковых квадратных матриц
{
 
 int aRows  = matrixA.length;
 
 double [] vectorC   = new double [aRows];
 double [][] matrixT = new double [aRows][aRows];
 double [] vectorB   = new double [aRows];
 double [] vectorX   = new double [aRows];
 
  for(int j2 = 0; j2 < aRows ; j2++)
        vectorB[j2] = vector[j2];
        
    
   for(int i1 = 0; i1 < aRows ; i1++)
     for(int j1 = 0; j1 < aRows ; j1++)
       matrixT[i1][j1] = matrixA[i1][j1];
   matrixT = transMat(matrixT);   
   Gauss(matrixT, vectorB, vectorX);

   
    for(int j2 = 0; j2 < aRows ; j2++)
        vectorC[j2] = vectorX[j2];
        

 return vectorC;
}

//double[][] solveMatM(double matrixB[][], double matrixA[][]) // Оператор MatLab /  для общего случая матрица B больше строк чем столбцов у A содержит
//{
 
// int aRows   = matrixA.length,
//     blocks  = int(matrixB.length / aRows + 0.5);

// double [][] matrixC  = new double [matrixB.length][aRows],
//             matrixT  = new double [aRows][aRows],
//             matrixT2 = new double [aRows][aRows];
 
// for(int h = 0; h < blocks ; h++)
// {
//   for(int i = 0; i < aRows ; i++)
//   for(int j = 0; j < aRows ; j++)
//     matrixT[i][j] = matrixB[i + h * aRows][j];
   
//   matrixT2 = solveMat(matrixT, matrixA);
   
//   for(int i = 0; i < aRows ; i++)
//   for(int j = 0; j < aRows ; j++)
//     matrixC[i + h * aRows][j] = matrixT2[i][j];
     
// }
 
// int currRows = blocks * aRows;
 
// //if((matrixB.length % 2) != 0)
// while(currRows < matrixB.length)
// {
//   double lastRow[] = Mat_to_VecT(matrixB, currRows);
//   double vectorX[] = new double[lastRow.length];
//   Gauss(matrixA, lastRow, vectorX);
//   for(int j = 0; j < vectorX.length ; j++)
//      matrixC[currRows][j] = vectorX[j];
//   currRows++;
// }
 
// return matrixC;
//}

//**************************************************
//    Вычисление детерминанта

double detMatrix(double[][] matrix)
{
    int              n1              = matrix.length;
    double [][]      m               = new double [n1][n1];
    double []        temp_Vector     = new double [n1];
    
    double           det;
    int              i, 
                     j,
                     k,
                     kK = 0;
    
    for(k = 0; k < n1; k++)
    {
      for(i = 0; i < n1; i++)
      {
        m[i][k] = matrix[i][k];
      }
    }
    
     // Прямой ход
    for(k = 1; k < n1; k++)
    {
      for(i = k; i < n1; i++)
      {
        
        if(m[k-1][k-1] == 0)
        {
          int i1 = k - 1;
          int j1;
          while((i1 < n1)&&(m[i1][k-1] == 0)) i1++;
          if(i1 >= n1)return 0;
      
          for(j1 = 0; j1 < n1 ; j1++)
           temp_Vector[j1] = m[i1][j1];
        
          for(j1 = 0; j1 < n1 ; j1++)
          m[i1][j1] = m[k - 1][j1];
        
          for(j1 = 0; j1 < n1 ; j1++)
          m[k - 1][j1] = temp_Vector[j1];  
          kK++;
        }
    
      for(j = k; j < n1; j++)
      {
        if(abs_db(m[k-1][k-1]) < 0.000001)return 0;
        m[i][j] = m[i][j] - m[i][k-1] * m[k-1][j] / (m[k-1][k-1]);
      }
      m[i][k-1] = 0;
     }
    }
   
    det = 1;
    for(i = 0; i < n1 ; i++)
    {
      if(abs_db(det) > 1e12) 
      {
          det = 1e12;
          i = n1;
      }
      else
      det *= m[i][i]; 
    }
    
    if((kK % 2) != 0)
    {
      det = -det;
    }
   
    return det;
}

double [][] diag(double v[])
{
  int n = v.length;
  
  double ret[][] = new double [n][n];
  
  for(int i = 0; i < n; i++)
    ret[i][i] = v[i];
    
   return ret;
  
}

double [] diag(double m[][])
{
  int n = m.length;
  
  double ret[] = new double [n];
  
  for(int i = 0; i < n; i++)
    ret[i] = m[i][i];
    
   return ret;
  
}

double [] getStr(double A[][], int k)
{
  int n = A.length;
  double ret[] = new double [n];
  for(int i = 0; i < n ; i++)
    ret[i] = A[i][k];
  return ret;
}

double [] getSlb(double A[][], int k)
{
  int n = A.length;
  double ret[] = new double [n];
  for(int i = 0; i < n ; i++)
    ret[i] = A[i][k];
  return ret;
}

double [][] bsxfun_plus(double A[], double B[])
{
  
  int n = A.length;
  double ret[][] = new double [n][n];
  
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
    ret[i][j] = A[i] + B[j];
  
  return ret;
}

double [][] bsxfun_minus(double A[], double B[])
{
  
  int n = A.length;
  double ret[][] = new double [n][n];
  
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
    ret[i][j] = A[i] - B[j];
  
  return ret;
}

double [][] bsxfun_times(double A[], double B[])
{
  
  int n = A.length,
      m = B.length;
  double ret[][] = new double [n][m];
  
  for(int i = 0; i < n; i++)
    for(int j = 0; j < m; j++)
    ret[i][j] = A[i] * B[j];
  
  return ret;
}

double [][] bsxfun_times(double A[], double B[][])
{
  
  int n = A.length,
      m = B[0].length;
  double ret[][] = new double [n][m];
  
  for(int i = 0; i < n; i++)
    for(int j = 0; j < m; j++)
    ret[i][j] = A[i] * B[i][j];
  
  return ret;
}

double[][] solveMatM(double matrixA[][], double matrixB[][])
{
  
  int aRows   = matrixB[0].length,
      addRows  = matrixA.length % aRows,
      blocks  = int((matrixA.length)/ aRows + 0.5);
 
 double [][] matrixAe,
             matrixC,
             matrixT  = new double [aRows][aRows],
             matrixT2 = new double [aRows][aRows];
 double [][] matrixC2  = new double [matrixA.length][aRows];
 
 if(addRows != 0)
 {
  blocks += 1;
  matrixAe  = new double [matrixA.length + aRows][aRows];
  matrixC  = new double [matrixA.length + aRows][aRows];
  
  for(int i = 0; i < matrixA.length; i++)
    for(int j = 0; j < aRows ; j++)
      matrixAe[i][j] = matrixA[i][j];
     
   for(int h = 0; h < blocks ; h++)
   {
     for(int i = 0; i < aRows ; i++)
     for(int j = 0; j < aRows ; j++)
       matrixT[i][j] = matrixAe[i + h * aRows][j];
     
     matrixT2 = transMat(operInvSlash(transMat(matrixB), transMat(matrixT)));
     
     for(int i = 0; i < aRows ; i++)
     for(int j = 0; j < aRows ; j++)
       matrixC[i + h * aRows][j] = matrixT2[i][j];
       
   }
 
   for(int i = 0; i < matrixA.length; i++)
     for(int j = 0; j < aRows ; j++)
       matrixC2[i][j] = matrixC[i][j];
  }
  else
  {
   matrixAe  = new double [matrixA.length][aRows];
   matrixC  = new double [matrixA.length][aRows];
  
  for(int i = 0; i < matrixA.length; i++)
    for(int j = 0; j < aRows ; j++)
      matrixAe[i][j] = matrixA[i][j];
     
   for(int h = 0; h < blocks ; h++)
   {
     for(int i = 0; i < aRows ; i++)
     for(int j = 0; j < aRows ; j++)
       matrixT[i][j] = matrixAe[i + h * aRows][j];
     
     matrixT2 = transMat(operInvSlash(transMat(matrixB), transMat(matrixT)));
     
     for(int i = 0; i < aRows ; i++)
     for(int j = 0; j < aRows ; j++)
       matrixC[i + h * aRows][j] = matrixT2[i][j];
       
   }
 
   for(int i = 0; i < matrixA.length; i++)
     for(int j = 0; j < aRows ; j++)
       matrixC2[i][j] = matrixC[i][j];
  }
  
  return(matrixC2);
}

// Реализация оператора R\s matlab.

double [] operInvSlash(double matrixA[][], double vectorB[])
{
 int aRows  = matrixA.length,
     aCols  = matrixA[0].length;

 double [] vectorX = new double [aRows];
 double [] vector = vectorB;//vectorCoef(vectorB, 1.0);
   Gauss(matrixA, vector, vectorX);

 return vectorX;
}

double [][] operInvSlash(double matrixA[][], double matrixB[][])
{
 int aRows  = matrixA.length,
     aCols  = matrixA[0].length;
 double [][] matrixC;
 
 if(aRows == aCols)
 {
 matrixC           = new double [aRows][aRows];
 double [] vectorB = new double [aRows];
 double [] vectorX = new double [aRows];
 
 int count = 0;
 
 for(count = 0; count < aRows; count++)
 {
   
   for(int j2 = 0; j2 < aRows ; j2++)
        vectorB[j2] = matrixB[j2][count];
          
   Gauss(matrixA, vectorB, vectorX);

   
    for(int j2 = 0; j2 < aRows ; j2++)
        matrixC[j2][count] = vectorX[j2];
        
 }
 }
 else
 {
   matrixC  = new double [aRows][aCols];
   matrixC = QR_Solver(matrixA, matrixB);
 }
 
 return matrixC;
}

double[] EXP_F(double A[])
{
 
  double ret[] = new double [A.length];
  
  for(int i = 0; i < A.length ; i++)
    ret[i] = exp_db(A[i]);
  return ret;
}

double[][] EXP_F(double A[][])
{
 
  double ret[][] = new double [A.length][A[0].length];
  
  for(int i = 0; i < A.length ; i++)
   for(int j = 0; j < A[0].length; j++)
     ret[i][j] = exp_db(A[i][j]);
  return ret;
}

double[][] COS_F(double A[][])
{
 
  double ret[][] = new double [A.length][A[0].length];
  
  for(int i = 0; i < A.length ; i++)
   for(int j = 0; j < A[0].length; j++)
     ret[i][j] = cos_db(A[i][j]);
  return ret;
}

double[][] SIN_F(double A[][])
{
 
  double ret[][] = new double [A.length][A[0].length];
  
  for(int i = 0; i < A.length ; i++)
   for(int j = 0; j < A[0].length; j++)
     ret[i][j] = sin_db(A[i][j]);
  return ret;
}


double[][] ForDig(double A[][])
{
  double ret[][] = new double [A.length][A[0].length];
  
  for(int i = 0; i < A.length ; i++)
   for(int j = 0; j < A[0].length; j++)
     ret[i][j] = (int)(A[i][j] * 10000) / 10000.0;
  return ret;
}

//*************************************

int [] createInd(int bg, int ed)
{
 
  int ret[] = new int [ed - bg],
      i     = 0;
  
  for(i = 0; i < ret.length ; i++)
    ret[i] = bg + i;
    
  return ret;
  
}

double [] getVec(double M[], int ind[])
{

  double ret[] = new double [ind.length];
  int   i     = 0;
  
  for(i = 0; i < ret.length ; i++)
    ret[i] = M[ind[i]];
    
  return ret;
 
}

double [][] getMat(double Ss[][], int indR[], int indC[])
{

  double ret[][] = new double [indR.length][indC.length];
  int   i       = 0,
        j       = 0;
  
  for(i = 0; i < indR.length ; i++)
  for(j = 0; j < indC.length ; j++)
    ret[i][j] = Ss[indR[i]][indC[j]];
    
  return ret;
 
}

double [][] getMat(double S[][], int indR[])
{

  double ret[][] = new double [indR.length][S[0].length];
  int   i       = 0,
        j       = 0;
  
  for(i = 0; i < indR.length ; i++)
  for(j = 0; j < S[0].length ; j++)
    ret[i][j] = S[indR[i]][j];
    
  return ret;
 
}

//zi(:,1:d)
double [][] getMat(double S[][], int indB, int indE)
{

  double ret[][] = new double [S.length][indE - indB + 1];
  int   i       = 0,
        j       = 0,
        n       = indE - indB + 1;
  
  for(i = 0; i < S.length ; i++)
  for(j = 0; j < n ; j++)
    ret[i][j] = S[i][j];
    
  return ret;
 
}

double[] setVec(double RCv[], int ind[], double TRm[])
{
 
  int i = 0;
  double[] ret = RCv;//vectorCoef(RCv, 1.0);
  for(i = 0; i < ind.length ; i++)
    ret[ind[i]] = TRm[i];
  return ret;
  
};

double[][] setMat(double RCv[][], int indR[], int indC[], double TRm[][])
{
 
  int i = 0,
      j = 0;
  double ret[][] = RCv;//coeffProdMat(1.0, RCv);    
  for(i = 0; i < indR.length ; i++)
  for(j = 0; j < indC.length ; j++)
    ret[indR[i]][indC[j]] = TRm[i][j];
  return ret;
  
};

double [][] setMat(double RCv[][], int indR[], double TRm[][])
{
 
  int i = 0,
      j = 0;
  
  double ret[][] = RCv;//coeffProdMat(1.0, RCv);  
  
  for(i = 0; i < indR.length ; i++)
  for(j = 0; j < TRm[0].length ; j++)
    ret[indR[i]][j] = TRm[i][j];
  
  return ret;
    
};

double [][][] set3DMass(double RCv[][][], int indR[], int indC[], int ind, double TRm[][])
{
 
  int i = 0,
      j = 0;
  
  double ret[][][] = coeffProd3DMass(1.0, RCv);  
  
  for(i = 0; i < indR.length ; i++)
  for(j = 0; j < TRm[0].length ; j++)
    ret[indR[i]][indC[j]][ind] = TRm[i][j];
  
  return ret;
    
};

double[][] addTo(double RCv[][], int ind_i[], int ind_j[], int ind_k[], double TRm[][])
{
  
 //q = S(j,i)*C; S(j,k) = q; S(k,j) = q'
 double tempR[][] = RCv,//coeffProdMat(1.0, RCv),
        tempT[][] = TRm;//coeffProdMat(1.0, TRm); 
 double q[][] = matrixMultiply(getMat(tempR, ind_j, ind_i),  tempT);
 tempR = setMat(tempR, ind_j, ind_k, q);
 tempR = setMat(tempR, ind_k, ind_j, transMat(q));
 return (tempR);
 
};

double[][] addToSimpl(double RCv[][], int ind_i[], int ind_k[], double TRm[][])
{
  
 //S(i,k) = S(i,i)*C; S(k,i) = S(i,k)'  
 double tempR[][] = RCv,//coeffProdMat(1.0, RCv),
        tempT[][] = TRm;//coeffProdMat(1.0, TRm); 
 double q[][] = matrixMultiply(getMat(tempR, ind_i, ind_i),  tempT);
 tempR = setMat(tempR, ind_i, ind_k, q);
 tempR = setMat(tempR, ind_k, ind_i, transMat(getMat(tempR, ind_i, ind_k)));
 return (tempR);
 
};

//X = reshape(1:D3*D3,[D3 D3]); XT = X'; Sds = (Sds + Sds(XT(:),:))/2
double [][] symmetrizeA(int D, double A[][])
{
  
  double ret[][] = new double [A.length][A[0].length];
  
  int X[][] = new int [D][D];
  
  int count = 0;
  for(int i1 = 0; i1 < D; i1++) 
  for(int i2 = 0; i2 < D; i2++) 
    X[i2][i1] = count++;
  
  int XT[][] = transMat(X);
  int fXT_ind = 0;
  int fXT[] = new int [D * D];
  for(int i1 = 0; i1 < D; i1++) 
  for(int i2 = 0; i2 < D; i2++) 
    fXT[fXT_ind ++] = XT[i2][i1];
  
  double A_XT[][] = new double [A.length][A[0].length];
  for(int i1 = 0; i1 < A.length; i1++) 
  for(int i2 = 0; i2 < A[0].length; i2++) 
    A_XT[i1][i2] = A[fXT[i1]][i2];
    
  A = matrixADD(A, A_XT);
  ret = coeffProdMat(0.5, A);
  
  return ret;
  
}

//X = reshape(1:D0*D0,[D0 D0]); XT = X'; Sds = (Sds + Sds(:,XT(:)))/2
double [][] symmetrizeB(int D, double A[][])
{
  
  double ret[][] = new double [A.length][A[0].length];
  
  int X[][] = new int [D][D];
  
  int count = 0;
  for(int i1 = 0; i1 < D; i1++) 
  for(int i2 = 0; i2 < D; i2++) 
    X[i2][i1] = count++;
  
  int XT[][] = transMat(X);
  int fXT_ind = 0;
  int fXT[] = new int [D * D];
  for(int i1 = 0; i1 < D; i1++) 
  for(int i2 = 0; i2 < D; i2++) 
    fXT[fXT_ind ++] = XT[i2][i1];
  
  double A_XT[][] = new double [A.length][A[0].length];
  for(int i1 = 0; i1 < A.length; i1++) 
  for(int i2 = 0; i2 < A[0].length; i2++) 
    A_XT[i1][i2] = A[i1][fXT[i2]];
    
  A = matrixADD(A, A_XT);
  ret = coeffProdMat(0.5, A);
  
  return ret;
  
}

//************************************************
//  Новая версия оператора \

double[][] matrix_minor(double A[][], int d)
{
  int    m = A.length,
         n = A[0].length;
  double ret[][] = new double [m][n];
  
  for (int i = 0; i < d; i++)
    ret[i][i] = 1;
  for (int i = d; i < m; i++)
    for (int j = d; j < n; j++)
      ret[i][j] = A[i][j];
  return ret;
}

/* c = a + b * s */
double [] vmadd(double a[], double b[], double s, int n)
{
  double ret[] = new double [n];
  for (int i = 0; i < n; i++)
    ret[i] = a[i] + s * b[i];
  return ret;
}

/* m = I - v v^T */
double [][] vmul(double v[], int n)
{
  double ret[][] = new double [n][n];
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      ret[i][j] = -2 *  v[i] * v[j];
  for (int i = 0; i < n; i++)
    ret[i][i] += 1;

  return ret;
}

/* ||x|| */
double vnorm(double A[], int n)
{
  double sum = 0;
  for (int i = 0; i < n; i++) sum += A[i] * A[i];
  return sqrt_db(sum);
}

/* y = x / d */
double[] vdiv(double A[], double d, int n)
{
  double ret[] = new double [n];
  for (int i = 0; i < n; i++) ret[i] = A[i] / d;
  return ret;
}

/* take c-th column of m, put in v */
double[] mcol(double A[][], int c)
{
  double ret[] = new double [A.length];
  for (int i = 0; i < A.length; i++)
    ret[i] = A[i][c];
  return ret;
}

double [][] matrix_mul(double x[][], double y[][])
{
  double r[][] = new double[x[0].length][y.length];
  for (int i = 0; i < x.length; i++)
    for (int j = 0; j < y[0].length; j++)
      for (int k = 0; k < x[0].length; k++)
        r[i][j] += x[i][k] * y[k][j];
  return r;
}

double [][] matrix_transpose(double A[][])
{
  double ret[][] = A;//coeffProdMat(1.0, A);
  for (int i = 0; i < ret.length; i++) 
  {
    for (int j = 0; j < i; j++) 
    {
      double t = ret[i][j];
      ret[i][j] = ret[j][i];
      ret[j][i] = t;
    }
  }
  return ret;
}

double [][] QR_Solver(double A[][], double B[][])
{
  
  double ret[][] = new double[1][1];
  //**********************************************
  // Построение QR разложения
  
  int m = A.length,
      n = A[0].length;
  double q[][][] = new double[m][0][0];
  double z[][], z1[][], Q[][], R[][];
  z = A;//coeffProdMat(1.0, A);
  
  for (int k = 0; (k < n) && (k < (m - 1)); k++) 
  {
    double e[] = new double[m],
           x[] = new double[m],
           a;
    z1 = matrix_minor(z, k);
    z = z1;//coeffProdMat(1.0, z1);

    x = mcol(z, k);
    a = vnorm(x, m);
    if (A[k][k] > 0) a = -a;

    for (int i = 0; i < m; i++)
      e[i] = (i == k) ? 1 : 0;

    e = vmadd(x, e, a, m);
    e = vdiv(e, vnorm(e, m), m);
    q[k] = vmul(e, m);
    z1 = matrixMultiply/*matrix_mul*/(q[k], z);
    z = z1;//coeffProdMat(1.0, z1);
  }
  z = null;
  Q = q[0];//coeffProdMat(1.0, q[0]);
  R = matrixMultiply/*matrix_mul*/(q[0], A);
  for (int i = 1; (i < n) && (i < (m - 1)); i++) 
  {
    z1 = matrixMultiply/*matrix_mul*/(q[i], Q);
    Q = z1;//coeffProdMat(1.0, z1);
  }
  z = matrixMultiply/*matrix_mul*/(Q, A);
  double Qt[][] = matrix_transpose(Q);
  Q = null;
  Q = new double[m][min(m, n)];
  for(int i = 0; i < m ; i++)
    for(int j = 0; j < min(m, n); j++)
     Q[i][j] = Qt[i][j];
  
  R = null;
  R = new double[min(m, n)][n];
  for(int i = 0; i < min(m, n) ; i++)
    for(int j = 0; j < n; j++)
     R[i][j] = z[i][j];
    
  //println("Матрица Q");
  //printMat(Q);
  //println("Матрица R");
  //printMat(R);
  //printMat(matrixMultiply(Q, R));
  //*************************************************
  
  //*************************************************
  // Решение уравнения Ax = B
  // Определение Y = Qt*B
  double Y[][] = matrixMultiply(transMat(Q), B);
  println("Матрица R");
  printMat(R);
  println("Матрица Y");
  printMat(Y);
  
   ret = new double[Y.length][Y.length];
  for(int i = 0; i < Y[0].length ; i++)
  {
    ret[i] = triMatSolve(R, Mat_to_Vec(Y, i));//vectorCoef(triMatSolve(R, Mat_to_Vec(Y, i)), 1.0);
  }
  
  //*************************************************
  return transMat(ret);
}

//**********************************************
// Вычисление решения Ax = B. Если A квадратная верхне триугольная матрица, B вектор
double[] triMatSolve(double matrixA[][], double vectorB[])
{
  
  int     aRows = matrixA.length,
          i,
          k;
  double  x[]   = new double [aRows];
  
  
  // Обратный ход
  
  for(k = (aRows-1); k >= 0; k--)
  {
    x[k] = vectorB[k];
    
    for(i = (aRows-1); i > k; i--)
    x[k] -= matrixA[k][i] * x[i];
    
    x[k] /= matrixA[k][k];
  }
    
  return x;
}

tp_p let_p(tp_p p)
{
  tp_p pl = new tp_p();
  if(mode_pol == 0)
  {
    pl.inputs = p.inputs;//coeffProdMat(1.0, p.inputs);
    pl.hyp = p.hyp;//vectorCoef(p.hyp, 1.0);
    pl.targets =p.targets; //vectorCoef(p.targets, 1.0);
  }
  else
  {
    pl.w = p.w;//vectorCoef(p.w, 1.0);
    pl.b = p.b;
  }
  return pl;
};

double[] unwrap(tp_p p)
{
  double ret_p[];
  if(mode_pol == 0)
  {
      ret_p = new double[ p.inputs.length * p.inputs[0].length + p.targets.length + p.hyp.length ];
      int p_index = 0;
    
      for (int i1 = 0; i1 < p.hyp.length; i1++)
        ret_p[p_index++] = p.hyp[i1];
    
      for (int i1 = 0; i1 < p.inputs[0].length; i1++)
        for (int j1 = 0; j1 < p.inputs.length; j1++)
          ret_p[p_index++] = p.inputs[j1][i1];
    
      for (int i1 = 0; i1 < p.targets.length; i1++)
        ret_p[p_index++] = p.targets[i1];
  }
  else
  {
      ret_p = new double[p.w.length + 1];
      int p_index = 0;
  
      ret_p[p_index++] = p.b;
      
      for (int i1 = 0; i1 < p.w.length; i1++)
        ret_p[p_index++] = p.w[i1];
  }
  
  return ret_p;

};

tp_p rewrap(tp_p _pl, double _dp[])
{
  tp_p ret = new tp_p();
  
  if(mode_pol == 0)
  {
      int  k = 0,
           n_hyp = _pl.hyp.length,
           n_inp_rows = _pl.inputs.length,
           n_inp_cols = _pl.inputs[0].length,
           n_targ = _pl.targets.length;
      ret.hyp = new double[n_hyp];
      ret.inputs = new double[n_inp_rows][n_inp_cols];
      ret.targets = new double[n_targ];
    
      for (int i = 0; i < n_hyp; i++)
      {
        ret.hyp[i] = _dp[k++];
      };
    
      for (int i = 0; i < n_inp_cols; i++)
        for (int j = 0; j < n_inp_rows; j++)
        {
          ret.inputs[j][i] = _dp[k++];
        };
    
      for (int i = 0; i < n_targ; i++)
      {
        ret.targets[i] = _dp[k++];
      };
  }
  else
  {
      int  k = 0;
      ret.w = new double[_pl.w.length];
      
      ret.b = _dp[k++];
      
      for (int i = 0; i < _pl.w.length; i++)
      {
        ret.w[i] = _dp[k++];
      }
  }
  
  return ret;

}



// Вектор-столбец → double[]
double[] csvToVector(double[][] m) {
  double[] v = new double[m.length];
  for (int i = 0; i < m.length; i++) v[i] = m[i][0];
  return v;
}
double[][] loadCSV(String filename) {
  String[] lines = loadStrings(filename);
  
  // считаем непустые строки
  int validCount = 0;
  for (String l : lines)
    if (l.trim().length() > 0) validCount++;

  double[][] result = new double[validCount][];
  int row = 0;
  for (String line : lines) {
    if (line.trim().length() == 0) continue;
    String[] parts = splitTokens(line, ",");
    result[row] = new double[parts.length];
    for (int j = 0; j < parts.length; j++)
      result[row][j] = Double.parseDouble(parts[j].trim());
    row++;
  }
  return result;
}
void loadMatlabData() {

  // --- mu0 и S0 (захардкодим, не нужен CSV) ---
 mu0L = new double[] {0, 0, Math.PI, Math.PI};  // [4]

  S0L = new double[][] {
    {0.01,  0,      0,       0      },
    {0,     0.01,   0,       0      },
    {0,     0,      0.0001,  0      },
    {0,     0,      0,       0.0001 }
  };  // [4][4]

  // --- policy.p ---
  double[][] pol_inputs = loadCSV("policy_inputs.csv");   // [100][6]
  double[][] pol_tg_2d  = loadCSV("policy_targets.csv");  // [100][1]
  double[][] pol_hyp_2d = loadCSV("policy_hyp.csv");      // [8][1]

  policy.p.inputs  = pol_inputs;

  policy.p.targets = new double[pol_tg_2d.length];
  for (int i = 0; i < pol_tg_2d.length; i++)
    policy.p.targets[i] = pol_tg_2d[i][0];

  policy.p.hyp = new double[pol_hyp_2d.length];
  for (int i = 0; i < pol_hyp_2d.length; i++)
    policy.p.hyp[i] = pol_hyp_2d[i][0];

  // --- dynmodel ---
  double[][] dyn_inputs  = loadCSV("dynmodel_inputs.csv");   // [1140][7]
  double[][] dyn_targets = loadCSV("dynmodel_targets.csv");  // [1140][4]
  double[][] dyn_hyp     = loadCSV("dynmodel_hyp.csv");      // [9][4]

  int ni = dyn_inputs.length;        // 1140
  int mi = dyn_inputs[0].length;     // 7
  int mt = dyn_targets[0].length;    // 4

  dynmodel = new dynmodel_t(ni, mi, mt);

  for (int i = 0; i < ni; i++) {
    for (int j = 0; j < mi; j++)
      dynmodel.inputs[i][j]  = dyn_inputs[i][j];
    for (int j = 0; j < mt; j++)
      dynmodel.targets[i][j] = dyn_targets[i][j];
  }

  dynmodel.hyp = new double[dyn_hyp.length][dyn_hyp[0].length];
  for (int i = 0; i < dyn_hyp.length; i++)
    for (int j = 0; j < dyn_hyp[0].length; j++)
      dynmodel.hyp[i][j] = dyn_hyp[i][j];

  dynmodel.n   = ni;
  dynmodel.mIn = mi;
  dynmodel.mTg = mt;

  println("loadMatlabData OK");
  println("  dynmodel: " + ni + "x" + mi + " -> " + ni + "x" + mt);
  println("  policy.p: " + pol_inputs.length + "x" + pol_inputs[0].length);
}
void draw_cost_graph(PGraphics pg, 
                     double[][] predCost,    // предсказанная стоимость
                     double[][] realCost,    // реальная стоимость
                     int iterCount,         // сколько итераций уже есть
                     int H)                 // горизонт
{
  // Размеры области графика
  int gx = 50,   gy = 50;   // левый верхний угол
  int gw = 700,  gh = 300;  // ширина и высота
  
  // Фон
  pg.fill(255);
  pg.noStroke();
  pg.rect(gx, gy, gw, gh);
  
  // Оси
  pg.stroke(0);
  pg.strokeWeight(1);
  pg.line(gx, gy + gh, gx + gw, gy + gh); // ось X
  pg.line(gx, gy, gx, gy + gh);            // ось Y
  
  // Подписи осей
  pg.fill(0);
  pg.textSize(14);
  pg.text("time in s", gx + gw / 2, gy + gh + 30);
  pg.text("immediate cost", gx - 40, gy + gh / 2);
  
  // Шкала Y: от 0 до 1
  for(int i = 0; i <= 5; i++) {
    float val = i / 5.0;
    int yPos = (int)(gy + gh - val * gh);
    pg.stroke(200);
    pg.line(gx, yPos, gx + gw, yPos);   // сетка
    pg.fill(0);
    pg.text(nf(val, 1, 1), gx - 30, yPos + 5);
  }
  
  // Шкала X: время в секундах
  float dt_val = 0.05;
  for(int i = 0; i <= 10; i++) {
    int xPos = (int)(gx + (float)i / 10.0 * gw);
    pg.stroke(200);
    pg.line(xPos, gy, xPos, gy + gh);   // сетка
    pg.fill(0);
    pg.text("" + i, xPos - 5, gy + gh + 15);
  }
  
  if(iterCount == 0) return;
  
  int lastIter = iterCount - 1;
  
  // --- Доверительный интервал предсказанной стоимости (синий) ---
  // Для простоты рисуем среднее по всем итерациям как полосу
  pg.noStroke();
  pg.fill(100, 150, 255, 80); // полупрозрачный синий
  
  // Верхняя граница (max по итерациям)
  // Нижняя граница (min по итерациям)
  for(int t = 0; t < H - 1; t++) {
    double tMax = 0, tMin = 1;
    for(int iter = 0; iter < iterCount; iter++) {
      if(predCost[iter][t] > tMax) tMax = predCost[iter][t];
      if(predCost[iter][t] < tMin) tMin = predCost[iter][t];
    }
    int x1 = (int)(gx + (float)t / H * gw);
    int x2 = (int)(gx + (float)(t+1) / H * gw);
    int yTop = (int)(gy + gh - tMax * gh);
    int yBot = (int)(gy + gh - tMin * gh);
    pg.rect(x1, yTop, x2 - x1, yBot - yTop);
  }
  
  // --- Среднее предсказанной стоимости (синяя линия) ---
  pg.stroke(0, 0, 255);
  pg.strokeWeight(2);
  pg.noFill();
  pg.beginShape();
  for(int t = 0; t < H; t++) {
    float mean = 0;
    for(int iter = 0; iter < iterCount; iter++)
      mean += predCost[iter][t];
    mean /= iterCount;
    int xPos = (int)(gx + (float)t / H * gw);
    int yPos = (int)(gy + gh - mean * gh);
    pg.vertex(xPos, yPos);
  }
  pg.endShape();
  
  // --- Реальная стоимость последнего прогона (красная линия) ---
  pg.stroke(255, 0, 0);
  pg.strokeWeight(2);
  pg.beginShape();
  for(int t = 0; t < H; t++) {
    int xPos = (int)(gx + (float)t / H * gw);
    int yPos = (int)(gy + gh - realCost[lastIter][t] * gh);
    pg.vertex(xPos, yPos);
  }
  pg.endShape();
  
  // --- Минимальная стоимость (чёрная пунктирная) ---
  pg.stroke(0);
  pg.strokeWeight(1);
  for(int t = 0; t < H - 1; t++) {
    if(t % 6 < 3) { // пунктир
      int x1 = (int)(gx + (float)t / H * gw);
      int x2 = (int)(gx + (float)(t+1) / H * gw);
      pg.line(x1, gy + gh, x2, gy + gh); // минимум = 0
    }
  }
  
  // Легенда
  pg.strokeWeight(2);
  pg.stroke(0, 0, 255); pg.line(gx + 20, gy + 20, gx + 50, gy + 20);
  pg.fill(0); pg.text("pred. cost", gx + 55, gy + 25);
  pg.stroke(255, 0, 0); pg.line(gx + 20, gy + 40, gx + 50, gy + 40);
  pg.fill(0); pg.text("cost of rollout", gx + 55, gy + 45);
  
  // Заголовок
  pg.textSize(16);
  pg.fill(0);
  pg.text("after " + iterCount + " policy searches", gx + gw/2 - 80, gy - 15);
}
