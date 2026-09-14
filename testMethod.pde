void TestPass(double truth, double result){
  if(Math.abs(truth-result)<1e-5){
    println("Test Pass");
  }
  else{
    println("Test Failed");
  }
}

void TestPassArr(double[] A, double[] B){
  if(A.length == B.length){
    for(int i = 0; i<A.length;i++){
      if(Math.abs(A[i]-B[i])<1e-5){
      }
      else{
        println("Test Failed");
        return;
      }
    }
    println("Test Pass");
  }
  else{
  println("Test Failed");
}
}
void TestPassMatrix(double[][] A, double[][] B){
  if(A.length == B.length && A[0].length == B[0].length){
    for(int i = 0; i<A.length;i++){
      for(int j = 0; j<A[0].length;j++){
        if(Math.abs(A[i][j]-B[i][j])<1e-5){
      }
      else{
        println("Test Failed");
        return;
      }
      }
    }
    println("Test Pass");
  }
  else{
  println("Test Failed");
}
}

void test_sin_db(){
  double x = PI;
  double result = sin_db(x);
  TestPass(result,0.00000); 
}
void test_cos_db(){
  double x = PI/2;
  double result = cos_db(x);
  TestPass(result,0.00000); 
}
void test_log_db(){
  double x = 1;
  double result = log_db(x);
  TestPass(result,0.00000); 
}
void test_exp_db(){
  double x = 1;
  double result = exp_db(x);
  TestPass(result, 2.71828); 
}
void test_sqrt_db(){
  double x = 100;
  double result = sqrt_db(x);
  TestPass(result,10);
}
void test_pow_db(){
  double x = 12;
  double y = 8;
  double result = pow_db(x,y);
  TestPass(result,429981696);
}
void test_abs_db(){
  double x = -1;
  double y = 1;
  double result = abs_db(x);
  TestPass(result,1);
  result = abs_db(y);
  TestPass(result,1);
}
void test_max_db(){
  double x = -1;
  double y = 1;
  double result = max_db(x,y);
  TestPass(result,1);
}
void test_min_db(){
  double x = -1;
  double y = 1;
  double result = min_db(x,y);
  TestPass(result,-1);
}
void test_constrain_db(){
  double x = 10;
  double a = 11;
  double b = 20;
  double result = constrain_db(x,a,b);
  TestPass(result,11);
}
void test_append_db_arr(){
  double x = 1;
  double[] A = {4,51};
  double[] result = append_db(A,x);
  double[] tr = {4,51,1};
  TestPassArr(result, tr);
}
void test_append_db_matrix(){
  double[][] A = {{1,2},{3,4}};
  double[][] B = {{1,2},{3,4},{0,0}};
  double[][] result = append_db(A);
  TestPassMatrix(result,B);
}
void test_transMat(){
  double[][] x = {{1,2},{3,4}};
  double[][] result = transMat(x);
  double[][] B = {{1,3},{2,4}};
  TestPassMatrix(result,B);
}
void test_matrixMultiply(){
  double[][] x = {{1,2},{3,4}};
  double[][] y = {{1,5},{-3,7}};
  double[][] result = matrixMultiply(x,y);
  double[][] B = {{-5,19},{-9,43}};
  TestPassMatrix(result,B);
}

void test_sparse_t(){
  println("==============");
  double[][] x = {{1,2},{0,4}};
  sparse_t y = new sparse_t(x);
  y.print();
  println("==============");
}
