
PrintWriter     output_X,
                output_Y,
                output_P;

void saveData()
{
  
  output_X = createWriter("dataX.txt");
  output_X.println(xModel.length + "\t" + xModel[0].length); 
  for(int i = 0 ; i < xModel.length ; i++)
  {
    for(int j = 0 ; j < xModel[0].length ; j++)
      output_X.print(xModel[i][j] + "\t"); 
    
    output_X.println();
  }   
  output_X.flush(); 
  output_X.close(); 
  
  output_Y = createWriter("dataY.txt");
  output_Y.println(yModel.length + "\t" + yModel[0].length); 
  for(int i = 0 ; i < yModel.length ; i++)
  {
    for(int j = 0 ; j < yModel[0].length ; j++)
      output_Y.print(yModel[i][j] + "\t"); 
    
    output_Y.println();
  }   
  output_Y.flush(); 
  output_Y.close(); 
 
  policyGP2 = new dynmodel_t(policy.p.inputs.length, policy.p.inputs[0].length, policy.p.targets.length);
  policyGP2.hyp    =   vectorCoef2(policy.p.hyp, 1.0);
  policyGP2.inputs =  policy.p.inputs; //coeffProdMat(1.0, policy.p.inputs);
  policyGP2.targets =  vectorCoef2(policy.p.targets, 1.0);
  
  printMatFor(policy.p.inputs);
  printVecFor(policy.p.targets);
  printVecFor(policy.p.hyp);
  
  output_P = createWriter("dataP.txt");
  output_P.println(policyGP2.inputs.length + "\t" + policyGP2.inputs[0].length); 
  for(int i = 0 ; i < policyGP2.inputs.length ; i++)
  {
    for(int j = 0 ; j < policyGP2.inputs[0].length ; j++)
      output_P.print(policyGP2.inputs[i][j] + "\t"); 
    
    output_P.println();
  }   
  
  output_P.println(policyGP2.hyp.length + "\t"); 
  for(int i = 0 ; i < policyGP2.hyp.length ; i++)
  {
    output_P.println(policyGP2.hyp[i][0] + "\t"); 
  }   
  
  output_P.println(policyGP2.targets.length + "\t"); 
  for(int i = 0 ; i < policyGP2.targets.length ; i++)
  {
    output_P.println(policyGP2.targets[i][0] + "\t"); 
  }   
  
  output_P.flush(); 
  output_P.close(); 
  
};


void readData() {
  int    m = 0, 
         n = 0;
  
  BufferedReader reader = createReader("dataX.txt");
  String line = null;
  try {
    if((line = reader.readLine()) != null)
    {
      String[] pieces = split(line, TAB); //<>// //<>//
      m = int(pieces[0]);
      n = int(pieces[1]);
      xx = new double [m][n];
 
    int i = 0;
    while ((line = reader.readLine()) != null) {
      pieces = split(line, TAB);
      for(int j = 0; j < n ; j++)
        xx[i][j] = (double)(float(pieces[j]));
     i++;
    }
    reader.close();
    
     printMat_f(xx);
     println();
    }
  } catch (IOException e) {
    e.printStackTrace();
  }
  
  reader = createReader("dataY.txt"); //<>// //<>//

  try {
    if((line = reader.readLine()) != null)
    {
      String[] pieces = split(line, TAB);
      m = int(pieces[0]);
      n = int(pieces[1]);
      yy = new double [m][n];
 
    int i = 0;
    while ((line = reader.readLine()) != null) {
      pieces = split(line, TAB);
      for(int j = 0; j < n ; j++)
        yy[i][j] = (double)(float(pieces[j]));
     i++;
    }
    reader.close();
    
     printMat_f(yy);
     println();
    }
  } catch (IOException e) {
    e.printStackTrace();
  }
  
  xModel =xx; //coeffProdMat(1.0, xx);
  yModel =yy;// coeffProdMat(1.0, yy);
 // mode = 1;
  
  reader = createReader("dataP.txt"); //<>// //<>//
  line = null;
  try {
    if((line = reader.readLine()) != null)
    {
      String[] pieces = split(line, TAB);
      m = int(pieces[0]);
      n = int(pieces[1]);
      policy.p.inputs = new double [m][n];
 
    int i = 0;
    while (i < m) {
      line = reader.readLine();
      pieces = split(line, TAB);
      for(int j = 0; j < n ; j++)
        policy.p.inputs[i][j] = (double)(float(pieces[j]));
     i++;
    }
   
   line = reader.readLine();
   pieces = split(line, TAB);
      m = int(pieces[0]);
      policy.p.hyp = new double [m];
 
    i = 0;
    while (i < m) {
      line = reader.readLine();
      pieces = split(line, TAB);
      policy.p.hyp[i] = (double)(float(pieces[0]));
     i++;
    }    
    
    line = reader.readLine();
    pieces = split(line, TAB);
      m = int(pieces[0]);
      policy.p.targets = new double [m];
 
    i = 0;
    while (i < m) {
      line = reader.readLine();
      pieces = split(line, TAB);
      policy.p.targets[i] = (double)(float(pieces[0]));
     i++;
    }    
        
    reader.close();
    
      printMatFor(policy.p.inputs); //<>// //<>//
      printVecFor(policy.p.targets);
      printVecFor(policy.p.hyp);

     println();
    }
  } catch (IOException e) {
    e.printStackTrace();
  }
  
} 
