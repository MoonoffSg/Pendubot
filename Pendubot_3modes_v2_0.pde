import com.reades.mapthing.*; //<>// //<>//
import net.divbyzero.gpx.*;
import net.divbyzero.gpx.parser.*;

//test
// Универсальная версия с 2 - я политиками

import org.quark.jasmine.*;
import processing.serial.*;
import controlP5.*;

void setup() {

  size(1200, 850, P2D);
  stroke(0, 0, 0);
  colorMode(RGB, 256);

  pgBR = createGraphics(800, 800, P3D);
  myScreen = new VidPort(vidHor, vidVer, vidHor / 2.0, vidVer / 2.0);

  currentTime = millis();

  mode = 100;
  stateW[0] = 0;
  stateW[1] = 0;

  indexH = 0;
  policy.p = let_p(policeInit(policy, 200));

  cp5 = new ControlP5(this);
  PFont font = createFont("Arial bold", 15, false);
  cp5.setFont(font);

  cp5.addButton("СИМУЛЯЦИЯ")
    .setPosition(100, 100)
    .setSize(150, 50)
    ;

  cp5.addButton("ПРОВЕРКА")
    .setPosition(100, 200)
    .setSize(150, 50)
    ;

  cp5.addButton("ОБУЧЕНИЕ")
    .setPosition(100, 300)
    .setSize(150, 50)
    ;

  cp5.addButton("ПРОЧИТАТЬ")
    .setPosition(100, 400)
    .setSize(150, 50)
    ;
  //cp5.addButton("Graph")
  //  .setPosition(100, 500)
  //  .setSize(150, 50)
  //  .setLabel("График")
  //  ;
  pgCost = createGraphics(800, 400, P2D);
  costHistory     = new double[N][H];
  costRealHistory = new double[N][H];
  test_sin_db();
  test_cos_db();
  test_log_db();
  test_exp_db();
  test_sqrt_db();
  test_pow_db();
  test_abs_db();
  test_max_db();
  test_min_db();
  test_constrain_db();
  test_append_db_arr();
  test_append_db_matrix();
  test_transMat();
  test_matrixMultiply();
  test_sparse_t();
  delay(1000);
}
public void Graph(int theValue)
{
  if (dynmodel == null) {
    loadMatlabData();
  }

  nargout_valueSh_t vRes = valueSh(
    let_p(policy.p), mu0L, S0L, dynmodel, policy, cost, H
  );

  // vRes.L — double[], оборачиваем в double[][]
  double[][] predCost = new double[1][vRes.L.length];
  for (int t = 0; t < vRes.L.length; t++)
    predCost[0][t] = vRes.L[t];

  // Graph.csv теперь 20×60
  double[][] realCostMatrix = loadCSV("Graph.csv");
  int iterCount = realCostMatrix.length;  // = 20

  pgCost.beginDraw();
  pgCost.background(255);
  draw_cost_graph(pgCost, predCost, realCostMatrix, iterCount, H);
  pgCost.endDraw();
  image(pgCost, 0, 450);
}



public void СИМУЛЯЦИЯ(int theValue) {

  mode_simulation = 1;
  mode = 0;
  stateW[0] = 0;
  stateW[1] = 0;
  stateW[2] = PI;
  stateW[3] = PI;

  // policy.p = new tp_p();
  //policy.p = let_p(policeInit(policy, 200));
  //policy.p.inputs = coeffProdMat(1.0, inpP);
  //policy.p.targets = vectorCoef(targP, 1.0);
  //policy.p.hyp = vectorCoef(hypP, 1.0);
}

public void ПРОВЕРКА(int theValue) {

  mode_simulation = 2;
  stateW[0] = 0;
  stateW[1] = 0;
  stateW[2] = -0.01;
  stateW[3] = 0;
}

public void ОБУЧЕНИЕ(int theValue) {

  println("Обучение");

  //stateW[0] = 0.977741787534506;
  //stateW[1] = -1.51097202374961;
  //stateW[2] = 3.12982827669948;
  //stateW[3] =  3.09640829252530;

  //int ind[] = {2, 3};
  //double vz[][] = {
  //                  {1.90184029997389,-2.79057025187072,0.0511770824422783,-0.0919869839039427},
  //                  {-2.79057025187072,4.35924250312016,-0.0683152835203531,0.133498013441253},
  //                  {0.0511770824422783,-0.0683152835203531,0.00226105603379251,-0.00294058781975162},
  //                  {-0.0919869839039427,0.133498013441253,-0.00294058781975162,0.00577907432012905}
  //                };

  //nargout_loss_cp_t rett = loss_pendubot(cost, stateW, vz);

  //println(rett.S2);
  //printVec(rett.dLds);

  //gTrig_Full_r test =  gTrigF(stateW, vz, ind, 1);
  //printVec(test.M);
  //printMat(test.dCdv);
  //printMat(test.dVdv, 16);//4
  //int D0 = 4,
  //    D1 = D0 + 2 * policy.angle.length,
  //    D2 = D1 + 1,
  //    D3 = D2 + D0;

  //double P[][] = new double [D0][D0 + D2];
  //for(int i1 = 0; i1 < D0; i1++)
  //{
  //   P[i1][i1 + D2] = 1.0;
  //}

  //for(int i1 = 0; i1 < difi.length; i1++)
  //{
  //  P[difi[i1]][difi[i1]] = 1.0;
  //}

  //printMat(P);

  //sparse_t P_sp = new sparse_t(P);

  //sparse_t P_sp_tr = transMat_SP(P_sp);

  //P_sp.print();
  //P_sp_tr.print();

  //double M[] = {0,
  //              0,
  //              3.14159265358979,
  //              3.14159265358979,
  //              1.22458556833818e-16,
  //              -0.999950001249979,
  //              1.22458556833818e-16,
  //              -0.999950001249979,
  //              0.0392646004583106,
  //              0.0199841132527624,
  //              -0.0210833078893056,
  //              -0.0119838892839487,
  //              -0.00266255526542893};

  //double S[][]  = {
  //                  {0.0100000000000000,  0,  0,  0, 0,  0,  0,  0,  0.00396398441456158,  0.000870810561497663,  -0.00169959311610240,  0.000516042564142335,  -0.000144426919017256},
  //                  {0,  0.0100000000000000,  0,  0,  0,  0,  0,  0,  -0.00159529303278272,  -0.00265238798005091,  -0.000269791083799460,  -3.67677463138513e-05,  0.000521672148730322},
  //                  {0,  0,  0.000100000000000000,  0,  -9.99950001249979e-05,  -1.22458556833818e-20,  0,  0,  0.000695237051140377,  0.000360647869139166,  -0.000338870318833434,  1.51424549046011e-06,  -1.30943570454497e-06},
  //                  {0,  0,  0,  0.000100000000000000,  0,  0,  -9.99950001249979e-05,  -1.22458556833818e-20,  -7.76064579396836e-05,  8.54821321556754e-05,  0.000232492152096411,  9.82006648391333e-07,  -5.22315335211565e-06},
  //                  {0,  0,  -9.99950001249979e-05,  0,  9.99900006666277e-05,  1.22446311641411e-20,  0,  0,  -0.000695202291315483,  -0.000360629837797541,  0.000338853376305810,  -1.51416978260189e-06,  1.30937023657872e-06},
  //                  {0,  0,  -1.22458556833818e-20,  0,  1.22446311641411e-20,  4.99950003618466e-09,  0,  0,  2.35523254111355e-08,  1.53305901852412e-08,  -2.41276499095706e-08,  2.00150693917960e-10,  -3.62423865655675e-10},
  //                  {0, 0,  0,  -9.99950001249979e-05,  0,  0,  9.99900006666277e-05,  1.22446311641411e-20,  7.76025778431263e-05,  -8.54778582983768e-05,  -0.000232480528166871,  -9.81957550922941e-07,  5.22289220968140e-06},
  //                  {0,  0,  0,  -1.22458556833818e-20,  0,  0,  1.22446311641411e-20,  4.99950003618466e-09,  -4.86106774643752e-08,  -3.18405692657174e-08,  4.99548838464105e-08,  -5.18155515301648e-10,  7.64703789847301e-10},
  //                  {0.00396398441456158,  -0.00159529303278272,  0.000695237051140377,  -7.76064579396836e-05,  -0.000695202291315483,  2.35523254111355e-08,  7.76025778431263e-05,  -4.86106774643752e-08,  0.0387302923748964,  0.0240626171502365,  -0.0358916719245596,  0.000534912279490187,  -0.000602148232016092},
  //                  {0.000870810561497663,  -0.00265238798005091,  0.000360647869139166,  8.54821321556754e-05,  -0.000360629837797541,  1.53305901852412e-08,  -8.54778582983768e-05,  -3.18405692657174e-08,  0.0240626171502365,  0.0160327883644741,  -0.0224190734281786,  0.000266013454945774,  -0.000457600797389232},
  //                  {-0.00169959311610240,  -0.000269791083799460,  -0.000338870318833434,  0.000232492152096411,  0.000338853376305810,  -2.41276499095706e-08,  -0.000232480528166871,  4.99548838464105e-08,  -0.0358916719245596,  -0.0224190734281786,  0.0395843030277243,  -0.000411299056830157,  0.000469579903891527},
  //                  {0.000516042564142335,  -3.67677463138513e-05,  1.51424549046011e-06,  9.82006648391333e-07,  -1.51416978260189e-06,  2.00150693917960e-10,  -9.81957550922941e-07,  -5.18155515301648e-10,  0.000534912279490187,  0.000266013454945774,  -0.000411299056830157,  9.14372966945713e-05,  -1.39316775998304e-05},
  //                  {-0.000144426919017256,  0.000521672148730322,  -1.30943570454497e-06,  -5.22315335211565e-06,  1.30937023657872e-06,  -3.62423865655675e-10,  5.22289220968140e-06,  7.64703789847301e-10,  -0.000602148232016092,  -0.000457600797389232,  0.000469579903891527,  -1.39316775998304e-05,  0.000167805742964902}
  //                };
  //printMat(S);
  //printVec(matrixMultiplyC_SP(P_sp, M));
  //printMat(matrixMultiply_SP(matrixMultiply_SP(P_sp, S), transMat_SP(P_sp)));

  //sparse_t PP_sp = new sparse_t(kron(P, P));
  //PP_sp.print();

  mode_simulation = 1; //<>//
  numbIter         = 0;
  int countIt = 0;
  get_rollout = rollout(gaussian(mu0, S0), H, plant, cost);
  println("Эпизод № " + countIt++);
  xModel =get_rollout.xLock; //coeffProdMat(1.0, get_rollout.xLock);
  yModel =get_rollout.yLock; //coeffProdMat(1.0, get_rollout.yLock);

  printMat(xModel);
  printMat(yModel);

  stateW[0] = 0;
  stateW[1] = 3.1;

  indexH = 0;

  policy.p = let_p(policeInit(policy, 200));

  //policy.p.hyp = vectorCoef(hypP, 1.0);
  //policy.p.targets = vectorCoef(targP, 1.0);
  //policy.p.inputs = coeffProdMat(1.0, inpP);

  mode = 1;
}

public void ПРОЧИТАТЬ(int theValue) {

  readData();
  println("Данные прочитаны");
}


void draw()
{

  if (mode_simulation == 1)
  {
    switch(mode)
    {
    case 0:
      deltaTime = millis() - currentTime;
      if (deltaTime > 100)
      {
        currentTime = millis();
        background(200);

        //u = PID(stateW);
        double stat_temp[] = new double [policy.poli.length];
        stat_temp[0] = stateW[0];
        stat_temp[1] = stateW[1];
        stat_temp[2] = sin_db(stateW[2]);
        stat_temp[3] = cos_db(stateW[2]);
        stat_temp[4] = sin_db(stateW[3]);
        stat_temp[5] = cos_db(stateW[3]);
        double zer[][] = new double[policy.poli.length][policy.poli.length];
        nargout_conpols_t con = conCat(policy, stat_temp, zer);
        u = con.M[0];
        stateW = simulate(stateW, u, plant);

        pgBR.beginDraw();
        pgBR.background(255);
        draw_Pendubot(pgBR, stateW, u, plant);
        pgBR.endDraw();
        image(pgBR, 350, 25);
      }
      break;
    case 1:

      deltaTime = millis() - currentTime;
      if (deltaTime > 100)
      {
        currentTime = millis();

        pgBR.beginDraw();
        pgBR.background(255);
        draw_rollout(pgBR, get_rollout.latentLock, plant, indexH);
        pgBR.endDraw();
        image(pgBR, 350, 25);

        indexH++;
        if (indexH >= H) mode = 2;
      }

      break;

    case 2:
      if (!isLearning)
      {
        isLearning = true;

        final double[][] xM = xModel;
        final double[][] yM = yModel;

        new Thread(new Runnable() {
          public void run() {
            try {
              numbIter++;
              if (numbIter >= N)
              {
                mode = 4;
              } else
              {
                println("=== Итерация " + numbIter + " ===");

                trainDynModel(xM, yM);
                println("trainDynModel - OK");

                dynmodel = gp0(dynmodel);
                println("gp0 - OK");

                learnPolicy();
                println("learnPolicy - OK");
                double[] mu0L = {0, 0, PI, PI};
                double[][] S0L = {
                  {0.01, 0, 0, 0},
                  {0, 0.01, 0, 0},
                  {0, 0, 0.0001, 0},
                  {0, 0, 0, 0.0001}
                };
                nargout_valueSh_t vRes = valueSh(
                  let_p(policy.p), mu0L, S0L, dynmodel, policy, cost, H
                  );
                for (int t = 0; t < H; t++)
                    costHistory[numbIter][t] = (float)vRes.L[t];

                mode = 3;
              }
            }
            catch(Exception e) {
              println("ОШИБКА: " + e.toString());
              e.printStackTrace();
              mode = 10;
            }
            finally {
              isLearning = false;
            }
          }
        }
        ).start();
      }
      break;

    case 3:
      try {
        println("=== applyController ===");
        applyController();
        println("applyController - OK");
        println("get_rollout.xLock: " + get_rollout.xLock.length + " x " + get_rollout.xLock[0].length);
        println("get_rollout.yLock: " + get_rollout.yLock.length + " x " + get_rollout.yLock[0].length);
        println("xModel before merge: " + xModel.length + " x " + xModel[0].length);
        println("yModel before merge: " + yModel.length + " x " + yModel[0].length);
        for (int t = 0; t < H; t++)
          costRealHistory[numbIter][t] = (float)get_rollout.L[t];
        indexH = 0;
        double xx2[][] = xModel;//coeffProdMat(1.0, xModel);
        double yy2[][] = yModel;//coeffProdMat(1.0, yModel);
        xModel = new double[xx2.length + get_rollout.xLock.length][get_rollout.xLock[0].length];

        for (int i = 0; i < xx2.length; i++)
          for (int j = 0; j < xx2[0].length; j++)
          {
            xModel[i][j] = xx2[i][j];
          }
        for (int i = 0; i < get_rollout.xLock.length; i++)
          for (int j = 0; j < get_rollout.xLock[0].length; j++)
          {
            xModel[i + xx2.length][j] = get_rollout.xLock[i][j];
          }

        yModel = new double[yy2.length + get_rollout.yLock.length][get_rollout.yLock[0].length];
        for (int i = 0; i < yy2.length; i++)
          for (int j = 0; j < yy2[0].length; j++)
          {
            yModel[i][j] = yy2[i][j];
          }
        for (int i = 0; i < get_rollout.yLock.length; i++)
          for (int j = 0; j < get_rollout.yLock[0].length; j++)
          {
            yModel[i + yy2.length][j] = get_rollout.yLock[i][j];
          }

        println("xModel after merge: " + xModel.length + " x " + xModel[0].length);
        println("yModel after merge: " + yModel.length + " x " + yModel[0].length);

        mode = 1;
      }
      catch(Exception e) {
        println("ОШИБКА в case 3: " + e.toString());
        e.printStackTrace();
        mode = 10;
      }
      break;

    case 4:
      printMatFor(policy.p.inputs);
      printVecFor(policy.p.targets);
      printVecFor(policy.p.hyp);
      printMatFor(policyGP.beta);

      mode = 10;
      break;
    }
  } else if (mode_simulation == 2)
  {
    deltaTime = millis() - currentTime;
    if (deltaTime > 100)
    {
      currentTime = millis();
      stateW = simulate(stateW, u, plant);

      pgBR.beginDraw();
      pgBR.background(255);
      draw_Pendubot(pgBR, stateW, 0, plant);
      pgBR.endDraw();
      image(pgBR, 350, 25);
      pgCost.beginDraw();
      pgCost.background(255);
      draw_cost_graph(pgCost, costHistory, costRealHistory, costIterCount, H);
      pgCost.endDraw();
      image(pgCost, 0, 450); // разместить под основным графиком
    }
  }
}
