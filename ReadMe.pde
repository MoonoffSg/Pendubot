НЕ РАБОТАЕТ/ DOESNT WORK

//Изменения.
// trig, trigF
//1 cost : p {0.5, 0.5}, target[], angle = {2, 3}
//1. Модель в definitions и отображение в models.
//dyni[]  = {0, 2, 3};
//2. Размеры массивов на 2 в variables 
//int          difi[] = {0, 1},
//             angle    = 1;
//3. rollout в basa
//4. Изменение в cost в definitions
//5. loss_pendulum в loss и в valueSh valueSh2 value basa, 
//6. trainDynModel 
//7. learnPolicy
//8. applyController
//9. policeInit в control
//10. Закоментировать   //policy.p.hyp = vectorCoef(hypP, 1.0);
  //policy.p.targets = vectorCoef(targP, 1.0);
  //policy.p.inputs = coeffProdMat(1.0, inpP);
//11. util pol 
//opt.lgth          = 75;
//        opt.MFEPLS        = 30;
//        opt.MSR           = 100;
//        opt.SIG           = 0.5;
//        opt.fh            = 1;
