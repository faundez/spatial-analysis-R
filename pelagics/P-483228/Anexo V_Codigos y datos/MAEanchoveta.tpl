//###########################################################################
// MODELO ANUAL EN EDADES ANCHOVETA V-X REGIONES
//###########################################################################

//###########################################################################

DATA_SECTION

//###########################################################################

//LEER DATOS "MAEanchoveta.dat" 
  init_int nanos  
  init_int nedades
  init_int ntallas
  init_vector edades(1,nedades)
  init_vector Tallas(1,ntallas)
  init_number M
  init_vector dt(1,4)
  init_vector msex(1,nedades)
  init_matrix matdat(1,nanos,1,9)
  init_matrix Ctot(1,nanos,1,nedades)
  init_matrix Ccru_a(1,nanos,1,nedades)
  init_matrix Ccru_pel(1,nanos,1,nedades)
  init_matrix Ccru_l(1,nanos,1,ntallas)
  init_matrix Wmed(1,nanos,1,nedades)
  init_matrix Win(1,nanos,1,nedades)
  init_matrix error_edad(1,nedades,1,nedades)
  int reporte_mcmc
//==============================================================================
// LEER controles y opciones
!! ad_comm::change_datafile_name("MAEanchoveta.ctl");
//==============================================================================
//1.Coeficientes de variación
  init_number sigmaR
  init_number cvpriorq_reclas
  init_number cvpriorq_pelaces
//2.Tamaños de muestra
  init_vector nmus(1,4)
//3.Fases de selectividad
  init_int    opt1_fase
  init_int    opt_Scru1
//4.Fases de capturabilidad
  init_int    opt_qrecl
  init_int    opt_qpela
  init_int    opt_qmph
//5.Parámetros de crecimiento 
  init_vector pars_Bio(1,5)
//6.Fases de estimación Lo y cv edad
  init_int    opt_VB
//7.Fase de estimación de M
  init_int    opt_M
//8.Fase de estimación condiciones iniciales
  init_int    opt_Ro
  init_int    opt_devR
  init_int    opt_devNo
//9.Fase de estimación de F
  init_int    opt_F
//10.Considera la matriz de asignación de error edad
  init_number erredad
// 11. PROYECCIÓN, CRITERIOS DE EXPLOTACIÓN Y ESCENARIOS DE RECLUTAMIENTO
// Número de puntos biológicos
  init_int    npbr
// Número de años de proyección
  init_int    nproy
// Define desde donde proyecto
  init_number opt_Proy
// Fpbr para el cálculo de CTP
  init_vector pbr(1,npbr) 
//%SSB
  init_number RPRmsy 
// Estrategia de explotación
  init_number Fmsy 
  init_int opt_Str 

//###########################################################################

INITIALIZATION_SECTION

//###########################################################################

// defino un valor inicial de log_reclutamiento promedio (factor de escala)
  log_Ro     12.54
  log_Lo      2
  log_cv_edad -2.52
  log_M       0
  log_qrecl   0
  log_qpela   0

//###########################################################################

PARAMETER_SECTION

//###########################################################################

// Selectividad paramétrica
 init_bounded_number A50f(-1,2,opt1_fase)  
 init_bounded_number log_rangof(-4,0,opt1_fase)

 init_bounded_number A50c(-1,2,opt_Scru1)  
 init_bounded_number log_rangoc(-4,0.6,opt_Scru1)

 init_bounded_number A50pela(-1,2,opt1_fase)  
 init_bounded_number log_rangopela(-4,0.6,opt1_fase)

// Parámetros reclutamientos y mortalidades
 init_bounded_number log_Ro(5,20,opt_Ro)
 init_bounded_vector log_desv_No(1,nedades-1,-10,10,opt_devNo)
 init_bounded_dev_vector log_desv_Rt(1,nanos,-10,10,opt_devR)
 init_bounded_vector log_Ft(1,nanos,-6,1.6,opt_F) // log  mortalidad por pesca por flota

// Capturabilidades
 init_number log_qrecl(opt_qrecl)
 init_number log_qpela(opt_qpela)
 init_number log_qmph(opt_qmph)

// Crecimiento y M
 init_bounded_number log_Lo(1,2.1,opt_VB)
 init_bounded_number log_cv_edad(-4,-1.8,opt_VB)
 init_bounded_number log_M(-0.3,0.4,opt_M)

//###########################################################################
//
// VARIABLES DE ESTADO
//
//###########################################################################  
 
  vector edad_rel(1,nedades);  //no se para que es?????
  // arreglos
  vector anos(1,nanos);
  vector Unos_edad(1,nedades);
  vector Unos_anos(1,nanos);
  vector Unos_tallas(1,ntallas);
  //==============================================
  // 1. SELECTIVIDAD
  //==============================================
  matrix Sel_f(1,nanos,1,nedades)
  matrix Scru(1,nanos,1,nedades)
  matrix Scru_pela(1,nanos,1,nedades)
  //==============================================
  // 2. PROBABILIDAD EDAD-TALLA
  //==============================================
  number Linf
  number k
  number Lo
  number cv_edad
  vector mu_edad(1,nedades);
  vector sigma_edad(1,nedades);
  matrix Prob_talla(1,nedades,1,ntallas)
  matrix P1(1,nedades,1,ntallas)
  matrix P2(1,nedades,1,ntallas)
  matrix P3(1,nedades,1,ntallas)
  //==============================================
  // 3. MORTALIDADES
  //==============================================
  number M
  matrix Ftot(1,nanos,1,nedades)
  matrix Z(1,nanos,1,nedades)
  matrix S(1,nanos,1,nedades)
  //==============================================
  // 4. ABUNDANCIA Y BIOMASAS
  //==============================================
  // 4.1. Abundancia inicial en equilibrio
  vector Neq(1,nedades);
  likeprof_number SSBo;
  // 4.2. Abundancia
  matrix N(1,nanos,1,nedades)
  matrix NM(1,nanos,1,nedades)
  // 4.3. Matrices y vectores abundancia derivadas
  matrix NVflo(1,nanos,1,nedades)
  matrix NVcru(1,nanos,1,nedades)
  matrix NVpel(1,nanos,1,nedades)
  matrix NVpel_l(1,nanos,1,ntallas)
  matrix NMD(1,nanos,1,nedades)

  sdreport_vector Reclutas(1,nanos);
  // 4.4. Vectores de biomasas derivadas
  sdreport_vector BD(1,nanos)
  sdreport_vector BT(1,nanos) 
  vector BMflo(1,nanos)
  vector BMpel(1,nanos)
  vector Bcru(1,nanos)
  //==============================================
  // 6. CAPTURAS observadas y predichas
  //==============================================
  // 6.1. Matrices de capturas predichas por edad y años
  matrix pred_Ctot(1,nanos,1,nedades)
  // 6.2. Matrices de proporción de capturas por edad y años
  matrix pobs_f(1,nanos,1,nedades)             
  matrix ppred_f(1,nanos,1,nedades)
  matrix pobs_crua(1,nanos,1,nedades)
  matrix ppred_crua(1,nanos,1,nedades)
  matrix pobs_pel(1,nanos,1,nedades)
  matrix ppred_pel(1,nanos,1,nedades)
  matrix pobs_crul(1,nanos,1,ntallas)
  matrix ppred_crul(1,nanos,1,ntallas)
  //=================================================
  // 7.INDICES DE ABUNDANCIA observadas y predichas
  //=================================================
  likeprof_number qrecl
  likeprof_number qpela
  vector Reclas(1,nanos);
  sdreport_vector Reclas_pred(1,nanos);
  vector Pelaces(1,nanos);
  sdreport_vector Pelaces_pred(1,nanos);
  vector MPH(1,nanos);
  sdreport_vector MPH_pred(1,nanos);
  vector Desemb(1,nanos);
  vector Desemb_pred(1,nanos);
  //=================================================
  // 8. REDUCCIÓN DE LA BIOMASA DESOVANTE
  //=================================================
  matrix Nv(1,nanos,1,nedades)
  matrix NDv(1,nanos,1,nedades)
  vector BDo(1,nanos);
  vector RPR(1,nanos)  
  sdreport_vector RPR2(1,nanos)
  //=================================================  
  // 9. LOGVEROSIMILITUD
  //=================================================
  matrix cvar(1,4,1,nanos)
  vector likeval(1,15);
  objective_function_value f
  //=================================================
  // 10. CAPTURA BIOLOGICAMENTE ACEPTABLE
  //=================================================
  number Fref
  vector Fpbr(1,nedades);
  vector Zpbr(1,nedades);
  vector CTP(1,nedades);
  //=================================================
  // 11.  PROYECCION
  //=================================================
  vector Sp(1,nedades);
  vector Np(1,nedades);
  vector Nvp(1,nedades)
  vector CTPp(1,nedades);
  matrix Yp(1,npbr,1,nproy)
  matrix Bp(1,npbr,1,nproy)
  matrix BDp(1,npbr,1,nproy)
  matrix BDvp(1,npbr,1,nproy)
  number RPRp
  matrix matFpbr(1,npbr,1,nproy)
  sdreport_vector YTP(1,npbr)
  sdreport_vector YTPp(1,npbr)

//###########################################################################  
// Estima nm y CV
//===========================================================================
  number suma1
  number suma2
  number suma3
  number suma4
  number nm1
  number nm2
  number nm3
  number nm4
  number cuenta1
  number cuenta2
  number cuenta3
  number cuenta4

//###########################################################################  

PRELIMINARY_CALCS_SECTION

//###########################################################################  

// LEE matriz de indices
  anos=column(matdat,1);    
  Reclas=column(matdat,2);
  cvar(1)=column(matdat,3);
  Pelaces=column(matdat,4);
  cvar(2)=column(matdat,5);
  MPH=column(matdat,6);
  cvar(3)=column(matdat,7);
  Desemb=column(matdat,8);
  cvar(4)=column(matdat,9);

// ARREGLOS para operaciones matriciales
  Unos_edad=1;       // con la edad
  Unos_anos=1;       // con el año
  Unos_tallas=1;     // con la talla
  reporte_mcmc=0;

//###########################################################################  

RUNTIME_SECTION

//###########################################################################  

  maximum_function_evaluations 200,1000,5000
  convergence_criteria  1e-3,1e-5,1e-6

//###########################################################################  

PROCEDURE_SECTION

//###########################################################################  

// FUNCIONES
  Eval_selectividad_logis();
  Eval_prob_talla_edad();
  Eval_mortalidades();
  Eval_abundancia();
  Eval_biomasas();
  Eval_capturas_predichas();
  Eval_indices();
  Eval_deinteres();
  Eval_logverosim();
  Eval_funcion_objetivo();
  Eval_CTP();
  Eval_mcmc();

//========================================================================== 

FUNCTION Eval_selectividad_logis

//==========================================================================
  int i;


    Sel_f=outer_prod(Unos_anos,(elem_div(Unos_edad,(1+exp(-1.0*log(19)*(edades-A50f)/exp(log_rangof))))));
    Scru_pela=outer_prod(Unos_anos,(elem_div(Unos_edad,(1+exp(-1.0*log(19)*(edades-A50pela)/exp(log_rangopela))))));    
   if (opt_Scru1>0){// evaluo si el indice es >0
    Scru=outer_prod(Unos_anos,(elem_div(Unos_edad,(1+exp(-1.0*log(19)*(edades-A50c)/exp(log_rangoc))))));}    
   else{
    Scru=1.;}


//==========================================================================

FUNCTION Eval_prob_talla_edad

//==========================================================================

// se supone proporcionalidad entre la integral de una pdf y la pdf

  Linf=pars_Bio(1);
  k=pars_Bio(2);
  if(opt_VB<0){
  Lo=pars_Bio(3);
  cv_edad=pars_Bio(4);}
  else{
  Lo=exp(log_Lo);
  cv_edad=exp(log_cv_edad);}


 int i, j;
 
// genero una clave edad-talla para otros calculos. Se modela desde L(1)
  mu_edad(1)=Lo;
  for (i=2;i<=nedades;i++){
   mu_edad(i)=Linf*(1-exp(-k))+exp(-k)*mu_edad(i-1);}

  sigma_edad=cv_edad*mu_edad;

  P1=elem_div(outer_prod(Unos_edad,Unos_tallas),sqrt(2*3.1416)*outer_prod(sigma_edad,Unos_tallas));
  P2=mfexp(elem_div(-square(outer_prod(Unos_edad,Tallas)-outer_prod(mu_edad,Unos_tallas)),2*square(outer_prod(sigma_edad,Unos_tallas))));
  P3=elem_prod(P1,P2);
  Prob_talla=elem_div(P3,outer_prod(rowsum(P3),Unos_tallas));// normalizo para que la suma sobre las edades sea 1.0

//==========================================================================

FUNCTION Eval_mortalidades

//==========================================================================

  if(opt_M>0){
  M=exp(log_M);}
  else{
  M=pars_Bio(5);}

  Ftot=elem_prod(Sel_f,outer_prod(mfexp(log_Ft),Unos_edad));
  Z=Ftot+M;
  S=mfexp(-1.0*Z);

//==========================================================================

FUNCTION Eval_abundancia

//==========================================================================

 int i, j;
  
//**************************************************************************
// RECLUTAS anuales a la edad 2
//**************************************************************************
  for (i=1;i<=nanos;i++)
  {
  N(i,1)=mfexp(log_Ro+log_desv_Rt(i)+0.5*square(sigmaR)); 
  }
 
//**************************************************************************
// ABUNDANCIA Y BIOMASA INICIAL EN EQUILIBRIO
//**************************************************************************
  Neq(1)=exp(log_Ro+0.5*square(sigmaR));
  for (i=2;i<=nedades;i++)
  {
  Neq(i)=Neq(i-1)*exp(-1*M);
  }
  Neq(nedades)=Neq(nedades)/(1-exp(-1*M));

  SSBo=sum(elem_prod(Neq*exp(-dt(4)*M),elem_prod(msex,colsum(Wmed)/nanos)));

//**************************************************************************
// ABUNDANCIA INICIAL
//**************************************************************************
  for (i=2;i<=nedades;i++)
  {
  N(1)(i)=Neq(i)*exp(log_desv_No(i-1)+0.5*square(sigmaR));
  }
//**************************************************************************
// SOBREVIVENCIA  por edad(a+1) y año(t+1)
//**************************************************************************
  for (i=2;i<=nanos;i++)
  {
      N(i)(2,nedades)=++elem_prod(N(i-1)(1,nedades-1),S(i-1)(1,nedades-1));
      N(i,nedades)+=N(i-1,nedades)*S(i-1,nedades); 
  }

//==========================================================================

FUNCTION Eval_biomasas

//==========================================================================

// matrices y vectores de abundancias derivadas
  NVcru=elem_prod(elem_prod(N,mfexp(-dt(1)*Z)),Scru);// Crucero Reclas
  NVpel=elem_prod(elem_prod(N,mfexp(-dt(2)*Z)),Scru_pela);// Pelaces

//--------------------------------------------------------------------------
// corrección por error de asignación de la edad
  if(erredad==1){
  NVcru=NVcru*error_edad;
  NVpel=NVpel*error_edad;}
//--------------------------------------------------------------------------
 
  NMD=elem_prod(elem_prod(N,mfexp(-dt(3)*Z)),outer_prod(Unos_anos,msex));// desovante y MPH
  NVflo=elem_prod(elem_prod(N,mfexp(-dt(4)*Z)),Sel_f);// explotable
  Reclutas=column(N,1);

// vectores de biomasas derivadas
  BD=rowsum(elem_prod(NMD,Wmed));//Desovante
  BT=rowsum(elem_prod(N,Win));//Total inicios de año
  BMflo=rowsum(elem_prod(NVflo,Wmed));// Biomasa explotable
  BMpel=rowsum(elem_prod(NVpel,Wmed));// pelaces
  Bcru=rowsum(elem_prod(NVcru,Win));// Reclas


//==========================================================================

FUNCTION Eval_capturas_predichas

//==========================================================================

// matrices de capturas predichas por edad y año
  pred_Ctot=(elem_prod(elem_div(Ftot,Z),elem_prod(1.-S,N)));

//--------------------------------------------------------------------------
// corrección por error de asignación de la edad
  if(erredad==1){
  pred_Ctot=pred_Ctot*error_edad;}
//--------------------------------------------------------------------------

// vectores de desembarques predichos por año
  Desemb_pred=rowsum(elem_prod(pred_Ctot,Wmed));

// matrices de proporcion de capturas por edad y año
  pobs_f=elem_div(Ctot,outer_prod(rowsum(Ctot+1e-10),Unos_edad));
 // pobs_f=Ctot;
  ppred_f=elem_div(pred_Ctot,outer_prod(rowsum(pred_Ctot),Unos_edad));

// matrices de capturas predichas por talla y año
// RECLAS
  pobs_crua=elem_div(Ccru_a,outer_prod(rowsum(Ccru_a+1e-10),Unos_edad));
//  pobs_crua=Ccru_a;
  ppred_crua=elem_div(NVcru,outer_prod(rowsum(NVcru),Unos_edad));

// PELACES

  // EN EDADES
  pobs_pel=elem_div(Ccru_pel,outer_prod(rowsum(Ccru_pel+1e-10),Unos_edad));
 // pobs_pel=Ccru_pel;
  ppred_pel=elem_div(NVpel,outer_prod(rowsum(NVpel),Unos_edad));


  // EN TALLAS
  pobs_crul=elem_div(Ccru_l,outer_prod(rowsum(Ccru_l+1e-10),Unos_tallas));
  ppred_crul=elem_div(NVpel*Prob_talla,outer_prod(rowsum(NVpel),Unos_tallas));


//==========================================================================

FUNCTION Eval_indices

//==========================================================================
 
 int i;

  Reclas_pred=exp(log_qrecl)*Bcru;
  Pelaces_pred=exp(log_qpela)*BMpel;
  MPH_pred=exp(log_qmph)*BD;

  qrecl=exp(log_qrecl);
  qpela=exp(log_qpela);

//==========================================================================

FUNCTION Eval_deinteres

//==========================================================================

// Rutina para calcular RPR
  Nv=N;// solo para empezar los calculos
  
 for (int i=2;i<=nanos;i++){
      Nv(i)(2,nedades)=++Nv(i-1)(1,nedades-1)*exp(-1.0*M);
      Nv(i)(nedades)=Nv(i-1)(nedades)*exp(-1.0*M);}

  NDv=elem_prod(Nv*exp(-dt(3)*M),outer_prod(Unos_anos,msex));
  BDo=rowsum(elem_prod(NDv,Wmed));
  RPR=elem_div(BD,BDo);
  RPR2=BD/SSBo;
 
//==========================================================================

FUNCTION Eval_logverosim

//==========================================================================

// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
  int i;

  suma1=0;  suma2=0;  suma3=0;   suma4=0;

  for (i=1;i<=nanos;i++)
  {
   if (Reclas(i)>0){
    suma1+=square((log(Reclas(i))-log(Reclas_pred(i)))/cvar(1,i));}
   if (Pelaces(i)>0){
    suma2+=square((log(Pelaces(i))-log(Pelaces_pred(i)))/cvar(2,i));}
   if (MPH(i)>0){
    suma3+=square((log(MPH(i))-log(MPH_pred(i)))/cvar(3,i));}
   if (Desemb(i)>0){
    suma4+=square((log(Desemb(i))-log(Desemb_pred(i)))/cvar(4,i));}
  }

//==========================================================================

FUNCTION Eval_funcion_objetivo

//==========================================================================

// se calcula la F.O. como la suma de las -logver
// lognormal

  likeval(1)=0.5*suma1;//Reclas
  likeval(2)=0.5*suma2;//pelaces
  likeval(3)=0.5*suma4;//MPH
  likeval(4)=0.5*suma3;//Desemb

// multinomial
  likeval(5)=-nmus(1)*sum(elem_prod(pobs_f,log(ppred_f)));
  likeval(6)=-nmus(2)*sum(rowsum(elem_prod(pobs_crua,log(ppred_crua))));
  likeval(7)=-nmus(3)*sum(rowsum(elem_prod(pobs_pel,log(ppred_pel))));

// tallas del pelaces
  likeval(8)=-nmus(4)*sum(rowsum(elem_prod(pobs_crul,log(ppred_crul))));

//  Reclutas
  likeval(9)=1./(2*square(sigmaR))*norm2(log_desv_Rt);

// q cruceros

  likeval(10)=1./(2*square(cvpriorq_reclas))*square(log_qrecl);
  likeval(11)=1./(2*square(cvpriorq_pelaces))*square(log_qpela);

// total
   f=sum(likeval);
 
  if(mceval_phase()){
    Eval_mcmc();}

//==========================================================================

FUNCTION  Eval_CTP

//==========================================================================

// re-calcula la CTP para el último año dado los PBR entregados.

  for (int i=1;i<=npbr;i++){
  Fpbr=Sel_f(nanos)*pbr(i);//
  Zpbr=Fpbr+M;
  CTP=elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),N(nanos)));
  YTP(i)=sum(elem_prod(CTP,Wmed(nanos)));}
//  RPRp=RPR(nanos);

  for (int i=1;i<=npbr;i++){ // ciclo de PBR

  Np=N(nanos);
  Sp=S(nanos);
  RPRp=RPR2(nanos);
  Nvp=Nv(nanos);

  for (int j=1;j<=nproy;j++){ // ciclo de años

// calcula la CTP para los PBR
  Np(2,nedades)=++elem_prod(Np(1,nedades-1),Sp(1,nedades-1));
  Np(nedades)+=Np(nedades)*Sp(nedades);
  //Np(1)=exp(log_Ro+0.5*square(sigmaR));
  Np(1)=mean(Reclutas(nanos-6,nanos));
  //Np(1)=mean(Reclutas(1,nanos-6));
  
  Nvp(2,nedades)=++Nvp(1,nedades-1)*exp(-1.0*M);
  Nvp(nedades)+=Nvp(nedades)*exp(-1.0*M);
  //Nvp(1)=exp(log_Ro+0.5*square(sigmaR));
  Nvp(1)=mean(Reclutas(nanos-6,nanos));
  //Nvp(1)=mean(Reclutas(1,nanos-6));
  
  Fref=pbr(i);
  if(opt_Str==0){

// Regla de decisión------------
  if(RPRp/RPRmsy<1){
  Fref=pbr(i)*RPRp/RPRmsy;}
//------------------------------
  }

  Fpbr=Sel_f(nanos)*Fref;
  matFpbr(i,j)=Fref;

  Zpbr=Fpbr+M;
  CTPp=elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Np));
  BDp(i,j)=  sum(elem_prod(elem_prod(elem_prod(Np,mfexp(-dt(3)*Zpbr)),msex),colsum(Wmed)/nanos));
  BDvp(i,j)= sum(elem_prod(elem_prod(Nvp*mfexp(-dt(3)*M),msex),colsum(Wmed)/nanos));
  
  Yp(i,j)=sum(elem_prod(CTPp,colsum(Wmed)/nanos)); // captura proyectada
  Bp(i,j)=sum(elem_prod(Np,colsum(Win)/nanos));  // biomasa proyectada
  RPRp=BDp(i,j);                            // RPR proyectado
  Sp=exp(-1.*Zpbr);                           

  }// FIN AÑOS

  if(opt_Proy==1){// proyectando desde RECLAS/PELACES para el mismo año
  YTPp(i)=YTP(i); // 
  }

  if(opt_Proy==2){// proyectando de fin de año pal otro
  YTPp(i)=Yp(i,1); // 
  }
 }

//###########################################################################  

REPORT_SECTION

//###########################################################################  

//-------------------------------
// INDICES DE ABUNDANCIA
//-------------------------------
  report << "reclasobs" << endl;
  report << Reclas << endl;
  report << "reclaspred" << endl;
  report << Reclas_pred << endl;
  report << "pelacesobs" << endl;
  report << Pelaces << endl;
  report << "pelacespred" << endl;
  report << Pelaces_pred << endl;
  report << "mphobs" << endl;
  report << MPH << endl;
  report << "mphpred" << endl;
  report << MPH_pred << endl;
  report << "desembarqueobs" << endl;
  report << Desemb << endl;
  report << "desembarquepred" << endl;
  report << Desemb_pred << endl;
//--------------------------------
// RECLUTAMIENTO
//--------------------------------
  report << "Reclutas" << endl;
  report << Reclutas << endl;
  report << "log_desv_Rt " << endl;
  report << log_desv_Rt << endl;
//--------------------------------
// BIOMASAS Y RPRs
//--------------------------------
  report << "SSB" << endl;
  report << BD << endl; 
  report << "BT" << endl;  
  report << BT << endl;
  report << "RPR_dinamico" << endl;
  report << RPR << endl;
  report << "RPR_equilibrio" << endl;
  report << RPR2 << endl;
//--------------------------------
// MORTALIDAD POR PESCA TOTAL
//--------------------------------
  report << "F" << endl;
  report << Ftot << endl;
//--------------------------------
// MATRIZ DE SOBREVIVENCIA
//--------------------------------
  report << "N " << endl;
  report << N << endl;
//--------------------------------
// SELECTIVIDADES
//--------------------------------
  report << "S_f" << endl;
  report << Sel_f << endl;
  report << "Scru_reclas" << endl;
  report << Scru << endl;
  report << "Scru_pelaces" << endl;
  report << Scru_pela << endl;
//--------------------------------
// PROPORCION DE LAS CAPTURAS
//--------------------------------
  report << "pf_obs " << endl;
  report << pobs_f << endl;
  report << "pf_pred " << endl;
  report << ppred_f << endl;
  report << "pobs_RECLAS" << endl;
  report << pobs_crua << endl;
  report << "ppred_RECLAS" << endl;
  report << ppred_crua << endl;
  report << "pobs_PELACES" << endl;
  report << pobs_pel << endl;
  report << "ppred_PELACES" << endl;
  report << ppred_pel << endl;
  report << "pobs_pel_tallas" << endl;
  report << pobs_crul << endl;
  report << "ppred_pel_tallas" << endl;
  report << ppred_crul << endl;
//--------------------------------
// VEROSIMILITUD
//--------------------------------
  report << "likeval Reclas Pelaces MPH Desemb pf preclas ppelaces ptallas desvR qrecl qpela" << endl;
  report << likeval << endl;
//--------------------------------
// CAPTURABILIDADES
//--------------------------------
  report << "q_reclas  q_pelaces" << endl;
  report << exp(log_qrecl)<<" "<<exp(log_qpela)<< endl;
  report << "M"<<endl;
  report << M << endl;
//--------------------------------
// CAPTURA BIOLOGICAMENTE ACEPTABLE
//--------------------------------
// Proyección de la captura y biomasa desovante
  report << "Capturas_proy" << endl;
  report << Yp<< endl;
  report << "BD_proy" << endl;
  report << BDp<< endl;
// Biomasa desovante sin pesca dinámica y equilibrio
  report << "BD0" << endl;
  report << BDo<< endl;
  report << "SSBo" << endl;
  report << SSBo<< endl;
// Reducción de la biomasa desovante en equilibrio
  report << "RPR_proy" << endl;
  report << BDp/SSBo<< endl;
  report << "Rprom"<<endl;
  report << Reclutas(nanos-6,nanos)<<endl;
// Criterio de explotación
  report << "Fpbr_proy" << endl;
  report << matFpbr<< endl;
// límites de biomasa y F relativos (diagrama fase)
  report << "B/Bmsy" << endl;
  report << RPR2/RPRmsy<< endl;
  report << "F/Fmsy" << endl;
  report << column(Ftot,max(edades))/Fmsy<< endl;
  report << "np"<<endl;
  report << Np <<endl;
 
//##########################################################################
//--------------------------------------------------------------------------

// ESTIMA nm y CV

//--------------------------------------------------------------------------
//##########################################################################

  suma1=0; suma2=0;nm1=1;cuenta1=0;cuenta2=0;

  for (int i=1;i<=nanos;i++){ //

   if (sum(pobs_f(i))>0){
      suma1=sum(elem_prod(ppred_f(i),1-ppred_f(i)));
      suma2=norm2(pobs_f(i)-ppred_f(i));
      nm1=nm1*suma1/suma2;
      cuenta1+=1;
   }}

  suma1=0; suma2=0;nm2=1;cuenta2=0;

  for (int i=1;i<=nanos;i++){ //

   if (sum(pobs_crua(i))>0){
      suma1=sum(elem_prod(ppred_crua(i),1-ppred_crua(i)));
      suma2=norm2(pobs_crua(i)-ppred_crua(i));
      nm2=nm2*suma1/suma2;
      cuenta2+=1;
   }}


  suma1=0; suma2=0;nm3=1;cuenta3=0;

  for (int i=1;i<=nanos;i++){ //

   if (sum(pobs_pel(i))>0){
      suma1=sum(elem_prod(ppred_pel(i),1-ppred_pel(i)));
      suma2=norm2(pobs_pel(i)-ppred_pel(i));
      nm3=nm3*suma1/suma2;
      cuenta3+=1;
   }}


  suma1=0; suma2=0;nm4=1;cuenta4=0;

  for (int i=1;i<=nanos;i++){ //

   if (sum(pobs_crul(i))>0){
      suma1=sum(elem_prod(ppred_crul(i),1-ppred_crul(i)));
      suma2=norm2(pobs_crul(i)-ppred_crul(i));
      nm4=nm4*suma1/suma2;
      cuenta4+=1;
   }}

  report << "nm_flota  nm_reclas  nm_pelaces    nm_pelaL" << endl;
  report<<pow(nm1,1/cuenta1)<<" "<<pow(nm2,1/cuenta2)<<" "<<pow(nm3,1/cuenta3)<<" "<<pow(nm3,1/cuenta4)<<endl;

  suma1=0;  suma2=0;  suma3=0;   suma4=0;  cuenta1=0;     cuenta2=0;   cuenta3=0;

  for (int i=1;i<=nanos;i++)
  {
   if (Reclas(i)>0){
    suma1+=square(log(Reclas(i))-log(Reclas_pred(i)));
    cuenta1+=1;}
   if (Pelaces(i)>0){
    suma2+=square(log(Pelaces(i))-log(Pelaces_pred(i)));
    cuenta2+=1;}
   if (MPH(i)>0){
    suma4+=square(log(MPH(i))-log(MPH_pred(i)));
   cuenta3+=1;}
  }


 report << "cv_recla  cv_pelaces  cv_mph" << endl;
 report<<sqrt(suma1/cuenta1)<<" "<<sqrt(suma2/cuenta2)<<" "<<sqrt(suma4/cuenta3)<<endl;
 

  if(erredad==1){
  report << "------------------------------------------------" << endl;
  report << " Matriz de error "<< endl;
  report << error_edad << endl;}

//==========================================================================

FUNCTION Eval_mcmc

//==========================================================================
  if(reporte_mcmc == 0)
  mcmc_report<<"f, RPR, CTP_66, CTP_60, CTP_50, Reclas, Pelaces"<<endl;
  mcmc_report<<f<<","<<RPR(nanos)<<","<<YTP(1)<<","<<YTP(2)<<","<<YTP(3)<<","<<Reclas_pred(nanos)<<","<<Pelaces_pred(nanos)<<endl;
  reporte_mcmc++;

//###########################################################################

GLOBALS_SECTION

//###########################################################################

  #include  <admodel.h>
  ofstream mcmc_report("mcmc.csv");


