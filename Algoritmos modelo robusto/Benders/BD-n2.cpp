/*********************************************
 * Concert Model
 * Author: Fabricio
 * Creation Date: 05-06-2020
 * B&B tree: one
 * Number of cut: root node 
 * Type of frac. cut: force
 *********************************************/


#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ilcplex/ilocplex.h>

using namespace std;


#define TimeLimit 21600 //36000
    #define EPSILON 1e-07
#define M 10e10 /// Parâmetro usado para representar um número suficientemente grande

ILOSTLBEGIN
double best_inc = 1.e30;
float alpha1;
float receita;

inline  int ij(int n, int i, int j){ /// Função usada para recuperar os índices das variáveis de decisão do suproblema vbeta_ij, vgama_ij
  return n*i + j;
}

inline  int kl(int n, int k, int l){ /// Função usada para recuperar os índices da variável de decisão do problema mestre z_kl, l!=k
  if(k<l){
    return (n-1)*k+l-1; 
  }else if(k>l){
    return (n-1)*k+l; 
  }
}

/// ==============================================
/// Benders decompositon class
/// ==============================================

class BendersDecomposition{
    public:
        /// Parâmetros das instâncias
        int n; // Quantidade de nós
        float alpha; // Fator de desconto no custo de transporte em um link entre hubs
        vector<double >  codx; // Coordenada x dos nós (AP)
        vector<double >  cody; // Coordenada y dos nós (AP)
        vector<vector<double > > w; // declara um vetor de vetores para representar uma matriz com n linhas e n colunas - Quantidade de demanda a ser enviada entre os nós i e j
        vector<vector<double > > r; // Receita obtida por enviar uma unidade de demanda entre os nós i e j
        vector<vector<double > > c; // Custos por enviar uma unidade de demanda entre os nós i e j
        vector<double > s; // Custos fixos de instalação de um hub
        vector<vector<double > > g; // Custos de operação nos links inter-hubs
        float soma; // Soma para obter a média geral dos custos fixos de instalação hub - Dados AP
		
		/// Parâmetros do Benders
        int NCuts; /// Guardar o número de cortes gerados
        int NNodes; /// Guardar quantos nós da árvore de B&C foram expandidos
        double ub; /// Limitante superior
        double lb; /// Limitante inferior
        double gap; /// Gap entre o UB e LB
        bool Activate_Lazy_Callback; /// Parâmetro passado por linha de comando que indica se a LazyCalback estará ativa ou não
        bool Activate_User_Callback; /// Parâmetro passado por linha de comando que indica se a BendersUserCutsCallback estará ativa ou não
        
        /// Variáveis do problema mestre 
        IloEnv env;
        IloModel mod;
        IloCplex cplex;
        IloObjective fo;
        IloNumVarArray veta; /// Variável de decisão do mestre
        IloNumVarArray vh; /// Variável de decisão do mestre
        IloNumArray h; /// variável auxiliar para variável de decisão do mestre
        IloNumVarArray vz; /// Variável de decisão do mestre
        IloNumArray z; /// variável auxiliar para variável de decisão do mestre

        /// SP parameter 
        IloEnv envsp;
        IloModel modsp;
        IloCplex cplexsp;
        IloObjective fosp;
        IloNumVar vbeta; /// Variável de decisão do subproblema dual
        IloNumVar vgama;  /// Variável de decisão do subproblema dual
        IloNumVarArray vtheta; /// Variável de decisão do subproblema dual
        IloNumVarArray vmu;  /// Variável de decisão do subproblema dual
        IloNumVarArray vrho;  /// Variável de decisão do subproblema dual
        IloNumVarArray vdelta;  /// Variável de decisão do subproblema dual
        IloNum beta; /// variável auxiliar para variável de decisão do subproblema dual
        IloNum gama;  /// variável auxiliar para variável de decisão do subproblema dual
        IloNumArray theta; /// variável auxiliar para variável de decisão do subproblema dual
        IloNumArray mu;  /// variável auxiliar para variável de decisão do subproblema dual
        IloNumArray rho;  /// variável auxiliar para variável de decisão do subproblema dual
        IloNumArray delta;  /// variável auxiliar para variável de decisão do subproblema dual
      IloRangeArray rd1; /// Vetor que armazena as restrições do dual
		  IloNumArray rd1UB; /// Vetor que armazena o limitante superior de cada restrição (+infinito)
		  IloNumArray rd1LB;  /// Vetor que armazena o limitante inferior de cada restrição (w[i][j]*(r[i][j]-c[i][k]))
		  IloRangeArray rd2; /// Vetor que armazena as restrições do dual
		  IloNumArray rd2UB; /// Vetor que armazena o limitante superior de cada restrição (+infinito)
		  IloNumArray rd2LB;  /// Vetor que armazena o limitante inferior de cada restrição (-w[i][j]*c[l][j]))
		  IloRangeArray rd3; /// Vetor que armazena as restrições do dual
		  IloNumArray rd3UB; /// Vetor que armazena o limitante superior de cada restrição (+infinito)
		  IloNumArray rd3LB;  /// Vetor que armazena o limitante inferior de cada restrição (-alpha*w[i][j]*c[k][l]))
	
    public:
		BendersDecomposition(): NCuts(0), NNodes(0), ub(M), lb(0.0),gap(1.0), Activate_Lazy_Callback(1), Activate_User_Callback(1){} /// Construtor da classe		
        void Read_data(char *arq);
        void Create_master(); /// Cria o problema mestre 
        void Create_SP(); /// Cria o subproblema de Benders (SPD)
        void GeraCortes(IloExprArray OptCut); /// Resolve o SPD e guarda as informações do corte na expressão OptCut
        void Display_configuration(); /// Imprime a solução corrente

};




/// ==============================================
/// Configuração da UserCutCallbackI
/// Adiciona corte de Benders a partir de soluções relaxadas do mestre
/// ==============================================

class BendersUserCutsCallback : public IloCplex::UserCutCallbackI, public BendersDecomposition {
	private:
        BendersDecomposition *bd;
    public:
        void main();
		BendersUserCutsCallback(IloEnv env, BendersDecomposition *xbd) :IloCplex::UserCutCallbackI(env),  bd(xbd){}
        IloCplex::CallbackI* duplicateCallback() const {
            return (new (getEnv()) BendersUserCutsCallback(*this));
        }
};
 
void BendersUserCutsCallback::main(){ 
	bd->NNodes = getNnodes();
    if(getNnodes() == 0){ /// Corte adicionados apenas no nó raiz
        double inc_value = getObjValue();
        if(best_inc != inc_value){ /// Evita gerar cortes para um mesma solução
			cout<<"."; /// Imprime na saída padrão um símbolo para mostrar que um corte foi gerado
            best_inc = inc_value;
   
			/// ========================
			/// Recupera o valor das variáveis do mestre que serão usados para resolver o SPD
			/// ========================
			
			getValues(bd->h, bd->vh);
			getValues(bd->z, bd->vz); 

			/// ========================
			/// Adiciona um corte de otimalidade
			/// ========================
			IloExprArray OptCut(getEnv()); /// Cria uma expressão chamada OptCut
			bd->GeraCortes(OptCut); /// Chama a função GeraCortes que resolve o SPD e guarda o corte gerado na expressão OptCut
			for(int i = 0; i < bd->n; i++){
			  for(int j = 0; j < bd->n; j++){
			    add(OptCut[ij(bd->n, i, j)] <= 0); /// Adiciona o corte de Benders relativo a expressão OptCut
			    bd->NCuts++; 
			  }
			}
			OptCut.endElements();
					
		}
    }
}

IloCplex::Callback getUserCutCallback(IloEnv env, BendersDecomposition *bd){// callback handle function
    return (IloCplex::Callback(new (env) BendersUserCutsCallback(env, bd)));
}
 

/// ==============================================
/// Configuração da LazyConstraintCallbackI
/// Adiciona corte de Benders a partir de soluções relaxadas do mestre
/// ==============================================

class BendersLazyCutCallback : public IloCplex::LazyConstraintCallbackI, public BendersDecomposition {
    private:
        BendersDecomposition *bd;
    public:
        void main();
		BendersLazyCutCallback(IloEnv env, BendersDecomposition *xbd) :IloCplex::LazyConstraintCallbackI(env),  bd(xbd){}
        IloCplex::CallbackI* duplicateCallback() const {
            return (new (getEnv()) BendersLazyCutCallback(*this));
        }
};

void BendersLazyCutCallback::main(){
	/// ========================
	/// Recupera as principais informações da árvore de B&C
	/// ========================	
	bd->NNodes = getNnodes(); /// Número de nós expandidos da árvore de B&C
  bd->ub = getBestObjValue(); /// Salva o valor da melhor limitante superior
  bd->lb = getIncumbentObjValue(); /// Guarda o valor da solução incumbent
  bd->gap = (bd->ub  - bd->lb )/bd->lb ; /// Atualiza o gap


	/// ========================
	/// Recupera o valor das variáveis do mestre que serão usados para resolver o SPD
	/// ========================
	
	getValues(bd->h, bd->vh);
	getValues(bd->z, bd->vz); 

	/// ========================
	/// Adiciona um corte de otimalidade
	/// ========================
	IloExprArray OptCut(getEnv()); /// Cria uma expressão chamada OptCut
	bd->GeraCortes(OptCut); /// Chama a função GeraCortes que resolve o SPD e guarda o corte gerado na expressão OptCut
	//bd->cplex.exportModel("mod1.lp");
	for(int i = 0; i < bd->n; i++){
	  for(int j = 0; j < bd->n; j++){
	    add(OptCut[ij(bd->n, i, j)] <= 0); /// Adiciona o corte de Benders relativo a expressão OptCut
	    bd->NCuts++; 
	  }
	}
	//bd->cplex.exportModel("mod2.lp");
	OptCut.endElements();
	cout<<"*"; /// Imprime na saída padrão um símbolo para mostrar que um corte foi gerado
}


IloCplex::Callback getLazyConstraintCallback(IloEnv env, BendersDecomposition *bd){// callback handle function
    return (IloCplex::Callback(new (env) BendersLazyCutCallback(env, bd)));
}

/// ==============================================
/// main program
/// ==============================================
void Read_main_arg(int argc, char *argv[],bool &Activate_Lazy_Callback, bool &Activate_User_Callback); /// Função que recebe os principais argumentos da entrada padrão
int main (int argc, char *argv[]){

    // ===============
    // bd variables
    // ===============

    BendersDecomposition *bd = new BendersDecomposition();

    try{

        Read_main_arg(argc, argv, bd->Activate_Lazy_Callback, bd->Activate_User_Callback);
        bd->Read_data(argv[1]);
        IloTimer crono(bd->env);

        //===========================
        // Creating Benders problems
        //========================== 
        bd->Create_master();		
        bd->Create_SP();

        // ==========================
        // cplex configurations
        // ==========================

       bd->cplex.setParam(IloCplex::Threads,1);
    //    bd->cplex.setParam(IloCplex::EpGap,0.00001);
   //     bd->cplex.setOut(bd->env.getNullStream());
   //     bd->cplexsp.setParam(IloCplex::Threads,1);
        bd->cplex.setParam(IloCplex::TiLim, TimeLimit);
        bd->cplex.setWarning(bd->env.getNullStream());
        bd->cplexsp.setWarning(bd->envsp.getNullStream());
        bd->cplexsp.setOut(bd->envsp.getNullStream());

        // 	bd->cplexsp.setParam(IloCplex::PreInd, false);
        // 	bd->cplex.setOut(cplexStream);
        //  bd->cplex.setParam(IloCplex::HeurFreq,-1);
        // 	bd->cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);
        // 	bd->cplexsp.setParam(IloCplex::RootAlg, IloCplex::Barrier);

        /// =======================
        /// Activating callback
        /// =======================

        
        if(bd->Activate_Lazy_Callback > 0){
            bd->cplex.use(getLazyConstraintCallback(bd->env, bd));
        }
        else{        
            bd->cplex.setOut(bd->env.getNullStream());
        }

        if(bd->Activate_User_Callback > 0){
            bd->cplex.use(getUserCutCallback(bd->env, bd));
        }
 

		int u = 0; /// Guardar o número de iterações
        crono.start();
 
        /// =======================
        ///  Loop principal do Benders 
        /// =======================
        while(bd->gap > EPSILON && bd->cplex.getTime() < TimeLimit){

			/// =======================
			///  Resolve-se o problema mestre 
			/// =======================
			bd->cplex.solve();
          
          if(bd->cplex.getTime() >= TimeLimit){
            if(bd->cplex.getStatus() == IloAlgorithm::Optimal){
              bd->ub = (double) bd->cplex.getObjValue(); /// Atualiza o UB
              bd->gap = (bd->ub  - bd->lb )/bd->lb; /// Atualiza o gap
              break;
            }
          }
			
			/// =======================
			///   Verifica se a Activate_Lazy_Callback está ativa ou não.
			///(i) Caso esteja ativa, haverá apenas uma iteração do método. 
			///Após o comando anterior já encontramos a solução ótima, pois
			///cortes de Benders são gerados dentro da árvore de B&C.
			///(ii) Caso não esteja ativada, haverá as iterações comuns de Benders.
			/// =======================
            if(bd->Activate_Lazy_Callback == 0){
			 	if(bd->cplex.getStatus() == IloAlgorithm::Optimal){
					bd->ub = (double) bd->cplex.getObjValue(); /// Atualiza o UB
			 	  bd->cplex.getValues(bd->vh,bd->h); /// Recupera a solução do mestre
			 	  bd->cplex.getValues(bd->vz,bd->z); /// Recupera a solução do mestre
					IloExprArray OptCut(bd->env); /// Cria uma expressão chamada OptCut
					bd->GeraCortes(OptCut); /// Chama a função GeraCortes que resolve o SPD e guarda o corte gerado na expressão OptCut
					for(int i = 0; i < bd->n; i++){
					  for(int j = 0; j < bd->n; j++){
					    bd->mod.add(OptCut[ij(bd->n, i, j)] <= 0); /// Adiciona o corte de Benders relativo a expressão OptCut
					  }
					}
					OptCut.endElements();
					bd->gap = (bd->ub  - bd->lb )/bd->lb ; /// Atualiza o gap
					u++;		      
					cout<<u<<"\t"<<bd->lb <<"\t"<<bd->ub <<"\t"<< (100*bd->gap)<<"\% \t"<<crono.getTime()<<"\t"<<bd->NCuts<<"\t"<<endl;
				}
				else{
					break;

				} 
			}
			else{
				bd->lb = (double) bd->cplex.getObjValue();
				bd->ub = (double) bd->cplex.getObjValue();
				bd->gap = 0.0;
				bd->cplex.getValues(bd->vh,bd->h);
			}
        }  
        //===========================


        //===========================
        // Print main output
        //===========================
        ofstream arq_saida("Resultados-bd2-n2.txt",std::ofstream::app);
		if (!arq_saida) {cerr << "Erro arquivo \n"; exit(0);}
		
        arq_saida<<argv[1]<<"\t"<<bd->alpha<<"\t"<<bd->r[0][0]<<"\t"<<u<<"\t"<<100 *  bd->gap<<"\t \t"<<bd->lb<<"\t"<<bd->ub<<"\t"<<(double) crono.getTime()<<"\t"<<bd->Activate_Lazy_Callback<<"\t"<<bd->Activate_User_Callback<<endl;
        cout<<argv[1]<<"\t"<<bd->alpha<<"\t"<<bd->r[0][0]<<"\t"<<u<<"\t"<<100 *  bd->gap<<"\t \t"<<bd->lb<<"\t"<<bd->ub<<"\t"<<(double) crono.getTime()<<"\t"<<bd->Activate_Lazy_Callback<<"\t"<<bd->Activate_User_Callback<<endl;

        bd->Display_configuration(); 
        delete bd;
    } 
    catch (IloException& exc) {
        cerr << "Error: " << exc << endl;
    }
    return 0;
}


void Read_main_arg(int argc, char *argv[],bool &Activate_Lazy_Callback, bool &Activate_User_Callback){
    Activate_Lazy_Callback = (argc > 2) ? atof(argv[2]) : 1; // Indica se as Lazy_callback estará ativa ou não
    Activate_User_Callback = (argc > 3) ? atof(argv[3]) : 1; // Indica se as User_cut_callback estará ativa ou não
    alpha1 = (argc > 4) ? atof(argv[4]) : 0.2;
    receita = (argc > 5) ? atof(argv[5]) : 20;
}


void  BendersDecomposition::Read_data(char name[]){
    ifstream arq(name);
	if (!arq) {cerr << "Erro arquivo \n"; exit(0);}
	
	arq >> n;
	
	codx = vector<double>(n); // Coordenada x dos nós (AP)
	cody = vector<double>(n); // Coordenada y dos nós (AP)
	w = vector<vector<double > > (n, vector<double>(n)); // declara um vetor de vetores para representar uma matriz com n linhas e n colunas - Quantidade de demanda a ser enviada entre os nós i e j
	r = vector<vector<double > > (n, vector<double>(n)); // Receita obtida por enviar uma unidade de demanda entre os nós i e j
	c = vector<vector<double > > (n, vector<double>(n)); // Custos por enviar uma unidade de demanda entre os nós i e j
	s = vector<double > (n); // Custos fixos de instalação de um hub
	g = vector<vector<double > > (n, vector<double>(n)); // Custos de operação nos links inter-hubs
	
	soma = 0; // Soma para obter a média geral dos custos fixos de instalação hub - Dados AP
	
	alpha = alpha1; //Fator de desconto entre nós inter-hubs AP
	
	for (int i = 0; i <  n; i++){ //Receita AP 20,30,50
	  for (int j = 0; j < n; j++){
	    r[i][j] = receita;
	  }
	}
	
	for (int i = 0; i <  n; i++){ //Coletar as coordenadas dos nós
	  arq >> codx[i];
	  arq >> cody[i];
	}
	
	for (int i = 0; i <  n; i++){ //Coletar custos AP
	  for (int j = 0; j < n; j++){
	    c[i][j] = 0.001 * sqrt((codx[i] - codx[j]) * (codx[i] - codx[j])  + (cody[i] - cody[j]) * (cody[i] - cody[j]));
	  }
	}
	
	for (int i = 0; i <  n; i++){ //Coletar demanda AP
	  for (int j = 0; j < n; j++){
	    arq >> w[i][j];
	  }
	}
	
	for (int i = 0; i <  n; i++){ //Coletar custos fixos instalação hub AP
	  arq >> s[i];
	  s[i] = s[i]/10;
	  soma = soma + s[i];
	}
	
	for (int i = 0; i <  n; i++){ //Custos de operação link entre hubs AP
	  for (int j = 0; j < n; j++){
	    g[i][j] = 0.1 * soma/n;
	  }
	}
	
    arq.close();
}

void BendersDecomposition::Create_master(){
	mod = IloModel(env); /// definição do modelo
	cplex = IloCplex(mod); /// definição do objetivo cplex
	
	///===================================================
	/// Definição das variáveis e parâmetros principais do mestre
	///===================================================
	float ls_veta = 0;
	for (int i = 0; i <  n; i++){ //Limitante superior para o veta
	  for (int j = 0; j < n; j++){
	    ls_veta += r[i][j]*w[i][j];
	  }
	}
	
	vh = IloNumVarArray(env, n, 0.0, 1.0, ILOINT);
	vz = IloNumVarArray(env, n * (n-1), 0.0, 1.0, ILOINT); /// Definição da variável vz_kl l!=k foi declarada como um vetorzão
	h = IloNumArray(env, n); /// Parâmetro que será usado para guardar o valor da variável vh e para guardar os coeficientes desta variável na função objetivo
	z = IloNumArray(env, n * (n-1) ); /// Definição do parâmetro que será usado para guardar o valor da variável vz_kl e para guardar os coeficientes desta variável na função objetivo
  veta = IloNumVarArray(env, n * n, 0.0, ls_veta, ILOFLOAT);
   

    /// ==============================
    /// Função objetivo - Mestre
    /// ==============================

    IloExpr expfo(env);
    for (int k = 0; k <  n; k++){
      expfo -= s[k] * vh[k];
      for (int l = 0; l <  n; l++){
        if (l!=k){
          expfo -= g[k][l] * vz[kl(n, k, l)];
        }
      }
    }
    for(int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        expfo += veta[ij(n, i, j)];
      }
    }
    IloAdd(mod, IloMaximize(env, expfo));
    expfo.end();
    
    ///=====================================================================================
    /// Restrição 1 do Problema Mestre -  forall(k in 1..n, l in 1..n, l!k)  z[k][l] <= h[k]
    ///=====================================================================================
    for (int k = 0; k < n; k++){
      for (int l = 0; l < n; l++){
        if (l!=k)
          mod.add(vz[kl(n, k, l)] <= vh[k]);
      }
    }  
    
    ///=====================================================================================
    /// Restrição 2 do Problema Mestre -  forall(k in 1..n, l in 1..n, l!k)  z[k][l] <= h[l]
    ///=====================================================================================
    for (int k = 0; k < n; k++){
      for (int l = 0; l < n; l++){
        if (l!=k)
          mod.add(vz[kl(n, k, l)] <= vh[l]);
      }
    }  
    
}

void BendersDecomposition::Create_SP(){
	modsp = IloModel(envsp); /// definição do modelo
	cplexsp = IloCplex(modsp); /// definição do objetivo cplex

	///===================================================
	/// Definição das variáveis e parâmetros principais do SPD
	///===================================================
	vbeta = IloNumVar(envsp, 0, IloInfinity, ILOFLOAT); /// Variável vbeta_ij
	vgama = IloNumVar(envsp, 0, IloInfinity, ILOFLOAT); /// Variável vgama_ij
	vtheta = IloNumVarArray(envsp, n, 0, IloInfinity, ILOFLOAT); /// Variável vtheta_ijk 
	vmu = IloNumVarArray(envsp, n, 0, IloInfinity, ILOFLOAT); /// Variável vmu_ijl
	vrho = IloNumVarArray(envsp, n * (n-1), 0, IloInfinity, ILOFLOAT); /// Variável vpho_ijkl l!=k
	vdelta = IloNumVarArray(envsp, n, -IloInfinity, IloInfinity, ILOFLOAT); /// Variável vdelta_ijk
	
	//beta = IloNum(envsp); ///  Não é necessário declarar parâmetro escalar
	//gama = IloNum(envsp); ///  Não é necessário declarar parâmetro escalar
	theta = IloNumArray(envsp, n); /// Definição do parâmetro que será usado para guardar o valor da variável theta_ijk
	mu = IloNumArray(envsp, n); /// Definição do parâmetro que será usado para guardar o valor da variável mu_ijl
	rho = IloNumArray(envsp, n * (n-1)); /// Definição do parâmetro que será usado para guardar o valor da variável rho_ijkl l!=k
	delta = IloNumArray(envsp, n); /// Definição do parâmetro que será usado para guardar o valor da variável delta_ijk
	
	
    /// ========================================
    /// Função objetivo simbólica do SPD
    /// Esta função objetivo é simbólica,
    /// pois a cada iteração ela irá se alterar
    /// =======================================

    IloExpr expfosp(envsp);
        expfosp += vbeta; 
        expfosp += vgama;
        for(int k = 0; k < n; k++){
          expfosp += vtheta[k]; /// Os coeficientes destas variáveis se alteram a cada iteração
          for(int l = 0; l < n; l++){
            if(l!=k){
              expfosp += vrho[kl(n, k, l)]; /// Os coeficientes destas variáveis se alteram a cada iteração
            }
          }
        }
        for(int l = 0; l < n; l++){
          expfosp += vmu[l]; /// Os coeficientes destas variáveis se alteram a cada iteração
        }
    fosp = IloAdd(modsp, IloMinimize(envsp,expfosp));
    expfosp.end();

    

  /// =========================
	/// Restrição 1 do subproblema
	/// ========================
	
	rd1 = IloRangeArray(envsp); /// Vetor que armazena as restrições do dual
	rd1UB = IloNumArray(envsp, n); /// Vetor que armazena o limitante superior de cada restrição (+infinito)
	rd1LB = IloNumArray(envsp, n);  /// Vetor que armazena o limitante inferior de cada restrição (w[i][j]*(r[i][j]-c[i][k]))
	
    for(int k = 0; k < n; k++){
	    IloExpr restsp1(envsp);
	    restsp1 += vbeta;
	    restsp1 += vdelta[k];
	    restsp1 += vtheta[k];
	    rd1.add(restsp1 <= 1); // w[i][j]*(r[i][j]-c[i][k])
	    restsp1.end();  
	    rd1UB[k] = +IloInfinity;
	   }
    modsp.add(rd1);	
    
    /// =========================
    /// Restrição 2 do subproblema
    /// ========================
    
    rd2 = IloRangeArray(envsp); /// Vetor que armazena as restrições do dual
    rd2UB = IloNumArray(envsp, n); /// Vetor que armazena o limitante superior de cada restrição (+infinito)
    rd2LB = IloNumArray(envsp, n);  /// Vetor que armazena o limitante inferior de cada restrição (-w[i][j]*c[l][j]))
    
   for(int l = 0; l < n; l++){
      IloExpr restsp2(envsp);
      restsp2 += vgama;
      restsp2 -= vdelta[l];
      restsp2 += vmu[l];
      rd2.add(restsp2 <= 1); // -w[i][j]*c[l][j]
      restsp2.end();
      rd2UB[l] = +IloInfinity;
    }
   modsp.add(rd2);
   
   /// =========================
   /// Restrição 3 do subproblema
   /// ========================
   
   rd3 = IloRangeArray(envsp); /// Vetor que armazena as restrições do dual
   rd3UB = IloNumArray(envsp, n * (n-1)); /// Vetor que armazena o limitante superior de cada restrição (+infinito)
   rd3LB = IloNumArray(envsp, n * (n-1));  /// Vetor que armazena o limitante inferior de cada restrição (-alpha*w[i][j]*c[k][l]))
   
  for(int k = 0; k < n; k++){
    for(int l = 0; l < n; l++){
      if(l!=k){
        IloExpr restsp3(envsp);
        restsp3 += vdelta[l];
        restsp3 -= vdelta[k];
        restsp3 += vrho[kl(n, k, l)];
        rd3.add(restsp3 <= 1); //-alpha*w[i][j]*c[k][l]
        restsp3.end();
        rd3UB[kl(n, k, l)] = +IloInfinity;
      }
    }
  }
  modsp.add(rd3);
  	
 
}

/// ==============================================
/// Método que irá montar um novo corte de Benders
/// a partir da resolução do SPD. A expressão para 
/// o corte será armazenada em OptCut
/// =============================================
void BendersDecomposition::GeraCortes(IloExprArray OptCut){
  
  ///=============================================================
  /// Evitar instabilidade numérica - resgate da solução do Mestre
  ///=============================================================
  for(int k = 0; k < n; k++){
    if(h[k] < 0.0001){
      h[k] = 0.0;
    }
    for (int l = 0; l <  n; l++){
      if (l!=k){
        if(z[kl(n,k,l)]< 0.0001){
          z[kl(n,k,l)]=0.0;
        } 
      }
    }
  }
  ///=============================================================  
  
  
  double sol_value = 0.0;
  for(int k = 0; k <  n; k++){
    sol_value -= s[k] * h[k];
    for(int l = 0; l <  n; l++){
      if (l!=k){
        sol_value -= g[k][l] * z[kl(n, k, l)];
      }
    }
  }  
  
	for(int i =0; i < n; i++){  //Para cada par ij será resolvido o subproblema dual e feito os cortes
	  for(int j =0; j < n; j++){
	    

		
		///=============================
		/// Atualiza FO do SPD ij
		///=============================	
	  
		    for(int k = 0; k < n; k++){
		      theta[k] = h[k];
		      for (int l = 0; l <  n; l++){
		        if (l!=k){
		          rho[kl(n, k, l)] = z[kl(n, k, l)];
		        }
		      }
		    }
		    for (int l = 0; l <  n; l++){
		      mu[l] = h[l];
		    }

			fosp.setLinearCoefs(vtheta, theta); /// Método que seta os coeficientes da variável vtheta_ijk na função objetivo
	    fosp.setLinearCoefs(vrho, rho); /// Método que seta os coeficientes da variável vrho_ijkl l!=k na função objetivo
	    fosp.setLinearCoefs(vmu, mu); /// Método que seta os coeficientes da variável vmu_ijl na função objetivo
		
		///=============================
		/// Atualiza as restrições do SPD ij
		///=============================	
		
		for(int k = 0; k < n; k++){
		  rd1LB[k] = w[i][j]*(r[i][j]-c[i][k]);
		}
		rd1.setBounds(rd1LB, rd1UB); /// Atualização dos limitantes das restrições duais 1
		
		for(int l = 0; l < n; l++){
		  rd2LB[l] = -w[i][j]*c[l][j];
		}
		rd2.setBounds(rd2LB, rd2UB); /// Atualização dos limitantes das restrições duais 2
		
		for(int k = 0; k < n; k++){
		  for(int l = 0; l < n; l++){
		    if(l!=k){
		      rd3LB[kl(n, k, l)] = -alpha*w[i][j]*c[k][l];
		    }
		  }
		}
		rd3.setBounds(rd3LB, rd3UB); /// Atualização dos limitantes das restrições duais 3
		
		
		///=============================
		/// Resolve o SPD ij
		///=============================
		cplexsp.solve(); /// Resolução do subproblema 
		
		if (cplexsp.getStatus() == IloAlgorithm::Optimal){ /// Verifica se uma solução ótima foi encontrada
    
			///================================
			/// Recupera-se a solução do SPD ij 
			///================================	
			  
			beta = cplexsp.getValue(vbeta);
		  gama = cplexsp.getValue(vgama);
		  cplexsp.getValues(vtheta, theta);
		  cplexsp.getValues(vmu, mu);
		  cplexsp.getValues(vrho, rho);
		  cplexsp.getValues(vdelta, delta);
		  
		  
		  sol_value += cplexsp.getObjValue(); /// Contabilização do valor da solução corrente, considerando também a contribuição das variáveis do mestre
		
		
			///=============================
			/// Monta-se o corte de Benders
			///=============================	
			IloNumExpr OptCutExpr(env);
			OptCutExpr += 0.0;
			OptCutExpr += veta[ij(n, i, j)];
		 	OptCutExpr -= beta;
		 	OptCutExpr -= gama;
		 	for(int k = 0; k < n; k++){
		 	  OptCutExpr -= vh[k] * theta[k];
		 	  for(int l = 0; l < n; l++){
		 	    if(l!=k){
		 	      OptCutExpr -= vz[kl(n, k, l)] * rho[kl(n, k, l)];
		 	    }
		 	  }
		 	}
		 	for(int l = 0; l < n; l++){
		 	  OptCutExpr -= vh[l] * mu[l];
		 	}
		 	OptCut.add(OptCutExpr);
			
		}
		else{ /// Como existe uma restrição no mestre que evita soluções inviáveis, o status do SPD tem que ser "Optimal", caso contrário existe algum erro na implementação ou na decomposição.
			cout<<"Erro: existe algum erro na implementação ou na decomposição. Corrija!!!"<<endl;
			exit(0);
		}
		
    }
	}
	
	lb = max(lb, sol_value); /// Atualiza-se o LB

	
}

void BendersDecomposition::Display_configuration(){
  cout<<"Concentradores: ";
  for(int k = 0; k < n; k++){
    if(h[k] >= 0.1){
      cout<<k+1<<"\t";
    }
  }
  cout<<endl;
}


