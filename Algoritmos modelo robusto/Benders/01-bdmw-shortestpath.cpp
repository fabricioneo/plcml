/*********************************************
 * Concert Model - Problema robusto completo
 * Author: Fabricio
 * Creation Date: 09/2022
 * B&B tree: one
 * Number of cut: root node 
 * Type of frac. cut: force
 *********************************************/


#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ilcplex/ilocplex.h>
#include <list> /// Novo include (Elisangela)

using namespace std;


#define TimeLimit 86400 //86400 //36000
#define EPSILON 1e-07 // 1e-07
#define M 10e10 /// Parâmetro usado para representar um número suficientemente grande

ILOSTLBEGIN
double best_inc = 1.e30;
double alpha1;
double receita;
double gamma1; // Define a porcentagem dos parâmetros de demanda que estão sob incerteza
int u = 0; /// Guardar o número de iterações

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
        double alpha; // Fator de desconto no custo de transporte em um link entre hubs
        double Gamma; // Parâmetro que define o tamanho do budget de incerteza - pertence ao intervalo [0, n*n]
        vector<double >  codx; // Coordenada x dos nós (AP)
        vector<double >  cody; // Coordenada y dos nós (AP)
        vector<vector<double > > w; // declara um vetor de vetores para representar uma matriz com n linhas e n colunas - Quantidade de demanda a ser enviada entre os nós i e j
        vector<vector<double > > w_chapeu; // Variação máxima de demanda sob incerteza a ser enviada entre os nós i e j
        vector<vector<double > > r; // Receita obtida por enviar uma unidade de demanda entre os nós i e j
        vector<vector<double > > c; // Custos por enviar uma unidade de demanda entre os nós i e j
        vector<double > s; // Custos fixos de instalação de um hub
        vector<double > s_chapeu; // Custos fixos de instalação de um hub com incerteza
        vector<vector<double > > g; // Custos de operação nos links inter-hubs
        vector<vector<double > > g_chapeu; // Custos de operação nos links inter-hubs com incerteza
        double soma; // Soma para obter a média geral dos custos fixos de instalação hub - Dados AP
		
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
        IloNumArray h2; /// variável auxiliar para variável de decisão do mestre
        IloNumArray z2; /// variável auxiliar para variável de decisão do mestre
        IloNumVar vq; // variável dual de decisão do mestre - inserida no modelo robusto
        IloNum q; // variável auxiliar dual de decisão do mestre - inserida no modelo robusto 
        IloNum q2; // variável auxiliar para variável de decisão do mestre
        IloNumVarArray vp; //  variáveis de decisão do mestre - inseridas no modelo robusto incerteza demanda
        IloNumVarArray vps; //  variáveis de decisão do mestre - inseridas no modelo robusto incerteza custos fixos instalacao
        IloNumVarArray vpg; //  variáveis de decisão do mestre - inseridas no modelo robusto incerteza custos fixos instalacao
        IloNumArray p; //  variáveis auxiliares para as variáveis inseridas no mestre do modelo robusto
        IloNumArray p2; /// variável auxiliar para variável de decisão do mestre
        IloNumArray ps; //  variáveis auxiliares para as variáveis inseridas no mestre do modelo robusto
        IloNumArray ps2; /// variável auxiliar para variável de decisão do mestre
        IloNumArray pg; //  variáveis auxiliares para as variáveis inseridas no mestre do modelo robusto
        IloNumArray pg2; /// variável auxiliar para variável de decisão do mestre

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
        void AtualizaLB(); /// Atualiza LB;
        void Display_configuration(); /// Imprime a solução corrente
        double Shortest_path_algorithm();

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
 
void BendersUserCutsCallback::main(){ /// Adiciona cortes a partir de soluções fracionárias
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
			bd->q=getValue(bd->vq);   //variável robusta (nova)
			getValues(bd->p, bd->vp); //variável robusta (nova)

			/// ========================
			/// Adiciona um corte de otimalidade
			/// ========================
			IloExprArray OptCut(getEnv()); /// Cria uma expressão chamada OptCut
			bd->GeraCortes(OptCut); /// Chama a função GeraCortes que resolve o SPD e guarda o corte gerado na expressão OptCut
			for(int i = 0; i < bd->n; i++){
			  for(int j = 0; j < bd->n; j++){
			    if(getValue(OptCut[ij(bd->n, i, j)]) >= 1e-07){
			      add(OptCut[ij(bd->n, i, j)] <= 0.000000); /// Adiciona o corte de Benders relativo a expressão OptCut
			      bd->NCuts++;
			    }
			   // add(OptCut[ij(bd->n, i, j)] <= 0.00); /// Adiciona o corte de Benders relativo a expressão OptCut
			   // bd->NCuts++; 
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
/// Adiciona corte de Benders a partir de soluções inteiras do mestre
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

void BendersLazyCutCallback::main(){ ///Adiciona cortes a partir de soluções inteiras (soluções incumbentes)
/*	/// ========================
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
	bd->q=getValue(bd->vq); //variável robusta (nova)
	getValues(bd->p, bd->vp); //variável robusta (nova)

	/// ========================
	/// Adiciona um corte de otimalidade
	/// ========================
	IloExprArray OptCut(getEnv()); /// Cria uma expressão chamada OptCut
	bd->GeraCortes(OptCut); /// Chama a função GeraCortes que resolve o SPD e guarda o corte gerado na expressão OptCut
	//bd->cplex.exportModel("mod1.lp");
	for(int i = 0; i < bd->n; i++){
	  for(int j = 0; j < bd->n; j++){
	    if(getValue(OptCut[ij(bd->n, i, j)]) >= EPSILON){
	      add(OptCut[ij(bd->n, i, j)] <= 0.000000);
	      bd->NCuts++; 
	    }
	    // add(OptCut[ij(bd->n, i, j)] <= 0); /// Adiciona o corte de Benders relativo a expressão OptCut
	    // bd->NCuts++; 
	  }
	}
	//bd->cplex.exportModel("mod2.lp");
	OptCut.endElements();
	cout<<"*"; /// Imprime na saída padrão um símbolo para mostrar que um corte foi gerado
 * */
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

       //bd->cplex.setParam(IloCplex::Threads,1);
       //bd->cplex.setParam(IloCplex::EpGap,0.00001);
       //bd->cplex.setOut(bd->env.getNullStream());
       //bd->cplexsp.setParam(IloCplex::Threads,1);
         bd->cplex.setParam(IloCplex::TiLim, TimeLimit);
         bd->cplex.setWarning(bd->env.getNullStream());
         bd->cplexsp.setWarning(bd->envsp.getNullStream());
         bd->cplexsp.setOut(bd->envsp.getNullStream());
         

        // 	bd->cplexsp.setParam(IloCplex::PreInd, false);
        // 	bd->cplex.setOut(cplexStream);
        //  bd->cplex.setParam(IloCplex::HeurFreq,-1);
        // 	bd->cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);
        // 	bd->cplexsp.setParam(IloCplex::RootAlg, IloCplex::Barrier);
       // bd->cplex.setParam(IloCplex::Param::Preprocessing::Reduce,1); //permite apenas reduções primais durante o pré-processamento - corrige problema usercut e lazy ativas juntas

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
        bd->cplex.setParam(IloCplex::Param::Preprocessing::Reduce,1); //permite apenas reduções primais durante o pré-processamento
 

		
        crono.start();
 
        /// =======================
        ///  Loop principal do Benders 
        /// =======================
        while(bd->gap > EPSILON){
          
          IloExprArray OptCut(bd->env); /// Cria uma expressão chamada OptCut
          bd->GeraCortes(OptCut); /// Chama a função GeraCortes que resolve o SPD e guarda o corte gerado na expressão OptCut
          for(int i = 0; i < bd->n; i++){
            for(int j = 0; j < bd->n; j++){
              bd->mod.add(OptCut[ij(bd->n, i, j)] <= 0.0000); /// Adiciona o corte de Benders relativo a expressão OptCut
            }
          }
          OptCut.endElements();
        
          

			/// =======================
			///  Resolve-se o problema mestre 
			/// =======================
			bd->cplex.solve();
			
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
			 	  bd->q = bd->cplex.getValue(bd->vq); /// Recupera a solução do mestre
			 	  bd->cplex.getValues(bd->vp,bd->p); /// Recupera a solução do mestre da variável robusta demanda
			 	  bd->cplex.getValues(bd->vps,bd->ps); /// Recupera a solução do mestre da variável robusta custos fixos instalacao
			 	  bd->cplex.getValues(bd->vpg,bd->pg); /// Recupera a solução do mestre da variável robusta custos fixos arcos
					bd->AtualizaLB();  
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
        ofstream arq_saida("Resultados-bdmw-robust-AP-correto.txt",std::ofstream::app);
		if (!arq_saida) {cerr << "Erro arquivo \n"; exit(0);}
		
        arq_saida<<argv[1]<<"\t"<<bd->alpha<<"\t"<<bd->r[0][0]<<"\t"<<u<<"\t"<<100 *  bd->gap<<"\t \t"<< fixed << setprecision(2)<<bd->lb<<"\t"<<bd->ub<<"\t"<<(double) crono.getTime()<<"\t"<<bd->Activate_Lazy_Callback<<"\t"<<bd->Activate_User_Callback<<"\t"<<gamma1<<endl;
		    cout<<argv[1]<<"\t"<<bd->alpha<<"\t"<<bd->r[0][0]<<"\t"<<u<<"\t"<<100 *  bd->gap<<"\t \t"<< fixed << setprecision(2)<<bd->lb<<"\t"<<bd->ub<<"\t"<<(double) crono.getTime()<<"\t"<<bd->Activate_Lazy_Callback<<"\t"<<bd->Activate_User_Callback<<"\t"<<gamma1<<endl;

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
    gamma1 = (argc > 6) ? atof(argv[6]) : 0;
}


void  BendersDecomposition::Read_data(char name[]){
    ifstream arq(name);
	if (!arq) {cerr << "Erro arquivo \n"; exit(0);}
	
	arq >> n;
	
	codx = vector<double>(n); // Coordenada x dos nós (AP)
	cody = vector<double>(n); // Coordenada y dos nós (AP)
	w = vector<vector<double > > (n, vector<double>(n)); // declara um vetor de vetores para representar uma matriz com n linhas e n colunas - Quantidade de demanda a ser enviada entre os nós i e j
	w_chapeu = vector<vector<double > > (n, vector<double>(n)); 
	r = vector<vector<double > > (n, vector<double>(n)); // Receita obtida por enviar uma unidade de demanda entre os nós i e j
	c = vector<vector<double > > (n, vector<double>(n)); // Custos por enviar uma unidade de demanda entre os nós i e j
	s = vector<double > (n); // Custos fixos de instalação de um hub
	s_chapeu = vector<double > (n);
	g = vector<vector<double > > (n, vector<double>(n)); // Custos de operação nos links inter-hubs
	g_chapeu = vector<vector<double > > (n, vector<double>(n));
	
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
	
	for (int i = 0; i <  n; i++){ //Determinar a variação máxima de incerteza na demanda de cada par (i,j) 
	  for (int j = 0; j < n; j++){
	    w_chapeu[i][j] = 0.3*w[i][j]; //(30% de w[i][j])
	  }
	}
	
	for (int i = 0; i <  n; i++){ //Coletar custos fixos instalação hub AP
	  arq >> s[i];
	  s[i] = s[i]/10;
	  soma = soma + s[i];
	}
	
	for (int i = 0; i <  n; i++){ //Variação máxima no custo fixo de instalação hub incerteza
	  s_chapeu[i] = 0.3*s[i];
	}
	
	for (int i = 0; i <  n; i++){ //Custos de operação link entre hubs AP
	  for (int j = 0; j < n; j++){
	    g[i][j] = 0.1 * soma/n;
	    //g[i][j] = 0.1 * ((s[i]+s[j])/2);  //alterei para a média dos custos entre os hubs - valor antigo: 0.1 * soma/n;
	  }
	}
	
	for (int i = 0; i <  n; i++){ //Determinar a variação máxima de incerteza nos custos fixos de instação de arcos hubs
	  for (int j = 0; j < n; j++){
	    g_chapeu[i][j] = 0.3*g[i][j]; //(30% de g[i][j])
	  }
	}
	
	Gamma = gamma1*((n*n)+n+(n*(n-1))); //Tamanho do budget de incerteza
	
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
	vq = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
	vp = IloNumVarArray(env, n * n, 0.0, IloInfinity, ILOFLOAT);
	vps = IloNumVarArray(env, n, 0.0, IloInfinity, ILOFLOAT);
	vpg = IloNumVarArray(env, n * (n-1), 0.0, IloInfinity, ILOFLOAT);
	h = IloNumArray(env, n); /// Parâmetro que será usado para guardar o valor da variável vh e para guardar os coeficientes desta variável na função objetivo
	z = IloNumArray(env, n * (n-1) ); /// Definição do parâmetro que será usado para guardar o valor da variável vz_kl e para guardar os coeficientes desta variável na função objetivo
	//q = IloNum(env);
	p = IloNumArray(env, n * n);
	ps = IloNumArray(env, n);
	pg = IloNumArray(env, n*(n-1));
	h2 = IloNumArray(env, n); /// Parâmetro que será usado para guardar o valor da variável vh e para guardar os coeficientes desta variável na função objetivo
	z2 = IloNumArray(env, n * (n-1) ); /// Definição do parâmetro que será usado para guardar o valor da variável vz_kl e para guardar os coeficientes desta variável na função objetivo
	//q2 = IloNum(env);
	p2 = IloNumArray(env, n * n);
	ps2 = IloNumArray(env, n);
	pg2 = IloNumArray(env, n*(n-1));
  veta = IloNumVarArray(env, n * n, 0.0, ls_veta, ILOFLOAT);
  
  ///====================================
  /// Cria ponto de Magnanti-Wong inicial
  ///====================================	
  
  for(int k = 0; k < n; k++){
      h2[k] = 1.00;
      h[k] = 0.00;
      for (int l = 0; l <  n; l++){
        if (l!=k){
          z2[kl(n,k,l)]=1.00;
          z[kl(n,k,l)]=0.00;
          pg2[kl(n,k,l)]=1.00;
          pg[kl(n,k,l)]=0.00;
        }
      }
    }
  
  q2 = 1.00;
  q = 0.00;
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      p2[ij(n,i,j)]=1.00;
      p[ij(n,i,j)]=0.00;
    }
  }
  
  for(int i = 0; i < n; i++){
    ps2[i]=1.00;
    ps[i]=0.00;
  }
   

    /// ==============================
    /// Função objetivo - Mestre
    /// ==============================

    IloExpr expfo(env);
    for (int k = 0; k <  n; k++){
      expfo += -s[k] * vh[k] - vps[k];
      for (int l = 0; l <  n; l++){
        if (l!=k){
          expfo += -g[k][l] * vz[kl(n, k, l)] - vpg[kl(n, k, l)];
        }
      }
    }
    for(int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
        expfo += veta[ij(n, i, j)]*w[i][j]-vp[ij(n, i, j)];
      }
    }
    expfo -= Gamma*vq;
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
    
    ///========================================================================================================
    /// Restrição 3 do Problema Mestre -  forall(i in 1..n, j in 1..n)  -q-p[i][j] <= -w_chapeu[i][j]veta[i][j]
    ///========================================================================================================
    for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
        mod.add(-vq-vp[ij(n, i, j)] <= -w_chapeu[i][j]*veta[ij(n,i,j)]);
      }
    }  
    
    ///========================================================================================================
    /// Restrição 4 do Problema Mestre -  forall(k in 1..n)  -q-ps[k] <= -s_chapeu[k]*vh[k]
    ///========================================================================================================
    for (int k = 0; k < n; k++){
      mod.add(-vq-vps[k] <= -s_chapeu[k]*vh[k]);
    }  
    
    ///==============================================================================================================
    /// Restrição 5 do Problema Mestre -  forall(k in 1..n, l in 1..n, l!k)  -q-pg[k][l] <= -g_chapeu[k][l]*vz[k][l]
    ///==============================================================================================================
    for (int k = 0; k < n; k++){
      for (int l = 0; l < n; l++){
        if (l!=k)
          mod.add(-vq-vpg[kl(n, k, l)] <= -g_chapeu[k][l]*vz[kl(n, k, l)]);
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
	rd1LB = IloNumArray(envsp, n);  /// Vetor que armazena o limitante inferior de cada restrição (r[i][j]-c[i][k])
	
    for(int k = 0; k < n; k++){
	    IloExpr restsp1(envsp);
	    restsp1 += vbeta;
	    restsp1 += vdelta[k];
	    restsp1 += vtheta[k];
	    rd1.add(restsp1 <= 1); // (r[i][j]-c[i][k])
	    restsp1.end();  
	    rd1UB[k] = +IloInfinity;
	   }
    modsp.add(rd1);	
    
    /// =========================
    /// Restrição 2 do subproblema
    /// ========================
    
    rd2 = IloRangeArray(envsp); /// Vetor que armazena as restrições do dual
    rd2UB = IloNumArray(envsp, n); /// Vetor que armazena o limitante superior de cada restrição (+infinito)
    rd2LB = IloNumArray(envsp, n);  /// Vetor que armazena o limitante inferior de cada restrição (-c[l][j]))
    
   for(int l = 0; l < n; l++){
      IloExpr restsp2(envsp);
      restsp2 += vgama;
      restsp2 -= vdelta[l];
      restsp2 += vmu[l];
      rd2.add(restsp2 <= 1); // -c[l][j]
      restsp2.end();
      rd2UB[l] = +IloInfinity;
    }
   modsp.add(rd2);
   
   /// =========================
   /// Restrição 3 do subproblema
   /// ========================
   
   rd3 = IloRangeArray(envsp); /// Vetor que armazena as restrições do dual
   rd3UB = IloNumArray(envsp, n * (n-1)); /// Vetor que armazena o limitante superior de cada restrição (+infinito)
   rd3LB = IloNumArray(envsp, n * (n-1));  /// Vetor que armazena o limitante inferior de cada restrição (-alpha*c[k][l]))
   
  for(int k = 0; k < n; k++){
    for(int l = 0; l < n; l++){
      if(l!=k){
        IloExpr restsp3(envsp);
        restsp3 += vdelta[l];
        restsp3 -= vdelta[k];
        restsp3 += vrho[kl(n, k, l)];
        rd3.add(restsp3 <= 1); //-alpha*c[k][l]
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

  ///=============================
  /// Cria ponto de Magnanti-Wong
  ///=============================	
  
  for(int k = 0; k < n; k++){
    h2[k] = 0.5*(h2[k]+h[k]);
    ps2[k] = 0.5*(ps2[k]+ps[k]);
    for (int l = 0; l <  n; l++){
      if (l!=k){
        z2[kl(n,k,l)]=0.5*(z2[kl(n,k,l)]+z[kl(n,k,l)]);
        pg2[kl(n,k,l)]=0.5*(pg2[kl(n,k,l)]+pg[kl(n,k,l)]);
      }
    }
  }
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      p2[ij(n,i,j)]=0.5*(p2[ij(n,i,j)]+p[ij(n,i,j)]);
    }
  }
    
  q2 = 0.5*(q2+q);

  
	for(int i =0; i < n; i++){  //Para cada par ij será resolvido o subproblema dual e feito os cortes
	  for(int j =0; j < n; j++){
	    

		
		///=============================
		/// Atualiza FO do SPD ij
		///=============================	
	  
		    for(int k = 0; k < n; k++){
		      theta[k] = h2[k];
		      for (int l = 0; l <  n; l++){
		        if (l!=k){
		          rho[kl(n, k, l)] = z2[kl(n, k, l)];
		        }
		      }
		    }
		    for (int l = 0; l <  n; l++){
		      mu[l] = h2[l];
		    }

			fosp.setLinearCoefs(vtheta, theta); /// Método que seta os coeficientes da variável vtheta_ijk na função objetivo
	    fosp.setLinearCoefs(vrho, rho); /// Método que seta os coeficientes da variável vrho_ijkl l!=k na função objetivo
	    fosp.setLinearCoefs(vmu, mu); /// Método que seta os coeficientes da variável vmu_ijl na função objetivo
		
		///=============================
		/// Atualiza as restrições do SPD ij
		///=============================	
		
		for(int k = 0; k < n; k++){
		  rd1LB[k] = (r[i][j]-c[i][k]);
		}
		rd1.setBounds(rd1LB, rd1UB); /// Atualização dos limitantes das restrições duais 1
		
		for(int l = 0; l < n; l++){
		  rd2LB[l] = -c[l][j];
		}
		rd2.setBounds(rd2LB, rd2UB); /// Atualização dos limitantes das restrições duais 2
		
		for(int k = 0; k < n; k++){
		  for(int l = 0; l < n; l++){
		    if(l!=k){
		      rd3LB[kl(n, k, l)] = -alpha*c[k][l];
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

	
}

void BendersDecomposition::AtualizaLB(){
  
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
  
  
  // double sol_value = 0.0;
  // for(int k = 0; k <  n; k++){
  //   sol_value -= s[k] * h[k];
  //   for(int l = 0; l <  n; l++){
  //     if (l!=k){
  //       sol_value -= g[k][l] * z[kl(n, k, l)];
  //     }
  //   }
  // }  
  // 
  // IloTimer cronosp(envsp); /// (Elisangela)
  // cronosp.start(); /// (Elisangela)
  // 
  // for(int i =0; i < n; i++){  //Para cada par ij será resolvido o subproblema dual e feito os cortes
  //   for(int j =0; j < n; j++){
  //     
  //     
  //     
  //     ///=============================
  //     /// Atualiza FO do SPD ij
  //     ///=============================	
  //     
  //     for(int k = 0; k < n; k++){
  //       theta[k] = h[k];
  //       for (int l = 0; l <  n; l++){
  //         if (l!=k){
  //           rho[kl(n, k, l)] = z[kl(n, k, l)];
  //         }
  //       }
  //     }
  //     for (int l = 0; l <  n; l++){
  //       mu[l] = h[l];
  //     }
  //     
  //     fosp.setLinearCoefs(vtheta, theta); /// Método que seta os coeficientes da variável vtheta_ijk na função objetivo
  //     fosp.setLinearCoefs(vrho, rho); /// Método que seta os coeficientes da variável vrho_ijkl l!=k na função objetivo
  //     fosp.setLinearCoefs(vmu, mu); /// Método que seta os coeficientes da variável vmu_ijl na função objetivo
  //     
  //     ///=============================
  //     /// Atualiza as restrições do SPD ij
  //     ///=============================	
  //     
  //     for(int k = 0; k < n; k++){
  //       rd1LB[k] = w[i][j]*(r[i][j]-c[i][k]);
  //     }
  //     rd1.setBounds(rd1LB, rd1UB); /// Atualização dos limitantes das restrições duais 1
  //     
  //     for(int l = 0; l < n; l++){
  //       rd2LB[l] = -w[i][j]*c[l][j];
  //     }
  //     rd2.setBounds(rd2LB, rd2UB); /// Atualização dos limitantes das restrições duais 2
  //     
  //     for(int k = 0; k < n; k++){
  //       for(int l = 0; l < n; l++){
  //         if(l!=k){
  //           rd3LB[kl(n, k, l)] = -alpha*w[i][j]*c[k][l];
  //         }
  //       }
  //     }
  //     rd3.setBounds(rd3LB, rd3UB); /// Atualização dos limitantes das restrições duais 3
  //     
  //     
  //     ///=============================
  //     /// Resolve o SPD ij
  //     ///=============================
  //     cplexsp.solve(); /// Resolução do subproblema 
  //     
  //     if (cplexsp.getStatus() == IloAlgorithm::Optimal){ /// Verifica se uma solução ótima foi encontrada
  //       
  //       sol_value += cplexsp.getObjValue(); /// Contabilização do valor da solução corrente, considerando também a contribuição das variáveis do mestre
  //       
  //       
  //     }
  //     else{ /// Como existe uma restrição no mestre que evita soluções inviáveis, o status do SPD tem que ser "Optimal", caso contrário existe algum erro na implementação ou na decomposição.
  //       cout<<"Erro: existe algum erro na implementação ou na decomposição. Corrija!!!"<<endl;
  //       exit(0);
  //     }
  //     
  //   }
  // }
  // 
  // double time = cronosp.getTime(); /// Tempo que foi gasto resolvendo o SP pelo solver (Elisangela)
  // 
  // double SPA = Shortest_path_algorithm();
  // cout<<"Comparando avaliação da soluções "<<sol_value<<"\t"<<SPA<<endl;
  // if(sol_value <= SPA - 0.1 or  sol_value >= SPA + 0.1){
  //   Display_configuration();
  //   exit(0);
  // }
  // cout<<"tempo "<<time<<"\t"<<cronosp.getTime()-time<<endl;
  // 
  // lb = max(lb, sol_value); /// Atualiza-se o LB
  
  double SPA = Shortest_path_algorithm();
  lb = max(lb, SPA); /// Atualiza-se o LB
  
}

double BendersDecomposition::Shortest_path_algorithm(){ /// (Elisangela)
  
  vector< vector<double> > hub_net_sp_matrix = vector< vector<double> > (n, vector<double> (n)); /// Shostest path considering only the hub network 
  vector< vector<double> > sp_matrix = vector< vector<double> > (n, vector<double> (n)); /// Shostest path for all nodes 
  
  ///===================================
  /// Calculate the shortest path from each pair of hubs considering only the hub network 
  ///===================================
  double sol_value = 0.0; /// Variável que computa o custo de instalação de hubs e de arcos entre hubs
  list<int> hub_list; /// Lista de hubs
  
  for(int k = 0; k < n; k++){ /// Monta uma lista de hubs
    if(h[k] >= (1 - EPSILON)){
      hub_list.insert(hub_list.end(), k);
      sol_value -= s[k];
    }
  }
  
  for(list<int>::iterator k = hub_list.begin(); k != hub_list.end(); k++){ /// Calculo do caminho mais curto na rede de hubs com apenas dois hubs
    for(list<int>::iterator m = hub_list.begin(); m != hub_list.end(); m++){
      if(k != m){
        if( z[kl(n, *k, *m)] >= (1 - EPSILON)){
          hub_net_sp_matrix[*k][*m] = alpha * c[*k][*m];
          sol_value -= g[*k][*m];
        }
        else {
          hub_net_sp_matrix[*k][*m] = M;
        }
      }
    }
  } 
  for(list<int>::iterator k = hub_list.begin(); k != hub_list.end(); k++){
    for(list<int>::iterator m = hub_list.begin(); m != hub_list.end(); m++){
      for(list<int>::iterator l = hub_list.begin(); l != hub_list.end(); l++){
        if(m != l and m != k and l != k){
          double cost_aux = hub_net_sp_matrix[*m][*k] + hub_net_sp_matrix[*k][*l];
          if(( cost_aux < hub_net_sp_matrix[*m][*l])){
            hub_net_sp_matrix[*m][*l] = cost_aux;
          }
        }
      }
    }
  }
  
  
  /// Compute shortest path for all pair i,j, also allowing the direct connections 
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      sp_matrix[i][j] = M; /// (Esse M não permite conexões diretas nesse caso. ???)
      for(list<int>::iterator k = hub_list.begin(); k != hub_list.end(); k++){
        sp_matrix[i][j] = min(c[i][*k]+c[*k][j], sp_matrix[i][j] );  
        for(list<int>::iterator m = hub_list.begin(); m != hub_list.end(); m++){
          if(k != m){
            double cost_aux = c[i][*k] + hub_net_sp_matrix[*k][*m] + c[*m][j];
            if(sp_matrix[i][j] > cost_aux){
              sp_matrix[i][j] = cost_aux;
            }
          }
        }
        
      }
      if( r[i][j] - sp_matrix[i][j] > EPSILON){
        sol_value += (r[i][j] - sp_matrix[i][j])* w[i][j]-p[ij(n, i, j)];
      }
    }
  } 
  
  sol_value -= Gamma*q;
  
  for(int k = 0; k < n; k++){
    sol_value -= ps[k];
    for(int l = 0; l < n; l++){
      if(l!=k){
        sol_value -=pg[kl(n, k, l)];
      }
    }
  }
  
  
  hub_net_sp_matrix.clear();
  hub_list.clear();
  return sol_value;
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


