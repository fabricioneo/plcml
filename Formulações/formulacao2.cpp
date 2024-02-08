/*********************************************
 * Concert Model
 * Autor: Fabricio
 * Data de criação: 15-03-2020 
 * Problem - Localização de hub com max do lucro - FORMULAÇÃO 2  
 *********************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ilcplex/ilocplex.h>

using namespace std;
ILOSTLBEGIN
  
  ILOMIPINFOCALLBACK3(Callback_mipinfo, double&, lb, double&, ub, double&, gap) {/// Callback nomeada Callback_mipinfo com 3  argumentos do tipo double
    ub = getBestObjValue(); // Melhor lower bound da árvore de B&C
    if(lb < getIncumbentObjValue()){
      lb = getIncumbentObjValue(); // Valor da melhor solução incumbent (solução incumbent corrente)  encontrada durante o B&C
    } 
    gap = getMIPRelativeGap();
  }

int main (int argc, char *argv[]){
    try{

		/**============================
		 *  Leitura dos dados
		 *=============================== */
		ifstream arq(argv[1]);
		if (!arq.is_open()){
			cout << "Error openning file: " << argv[1] << endl;
			arq.close();
			exit(EXIT_FAILURE);
		}
		
		int n; // Quantidade de nós
		arq >> n;
		float alpha; // Fator de desconto no custo de transporte em um link entre hubs
		vector<double >  codx(n); // Coordenada x dos nós (AP)
		vector<double >  cody(n); // Coordenada y dos nós (AP)
  	vector<vector<double > > w(n, vector<double>(n)); // declara um vetor de vetores para representar uma matriz com n linhas e n colunas - Quantidade de demanda a ser enviada entre os nós i e j
  	vector<vector<double > > r(n, vector<double>(n)); // Receita obtida por enviar uma unidade de demanda entre os nós i e j
  	vector<vector<double > > c(n, vector<double>(n)); // Custos por enviar uma unidade de demanda entre os nós i e j
  	vector<double > s(n); // Custos fixos de instalação de um hub
  	vector<vector<double > > g(n, vector<double>(n)); // Custos de operação nos links inter-hubs
  	float soma = 0; // Soma para obter a média geral dos custos fixos de instalação hub - Dados AP
  	
  	// Sintaxe para o executável: make comp-form2
		
		// ========================================================================================
		// COLETAR DADOS AP  (Exemplo de sintaxe para rodar: ./form2 10LTm.txt)
		// ========================================================================================

		// if (argc >= 3) //Coletar fator de desconto nos links inter-hubs AP
		//   alpha = atof(argv[2]);
		// else alpha = 0.2;
		// 
		// for (int i = 0; i <  n; i++){ //Receita AP 20,30,50
		//   for (int j = 0; j < n; j++){
		//     if (argc >= 4) 
		//       r[i][j] = atof(argv[3]);
		//     else  r[i][j] = 20;
		//   }
		// }
		// 
		// for (int i = 0; i <  n; i++){ //Coletar as coordenadas dos nós
		//       arq >> codx[i];
		//       arq >> cody[i];
		// }
		// 
		// for (int i = 0; i <  n; i++){ //Coletar custos AP
		//   for (int j = 0; j < n; j++){
		//     c[i][j] = 0.001 * sqrt((codx[i] - codx[j]) * (codx[i] - codx[j])  + (cody[i] - cody[j]) * (cody[i] - cody[j]));
		//   }
		// }
		// 
		// for (int i = 0; i <  n; i++){ //Coletar demanda AP
		//   for (int j = 0; j < n; j++){
		//     arq >> w[i][j];
		//   }
		// }
		// 
		// for (int i = 0; i <  n; i++){ //Coletar custos fixos instalação hub AP
		//   arq >> s[i];
		//   s[i] = s[i]/10;
		//   soma = soma + s[i];
		// }
		// 
		// for (int i = 0; i <  n; i++){ //Custos de operação link entre hubs AP
		//   for (int j = 0; j < n; j++){
		//     g[i][j] = 0.1 * soma/n;
		//   }
		// }




		// ========================================================================================
		// COLETAR DADOS CAB (Exemplo de sintaxe para rodar: ./form2 cab.txt 0.8 1000 150 15)
		// ========================================================================================
		
		
		for (int i = 0; i <  n; i++){ //Coletar demanda e custos CAB
		  for (int j = 0; j < n; j++){
		    arq >> w[i][j];
		    arq >> c[i][j];
		    soma = soma + w[i][j];
		  }
		}

		for (int i = 0; i <  n; i++){ //Escalonando demanda CAB
		  for (int j = 0; j < n; j++){
		    w[i][j] = w[i][j]/soma;
		  }
		}

		if (argc >= 3) //Coletar fator de desconto nos links inter-hubs CAB
		  alpha = atof(argv[2]);
		else alpha = 0.2;

		for (int i = 0; i <  n; i++){ //Coletar receita CAB
		  for (int j = 0; j < n; j++){
		    //r[i][j] = 100*(i+j+1);
		   r[i][j] = atoi(argv[3]);
		  }
		}

		for (int i = 0; i <  n; i++){ //Coletar custo fixo de instalação CAB
		  s[i] = atoi(argv[4]);
		  }


		for (int i = 0; i <  n; i++){ //Coletar custo de operar links inter-hubs CAB
		  for (int j = 0; j < n; j++){
		    g[i][j] = atoi(argv[5]);
		  }
		}



		/**============================
		*  Declaração do modelo
		*=============================== */
		
		IloEnv env;
		IloModel mod(env);
		IloCplex cplex(mod);
		
		IloNumVarArray h(env, n, 0, 1, ILOBOOL); // ILOBOOL h[k] indica se um hub é localizado no nó k 
		
		IloArray<IloNumVarArray> z(env, n); //  ILOBOOL z[k][l] indica se um link inter-hub é operado entre os hubs l e k
		for(int k = 0; k < n; k++)   
			z[k] = IloNumVarArray(env, n, 0, 1, ILOBOOL); 
		
		IloArray<IloArray<IloArray<IloNumVarArray>>> x(env, n); // x[i][j][k][l] representa a fração da demanda entre os nós i e j que é roteada entre os hubs k e l
		for(int i = 0; i < n; i++){
		  x[i] = IloArray<IloArray<IloNumVarArray>>(env, n);
		  for(int j = 0; j < n; j++){
		    x[i][j] = IloArray<IloNumVarArray>(env, n);
		    for(int k=0; k<n; k++){
		      x[i][j][k] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT); 
		    }
		  }
		}
		
		IloArray<IloArray<IloNumVarArray>> a(env, n); // a[i][j][k] é a fração da demanda entre os nós i e j que é roteada através de um caminho no qual o primeiro hub é k
		for(int i = 0; i < n; i++){
		  a[i] = IloArray<IloNumVarArray>(env, n);
		  for(int j = 0; j < n; j++){
		    a[i][j] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
		  }
		}
		
		IloArray<IloArray<IloNumVarArray>> b(env, n); // b[i][j][l] é a fração da demanda entre os nós i e j que é roteada através de um caminho no qual o último hub é l
		for(int i = 0; i < n; i++){
		  b[i] = IloArray<IloNumVarArray>(env, n);
		  for(int j = 0; j < n; j++){
		    b[i][j] = IloNumVarArray(env, n, 0, IloInfinity, ILOFLOAT);
		  }
		}
		
		
		// ====================================Formulação 2=============================================================== 
		// maximize sum(i in 1..n, j in 1..n, k in 1..n, l in 1..n) r[i][j] * w[i][j] * a[i][j][k] 
		// - [ sum(i in 1..n, j in 1..n, k in 1..n) c[i][k] * w[i][j] * a[i][j][k]
		// + sum(i in 1..n, j in 1..n, l in 1..n) c[l][j] * w[i][j] * b[i][j][l]
		// + sum(i in 1..n, j in 1..n, k in 1..n, l in 1..n) alpha * c[k][l] * x[i][j][k][l] + sum(k in 1..n) s[k] * h[k] 
		// + sum(k in 1..n, l in 1..n) g[k][l] * z[k][l] ];
		// ===============================================================================================================

		IloExpr expfo(env);
		for (int i = 0; i < n; i++){
		  for (int j = 0; j < n; j++){
		    for (int k = 0; k < n; k++){
		      expfo += (r[i][j] - c[i][k]) * w[i][j] * a[i][j][k];
		    }
		    for (int l = 0; l < n; l++){
		      expfo -= c[l][j] * w[i][j] * b[i][j][l];
		    }
		    for (int k = 0; k < n; k++){
		      for (int l = 0; l < n; l++){
		        expfo -= alpha * c[k][l] * w[i][j] * x[i][j][k][l];
		      }
		    }
		  }
		}
		for (int k = 0; k <  n; k++){
			expfo -= s[k] * h[k];
			for(int l = 0; l < n; l++){
				expfo -=  g[k][l] * z[k][l];
			}
		}
		  
		IloAdd(mod, IloMaximize(env, expfo));
		expfo.end();

		//===========================================================================
		// forall(i in 1..n, j in 1..n)  sum(k in 1..n) a[i][j][k] <= 1 
		//===========================================================================
		for (int i = 0; i < n; i++){
		  for (int j =0; j < n; j++){
		  IloExpr r1(env);
		    for (int k = 0; k < n; k++){
		      r1 += a[i][j][k];
		    }
		      mod.add(r1 <= 1);
		      r1.end();
		  }
		}
		
		//===========================================================================
		// forall(i in 1..n, j in 1..n)  sum(l in 1..n) b[i][j][l] <= 1 
		//===========================================================================
		for (int i = 0; i < n; i++){
		  for (int j =0; j < n; j++){
		    IloExpr r2(env);
		    for (int l = 0; l < n; l++){
		      r2 += b[i][j][l];
		    }
		    mod.add(r2 <= 1);
		    r2.end();
		  }
		}
		
		//==========================================================================================================================================
		// forall(i in 1..n, j in 1..n, k in 1..n)  a[i][j][k] + sum(l in 1..n, l!k) a[i][j][k] + x[i][j][l][k] == b[i][j][k] + sum(l in 1..n, l!k) x[i][j][k][l]
		//==========================================================================================================================================
		for (int i = 0; i < n; i++){
		  for (int j =0; j < n; j++){
		    for (int k =0; k < n; k++){
		    IloExpr r3(env);
		    r3 += a[i][j][k] - b[i][j][k];
		      for (int l = 0; l < n; l++){
		        if (l!=k){
		          r3 += x[i][j][l][k] - x[i][j][k][l];
		        }
		      }
		    mod.add(r3 == 0);
		    r3.end();
		    }
		  }
		}
		
		//=============================================================
		// forall(i in 1..n, j in 1..n, k in 1..n) a[i][j][k] <= h[k], i != k
		//=============================================================
		for (int i = 0; i < n; i++){
		  for (int j = 0; j < n; j++){
		    for (int k = 0; k < n; k++){
		      mod.add(a[i][j][k] <= h[k]);
		    }
		  }
		}
		
		//=============================================================
		// forall(i in 1..n, j in 1..n, l in 1..n) b[i][j][l] <= h[l], j != l
		//=============================================================
		for (int i = 0; i < n; i++){
		  for (int j = 0; j < n; j++){
		    for (int l = 0; l < n; l++){
		      mod.add(b[i][j][l] <= h[l]);
		    }
		  }
		}
		
		//==================================================================================================
		// forall(i in 1..n, j in 1..n, k in 1..n, l in 1..n, l!k) x[i][j][k][l] <= z[k][l]
		//==================================================================================================
		for (int i = 0; i < n; i++){
		  for (int j = 0; j < n; j++){
		    for (int k = 0; k < n; k++){
		      for (int l = 0; l < n; l++){
		        if (l!=k)
		        mod.add(x[i][j][k][l] <= z[k][l]);
		      }
		    }
		  }
		}
		
		//===================================================
		// forall(k in 1..n, l in 1..n, l!k)  z[k][l] <= h[k]
		//===================================================
		for (int k = 0; k < n; k++){
		  for (int l = 0; l < n; l++){
		    if (l!=k)
		    mod.add(z[k][l] <= h[k]);
		  }
		}
		
		//===================================================
		// forall(k in 1..n, l in 1..n, l!k)  z[k][l] <= h[l]
		//===================================================
		for (int k = 0; k < n; k++){
		  for (int l = 0; l < n; l++){
		    if (l != k)
		    mod.add(z[k][l] <= h[l]);
		  }
		}
 

        /// ==========================
        /// configurações do cplex
        /// ==========================

        //cplex.setParam(IloCplex::EpGap, 0.0000001); // Definindo uma tolerancia
		    cplex.setParam(IloCplex::TiLim, 86400); // Tempo limite de resolução
		    cplex.setWarning(env.getNullStream()); // Eliminar warnings
		//  cplex.setOut(mono->env.getNullStream()); // Eliminar os logs do solver
    //  cplex.setParam(IloCplex::Threads, 4); // Definir a quantidade de threads
     cplex.setParam(IloCplex::Param::Benders::Strategy, 3); // Ativar Benders do solver


        ///==============================
        /// Resolvendo o problema
        ///==============================

        IloTimer crono(env);// Variável para coletar o tempo
        double lb = 0; /// Callback_mipinfo
        double ub = 10e10; /// Callback_mipinfo
        double gap = 1; /// Callback_mipinfo
        cplex.use(Callback_mipinfo(env, lb, ub, gap));/// Callback_mipinfo
		    
		    for(int i = 0; i< n; i++){   //Pré-fixando variáveis
		      for(int j = 0; j < n; j++){
		        for(int k = 0; k < n; k++){
		          if (r[i][j] < c[i][k]){
		              a[i][j][k].setBounds(0, 0);
		          }
		        }
		        for(int l = 0; l < n; l++){
		            if (r[i][j] < c[l][j]){
		              b[i][j][l].setBounds(0, 0);
		            }
		        }
		      }
		    }
		    
          
        crono.start();
        cplex.solve();
        crono.stop();
        
       
        if(cplex.getStatus() == IloAlgorithm::Optimal){
            lb = cplex.getObjValue();
            ub = cplex.getObjValue();
            gap = 0.0;
        }


        ///=====================================
        /// Salvando os resultados - CAB
        ///=====================================

	      printf("%s\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.6f\n",  argv[1], n, alpha,  100 *  gap,  lb, ub, (double) crono.getTime());

        cout << " Valor função objetivo: " << cplex.getObjValue () << endl ;

        printf("Hubs Instalados:");
        for(int j = 0; j < n; j++){
          if( cplex.getValue(h[j]) >= 0.1){
            printf(" %d\t \n", j+1);
          }
        }
        cout<<endl;

        /**=====================================
         *  Apresenta a configuração final
         * ====================================*/




		FILE *re;
		re = fopen("ResultadosF2.txt","aw+");
		// fprintf(re, "\n Informações Gerais: " "%s\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f \n",  argv[1], n, alpha,  r[0][0],  s[0], g[0][0] );
		// fprintf(re, "\n Valor função objetivo: " "%f\t \n", (double) cplex.getObjValue ());
		// fprintf(re, "\n Tempo de CPU: " "%f\t \n", (double) crono.getTime());
		// fprintf(re, "\n lb: " "%1.2f\t \t" "ub: " "%1.2f\t \t" "gap: " "%1.2f \n", lb, ub, 100 *  gap);

		float cont = 0;
		float demanda = 0;
		for(int i = 0; i < n; i++){
		  for(int j = 0; j < n; j++){
		    for(int k = 0; k < n; k++){
		      if(cplex.getValue(a[i][j][k])>0.001){
		        cont = cont + 1;
		        demanda = demanda + w[i][j];
		      }
		    }
		  }
		}


		// fprintf(re, "\n Quantidade de pares atendidos: " "%f \t \t"  "Porcentagem de pares atendidos: " "%f \t \n", cont, cont/600);
		// fprintf(re, "\n Porcentagem demanda satisfeita: " "%1.6f \t", demanda);
		
		fprintf(re, "%s \t %d \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \n",  argv[1], n, alpha,  r[0][0],  s[0], g[0][0], (double) cplex.getObjValue (), (double) crono.getTime(), cont/600, demanda);
		

		// fprintf(re, "\n \n Hubs: " );
		// for(int j = 0; j < n; j++){
		//   if( cplex.getValue(h[j]) >= 0.1){
		//     fprintf(re, "h[%d] = %f \t ", j+1, cplex.getValue(h[j]));
		//   }
		// }
		// 
		// fprintf(re, "\n \n Variável x: \n \n");
		// for(int i = 0; i < n; i++){
		//   for(int j = 0; j < n; j++){
		//     for(int k = 0; k < n; k++){
		//       for(int l =0; l< n; l++){
		//         if(cplex.getValue(z[k][l])>0.001)
		//           fprintf(re, "x[%d][%d][%d][%d]= %f\t \t", i+1, j+1, k+1, l+1, cplex.getValue(z[k][l]));
		//       }
		//     }
		//   }
		// }
		// 
		// 
		// fprintf(re, "\n \n Variável a: \n \n");
		// for(int i = 0; i < n; i++){
		//   for(int j = 0; j < n; j++){
		//     for(int k = 0; k < n; k++){
		//       if(cplex.getValue(a[i][j][k])>0.001)
		//         fprintf(re, "a[%d][%d][%d]=%f \t \t", i+1, j+1, k+1, cplex.getValue(a[i][j][k]) );
		//     }
		//   }
		// }
		// 
		// fprintf(re, "\n \n Variável b: \n \n");
		// for(int i = 0; i < n; i++){
		//   for(int j = 0; j < n; j++){
		//     for(int k = 0; k < n; k++){
		//       if(cplex.getValue(b[i][j][k])>0.001)
		//         fprintf(re, "b[%d][%d][%d]=%f \t \t", i+1, j+1, k+1, cplex.getValue(b[i][j][k]));
		//     }
		//   }
		// }
		// 
		// fprintf(re, "\n \n Variável z: \n \n");
		// for(int i = 0; i < n; i++){
		//   for(int j = 0; j < n; j++){
		//         if(cplex.getValue(z[i][j])>0.001)
		//         fprintf(re, "z[%d][%d]=%f \t \t", i+1, j+1, cplex.getValue(z[i][j]));
		//     }
		//   }

	//	fprintf(re, "\n \n ====================================================================== \n \n");
		
		// FILE *re;
		// re = fopen("ResultadosF2_LP.txt","aw+"); 
		// fprintf(re, "%1.2f\t%1.2f\t%1.2f\t%1.2f\t%f\t%f \n", alpha,  r[0][0],  s[0], g[0][0], (double) cplex.getObjValue (), (double) crono.getTime());
		
		
    
    ///=====================================
    /// Salvando os resultados - AP
    ///=====================================
    
    // printf("%s\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.6f\n",  argv[1], n, alpha,  100 *  gap,  lb, ub, (double) crono.getTime());
    // 
    // cout << " Valor função objetivo: " << cplex.getObjValue () << endl ;
    // 
    // printf("Hubs Instalados:");
    // for(int j = 0; j < n; j++){
    //   if( cplex.getValue(h[j]) >= 0.1){
    //     printf(" %d\t \n", j+1);
    //   }
    // }
    // cout<<endl;
    // 
    // /**=====================================
    //  *  Apresenta a configuração final
    //  * ====================================*/
    // 
    // FILE *re;
    // re = fopen("ResultadosF2_AP.txt","aw+");
    // fprintf(re, "\n Informações Gerais: " "%s\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f \n",  argv[1], n, alpha,  r[0][0], (double) cplex.getObjValue (), (double) crono.getTime() );
    // fprintf(re, "\n Valor função objetivo: " "%f\t \n", (double) cplex.getObjValue ());
    // fprintf(re, "\n Tempo de CPU: " "%f\t \n", (double) crono.getTime());
    // 
    // 
    // float cont = 0;
    // for(int i = 0; i < n; i++){
    //   for(int j = 0; j < n; j++){
    //     for(int k = 0; k < n; k++){
    //       if(cplex.getValue(a[i][j][k])>0.001){
    //         cont = cont + 1;
    //       }
    //     }
    //   }
    // }
    // 
    // 
    // fprintf(re, "\n Quantidade de pares atendidos: " "%f \t"  "Porcentagem de pares atendidos: " "%f \t \n", cont, cont/1600);
    // 
    // fprintf(re, "\n Hubs: \n" );
    // for(int j = 0; j < n; j++){
    //   if( cplex.getValue(h[j]) >= 0.1){
    //     fprintf(re, "%d\t \n", j+1);
    //   }
    // }

    // fprintf(re, "\n Variável x: \n");
    // for(int i = 0; i < n; i++){
    //   for(int j = 0; j < n; j++){
    //     for(int k = 0; k < n; k++){
    //       for(int l =0; l< n; l++){
    //         if(cplex.getValue(z[k][l])>0.001)
    //           fprintf(re, "x[%d][%d][%d][%d] = %f \t \n", i, j, k, l, cplex.getValue(x[i][j][k][l]) );
    //       }
    //     }
    //   }
    // }
    // 
    // 
    // fprintf(re, "\n Variável a: \n");
    // for(int i = 0; i < n; i++){
    //   for(int j = 0; j < n; j++){
    //     for(int k = 0; k < n; k++){
    //       if(cplex.getValue(a[i][j][k])>0.001)
    //         fprintf(re, "a[%d][%d][%d] = %f \t \n", i, j, k, cplex.getValue(a[i][j][k]) );
    //     }
    //   }
    // }
    // 
    // fprintf(re, "\n Variável b: \n");
    // for(int i = 0; i < n; i++){
    //   for(int j = 0; j < n; j++){
    //     for(int k = 0; k < n; k++){
    //       if(cplex.getValue(b[i][j][k])>0.001)
    //         fprintf(re, "b[%d][%d][%d] = %f \t \n", i, j, k, cplex.getValue(b[i][j][k]) );
    //     }
    //   }
    // }
    // 
    // fprintf(re, "\n Variável z: \n");
    // for(int i = 0; i < n; i++){
    //   for(int j = 0; j < n; j++){
    //     if(cplex.getValue(z[i][j])>0.001)
    //       fprintf(re, "z[%d][%d] = %f \t \n", i, j, cplex.getValue(z[i][j]) );
    //   }
    // }

  //  fprintf(re, "======================================================================");
    
 
     
    } 
    catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    }
    return 0;

}



