// ************************************************************
//   Exemplo de leitura de uma instância de problemas especificos para job shop scheduling
//   Este programa lê um arquivo texto e imprime o seu
//   conteudo na tela.
// ************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <dirent.h>

using namespace std;

int main()
{
	FILE *arq;
  	vector<int> maq;
  	vector<int> T;
  	char Linha[100];
  	char *result;
  	char *pch;
  
   
//////////////////////////////////////////////////////////////////////////////////////
//***********************Fim tratamento de excessões*********************************/
/////////////////////////////////////////////////////////////////////////////////////

	//arq = fopen(diretorio[p-1].c_str(), "rt");
	
	arq = fopen("orb03", "rt");

	//i = 1;
  
  	int cont = 0; //Conta o número de linhas
  
  	while (!feof(arq))
  	{  	
		// Lê uma linha (inclusive com o '\n')
    	result = fgets(Linha, 100, arq);  // o 'fgets' lê até 99 caracteres ou até o '\n'
      	
      	///Lê-se a partir da linha 3, por causa dos formatos
      	if(cont > 3){
			if (result){
				pch = strtok(Linha, " ");//Armazena o conteúdo da linha
			
				while(pch != NULL)
				{
					maq.push_back(atoi(pch) + 1);//Armazena o número da máquina +1
					pch = strtok(NULL, " ");
					
					T.push_back(atoi(pch)); //Armazena o número do tempo
					pch = strtok(NULL, " ");			
				}
		  	}		  
	  	}	  
		cont++;
  	}
  	
	fclose(arq); //Fecha arquivo da instância do problema

  	int m = maq[0];
  	int t = T[0];
      
  	cout << m << " x " << t; //Demostra tamanho da instância do problema
      
  	int JT[m][t];//Matriz onde será armazenado o tempo, separado por máquina

  	for(int i = 0; i < m; i++)
  		for(int j = 0; j < t; j++){
  			JT[i][j] = T[(j + 1) + (i * t) ];
  		}
  		
  	//////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	  	
	  	/*
	  	** A partir daqui é somente exibição do conteúdo
	  	*/
	  	
  	cout << "\n\n";
	
	cout << "Maquina: ";
		  
	maq.erase(maq.begin());
	T.erase(T.begin());
		  
	for(std::vector<int>::iterator it = maq.begin(); it != maq.end(); ++it)
		cout << ' ' << *it;
		
	cout << "\n\n";
	
	cout << "Job: ";
			
	for(std::vector<int>::iterator it = T.begin(); it != T.end(); ++it)
		cout << ' ' << *it;
  
  	cout << "\n";
  
    for(int i = 0; i < m; i++){
    	for(int j = 0; j < t; j++){
	  		cout << JT[i][j];
	  		cout << " ";
	  	}
	  	cout << "\n";
    }
	  	
  	system("pause");
  
  	return 0;
}
