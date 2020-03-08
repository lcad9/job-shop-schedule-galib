// ************************************************************
//   Exemplo de leitura de uma inst�ncia de problemas especificos para job shop scheduling
//   Este programa l� um arquivo texto e imprime o seu
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
//***********************Fim tratamento de excess�es*********************************/
/////////////////////////////////////////////////////////////////////////////////////

	//arq = fopen(diretorio[p-1].c_str(), "rt");
	
	arq = fopen("orb03", "rt");

	//i = 1;
  
  	int cont = 0; //Conta o n�mero de linhas
  
  	while (!feof(arq))
  	{  	
		// L� uma linha (inclusive com o '\n')
    	result = fgets(Linha, 100, arq);  // o 'fgets' l� at� 99 caracteres ou at� o '\n'
      	
      	///L�-se a partir da linha 3, por causa dos formatos
      	if(cont > 3){
			if (result){
				pch = strtok(Linha, " ");//Armazena o conte�do da linha
			
				while(pch != NULL)
				{
					maq.push_back(atoi(pch) + 1);//Armazena o n�mero da m�quina +1
					pch = strtok(NULL, " ");
					
					T.push_back(atoi(pch)); //Armazena o n�mero do tempo
					pch = strtok(NULL, " ");			
				}
		  	}		  
	  	}	  
		cont++;
  	}
  	
	fclose(arq); //Fecha arquivo da inst�ncia do problema

  	int m = maq[0];
  	int t = T[0];
      
  	cout << m << " x " << t; //Demostra tamanho da inst�ncia do problema
      
  	int JT[m][t];//Matriz onde ser� armazenado o tempo, separado por m�quina

  	for(int i = 0; i < m; i++)
  		for(int j = 0; j < t; j++){
  			JT[i][j] = T[(j + 1) + (i * t) ];
  		}
  		
  	//////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	  	
	  	/*
	  	** A partir daqui � somente exibi��o do conte�do
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
