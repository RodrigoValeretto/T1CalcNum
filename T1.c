#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define TAM 500
#define ERRO pow(10, -8)

int calculaErro(double *K1, double* K, double erro)
{
    double vet[TAM];
    double maior;
    double moduloV;

    for (int i = 0; i < TAM; i++)
        vet[i] = (K1[i] - K[i]);

    maior = fabs(vet[0]);
    for(int i = 1; i < TAM; i++)
    {
        moduloV = fabs(vet[i]);
        if(moduloV > maior)
            maior = moduloV;
    }

    if(maior < erro)
        return 1;
    else
        return 0;
    
    
}

double ** criaMatriz()
{
    double **matriz = (double**)calloc(TAM,sizeof(double));
    
    for(int i = 0; i < TAM; i++)
        matriz[i] = (double*)calloc(TAM,sizeof(double));

    for(int i = 0; i < TAM; i++)
        matriz[i][i] = 5;

    for(int i = 0; i < TAM-1; i++)
        {
            matriz[i][i+1] = -1;
            matriz[i+1][i] = -1;
        }

    for(int i = 0; i < TAM-3; i++)
    {
        matriz[i+3][i] = -1;
        matriz[i][i+3] = -1;
    }

    return matriz;

}

double * gseidel(double ** matriz,double * b, double *K, double * K1)
{
    int i = 0;
    double somatorio;
    double somatorio2;

    for(int kmax = 0; kmax < 5*TAM; kmax++)
    {
        somatorio = 0;
        i = 0;
        for(int j = i+1; j < TAM; j++)
        {
            somatorio = somatorio + matriz[i][j]*K[j];
        }
        
        K1[i] = (b[i] - somatorio)/matriz[i][i];
        for(i = 1; i < TAM ; i++)
        {
            somatorio = 0;
            somatorio2 = 0;
            for(int j = i+1; j < TAM; j++)
            {
                somatorio = somatorio + matriz[i][j]*K[j];
            }

            for(int j = 0; j <= i-1; j++)
            {
                somatorio2 = somatorio2 + matriz[i][j]*K1[i];
            }

            K1[i] = (b[i] - somatorio2 - somatorio)/matriz[i][i];
        }

        if(calculaErro(K1, K, ERRO) == 1)
            return K1;
        else
            if(kmax > 0)
                for(int j = 0; j < TAM; j++)
                    K[j] = K1[j];
    }

    return NULL;
}

void printmatriz(double **matriz)
{
    for(int i = 0; i < TAM; i++)
    {
        for(int j = 0; j < TAM; j++)
            printf(" %lf ", matriz[i][j]);
        printf("\n");
    }
}

int main()
{
    double ** matriz = criaMatriz();
    double * b = (double*)calloc(TAM, sizeof(double));
    double * K = (double*)calloc(TAM, sizeof(double));
    double * K1 = (double*)calloc(TAM, sizeof(double));
    double somatorio = 0;

/*
    for(int i = 0; i <  TAM; i++)
    {
        somatorio = 0;
        for(int j = 0; j < TAM; j++)
            somatorio = somatorio + matriz[i][j];
        
        b[i] = somatorio;
    }
*/
    for(int i = 0; i < TAM; i++)
        b[i] = (1.0/(i+1));

    K1 = gseidel(matriz, b, K, K1);

    for(int i = 0; i < TAM; i++)
        printf(" %lf ", K1[i]);

}