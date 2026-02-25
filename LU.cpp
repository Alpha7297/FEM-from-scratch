#include<cmath>
#include<cstdio>
#include<omp.h>

int np=16;

void swap(double* a,double* b){
    double temp=*a;
    *a=*b;
    *b=temp;
}

int mydgetrf(double *A,int *ipiv,int n,int on,int s){
    for(int sb=0;sb<n;sb+=np){
        int limit=(sb+np<n)?sb+np:n;
        for(int k=sb;k<limit;k++){
            int pivot_row=k;
            double max_val=fabs(A[k*on+k]);
            for(int i=k+1;i<n;i++){
                double val=fabs(A[i*on+k]);
                if(val>max_val){
                    max_val=val;
                    pivot_row=i;
                }
            }
            ipiv[k]=pivot_row;

            if(pivot_row!=k){
                for(int j=sb;j<limit;j++){
                    swap(&A[k*on+j],&A[pivot_row*on+j]);
                }
            }

            double diag=A[k*on+k];
            if(fabs(diag)>1e-18){
                double inv_diag=1.0/diag;
                for(int i=k+1;i<n;i++){
                    A[i*on+k]*=inv_diag;
                }
            }

            for(int i=k+1;i<n;i++){
                double aik=A[i*on+k];
                for(int j=k+1;j<limit;j++){
                    A[i*on+j]-=aik*A[k*on+j];
                }
            }
        }

        #pragma omp parallel for schedule(static)
        for(int b=0;b<2;b++){
            int j_start,j_end;
            if(b==0){ j_start=0; j_end=sb; }
            else{ j_start=limit; j_end=n; }

            for(int k=sb;k<limit;k++){
                if(ipiv[k]!=k){
                    double *row1=A+k*on;
                    double *row2=A+ipiv[k]*on;
                    for(int j=j_start;j<j_end;j++){
                        swap(&row1[j],&row2[j]);
                    }
                }
            }
        }

        for(int i=sb;i<limit;i++){
            #pragma omp parallel for schedule(static)
            for(int k=i+1;k<limit;k++){
                double aik=A[k*on+i];
                for(int j=limit;j<n;j++){
                    A[k*on+j]-=aik*A[i*on+j];
                }
            }
        }

        #pragma omp parallel for schedule(static)
        for(int i=limit;i<n;i++){
            for(int k=sb;k<limit;k++){
                double aik=A[i*on+k];
                for(int j=limit;j<n;j++){
                    A[i*on+j]-=aik*A[k*on+j];
                }
            }
        }
    }
    return 1;
}

void mydtrsv(char UPLO,double *A,double *B,int n,int on,int *ipiv){
    if(UPLO=='L'||UPLO=='l'){
        for(int i=0;i<n;i++){
            double sum=B[i];
            if(i>100){
                #pragma omp parallel for reduction(+:sum) schedule(static)
                for(int j=0;j<i;j++){
                    sum-=A[i*on+j]*B[j];
                }
            }else{
                for(int j=0;j<i;j++){
                    sum-=A[i*on+j]*B[j];
                }
            }
            B[i]=sum;
        }
    }else if(UPLO=='U'||UPLO=='u'){
        for(int i=n-1;i>=0;i--){
            double sum=B[i];
            if(n-i>100){
                #pragma omp parallel for reduction(+:sum) schedule(static)
                for(int j=i+1;j<n;j++){
                    sum-=A[i*on+j]*B[j];
                }
            }else{
                for(int j=i+1;j<n;j++){
                    sum-=A[i*on+j]*B[j];
                }
            }
            B[i]=sum/A[i*on+i];
        }
    }
}

void my_solver(int n,double *A,double *b){
    int *ipiv=(int*)malloc(n*sizeof(int));
    for(int i=0;i<n;i++) ipiv[i]=i;

    if(mydgetrf(A,ipiv,n,n,0)==0){
        free(ipiv);
        printf("error矩阵奇异\n");
        return;
    }

    for(int i=0;i<n;i++){
        if(ipiv[i]!=i){
            swap(&b[i],&b[ipiv[i]]);
        }
    }

    mydtrsv('L',A,b,n,n,ipiv);
    mydtrsv('U',A,b,n,n,ipiv);

    free(ipiv);
}