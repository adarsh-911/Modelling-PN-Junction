#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N 90
#define iter 5000
#define L 4e-7
#define epsi 8.854e-12
#define er 11.7
#define ev 1.602e-19
#define Na 1e21
#define Nd 1e21
#define Vapp -1
#define Kb 1.38e-23
#define T 300
#define Nc 2.8e25
#define Nv 1.04e25
#define Ni 1.5e16
#define muN 0.135 //0.02835
#define muP 0.048 //0.1339
#define tau_n 1e-5
#define tau_p 1e-5
#define delta_x L/(N-1)
#define delta_t 1e-9
#define Eg 1.12
#define smol 1e-9

long double phi(long double v){  
    return ev*v/(Kb*T);
}

void Mesh(long double m[], long double position[],int n){
      position[0]=0;
      for(int i=0;i<n-1;i++){
            m[i] = L*(1-exp(-(i+1)*log(2)/(n-1)));
            position[i] = position[i-1] + m[i-1];
            position[n-1]+=m[i];
        }       
}

long double Bern(long double pi){
    return pi/(exp(pi)-1);
  }

long double F(long double v_prev, long double v_pres, long double v_next, long double n, long double p, long double mes0, long double mes1, int side){
      long double func;
      if (side == 0) func = 2*((v_prev*mes1 - v_pres*(mes0+mes1) + v_next*mes0)/(mes1*mes0*(mes0+mes1))) + ev*(-Na + p - n)/(er*epsi);
      else if(side == 1) func = 2*((v_prev*mes1 - v_pres*(mes0+mes1) + v_next*mes0)/(mes1*mes0*(mes0+mes1))) + ev*(Nd + p - n)/(er*epsi);
    return func;  
}

int main(){
      
      FILE *n_file, *p_file;
      long double n_eq[N], p_eq[N];

      n_file = fopen("n.txt", "r");
      p_file = fopen("p.txt", "r");
      for (int i = 0; i < N; i++) {
          if (fscanf(n_file, "%Lf", &n_eq[i]) != 1) {
              fclose(n_file);
          }
          if (fscanf(p_file, "%Lf", &p_eq[i]) != 1) {
              fclose(p_file);
          }
      }
      
      fclose(n_file);
      fclose(p_file);
      
      long double n[N] = {0};
      long double p[N] = {0};
      long double V[N] = {0};
      
      V[0] = Vapp - Eg - (Kb*T/ev)*log(Na/Nv);
      V[N-1] = (Kb*T/ev)*log(Nd/Nc);
      
      long double jacob_matrix[N][N] = {0};
      long double deltaV_matrix[N] = {0};
      long double F_matrix[N] = {0};
      long double p_func[N] = {0};
      long double n_func[N] = {0};
      long double Jp_matrix[N][N] = {0};
      long double Jn_matrix[N][N] = {0}; 
      long double deltaP_matrix[N] = {0};
      long double deltaN_matrix[N] = {0};
      long double mes[N] = {0};
      long double x[N] = {0};
      
      Mesh(mes,x,N);
      
      for(int a=0;a<iter;a++){
      
            for(int b=1;b<N-1;b++){
                  if(b<=(int)((N+1)/2)) F_matrix[b] = (-1)*(F(V[b-1],V[b],V[b+1],n[b],p[b],mes[b-1],mes[b],0));
                  else F_matrix[b] = (-1)*(F(V[b-1],V[b],V[b+1],n[b],p[b],mes[b-1],mes[b],1));
                  p_func[b] = p_eq[b];
                  n_func[b] = -n_eq[b];
                }
                
            F_matrix[0] = 0;
            F_matrix[N-1] = 0;
            
            p_func[0] = Na;
            p_func[N-1] = pow(Ni,2)/Nd;
            
            n_func[0] = pow(Ni,2)/Na;
            n_func[N-1] = Nd;
                
            jacob_matrix[0][0] = 1;
            jacob_matrix[N-1][N-1] = 1;          
                  
            Jp_matrix[0][0] = 1;
            Jp_matrix[N-1][N-1] = 1;
            
            Jn_matrix[0][0] = 1;
            Jn_matrix[N-1][N-1] = 1; 
            
            for(int i=1;i<N-1;i++){
                  for(int j=0;j<N;j++){
                        if(i-j == 1) jacob_matrix[i][j] = 2/(mes[i-1]*(mes[i]+mes[i-1]));                          
                        else if(j-i == 1) jacob_matrix[i][j] = 2/(mes[i]*(mes[i]+mes[i-1]));
                        else if(i == j) jacob_matrix[i][j] = -2/(mes[i]*mes[i-1]) - pow(ev,2)*(n[i] + p[i])/(er*epsi*Kb*T);                                                   
                        }
                    }
            
            /*for(int i=0;i<N;i++){
                  for(int j=0;j<N;j++){
                        printf("%Le ",jacob_matrix[i][j]);
                    }printf("\n");
              }*/
            
            for(int i=0;i<N-1;i++){      
                  for(int j=i+1;j<N;j++){              
                        long double factor = jacob_matrix[j][i]/jacob_matrix[i][i];
                        for(int k=i;k<N;k++){
                              jacob_matrix[j][k] -= factor*jacob_matrix[i][k];
                          }
                        F_matrix[j] -= factor*F_matrix[i];
                    }
              }
                
            for(int p=N-1;p>=0;p--){
                  float sum = 0;
                  for(int q=p+1;q<N;q++){
                        sum += jacob_matrix[p][q]*deltaV_matrix[q];     
                    }
                  deltaV_matrix[p] = (F_matrix[p] - sum)/jacob_matrix[p][p];
              }
            
            deltaV_matrix[0] = 0;
            deltaV_matrix[N-1] = 0;
            
            for(int i=0;i<N;i++){
                  V[i] += deltaV_matrix[i];
              }                      
            
            for(int i=1;i<N-1;i++){
                  for(int j=0;j<N;j++){
                              if(i-j==1){
                                    Jp_matrix[i][j] = (-2)*Kb*T*muP*tau_p*Bern(phi(V[i])-phi(V[i-1]))/(ev*(mes[i]+mes[i-1])*mes[i-1]);
                                    Jn_matrix[i][j] = 2*Kb*T*muN*tau_n*Bern(phi(V[i-1])-phi(V[i]))/(ev*mes[i-1]*(mes[i]+mes[i-1]));
                                }
                              else if(j-i==1){
                                    Jp_matrix[i][j] = (-2)*Kb*T*muP*tau_p*Bern(phi(V[i])-phi(V[i+1]))/(ev*(mes[i]+mes[i-1])*mes[i]);
                                    Jn_matrix[i][j] = 2*Kb*T*muN*tau_n*Bern(phi(V[i+1])-phi(V[i]))/(ev*mes[i]*(mes[i]+mes[i-1]));
                                }
                              else if(i==j){
                                    Jp_matrix[i][j] = 2*Kb*T*muP*tau_p*(Bern(phi(V[i+1])-phi(V[i]))/mes[i] + Bern(phi(V[i-1])-phi(V[i]))/mes[i-1])/(ev*(mes[i]+mes[i-1])) + 1;
                                    Jn_matrix[i][j] = (-2)*Kb*T*muN*tau_n*(Bern(phi(V[i])-phi(V[i+1]))/mes[i] + Bern(phi(V[i])-phi(V[i-1]))/mes[i-1])/((mes[i]+mes[i-1])*ev) - 1;
                                }                       
                        }
                    }
            
            for(int i=0;i<N-1;i++){      
                  for(int j=i+1;j<N;j++){
                        long double factor1 = Jp_matrix[j][i]/Jp_matrix[i][i];
                        long double factor2 = Jn_matrix[j][i]/Jn_matrix[i][i];
                        for(int k=i;k<N;k++){
                              Jp_matrix[j][k] -= factor1*Jp_matrix[i][k];
                              Jn_matrix[j][k] -= factor2*Jn_matrix[i][k];
                          }
                        p_func[j] -= factor1*p_func[i];
                        n_func[j] -= factor2*n_func[i];
                    }
              }
              
            for(int m=N-1;m>=0;m--){
                  long double sum1 = 0;
                  long double sum2 = 0;
                  for(int o=m+1;o<N;o++){
                        sum1 += Jp_matrix[m][o]*deltaP_matrix[o];
                        sum2 += Jn_matrix[m][o]*deltaN_matrix[o];
                    }
                  deltaP_matrix[m] = (p_func[m] - sum1)/Jp_matrix[m][m];
                  deltaN_matrix[m] = (n_func[m] - sum2)/Jn_matrix[m][m];
              }
            
            for(int i=0;i<N;i++){
                  p[i] = deltaP_matrix[i];
                  n[i] = deltaN_matrix[i];
                  printf("%Le\t%Le\t%Le\t%Le\n",n[i],p[i],V[i],x[i]);
            }               
            
            for(int i=0;i<N;i++){
                  for(int j=0;j<N;j++){
                        jacob_matrix[i][j]=0;
                        Jp_matrix[i][j]=0;
                        Jn_matrix[i][j]=0;
                    }
              }
        }
      
      for(int i=0;i<N;i++){
            printf("%Le\t%Le\t%Le\n",n[i],p[i],V[i]);
      }
             
      FILE *nfile = fopen("n_ss.txt", "w");
      FILE *pfile = fopen("p_ss.txt", "w");
      FILE *xfile = fopen("x_ss.txt", "w");
      FILE *vfile = fopen("v_ss.txt", "w");

      for(int i = 0; i < N; i++) {
            fprintf(nfile, "%Le\n",n[i]);
            fprintf(pfile, "%Le\n",p[i]);
            fprintf(xfile, "%Le\n",x[i]);
            fprintf(vfile, "%Le\n",V[i]);
            
            /*fprintf(nfile, "%f\n",log(n[i]));
            fprintf(pfile, "%f\n",log(p[i]));
            fprintf(xfile, "%Le\n",x[i]);
            fprintf(vfile, "%Le\n",V[i]);*/
        }

      fclose(nfile);
      fclose(pfile);
      fclose(xfile);
      fclose(vfile);
      
  return 0;
}
