#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <climits>
#include <limits.h>
#include <stack>
using namespace std;


const double a=-3.4;
const double b=0;
const double c=0.1;
double loop_energy[31][5];
double V[2000][2000],WM[2000][2000];
double WW[2000][2000];
char CC[2000];
int backtrack[2000][2000]; //backtracking for V
int backtrack_WM[2000][2000]; //backtracking for WM
int WM_backtrack_k[2000][2000];


double AU[4][4]={{INT_MAX,INT_MAX,INT_MAX,-0.9}, {INT_MAX,INT_MAX,-2.2,INT_MAX}, {INT_MAX,-2.1,INT_MAX,-0.6},{-1.1,INT_MAX,-1.4,INT_MAX}};
double UA[4][4]={{INT_MAX,INT_MAX,INT_MAX,-1.3}, {INT_MAX,INT_MAX,-2.4,INT_MAX}, {INT_MAX,-2.1,INT_MAX,-1},{-0.9,INT_MAX,-1.3,INT_MAX}};

double CG[4][4]={{INT_MAX,INT_MAX,INT_MAX,-2.1}, {INT_MAX,INT_MAX,-3.3,INT_MAX}, {INT_MAX,-2.4,INT_MAX,-1.4},{-2.1,INT_MAX,-2.1,INT_MAX}};
double GC[4][4]={{INT_MAX,INT_MAX,INT_MAX,-2.4}, {INT_MAX,INT_MAX,-3.4,INT_MAX}, {INT_MAX,-3.3,INT_MAX,-1.5},{-2.2,INT_MAX,-2.5,INT_MAX}};

double GU[4][4]={{INT_MAX,INT_MAX,INT_MAX,-1.3}, {INT_MAX,INT_MAX,-2.5,INT_MAX}, {INT_MAX,-2.1,INT_MAX,-0.5},{-1.4,INT_MAX,1.3,INT_MAX}};
double UG[4][4]={{INT_MAX,INT_MAX,INT_MAX,-1}, {INT_MAX,INT_MAX,-1.5,INT_MAX}, {INT_MAX,-1.4,INT_MAX,.3},{-0.6,INT_MAX,-0.5,INT_MAX}};


string s;

void array_initialization();
void read_loop_energy();
void print_array(int n);
double eH(int a,int b);
double eL(int i,int j, int ii, int jj);
double eS(int i,int j, int x, int y);
void folding(int i, int j);
void foldingV(int i,int j);
void foldingWM(int i,int j);



int main(){

    ifstream fin3;
    fin3.open("sequence.txt");
    getline(fin3,s);
    cout<<s<<endl;
    cout<<s.length()<<endl;

    array_initialization();

    read_loop_energy();


    for(int d=4;d<s.length();d++){

        for(int i=1;i<=s.length()-d;i++){

            int j=i+d;

            //**********************************Recursive function for V array**********************************
            if((s[i-1]=='C' && s[j-1]=='G') || (s[i-1]=='G' && s[j-1]=='C') || (s[i-1]=='A' && s[j-1]=='U') || (s[i-1]=='U' && s[j-1]=='A') || (s[i-1]=='G' && s[j-1]=='U') || (s[i-1]=='U' && s[j-1]=='G')){
                V[i][j]=min(eH(i,j),eS(i,j,i+1,j-1)+V[i+1][j-1]);

                int i_bar;

                double possible_loop_min_val=INT_MAX;
                for(int ii=i+1;ii<j;ii++){
                    for(int jj=ii+1;jj<j;jj++){

                         double temp1=eL(i,j,ii,jj)+V[ii][jj];

                         if(temp1<possible_loop_min_val) {
                            possible_loop_min_val=temp1;
                            i_bar=ii;
                         }
                    }
                }

                V[i][j]=min(V[i][j],possible_loop_min_val);

                int k_bar=-7;

                double possible_multi_loop_min_val=INT_MAX;
                for(int k=i+2;k<j;k++){
                    double temp2=WM[i+1][k-1]+WM[k][j-1]+a;
                    if(temp2<possible_multi_loop_min_val) {
                        possible_multi_loop_min_val=temp2;
                        k_bar=k;
                    }
                }


                V[i][j]=min(V[i][j],possible_multi_loop_min_val);


                 //***************Keeping track of V values in backtrack 2D array**********************************

                if(eH(i,j)<=eS(i,j,i+1,j-1)+V[i+1][j-1] && eH(i,j)<=possible_loop_min_val && eH(i,j)<=possible_multi_loop_min_val) backtrack[i][j]=-1;  //-1 means Hairpin loop
                else if(eS(i,j,i+1,j-1)+V[i+1][j-1]<=eH(i,j) && eS(i,j,i+1,j-1)+V[i+1][j-1]<=possible_loop_min_val && eS(i,j,i+1,j-1)+V[i+1][j-1]<=possible_multi_loop_min_val) backtrack[i][j]=-2;  //-2 means stack
                else if (possible_loop_min_val<=eH(i,j) && possible_loop_min_val<=eS(i,j,i+1,j-1)+V[i+1][j-1] && possible_loop_min_val<=possible_multi_loop_min_val) backtrack[i][j]=i_bar;
                else if(possible_multi_loop_min_val<=eH(i,j) && possible_multi_loop_min_val<=eS(i,j,i+1,j-1)+V[i+1][j-1] && possible_multi_loop_min_val<=possible_loop_min_val) {

                        backtrack[i][j]=-6;  //means minimum value came from multi-loop
                        WM_backtrack_k[i][j]=k_bar;
                }

            }
            else{
                V[i][j]=INT_MAX;

            }


            //**********************************Recursive function for WM array**********************************


            WM[i][j]=min(V[i][j]+b,WM[i][j-1]+c);
            WM[i][j]=min(WM[i][j],WM[i+1][j]+c);

            double possible_WM_min_val=INT_MAX;
            int kk_bar=-9;
            for(int k=i+1;k<=j;k++){
                double temp3=WM[i][k-1]+WM[k][j];
                if(temp3<possible_WM_min_val) {
                    kk_bar==k;
                    possible_WM_min_val=temp3;
                }
            }

            WM[i][j]=min(WM[i][j],possible_WM_min_val);


            //***************Keeping track of WM values in WM_backtrack 2D array**********************************

                if(V[i][j]+b<=WM[i][j-1]+c && V[i][j]+b<=WM[i+1][j]+c && V[i][j]+b<=possible_WM_min_val) backtrack_WM[i][j]=-1;  //-1 means V[i][j]+b<
                else if(WM[i][j-1]+c<=V[i][j]+b && WM[i][j-1]+c<=WM[i+1][j]+c && WM[i][j-1]+c<=possible_WM_min_val) backtrack_WM[i][j]=-2;  //-2 means WM[i][j-1]+c
                else if (WM[i+1][j]+c<=V[i][j]+b && WM[i+1][j]+c<=WM[i][j-1]+c && WM[i+1][j]+c<=possible_WM_min_val) backtrack_WM[i][j]=-3; //-3 means WM[i+1][j]+c
                else if(possible_WM_min_val<=V[i][j] && possible_WM_min_val<=WM[i+1][j]+c && possible_WM_min_val<=WM[i][j-1]+c) {

                        backtrack_WM[i][j]=kk_bar;  //means minimum value came from multi-multi-loop

                }

               //  if(i==27 && j==45) printf("%lf %lf %lf %lf %d\n",V[i][j]+b,WM[i][j-1]+c, WM[i+1][j]+c, possible_WM_min_val, backtrack_WM[i][j]);


            //********************************Recursive function For WW********************************

            WW[i][j]=min(WW[i+1][j],WW[i][j-1]);

            WW[i][j]=min(WW[i][j],V[i][j]);

            double possible_WW_min_val=INT_MAX;
            double kk=-5;
            for(int k=i+1;k<j;k++){
                double temp4=WW[i][k]+WW[k+1][j];

                if(temp4<possible_WW_min_val){
                        possible_WW_min_val=temp4;
                        kk=k;
                }
            }


            WW[i][j]=min(WW[i][j],possible_WW_min_val);

        }


    }


    //********************************Recursive function For W********************************

    cout<<endl;

    cout<<"Final Result: "<<WW[1][s.length()]<<endl;

    folding(1,s.length());

    for(int x=1;x<=s.length();x++){
        cout<<CC[x]<<"   ";
     }
    cout<<endl;

    ofstream fout;
    fout.open ("suboptstructure.txt");
    fout<<s<<endl;
    for(int x=1;x<=s.length();x++){
        fout<<CC[x];
     }
    fout<<endl;
    fout.close();
    fin3.close();

    // ***** Generate Plot using ViennaRNA Package ***** //
    //system("RNAplot < suboptstructure.txt");

    return 0;

}

void folding(int i, int j){

    if(i==j) {
        CC[i]='.';
        return;
    }
    if(V[i][j]==WW[i][j]){

            foldingV(i,j);
    }

   else if(WW[i][j-1]==WW[i][j]){
            CC[j]='.';
            j=j-1;
            folding(i,j);

    }
    else if(WW[i+1][j]==WW[i][j]){
            CC[i]='.';
            i=i+1;
            folding(i,j);
    }

    else{
        int index=-1;
        for(int k=i+1;k<j;k++){

            double temp6=WW[i][k]+WW[k+1][j];

            if(temp6 == WW[i][j]){
                index=k;

                break;
            }
        }

        folding(i,index);
        folding(index+1,j);
    }

}


void foldingV(int i,int j){

    if(i==j) {
        CC[i]='.';
        return;
    }
    else if(backtrack[i][j]==-1){
        CC[j]=')';
        CC[i]='(';
        for(int k=i+1;k<j;k++){
              CC[k]='.';
        }
        foldingV(i+1,j-1);
    }
    else if(backtrack[i][j]==-2){
        CC[j]=')';
        CC[i]='(';
        foldingV(i+1,j-1);
    }
    else if(backtrack[i][j]>0 && backtrack[i][j]<INT_MAX){
        CC[j]=')';
        CC[i]='(';

        int possible_loop_min_val=INT_MAX;
        int ii=backtrack[i][j];
        int jbar=-5;

        for(int jj=ii+1;jj<j;jj++){
            double temp1=eL(i,j,ii,jj)+V[ii][jj];
            if(temp1<possible_loop_min_val) {
                possible_loop_min_val=temp1;
                jbar=jj;
            }
        }

        for(int k=i+1;k<ii;k++){
              CC[k]='.';
        }


        for(int k=jbar+1;k<j;k++){
              CC[k]='.';
        }

        foldingV(ii,jbar);

    }
    else if(backtrack[i][j]==-6){
        CC[j]=')';
        CC[i]='(';
        int kk=WM_backtrack_k[i][j];
        foldingWM(i+1,kk-1);
        foldingWM(kk,j-1);
    }

}


void foldingWM(int i,int j){

    if(i==j) {
        CC[i]='.';
        return;
    }
    else if(backtrack_WM[i][j]==-1){
        foldingV(i,j);

    }
    else if(backtrack_WM[i][j]==-2){
        CC[j]='.';
        foldingWM(i,j-1);
    }
    else if(backtrack_WM[i][j]==-3){
        CC[i]='.';
        foldingWM(i+1,j);

    }
    else if(backtrack_WM[i][j]>0 && backtrack_WM[i][j]<INT_MAX){

        int k_k_bar=backtrack_WM[i][j];
        foldingWM(i,k_k_bar-1);
        foldingWM(k_k_bar,j);

    }

}



void array_initialization(){

    for(int i=1;i<=s.length();i++){
        for(int j=1;j<=s.length();j++){
            if(j-i<4){
                V[i][j]=INT_MAX;
                WM[i][j]=INT_MAX;
                WW[i][j]=INT_MAX;
            }
            backtrack[i][j]=0;
            backtrack_WM[i][j]=0;

        }

    }

}

void read_loop_energy(){
    ifstream fin1;
    fin1.open("loop_energy.txt");
    for(int i=1;i<=30;i++){
        for(int j=1;j<=4;j++){
            fin1>>loop_energy[i][j];
        }
    }
}

double eH(int a,int b){

    int dif=b-a-1;
    if(dif==0) return INT_MAX;
    else return loop_energy[dif][4];

}

double eL(int i,int j, int ii, int jj){


    if((s[ii-1]=='C' && s[jj-1]=='G') || (s[ii-1]=='G' && s[jj-1]=='C') || (s[ii-1]=='A' && s[jj-1]=='U') || (s[ii-1]=='U' && s[jj-1]=='A') || (s[ii-1]=='G' && s[jj-1]=='U') || (s[ii-1]=='U' && s[jj-1]=='G')){
        int n1=ii-i-1;
        int n2=j-jj-1;

        if(n1+n2>=1){       //***********loop exists************

            if((n1==0 && n2>=1) || (n1>=1 && n2==0)){        //***********Bulge exists************
                return loop_energy[n1+n2][3];
            }
            else{                  //***********Internal loop exists************
                return loop_energy[n1+n2][2];
            }
        }
        else return INT_MAX;

    }
    else return INT_MAX;

}

double eS(int i,int j,int x,int y){

    if(j-i>2){
        if(s[i-1]=='A' && s[j-1]=='U'){


            if(s[x-1]=='A' && s[y-1]=='U'){
                 return AU[0][3];
            }
            else if(s[x-1]=='U' && s[y-1]=='A'){
                 return AU[3][0];
            }
            else if(s[x-1]=='C' && s[y-1]=='G'){
                 return AU[1][2];
            }
            else if(s[x-1]=='G' && s[y-1]=='C'){
                 return AU[2][1];
            }
            else if(s[x-1]=='G' && s[y-1]=='U'){
                 return AU[2][3];
            }
            else if(s[x-1]=='U' && s[y-1]=='G'){
                 return AU[3][2];
            }
            else return INT_MAX;


        }
        else if(s[i-1]=='U' && s[j-1]=='A'){
            if(s[x-1]=='A' && s[y-1]=='U'){
                 return UA[0][3];
            }
            else if(s[x-1]=='U' && s[y-1]=='A'){
                 return UA[3][0];
            }
            else if(s[x-1]=='C' && s[y-1]=='G'){
                 return UA[1][2];
            }
            else if(s[x-1]=='G' && s[y-1]=='C'){
                 return UA[2][1];
            }
            else if(s[x-1]=='G' && s[y-1]=='U'){
                 return UA[2][3];
            }
            else if(s[x-1]=='U' && s[y-1]=='G'){
                 return UA[3][2];
            }
            else return INT_MAX;
        }
        else if(s[i-1]=='C' && s[j-1]=='G'){
            if(s[x-1]=='A' && s[y-1]=='U'){
                 return CG[0][3];
            }
            else if(s[x-1]=='U' && s[y-1]=='A'){
                 return CG[3][0];
            }
            else if(s[x-1]=='C' && s[y-1]=='G'){
                 return CG[1][2];
            }
            else if(s[x-1]=='G' && s[y-1]=='C'){
                 return CG[2][1];
            }
            else if(s[x-1]=='G' && s[y-1]=='U'){
                 return CG[2][3];
            }
            else if(s[x-1]=='U' && s[y-1]=='G'){
                 return CG[3][2];
            }
            else return INT_MAX;

        }
        else if(s[i-1]=='G' && s[j-1]=='C'){
            if(s[x-1]=='A' && s[y-1]=='U'){
                 return GC[0][3];
            }
            else if(s[x-1]=='U' && s[y-1]=='A'){
                 return GC[3][0];
            }
            else if(s[x-1]=='C' && s[y-1]=='G'){
                 return GC[1][2];
            }
            else if(s[x-1]=='G' && s[y-1]=='C'){
                 return GC[2][1];
            }
            else if(s[x-1]=='G' && s[y-1]=='U'){
                 return GC[2][3];
            }
            else if(s[x-1]=='U' && s[y-1]=='G'){
                 return GC[3][2];
            }
            else return INT_MAX;

        }
        else if(s[i-1]=='G' && s[j-1]=='U'){
            if(s[x-1]=='A' && s[y-1]=='U'){
                 return GU[0][3];
            }
            else if(s[x-1]=='U' && s[y-1]=='A'){
                 return GU[3][0];
            }
            else if(s[x-1]=='C' && s[y-1]=='G'){
                 return GU[1][2];
            }
            else if(s[x-1]=='G' && s[y-1]=='C'){
                 return GU[2][1];
            }
            else if(s[x-1]=='G' && s[y-1]=='U'){
                 return GU[2][3];
            }
            else if(s[x-1]=='U' && s[y-1]=='G'){
                 return GU[3][2];
            }
            else return INT_MAX;

        }
        else if(s[i-1]=='U' && s[j-1]=='G'){
            if(s[x-1]=='A' && s[y-1]=='U'){
                 return UG[0][3];
            }
            else if(s[x-1]=='U' && s[y-1]=='A'){
                 return UG[3][0];
            }
            else if(s[x-1]=='C' && s[y-1]=='G'){
                 return UG[1][2];
            }
            else if(s[x-1]=='G' && s[y-1]=='C'){
                 return UG[2][1];
            }
            else if(s[x-1]=='G' && s[y-1]=='U'){
                 return UG[2][3];
            }
            else if(s[x-1]=='U' && s[y-1]=='G'){
                 return UG[3][2];
            }
            else return INT_MAX;

        }


    }
    else return INT_MAX;

}

void print_array(int n){

    for(int i=1;i<=n;i++){
        for(int j=5;j<=n;j++){
            if(backtrack_WM[i][j]==INT_MAX)
                cout<<"I ";
            else
                cout<<backtrack_WM[i][j]<<"   ";
        }
        cout<<endl<<endl;
    }

     cout<<endl<<endl<<endl<<endl<<endl<<endl;

}



