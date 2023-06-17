#include "mex.h" 
#include <vector>
#include<algorithm>
#include<cmath>
#include<numeric>
#include<iostream>
using namespace std;


void isNeighbor(int SET1,vector<int> SET2,const mxArray *prhs[],double tolerance,bool &Flag)

{
    double *NeighborMatrix;  
    NeighborMatrix = mxGetPr(prhs[0]);
    int rNeighborMatrix = (int)mxGetM(prhs[0]);   //Get row
    int cNeighborMatrix = (int)mxGetN(prhs[0]);   //Get col 
    
    double *BaseRegionCenter; 
    BaseRegionCenter = mxGetPr(prhs[3]); 
    int rBaseRegionCenter = (int)mxGetM(prhs[3]);   
    int cBaseRegionCenter = (int)mxGetN(prhs[3]);  
    
    vector<double> X1;
    X1.resize(cBaseRegionCenter);    //initialize X1
    int TmpRow = SET1-1;
    for (int i=0;i<cBaseRegionCenter;i++)   X1[i] = (double) BaseRegionCenter[ i* rBaseRegionCenter + TmpRow];  
   
    vector<vector<double> > X2;              //initialize X2
    X2.resize(SET2.size());
    for (int i=0;i<SET2.size();i++)  X2[i].resize(cBaseRegionCenter);
    
    
   //---------------
    for(int i=0;i<X2.size();i++)       //X2 = X2-X1
    {
        TmpRow = SET2[i]-1;
        for(int j=0;j<X2[0].size();j++)
        { 
            X2[i][j] = (double)BaseRegionCenter[j* rBaseRegionCenter + TmpRow]-X1[j];
        };
    }  
    vector<double> Dist;   //distance between X1 and X2;
    Dist.resize(SET2.size()); 
    for(int i=0;i<X2.size();i++){
        Dist[i] = 0.0;
        for (int j=0;j<X2[0].size();j++)     Dist[i]+=X2[i][j]*X2[i][j];
        Dist[i] = sqrt(Dist[i]);        
    }

    vector<int> Pos;
    for(int i=0;i<SET2.size();i++)
      if (Dist[i] < tolerance) Pos.push_back(i);
    
    if (Pos.size()>0)  
    {
        int sum = 0;
        for (int i = 0; i<Pos.size();i++) 
        {   int TmpRow = SET1-1,TmpCol = SET2[Pos[i]]-1;
            sum+=(int)NeighborMatrix[ TmpCol*rNeighborMatrix + TmpRow ];
        }
        if (sum>0) Flag = true;        
     }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//nlhs: number of output parameters
	//plhs: output pointer arrays
	//nrhs: number of input parameters
	//prhs: output pointer arrays 
	double tolerance = mxGetScalar(prhs[4]); 
       
    double *NeighborMatrix;    
// ----------- Get NeighborMatrix ----------------     
    NeighborMatrix = mxGetPr(prhs[0]); 
    int rNeighborMatrix = (int)mxGetM(prhs[0]);   //row
    int cNeighborMatrix = (int)mxGetN(prhs[0]);   //col 
    

    
// ------------ Get ConvexRegion ----------------
    double *ConvexRegionx;    
    ConvexRegionx = mxGetPr(prhs[1]); 
    int rConvexRegionx = (int)mxGetM(prhs[1]);            

  
// ------------ Get LastConvexRegion -------------
    double *LastConvexRegion; 
    LastConvexRegion = mxGetPr(prhs[2]); 
    int rLastConvexRegion = (int)mxGetM(prhs[2]);  
    int cLastConvexRegion = (int)mxGetN(prhs[2]); 
    
// ------------ Get BaseRegionCenter -------------
    double *BaseRegionCenter; 
    BaseRegionCenter = mxGetPr(prhs[3]); 
    int rBaseRegionCenter = (int)mxGetM(prhs[3]);   
    int cBaseRegionCenter = (int)mxGetN(prhs[3]);         
       

// ---------------------- Initialization finished --------------------------
     int k = 0,i=0,j=0; 
     int ii = 0, jj = 0;
     vector<vector<double> > SubNeighborMatrix;
     int rSubNeighborMatrix = cLastConvexRegion;
     int cSubNeighborMatrix = cNeighborMatrix;
     SubNeighborMatrix.resize(rSubNeighborMatrix);
     for (i=0;i<rSubNeighborMatrix;i++) SubNeighborMatrix[i].resize(cSubNeighborMatrix);

     vector<vector<int> > TMP; 
     for (i=0;i<rLastConvexRegion;i++)    
     {
         vector<int> set2;
         set2.resize(cLastConvexRegion);
         for (j=0;j<cLastConvexRegion;j++) set2[j] = (int)LastConvexRegion[j * rLastConvexRegion + i];
         vector<int> SET2(set2); 

         
         for(ii=0;ii<rSubNeighborMatrix;ii++)
             for(jj=0;jj<cSubNeighborMatrix;jj++)  
             { int TmpRow = set2[ii]-1; SubNeighborMatrix[ii][jj] = NeighborMatrix[jj * rNeighborMatrix + TmpRow];}
         
         vector<int> GetNeighbors;
         for(ii=0;ii<rSubNeighborMatrix;ii++)
             for(jj=0;jj<cSubNeighborMatrix;jj++)  
                 if (SubNeighborMatrix[ii][jj])   GetNeighbors.push_back(jj+1);
         
         sort(GetNeighbors.begin(),GetNeighbors.end());
         vector<int>::iterator pos = unique(GetNeighbors.begin(),GetNeighbors.end());
         GetNeighbors.erase(pos,GetNeighbors.end());
         
         sort(SET2.begin(),SET2.end());
         vector<int>::iterator iter = set_difference(GetNeighbors.begin(),GetNeighbors.end(),SET2.begin(),SET2.end(),GetNeighbors.begin());
         GetNeighbors.resize(iter-GetNeighbors.begin()); 

        int Len = (int)GetNeighbors.size();      
         for (j=0;j<Len;j++)
         {
              int SET1 = (int)ConvexRegionx[GetNeighbors[j]-1];             
              vector<int>::iterator iter = find(SET2.begin(),SET2.end(),SET1);
              if (iter != SET2.end())   continue;
            
              bool Flag = false;  
              
           isNeighbor(SET1,SET2,prhs,tolerance,Flag);              
           if(Flag)
           { 
             k = k+1;
             TMP.resize(k);
             TMP[k-1].resize(SET2.size()+1);
             TMP[k-1][0] = SET1;
             for (jj = 1;jj<SET2.size()+1;jj++)   TMP[k-1][jj] = SET2[jj-1];
             sort(TMP[k-1].begin(),TMP[k-1].end());
            }
         
     }
     
    }  

 int row = (int)TMP.size(); 
     if (row<1)
     {
         plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
         double *y = mxGetPr(plhs[0]);
         y[0] = -1;
     }
     else
     {
        int col = (int)TMP[0].size();
         plhs[0] = mxCreateDoubleMatrix(row, col, mxREAL);
         double *y = mxGetPr(plhs[0]); //获得矩阵的第一个元素的指针
         for(int i = 0; i < row; i++) 
             for(int j = 0; j < col; j++)
                 y[j*row+i] = TMP[i][j]; 
       }
 }

