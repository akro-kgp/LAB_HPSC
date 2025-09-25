#include<bits/stdc++.h>
using namespace std;

int main(){
    int m,n;
    cin>>m>>n;
    vector<vector<int>> A(m,vector<int>(n,0));

    //input the matrix
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            cin>>A[i][j];
        }
    }

    //the input vector
    vector<int>x;
    for(int i=0;i<n;i++){
        int l;
        cin>>l;
        x.push_back(l);
    }

    //product Ax
    vector<int>prod(m,0);

    for(int i=0;i<m;i++){
        
    }
    return 0;
}
