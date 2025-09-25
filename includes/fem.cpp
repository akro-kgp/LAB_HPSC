#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib> // For exit()

using namespace std;

// --- Function Prototypes ---
void quad1 (double *p, double *w,int ngl);
void shapefn (double rvalue, double svalue, double *shape, double * dhdr, double *dhds);
void jacobi(double * jacobian, int nnel,int n, double *dhdr, double*dhds,double*xcoord, double *ycoord );
double det(double *a, int n);
void invert(double*a, int n,double det, double *ai);
void derivshape(double *dhdx,double *dhdy,int n,int nnel,double *dhdr,double *dhds,double *jinv);
void bmat(double *be,int nnel,int n, double *dhdx,double *dhdy);
void btdb(double*kb, double *be, double wtx, double wty, double detj,int n, int nnel);
void quad2(int nglx,int ngly,int mnode,int ngl, double * gqp, double *gqw );
void stiffmatgen(int n,int N, int iel, int nel, int *nd, double* k, double * K );
void applycond (double *k, double *f,int n,int m, int * bcdof, double * bcval );


// --- Function Implementations ---

void quad1 (double *p, double *w,int ngl)
{
  if(ngl==1)
  {
    p[ngl-1]=0;w[ngl-1]=2;
  }
if(ngl==2)
{
  p[0]=-0.577350269189625765;
  p[1]=-p[0];
  w[0]=1.0;
  w[1]=w[0];
}
if(ngl==3)
{
  p[0]=0.774596669241483;
  p[1]=0;
  w[0]=0.5555555555555556;
  p[2]=-p[0];
  w[1]=0.8888888888888889;
  w[2]=w[1];
}
if (ngl==4)
{
  p[0]=-0.861136311594053;
    p[1]=-0.339981043584856;
    p[2]=-p[1];
    p[3]=-p[0];
    w[0]=0.347854845137454;
    w[1]=0.652145154862546;
    w[2]=w[1];
    w[3]=w[0];

}
}


void shapefn (double rvalue, double svalue, double *shape, double * dhdr, double *dhds)
{
  shape[0]=0.25*(1-rvalue)*(1-svalue);
  shape[1]=0.25*(1+rvalue)*(1-svalue);
  shape[2]=0.25*(1+rvalue)*(1+svalue);
  shape[3]=0.25*(1-rvalue)*(1+svalue);

  dhdr[0]=  -0.25*(1-svalue);
  dhdr[1]=   0.25*(1-svalue);
  dhdr[2]=   0.25*(1+svalue);
  dhdr[3]=  -0.25*(1+svalue);

  dhds[0]=-0.25*(1-rvalue);
  dhds[1]=-0.25*(1+rvalue);
  dhds[2]=0.25*(1+rvalue);
  dhds[3]=0.25*(1-rvalue);
}

void jacobi(double * jacobian, int nnel,int n, double *dhdr, double*dhds,double*xcoord, double *ycoord )
{
  int i,j;
  for(i=0;i<n;i++)
  for(j=0;j<n;j++)
  jacobian[n*i+j]=0;

for(i=0;i<nnel;i++)
{
  jacobian[0]+=(dhdr[i]*xcoord[i]);
  jacobian[1]+=(dhdr[i]*ycoord[i]);
  jacobian[2]+=(dhds[i]*xcoord[i]);
  jacobian[3]+=(dhds[i]*ycoord[i]);
}
}

double det(double *a, int n)
{
  double detj;
  detj=(a[0]*a[3])-(a[1]*a[2]);
  return detj;
}


void invert(double*a, int n,double det, double *ai)
{
  double adj[2][2];
  adj[0][0]=a[3];
  adj[0][1]=-a[1];
  adj[1][0]=-a[2];
  adj[1][1]=a[0];

  for(int i=0;i<n;i++)
  {
      for(int j=0;j<n;j++)
      {
         ai[n*i+j]=(adj[i][j]/det);
      }
  }
}


void derivshape(double *dhdx,double *dhdy,int n,int nnel,double *dhdr,double *dhds,double *jinv)
{
  int i;
  for(i=0;i<nnel;i++)
  {
    dhdx[i]=jinv[n*0+0]*dhdr[i]+jinv[n*0+1]*dhds[i];
    dhdy[i]=jinv[n*1+0]*dhdr[i]+jinv[n*1+1]*dhds[i];
  }
}

void bmat(double *be,int nnel,int n, double *dhdx,double *dhdy)
{
  // here not parametrized!!
  int i,j;
  for(i=0;i<n;i++)
  {
    for(j=0;j<nnel;j++)
    {
      if(i==0)
      {
        be[nnel*i+j]=dhdx[j];
      }
      if(i==1)
      {
        be[nnel*i+j]=dhdy[j];
      }
    }
  }
}

void btdb(double*kb, double *be, double wtx, double wty, double detj,int n, int nnel)
{
  int i,j,k;
  double res[nnel][nnel];
  double bet[nnel][n];

  for(i=0;i<nnel;i++)
  {
    for(j=0;j<n;j++)
    {
      bet[i][j]=be[nnel*j+i];
    }
  }

  for(i=0;i<nnel;i++)
    for(j=0;j<nnel;j++)
    res[i][j]=0;

  for(i=0;i<nnel;i++)
  {
    for(j=0;j<nnel;j++)
    {
      for(k=0;k<n;k++)
      {
        res[i][j]+=(bet[i][k]*be[nnel*k+j]);
      }
    }
  }

  for(i=0;i<nnel;i++)
  {
    for(j=0;j<nnel;j++)
    {
      kb[nnel*i+j]+=((res[i][j])*wtx*wty*detj);
    }
  }
}

void quad2(int nglx,int ngly,int mnode,int ngl, double * gqp, double *gqw )
{
  double px[nglx];double py[ngly];
  double wx[nglx];double wy[ngly];int i;

  for(i=0;i<nglx;i++)
  {
      wx[i]=0;px[i]=0;
  }

  for(i=0;i<ngly;i++)
  {
      wy[i]=0;py[i]=0;
  }

  quad1(px,wx,nglx);
  quad1(py,wy,ngly);

  for(i=0;i<nglx;i++)
  {
    gqp[mnode*i+0]=px[i];
    gqw[mnode*i+0]=wx[i];
  }

  for(i=0;i<ngly;i++)
  {
    gqp[mnode*i+1]=py[i];
    gqw[mnode*i+1]=wy[i];
  }
}

void stiffmatgen(int n,int N, int iel, int nel, int *nd, double* k, double * K )
{
  int i,j;
  for(i=0;i<n;i++)
  {
    for(j=0;j<n;j++)
    {
      K[N*(nd[i]-1)+(nd[j]-1)]+=k[n*i+j];
    }
  }
}

void applycond (double *k, double *f,int n,int m, int * bcdof, double * bcval )
{
  int i,j,a;
  for(i=0;i<n;i++)
  {
    a=bcdof[i]-1;
    for(j=0;j<m;j++)
    {
      k[a*m+j]=0;
    }
    k[a*m+a]=1;
    f[a]=bcval[i];
  }
}

int main()
{
    int i,j,k,nnode,mnode,nel,mel;
    double *coord;
    int *elem;
    k=0;i=0;j=0;

    ifstream f1("nodeinfo.txt");
    if(!f1.is_open())
    {
        cout << "failed to open nodeinfo matrix" << endl;
        return 1;
    }
    cout << "succesfully opened nodeinfo matrix" << endl;
    f1 >> nnode >> mnode;
    f1.close();

    ifstream f2("eleminfo.txt");
    if(!f2.is_open())
    {
        cout << "failed to open eleminfo matrix" << endl;
        return 1;
    }
    cout << "succesfully opened eleminfo matrix" << endl;
    f2 >> nel >> mel;
    f2.close();

    coord = new double[nnode * mnode];
    elem = new int[nel * mel];

    ifstream f3("unstrucquad_nodes.txt");
    if(!f3.is_open())
    {
        cout << "failed to open unstrucquad_nodes" << endl;
        return 1;
    }
    cout << "succesfully opened unstrucquad_nodes" << endl;
    for(i=0;i<nnode;i++)
    {
        for(j=0;j<mnode;j++)
        {
            f3 >> coord[mnode*i+j];
        }
    }
    f3.close();

    cout << "--- Coordinates ---" << endl;
    for(i=0;i<nnode;i++)
    {
        for(j=0;j<mnode;j++)
        {
            cout << coord[mnode*i+j] << " ";
        }
        cout << endl;
    }

    ifstream f4("unstrucquad_elements.txt");
    if(!f4.is_open())
    {
        cout << "failed to open unstrucquad_elements" << endl;
        return 1;
    }
    cout << "succesfully opened unstrucquad_elements" << endl;
    for(i=0;i<nel;i++)
    {
        for(j=0;j<mel;j++)
        {
            f4 >> elem[mel*i+j];
        }
    }
    f4.close();

    cout << "--- Elements ---" << endl;
    for(i=0;i<nel;i++)
    {
        for(j=0;j<mel;j++)
        {
            cout << elem[mel*i+j] << " ";
        }
        cout << endl;
    }

    int nnel=4;
    int ngqx,ngqy; int ngl;
    ngqx=2;ngqy=2;
    if(ngqx>ngqy)
        ngl=ngqx;
    else
        ngl=ngqy;
    int dofpernode=1;int nbc;
    int tst=0;

    for(i=0;i<nnode;i++)
    {
      if(coord[mnode*i+0]==0&&coord[mnode*i+1]!=0&&coord[mnode*i+1]!=1)
      {
        tst++;
      }
      if(coord[mnode*i+0]==1&&coord[mnode*i+1]!=0&&coord[mnode*i+1]!=1)
      {
        tst++;
      }
      if(coord[mnode*i+1]==0)
      {
        tst++;
      }
      if(coord[mnode*i+1]==1)
      {
        tst++;
      }
    }
    nbc=tst;

    double*gqw; double*gqp;double*ycoord;double*xcoord;
    int totdof;
    int *nd; double *be;double *shape; double *dhdr; double * dhds;
    double *jacobian; double * jinv; double *temp;double *dhdx;double *dhdy;
    totdof=nnel*1; double *p; double *K;int *bcnode; double *bcval;
    double *F;int iel;double qx,qy,wtx,wty,detj;
    int sysdof=nnode*dofpernode;

    gqw = new double[ngl * mnode];
    gqp = new double[ngl * mnode];
    nd = new int[nnel];
    xcoord = new double[nnel];
    ycoord = new double[nnel];
    shape = new double[nnel];
    dhdr = new double[nnel];
    dhds = new double[nnel];
    jacobian = new double[mnode*mnode];
    jinv = new double[mnode*mnode];
    temp = new double[mnode*mnode];
    dhdx = new double[nnel];
    dhdy = new double[nnel];
    be = new double[4*2]; // Not parameterized
    p = new double[nnel*nnel];
    K = new double[sysdof*sysdof];
    bcnode = new int[nbc];
    bcval = new double[nbc];
    F = new double[sysdof];

    //boundary cond<--------------------------------
    const double PI = 3.141592653589793;
    k=0;
    for(i=0;i<nnode;i++)
    {
      if(coord[mnode*i+0]==0&&coord[mnode*i+1]!=0&&coord[mnode*i+1]!=1)
      {
        bcnode[k]=i+1;
        bcval[k]=0;
        k++;
      }
      if(coord[mnode*i+0]==1 && coord[mnode*i+1]!=0 && coord[mnode*i+1]!=1)
    {
        bcnode[k] = i+1;
        bcval[k] = sin(PI * coord[mnode*i+1]); // Use the node's y-coordinate
        k++;
    }
      if(coord[mnode*i+1]==0)
      {
        bcnode[k]=i+1;
        bcval[k]=0;
        k++;
      }
      if(coord[mnode*i+1]==1)
      {
        bcnode[k]=i+1;
        bcval[k]=0;
        k++;
      }
    }

    k=0;i=0;iel=0;
    for(i=0;i<sysdof;i++)
    {
        F[i]=0;
        for(j=0;j<sysdof;j++)
        {
            K[sysdof*i+j]=0;
        }
    }
    for(i=0;i<nnel;i++)
        for(j=0;j<mnode;j++)
            be[mnode*i+j]=0;

    quad2(ngqx,ngqy,mnode,ngl,gqp,gqw);

    for(iel=0;iel<nel;iel++)
    {
        for(i=0;i<nnel;i++)
        {
          nd[i]=elem[mel*iel+i];
          xcoord[i]=coord[(mnode*(nd[i]-1))+0];
          ycoord[i]=coord[(mnode*(nd[i]-1))+1];
        }

        for(i=0;i<nnel;i++)
            for(j=0;j<nnel;j++)
                p[nnel*i+j]=0;

        for(i=0;i<ngqx;i++)
        {
          qx=gqp[mnode*i+0];
          wtx=gqw[mnode*i+0];
          for(j=0;j<ngqy;j++)
          {
              qy=gqp[mnode*j+1];
              wty=gqw[mnode*j+1];
              shapefn(qx,qy,shape,dhdr,dhds);
              jacobi(jacobian,nnel,mnode,dhdr,dhds,xcoord,ycoord);
              detj=det(jacobian,mnode);
              invert(jacobian,mnode,detj,jinv);
              derivshape(dhdx,dhdy,mnode,nnel,dhdr,dhds,jinv);
              bmat(be,nnel,mnode,dhdx,dhdy);
              btdb(p,be,wtx,wty,detj,mnode,nnel);
          }
        }
        stiffmatgen(nnel,sysdof,iel,nel,nd,p,K);
    }
    applycond(K,F,nbc,sysdof,bcnode,bcval);

    cout << "\n n = " << sysdof << endl;

    ofstream o1("kinfo.txt");
    o1 << sysdof << endl;
    o1.close();

    ofstream o2("Kmat.txt");
    for(k=0;k<sysdof;k++)
    {
        for(int m=0;m<sysdof;m++)
        {
            o2 << K[sysdof*k+m] << endl;
        }
    }
    o2.close();

    ofstream o3("Fvec.txt");
    for(k=0;k<sysdof;k++)
    {
        o3 << F[k] << endl;
    }
    o3.close();

    // --- Cleanup ---
    delete[] coord;
    delete[] xcoord;
    delete[] ycoord;
    delete[] nd;
    delete[] elem;
    delete[] gqp;
    delete[] gqw;
    delete[] dhds;
    delete[] dhdr;
    delete[] shape;
    delete[] jacobian;
    delete[] jinv;
    delete[] temp;
    delete[] dhdx;
    delete[] dhdy;
    delete[] K;
    delete[] F;
    delete[] bcnode;
    delete[] bcval;
    delete[] p;
    delete[] be;

    return 0;
}

