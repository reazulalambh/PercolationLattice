
#include<bits/stdc++.h>
#include<chrono>

using namespace std;
using namespace std::chrono;
 int filenum =0;
 int E= 0;
int L= 0; /* Linear dimension */
int N= 0;
int M= 0;
int mm= 2; //bondnum
vector<int>abs_x;                       //int abs_x[N];
vector<int>abs_y;                                   //int abs_y[N];
vector<int>relative_x;                                       //int relative_x[N];
vector<int> relative_y;
vector<int>clust_id;           // int clust_id[N];
vector<int>ptr;
long double pc =0;
vector<int>bond1;
vector<int>bond2;
vector<int>bond_id;
vector<int>sbond_id;
vector< vector <int> >clusters;
std::random_device r;
std::seed_seq seed{r(), r(), r()};

//auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937 gen(seed);



 string wrapstr, timestr;
// string opstr = "orderpm_"+ to_string(mm)+"_L"+ to_string(L) +"_E"+ to_string(E)+".dat";
 //string entrstr = "entropy_" + to_string(mm)+"_L"+ to_string(L) +"_E"+ to_string(E)+".dat";






int findroot(int i)
{
    if (ptr[i]<0) return i;
    return ptr[i] = findroot(ptr[i]);
}

int delta_x(int unchange, int change)
{

   int delxr = relative_x[unchange] - relative_x[change];
        int delxc = abs_x[change] - abs_x[unchange];

        if(abs(delxc)>1)
        {
            if(delxc>0)
            {
                delxc=-1;
            }
            else{delxc= 1;}
        }
        int delx = delxr+ delxc;

        return delx;
}

int delta_y(int i, int j)
{

   int delyr = relative_y[i] - relative_y[j];
    int delyc = abs_y[j] - abs_y[i];

    if(abs(delyc)>1)
        {
            if(delyc>0)
            {
                delyc=-1;
            }
            else{delyc= 1;}
        }


    int dely = delyr+ delyc;
    return dely;

}


void changeco (int fixed, int change)
{
                int tempclassids = clust_id[change];
                int delx = delta_x(fixed,change);
                int dely = delta_y(fixed,change);
                int tempclassidl = clust_id[fixed];
                for(int i =0; i<clusters[tempclassids].size(); i++)
                {
                   int j = clusters[tempclassids][i];
                   clusters[tempclassidl].push_back(j);
                   clust_id[j]= clust_id[fixed];
                   relative_x[j]= relative_x[j]+ delx;
                   relative_y[j]= relative_y[j] +dely;
                }
               // clusters[tempclassids].erase(clusters[tempclassids].begin(),clusters[tempclassids].end());
                vector<int>().swap(clusters[tempclassids]);

                //clusters[tempclassids].clear();
}


void abs_coordinates()
{

   //initialization of abs x and abs y
 //cout<<"abs positions :"<<endl;
     for(int i = 0; i<N; i++)
    {
        int ypos = i/L;
        abs_y.push_back(ypos);
        int xpos = i- L*ypos;
        abs_x.push_back(xpos);

    //cout<<abs_x[i]<<" , "<<abs_y[i]<<endl;
    }

}

void bondcreate()
{

  for(int a=0; a<M; a++)
  {
      bond_id.push_back(a);
  }


   for(int i=0; i<2*L-1; i++)
   {
       if(i%2==0)
       {
           for(int j=(i/2)*L; j<(i/2)*L + L-1;j++)
           {
               bond1.push_back(j);
               bond2.push_back(j+1);
           }
       }
       else
       {
           for(int j=((i-1)/2)*L; j<((i-1)/2)*L +L; j++)
           {
               bond1.push_back(j);
               bond2.push_back(j+L);
           }
       }
   }

    for(int k =0; k<N; k+=L)  // boundary condition
          {
               bond1.push_back(k);
               bond2.push_back(k+L-1);
          }

          for(int l=0; l<=L-1; l++) //boundary condition
           {
               bond1.push_back(l);
               bond2.push_back(l+L*(L-1));
           }

}


void shufflebonds()
{
    for(int h=0; h<bond_id.size(); h++)
    {
        sbond_id.push_back(bond_id[h]);
    }




  shuffle(sbond_id.begin(),sbond_id.end(), gen);


    for(int i =0; i<N; i++)
    {
        relative_x.push_back(0);
        relative_y.push_back(0);
        clust_id.push_back(i);

        vector<int>clust;
        clust.push_back(i);
        clusters.push_back(clust);
        //clust.erase(clust.begin(), clust.end());
        vector<int>().swap(clust);
    }

}

 vector<unsigned long long int>product_of_clustersize;
int bondselection(int A)
{


   for(int i=0; i<mm; i++)
   {
      int e1 = sbond_id[A+ i];
      int s1 = bond1[e1];
      int s2 = bond2[e1];
      int s1_root = findroot(s1);
      int s2_root = findroot(s2);
      int s1_size= -ptr[s1_root];
      int s2_size = -ptr[s2_root];
     unsigned long long int product = (unsigned long long int) s1_size * (unsigned long long int) s2_size;
      product_of_clustersize.push_back(product);
   }

    unsigned long long int selected_bond_size = product_of_clustersize[0];
   int selected_bond = A;
   for(int j=1; j<mm; j++)
   {
       if(product_of_clustersize[j]<selected_bond_size)
       {
           selected_bond_size= product_of_clustersize[j];
           selected_bond= A+j;

       }
   }


    int original_selection = sbond_id[selected_bond];

    swap(sbond_id[A],sbond_id[selected_bond]);

    uniform_int_distribution<int> dist(A+1,sbond_id.size()-1);
    for(int k=1; k<mm; k++)
        {


            int f = dist(gen);
            swap(sbond_id[A+k], sbond_id[f]);
        }



    //product_of_clustersize.erase(product_of_clustersize.begin(),product_of_clustersize.end());
    vector<unsigned long long int>().swap(product_of_clustersize);

    return original_selection;

}





void percolation()
{
     for (int i=0; i<N; i++){ptr.push_back(-1);}
    long double p = 0;
    int pp=0;
        for(int a=0; a<sbond_id.size(); a++)
    {
        pp++;
        p= (long double)pp/(long double)M;
        int x ,y, a1;


       // a1 = sbond_id[a];

       if(a<=M-mm)
       {a1 = bondselection(a);}
       else
       {
           a1 = sbond_id[a];
       }
        x= bond1[a1];

        y= bond2[a1];

        int x1= findroot(x);
       int y1 = findroot(y);

        if( ptr[x1]==-1 && ptr[y1]==-1 )
        {
                ptr[x1]+= -1 ;
                ptr[y1]= x1 ;

                changeco(x,y);

        }

       else
     {
       if(x1!= y1)
       {
           if(-ptr[x1]>=-ptr[y1])
           {

                ptr[x1] += ptr[y1];
              ptr[y1]= x1;
                changeco(x,y);

           }
           else //(-ptr[x1]<-ptr[y1])
           {

            ptr[y1] += ptr[x1];
               ptr[x1]= y1;
             changeco(y,x);

           }
       }
           else if(p>0.4)
           {
              int xcomp = relative_x[x] - relative_x[y];
              int ycomp = relative_y[x] - relative_y[y];

              if(abs(xcomp)>1|| abs(ycomp)>1)
               {
               // cout<<"percolation done"<<endl;
                   pc=p;//cout<<"pc = "<<pc<<endl;
                   //cout<<"a+1 = "<<a+1<<endl;
                    break;
               }

           }

       }
    }
}

vector<long double>spanprob;
vector<long double>occprob;

void span()
{

    for(int u =0; u<=M; u+=1)
   {
      long double ss = (long double)u/(long double)M;
      occprob.push_back(ss);
   }
   vector<int>countpc(occprob.size(),0);


        abs_coordinates();
        bondcreate();

 auto start = high_resolution_clock::now();
    for(int e =0; e<E; e++)
    {

        shufflebonds();


        percolation();


        for(int v =0; v<occprob.size(); v+=1)
        {

            if(pc<= occprob[v]) {countpc[v]+=1;}

        }
       // cout<<"L = "<<L<<" ensemble "<<e<<" done"<<endl;

        //sbond_id.erase(sbond_id.begin(),sbond_id.end());
        vector<int>().swap(sbond_id);
        // ptr.erase(ptr.begin(), ptr.end());
         vector<int>().swap(ptr);
          pc =0;
           //clusters.erase(clusters.begin(), clusters.end());
           vector< vector <int> >().swap(clusters);

         //relative_x.erase(relative_x.begin(), relative_x.end());
         vector<int>().swap(relative_x);
          //relative_y.erase(relative_y.begin(), relative_y.end());
          vector<int>().swap(relative_y);

        //clust_id.erase(clust_id.begin(), clust_id.end());
        vector<int>().swap(clust_id);



    }

   int h =0;
    for(int v =0; v<countpc.size(); v+=1)
    {
        spanprob.push_back((long double)countpc[v]/(long double)E);
        h++;
    }
 auto stop = high_resolution_clock::now();

 auto duration = duration_cast<seconds>(stop-start);

 ofstream fout (timestr);
	fout<<"Time taken for Lattice Size "<<L<<":"<<duration.count()/(double)E<<"seconds"<<endl;
	fout.close();

    ofstream fout1(wrapstr);
for(int y=0; y<occprob.size(); y++)
{
    //cout<<occprob[y]<<" "<<spanprob[y]<<endl;
    fout1<<fixed<<setprecision(20)<<occprob[y]<<" "<<setprecision(20)<<spanprob[y]<<endl;

}
fout1.close();

countpc.erase(countpc.begin(),countpc.end());
bond1.erase(bond1.begin(),bond1.end());
bond2.erase(bond2.begin(),bond2.end());
spanprob.erase(spanprob.begin(),spanprob.end());
occprob.erase(occprob.begin(),occprob.end());
bond_id.erase(bond_id.begin(), bond_id.end());
abs_x.erase(abs_x.begin(), abs_x.end());
abs_y.erase(abs_y.begin(),abs_y.end());
}

int main(int argc,char* argv[])
{

    mm= atoi(argv[1]);
    L=atoi(argv[2]);
    E= atoi(argv[3]);
    filenum = atoi (argv[4]);

    N= (L*L);
    M= (2*L*L);
    wrapstr = "wrap_" + to_string(mm)+"_L"+ to_string(L) +"_E"+ to_string(E)+ "_n" + to_string(filenum) + ".dat";
    timestr = "wraptime_" + to_string(mm) +"_L"+ to_string(L)+"_E"+ to_string(E)+".dat";






    span();



    return 0;
}


