#include<bits/stdc++.h>
#include<chrono>

using namespace std;

using namespace std::chrono;

int E= 0;
int L= 0; /* Linear dimension */
int N = 0;
int M = 0;
int mm = 0; //bondnum
int filenum =0;

vector<int> ptr;

vector<int>bond1;
vector<int>bond2;
vector<int>bond_id;
vector<int>sbond_id;
vector<long double>occprob;
vector<long double> percstrength;
vector<long double>entropy;
vector<long double>tempentropy;

std::random_device r;
std::seed_seq seed{r(), r(), r()};
std::mt19937 gen(seed);



 //string wrapstr = "wrap_" + to_string(mm)+"_L"+ to_string(L) +"_E"+ to_string(E)+".dat";
 string opstr;
 string entrstr;
 string timestr;



int findroot(int i)
{
    if (ptr[i]<0) return i;
    return ptr[i] = findroot(ptr[i]);
}



double entropy_term(int i)
{
    int clustsize = -ptr[i];
          long double meu = (long double)clustsize/(long double)N;
          long double term = -1*meu*logl(meu);
           return term;
}

void ptr_entropy_change(int b_num,int rootfix, int rootchange)
{
     tempentropy[b_num]= tempentropy[b_num-1] - entropy_term(rootfix) -entropy_term(rootchange);

              ptr[rootfix] += ptr[rootchange];

              ptr[rootchange]= rootfix;

              tempentropy[b_num]+= entropy_term(rootfix);
}



void bondcreate()
{


 for (int i=0; i<=M; i++)
  {
      percstrength.push_back(0);
      entropy.push_back(0);
  }
   percstrength[0]= (long double)1/(long double)N;
 entropy[0] = (long double)logl(N);





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



   for (int i=0; i<=M; i++)
  {

      tempentropy.push_back(0);
  }

   tempentropy[0] = (long double)logl(N);




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
     unsigned long long int product = (unsigned long long int) s1_size* (unsigned long long int) s2_size;
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



   // product_of_clustersize.erase(product_of_clustersize.begin(), product_of_clustersize.end());
    vector<unsigned long long int>().swap(product_of_clustersize);

    return original_selection;

}





void percolation()
{
     for (int i=0; i<N; i++){ptr.push_back(-1);}

       int big = 0;
        for(int a=0; a<sbond_id.size(); a++)
    {


        int x ,y, a1;


       // a1 = sbond_id[a];

       if(a<=M-mm)
        {
           a1 = bondselection(a);
        }
       else
        {
            a1= sbond_id[a];
        }


        x= bond1[a1];

        y= bond2[a1];

        int x1= findroot(x);
       int y1 = findroot(y);

        if( ptr[x1]==-1 && ptr[y1]==-1 )
        {
                ptr_entropy_change(a+1,x1,y1);



        }

       else
     {
       if(x1!= y1)
       {
           if(-ptr[x1]>=-ptr[y1])
           {

                ptr_entropy_change(a+1,x1,y1);


           }
           else //(-ptr[x1]<-ptr[y1])
           {

            ptr_entropy_change(a+1,y1,x1);

           }
       }
           else
           {
              tempentropy[a+1]= tempentropy[a];

           }

       }

        if(tempentropy[a+1]<0){tempentropy[a+1]=0;}

        int ll = findroot(x);
         if(- ptr[ll]>big){big = - ptr[ll];}
    percstrength[a+1]+= (long double)big/(long double)N;
    entropy[a+1]+= tempentropy[a+1];



    }
}



void span()
{

    for(int u =0; u<=M; u+=1)
   {
     long double ss = (long double)u/(long double)M;
      occprob.push_back(ss);
   }

        bondcreate();

     auto start = high_resolution_clock::now();
    for(int e =0; e<E; e++)
    {

        shufflebonds();


        percolation();



        //cout<<"ensemble "<<e<<" done"<<endl;

        //sbond_id.erase(sbond_id.begin(), sbond_id.end());
        vector<int>().swap(sbond_id);
        //ptr.erase(ptr.begin(), ptr.end());
        vector<int>().swap(ptr);
        //tempentropy.erase(tempentropy.begin(), tempentropy.end());
        vector<long double>().swap(tempentropy);




    }

    auto stop = high_resolution_clock::now();

 auto duration = duration_cast<seconds>(stop-start);

 ofstream fout (timestr);
	fout<<"Time taken for Lattice Size "<<L<<":"<<duration.count()/(double)E<<"seconds"<<endl;
	fout.close();
    for(int i=1; i<=M; i++)
  {
      entropy[i]= entropy[i]/(long double)E;
      percstrength[i]= percstrength[i]/(long double)E;
  }
    ofstream file1 (opstr);
    ofstream file2(entrstr);

    for(int y=0; y<occprob.size(); y++)
{

    file1<<fixed<<setprecision(20)<<occprob[y]<<" "<<setprecision(20)<<percstrength[y]<<endl;
    file2<<fixed<<setprecision(20)<<occprob[y]<<" "<<setprecision(20)<<entropy[y]<<endl;

}





 //bond1.erase(bond1.begin(), bond1.end());
 vector<int>().swap(bond1);

 //bond2.erase(bond2.begin(), bond2.end());
 vector<int>().swap(bond2);

 //occprob.erase(occprob.begin(), occprob.end());
 vector<long double>().swap(occprob);

 // bond_id.erase(bond_id.begin(), bond_id.end());
  vector<int>().swap(bond_id);

 //entropy.erase(entropy.begin(), entropy.end());
 vector<long double>().swap(entropy);
 //percstrength.erase(percstrength.begin(), percstrength.end());
 vector<long double>().swap(percstrength);

}

int main(int argc,char* argv[])
{


    mm= atoi(argv[1]);
    L=atoi(argv[2]);
    E= atoi(argv[3]);
    filenum = atoi(argv[4]);

    N= (L*L);
    M= (2*L*L);

    opstr = "orderpm_"+ to_string(mm)+"_L"+ to_string(L) +"_E"+ to_string(E)+ "_n"+ to_string(filenum) + ".dat";
    entrstr = "entropy_" + to_string(mm)+"_L"+ to_string(L) +"_E"+ to_string(E)+ "_n"+ to_string(filenum) + ".dat";
    timestr = "opentime_" + to_string(mm)+"_L"+ to_string(L) +"_E"+ to_string(E)+ "_n"+ to_string(filenum) + ".dat";

   //auto start = high_resolution_clock::now();
    span();
    //auto stop = high_resolution_clock::now();

	//auto duration = duration_cast<seconds>(stop-start);

	//fout<<"Time taken for Lattice Size "<<L<<":"<<duration.count()/(double)E<<"seconds"<<endl;
    return 0;
}



