#include<bits/stdc++.h>
#include<chrono>
using namespace std;
using namespace std::chrono;


int filenum =0;
int E= 1;
int L= 0; /* Linear dimension */
int N= 0;
int M= 0;
int mm= 1;
string opstr;
string entrstr;
string timestr;
std::random_device r;
std::seed_seq seed{r(), r()};
std::mt19937 gen(seed);


vector<int>shuffled_sites;
vector<int>ptr;
vector< vector <int> >nb; //neighbor bonds, right left up down
vector<long double>occprob;
vector<long double> percstrength;
vector<long double>entropy;
vector<long double>tempentropy;

int findroot(int i)
{
    if (ptr[i]<0) {return i;}

    else{return ptr[i] = findroot(ptr[i]);}
}

void static_initialization()
{



    for (int i=0; i<=N; i++)
  {
      percstrength.push_back(0);
      entropy.push_back(0);
  }
    percstrength[0]= (long double)1/(long double)M;
    entropy[0] = logl((long double)M);




//bond up down right left initialization

   nb.resize(N,vector<int>());

    for(int i=0; i<N; i++)
    {
        nb[i].push_back(i); //right bond


        //left bond

        if (i%L==0)
        {
            nb[i].push_back(i+L-1) ;
        }

        else
        {
            nb[i].push_back(i-1) ;
        }

        //up bond
        nb[i].push_back(i+N);

        //down bond
        if(i<L)
        {
            nb[i].push_back(M-L+i);
        }

        else
        {
            nb[i].push_back(L*(L-1) + i);
        }



    }

}


void variable_initialization()
{


    for (int i=0; i<=N; i++)
  {

      tempentropy.push_back(0);
  }

   tempentropy[0] = logl((long double)M);



     //shuffled_site_initialization

    for(int i=0; i<N; i++)
    {
        shuffled_sites.push_back(i);
    }

     //shuffle

    shuffle(shuffled_sites.begin(),shuffled_sites.end(), gen);



    //ptr initiallization

    for(int i=0; i<M; i++)
    {
        ptr.push_back(-1);
    }



}


void static_cleaner()
{
    vector< vector <int> >().swap(nb);
    vector<long double>().swap(occprob);
    vector<long double>().swap(entropy);
    vector<long double>().swap(percstrength);

}

void variable_cleaner()
{
    vector<int>().swap(shuffled_sites);
    vector<int>().swap(ptr);
    vector<long double>().swap(tempentropy);

}



int site_selection(int A)
{

    long double product_of_clustersize = 0;
    long double min_size = 4 * logl((long double)M);
    int site = 0;
    int chosen_site;
    int rb, lb, ub, db;
    long double cr, cl, cu, cd;
    int position = 0;

    for(int i=0; i<mm; i++)
    {
        site = shuffled_sites[A+i];
        rb = nb[site][0];
        lb = nb[site][1];
        ub = nb[site][2];
        db = nb[site][3];

        int root1 = findroot(rb);
        int root2 = findroot(lb);
        int root3 = findroot(ub);
        int root4 = findroot(db);


        cr = (long double) -ptr[root1];
        cl = (long double) -ptr[root2];
        cu = (long double) -ptr[root3];
        cd = (long double) -ptr[root4];
        long double lr = logl(cr);
        long double ll = logl(cl);
        long double lu = logl(cu);
        long double ld = logl(cd);

        product_of_clustersize = lr + ll + lu + ld;

        if(product_of_clustersize<min_size)
        {
            min_size = product_of_clustersize;
            chosen_site = site;
            position = A + i;

        }

    }
    swap(shuffled_sites[A],shuffled_sites[position]);

    uniform_int_distribution<int> dist(A+1,N-1);

    for(int k=1; k<mm; k++)
        {


            int f = dist(gen);
            swap(shuffled_sites[A+k], shuffled_sites[f]);
        }

        return chosen_site;

}


long double entropy_term(int i)
{
    int clustsize = -ptr[i];
          long double meu = (long double)clustsize/(long double)M;
          long double term = -1*meu*logl(meu);
           return term;
}

void percolation()
{
    int big = 0;

     for(int a=0; a<shuffled_sites.size(); a++)
     {
         int a1, r1;
         //cout<<a<<endl;



         if(mm>1 && a<=N-mm)
        {
            a1= site_selection(a);
        }
       else
        {
            a1= shuffled_sites[a];
        }

         r1 = findroot(nb[a1][0]);

         tempentropy[a+1]= tempentropy[a];

          for(int i=0; i<4; i++)
          {
               int bond = nb[a1][i];
               int r2 = findroot(bond);

               if(i!=0 && r1 != r2)
               {

                       if(-ptr[r1]>-ptr[r2])
                       {
                          tempentropy[a+1]-= (entropy_term(r1) +entropy_term(r2));

                           ptr[r1] += ptr[r2];
                           ptr[r2] = r1;

                           tempentropy[a+1]+= entropy_term(r1);
                       }

                       else
                       {
                           tempentropy[a+1]-= (entropy_term(r1) +entropy_term(r2));

                           ptr[r2]+= ptr[r1];
                           ptr[r1] = r2;
                           r1 = r2;

                           tempentropy[a+1]+= entropy_term(r1);


                       }

               }


          }

          if(tempentropy[a+1]<0){tempentropy[a+1]=0;}
          //cout<<tempentropy[a+1]<<endl;
            r1 = findroot(nb[a1][0]);
          if(-ptr[r1]>big){big = -ptr[r1];}
         // cout<<big<<endl;
          percstrength[a+1]+= (long double)big/(long double)M;
          entropy[a+1]+= tempentropy[a+1];
     }
}




void process()
{

    for(int u =0; u<=N; u+=1)
   {
     long double ss = (long double)u/(long double)N;
      occprob.push_back(ss);
   }

 static_initialization();

     auto start = high_resolution_clock::now();
    for(int e =0; e<E; e++)
    {

       variable_initialization();




        percolation();




        cout<<"ensemble "<<e<<" done"<<endl;

       variable_cleaner();



    }

    auto stop = high_resolution_clock::now();

 auto duration = duration_cast<seconds>(stop-start);

 ofstream fout (timestr);
	fout<<"Time taken for Lattice Size "<<L<<":"<<duration.count()/(double)E<<"seconds"<<endl;
	fout.close();
    for(int i=1; i<=N; i++)
  {
      entropy[i]= entropy[i]/(long double)E;
      percstrength[i]= percstrength[i]/(long double)E;
  }
    ofstream file1(opstr);

    ofstream file2(entrstr);




    for(int y=0; y<occprob.size(); y++)
{

    file1<<fixed<<setprecision(20)<<occprob[y]<<" "<<setprecision(20)<<percstrength[y]<<endl;

    file2<<fixed<<setprecision(20)<<occprob[y]<<" "<<setprecision(20)<<entropy[y]<<endl;


}


static_cleaner();

}

int main(int argc,char* argv[])
{
   //ofstream fout ("timeop.dat");

    mm= atoi(argv[1]);
    L=atoi(argv[2]);
    E= atoi(argv[3]);
    filenum = atoi(argv[4]);

    N= (L*L);
    M= (2*L*L);

    opstr = "orderpmsite_"+ to_string(mm)+"_L"+ to_string(L) +"_E"+ to_string(E)+ "_n"+ to_string(filenum) + ".dat";
    entrstr = "entropysite_" + to_string(mm)+"_L"+ to_string(L) +"_E"+ to_string(E)+ "_n"+ to_string(filenum) + ".dat";
    timestr = "opentimesite_" + to_string(mm)+"_L"+ to_string(L) +"_E"+ to_string(E)+ "_n"+ to_string(filenum) + ".dat";

   //auto start = high_resolution_clock::now();
    process();
    //auto stop = high_resolution_clock::now();

	//auto duration = duration_cast<seconds>(stop-start);

	//fout<<"Time taken for Lattice Size "<<L<<":"<<duration.count()/(double)E<<"seconds"<<endl;
    return 0;
}

