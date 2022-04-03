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
int EMPTY=0;
long double mass=0;
double pc =0;
int ccheck =0;
string fractalstr, timestr, distrstr;
std::random_device r;
std::seed_seq seed{r(), r(), r()};
std::mt19937 gen(seed);

//***********************************************************************************************************************************************

class cluster
{
public:
    vector<int>site;
    vector<int>bond;
};
//***********************************************************************************************************************************************
vector<cluster>clusterBS;
vector<int>shuffled_sites;
vector<int>ptr;
vector<int>abs_x;
vector<int>abs_y;
vector<int>relative_x;
vector<int>relative_y;
vector< vector <int> >nn;  //serial right, left, up, down
vector< vector <int> >nb; //neighbor bonds, right left up down
vector<int>cluster_id_site;
vector<int>cluster_id_bond;
//***********************************************************************************************************************************************
int findroot(int i)
{
    if (ptr[i]<0) {return i;}

    else{return ptr[i] = findroot(ptr[i]);}
}
//***********************************************************************************************************************************************
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
//************************************************************************************************************************************************
int delta_y(int fixed, int change)
{

   int delyr = relative_y[fixed] - relative_y[change];
    int delyc = abs_y[change] - abs_y[fixed];

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
//***********************************************************************************************************************************************
void static_initialization()
{


     //absolute coordinate initialization
      for(int i = 0; i<N; i++)
    {
        int ypos = i/L;
        abs_y.push_back(ypos);
        int xpos = i- L*ypos;
        abs_x.push_back(xpos);
    }

    //nn initialization

    for(int i=0; i<N; i++)
    {
        vector<int>neighbours;
        for(int j=0; j<4; j++)
        {
            neighbours.push_back(0);
        }
        nn.push_back(neighbours);
        vector<int>().swap(neighbours);
    }

for (int i=0; i<N; i++) {
nn[i][0] = (i+1)%N;
nn[i][1] = (i+N-1)%N;
nn[i][2] = (i+L)%N;
nn[i][3] = (i+N-L)%N;
if (i%L==0) nn[i][1] = i+L-1;
if ((i+1)%L==0) nn[i][0] = i-L+1;
}

/*
//nearest neighbers printing

for(int i=0; i<nn.size(); i++)
{
    cout<<"neighbours of site "<<i<<" :"<<endl;

    for(int j=0; j<nn[i].size(); j++)
    {
        cout<<nn[i][j]<<"\t";
    }
    cout<<endl;
}

*/





       // cout<<"bond initialization starts" <<endl;
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

// cout<<"bond initialization check"<<endl;

    }

   /* cout<<"nb check :"<<endl;
    for(int i =0; i<nb.size(); i++)
    {
        cout<<"neighbouring bonds of site "<<i<<" : "<<endl;

        for(int j=0; j<nb[i].size(); j++)
        {
            cout<<nb[i][j]<<"\t";
        }
        cout<<endl;
    }

    */


}

//***********************************************************************************************************************************************

void variable_initialization()
{


     //shuffled_site_initialization

    for(int i=0; i<N; i++)
    {
        shuffled_sites.push_back(i);
    }

     //shuffle

    shuffle(shuffled_sites.begin(),shuffled_sites.end(), gen);



    //clusterBS initialization

    cluster mycluster;
    for(int i =0; i<M; i++)
    {

        mycluster.bond.push_back(i);
        clusterBS.push_back(mycluster);
        vector<int>().swap(mycluster.bond);
    }



    //ptr initiallization

    for(int i=0; i<N; i++)
    {
        ptr.push_back(EMPTY);
    }

    //relative co-ordinates initialization

    for(int i =0; i<N; i++)
    {
        relative_x.push_back(0);
        relative_y.push_back(0);
    }


    //cluster_id_site initialization

    for(int i=0; i<N; i++)
    {
        cluster_id_site.push_back(EMPTY);
    }

    //cluster_iD_bond initialization

    for(int i=0; i<M; i++)
    {
        cluster_id_bond.push_back(i);
    }

}

//*****************************************************************************************************************************************

void static_cleaner()
{

    vector<int>().swap(abs_x);
    vector<int>().swap(abs_y);
    vector< vector <int> >().swap(nn);
    vector< vector <int> >().swap(nb);
  
    

}


//*****************************************************************************************************************************************


void variable_cleaner()
{
    vector<cluster>().swap(clusterBS);
    vector<int>().swap(shuffled_sites);
    vector<int>().swap(ptr);
    vector<int>().swap(relative_x);
    vector<int>().swap(relative_y);
    vector<int>().swap(cluster_id_bond);
    vector<int>().swap(cluster_id_site);
    pc =0;

}

//******************************************************************************************************************************************


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

        cr = (long double) (clusterBS[cluster_id_bond[rb]].bond.size());
        cl = (long double) (clusterBS[cluster_id_bond[lb]].bond.size());
        cu = (long double) (clusterBS[cluster_id_bond[ub]].bond.size());
        cd = (long double) (clusterBS[cluster_id_bond[db]].bond.size());
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

//******************************************************************************************************************************************

void percolation()
{
    double p = 0;
    for(int a=0; a<shuffled_sites.size(); a++)
    {
        p = (double)(a+1)/(double)N;
        //cout<<p<<endl;

        int a1, a2, r1, r2;

        a1= site_selection(a);
        //cout<<a1<<endl;

        r1= a1;

        ptr[r1]=-1;


         for(int i=0; i<4; i++)
         {

             a2 = nn[a1][i];

             int bond = nb[a1][i];

             if(ptr[a2]==EMPTY)
             {
                 if(i==0)
                {
                    cluster_id_site[a1] = cluster_id_bond[bond];
                    clusterBS[cluster_id_site[a1]].site.push_back(a1);
                }
                else
                {
                    cluster_id_bond[bond] = cluster_id_site[a1];
                    clusterBS[cluster_id_site[a1]].bond.push_back(bond);
                }


             }
             else //a2 not empty
             {


                 r2 = findroot(a2);
                 if(r2 != r1)
                 {
                     if (-ptr[r1]>-ptr[r2])
                     {
                         ptr[r1] += ptr[r2];
                         ptr[r2] = r1;

                         int tempclustids = cluster_id_site[a2];
                         int tempclustidl = cluster_id_site[a1];


                         int delx = delta_x(a1,a2);
                         int dely = delta_y(a1,a2);

                         for(int k=0; k<clusterBS[tempclustids].site.size(); k++)
                         {
                            int j = clusterBS[tempclustids].site[k];

                            relative_x[j]= relative_x[j]+ delx;
                            relative_y[j]= relative_y[j] +dely;

                            cluster_id_site[j] = tempclustidl;

                            clusterBS[tempclustidl].site.push_back(j);

                         }
                         vector<int>().swap(clusterBS[tempclustids].site);

                         for(int k=0; k<clusterBS[tempclustids].bond.size(); k++)
                         {
                            int j = clusterBS[tempclustids].bond[k];

                            cluster_id_bond[j] = tempclustidl;
                            clusterBS[tempclustidl].bond.push_back(j);

                         }
                         vector<int>().swap(clusterBS[tempclustids].bond);



                     }

                     else  //cluster of neighbour large
                     {
                         ptr[r2] += ptr[r1];
                         ptr[r1] = r2;
                         r1 = r2;

                         if(i==0)
                         {
                             cluster_id_site[a1] = cluster_id_site[a2];
                             clusterBS[cluster_id_site[a1]].site.push_back(a1);

                             int delx = delta_x(a2,a1);
                             int dely = delta_y(a2,a1);
                             relative_x[a1]= relative_x[a1]+ delx;
                             relative_y[a1]= relative_y[a1] +dely;

                         }
                         else
                         {
                             int tempclustids = cluster_id_site[a1];
                             int tempclustidl = cluster_id_site[a2];


                         int delx = delta_x(a2,a1);
                         int dely = delta_y(a2,a1);

                         for(int k=0; k<clusterBS[tempclustids].site.size(); k++)
                         {
                            int j = clusterBS[tempclustids].site[k];

                                        relative_x[j]= relative_x[j]+ delx;
                                        relative_y[j]= relative_y[j] +dely;

                                        cluster_id_site[j] = tempclustidl;

                                        clusterBS[tempclustidl].site.push_back(j);

                                    }
                                    vector<int>().swap(clusterBS[tempclustids].site);

                                    for(int k=0; k<clusterBS[tempclustids].bond.size(); k++)
                                    {
                                        int j = clusterBS[tempclustids].bond[k];

                                        cluster_id_bond[j] = tempclustidl;
                                        clusterBS[tempclustidl].bond.push_back(j);

                                    }
                                    vector<int>().swap(clusterBS[tempclustids].bond);


                         }
                     }


                 }
                 else
                 {
                     int xcomp = relative_x[a1] - relative_x[a2];
                                int ycomp = relative_y[a1] - relative_y[a2];

                                if(abs(xcomp)>1|| abs(ycomp)>1)
                                {
                                    // cout<<"percolation done"<<endl;
                                    pc=p;
									ccheck += 1;
                                    //cout<<"a+1 = "<<a+1<<endl;
                                   
                                }
                 }
             }



         }
		 
		 if(ccheck>=1)
		 {
			int rightbond = nb[a1][0];
			mass +=  (clusterBS[cluster_id_bond[rightbond]].bond.size());
		 }

         if(pc>0){break;}
    }

}


//**************************************************************************************************************************
void wrap()
{


	 vector<long double>sizedistr(M,0.0);


    static_initialization();


    auto start = high_resolution_clock::now();

    for(int e =0; e<E; e++)
    {

        variable_initialization();

        percolation();
		ccheck =0;
		
		for(int aa=0; aa<ptr.size(); aa++)
				   {
					   int ch = ptr[aa];
					   if (ch<0)
					   {
						  int RB = nb[aa][0];

						  int clustsize = (long double) (clusterBS[cluster_id_bond[RB]].bond.size());
						  
						  sizedistr[clustsize-1]+= 1; 
					   }
				   }

        cout<<"L = "<<L<<" ensemble "<<e<<" done"<<endl;
        variable_cleaner();
    }
    auto stop = high_resolution_clock::now();

 auto duration = duration_cast<seconds>(stop-start);

 ofstream fout (timestr);
	fout<<"Time taken for Lattice Size "<<L<<":"<<duration.count()/(double)E<<"seconds"<<endl;
	fout.close();
	ofstream file (distrstr);
    for(int v =0; v<M; v+=1)
    {
       long double ccc = sizedistr[v];
	   if(ccc!=0)
	   {file<<v+1<<"	"<<fixed<<setprecision(20)<<ccc/(long double)(E*M)<<endl;}
    }
file.close();
 
ofstream fout1(fractalstr);
 fout1<<L<<"	"<<fixed<<setprecision(20)<<mass/(long double)E<<endl;
fout1.close();


static_cleaner();
}


    int main(int argc,char* argv[])
{

    mm= atoi(argv[1]);
    L=atoi(argv[2]);
    E= atoi(argv[3]);
    filenum = atoi (argv[4]);

    N= (L*L);
    M= (2*L*L);
    EMPTY = -M-1;
    fractalstr = "fractalsite_" + to_string(mm)+"_L"+ to_string(L) +"_E"+ to_string(E)+ "_n" + to_string(filenum) + ".dat";
    timestr = "nsdftimesite_" + to_string(mm) +"_L"+ to_string(L)+"_E"+ to_string(E)+".dat";
	distrstr = "nspcsite_" + to_string(mm)+"_L"+ to_string(L) +"_E"+ to_string(E)+ "_n" + to_string(filenum) + ".dat";
    wrap();




    return 0;
}



