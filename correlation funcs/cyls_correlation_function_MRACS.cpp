// Second Order STAtistics (SOSTA)of MRACS code, which include:
// (1) two-point correlation (with or without volume-average) function
// (2) density variance with top-hat filter
#include"mracs.h"

#define R0 0.5           // Mpc/h
#define R1 50.           // Mpc/h
#define NUMTEST 10
#define NUMTEST1 20



int main()
{
    
    read_parameter();
    std::string ifname{"/data1/quijote/lzy_pos/dm/halo/real_space/qjt_halo_0.bin"};
    std::string outname{"/data1/quijote/lzy_pos/result/qjt_halo_0.txt"};


    std::vector<double> theta_log,r_log;
    std::vector<double> rofh, xi_r, var_r,var_r1,der_r;
        
    for(int i = 0; i < NUMTEST; ++i) 
    {
        r_log.push_back(50+12*i);
    }

    for(int i = 0; i < NUMTEST1; ++i) 
    {
        theta_log.push_back((-0.98+0.103*i));
    }

    



    for(int ih=799; ih<1000; ++ih)
    {
        ifname=("/data1/quijote/lzy_pos/dm/halo/red_space/qjthalorsd"+std::to_string(ih)+".bin");
        std::cout << ifname << ", "; std::cout << std::endl;
        std::vector<Particle> p   =  read_in_Halo_4vector(ifname);






        
        
        double tot {0};
        #ifdef IN_PARALLEL
        #pragma omp parallel for reduction (+:tot)
        #endif
        for(auto x : p) tot += x.weight;
        tot /= p.size();
        std::cout << "total = " << tot << std::endl;



        auto s = sfc(p);


        double psize=p.size();

        double gridnum=GridNum;

        
        for(int i=0; i<GridNum; ++i)
        {
            s[i]=s[i]*gridnum/(psize)-1;
            

        }




        double itheta;


        auto sc  = sfc_r2c(s);

        double twoofthree=2/3;



        force_kernel_type(5);   //cylinder kernal
        
        int count1 = 0;
        for(int i=0; i<NUMTEST; ++i)
        {
            for(int j=0; j<NUMTEST1; ++j)
            {

             
                
            
                itheta=acos(theta_log[j]);

                double ir=r_log[i]*sin(itheta);
                double iih=r_log[i]*cos(itheta);

                if(iih<0)
                {
                    iih=-iih;
                }



                






                auto w0 = wfc(ir,0,iih);

                auto c0 = convol_c2r(sc, w0);














                double point1_value = inner_product(s,c0,GridNum)/(GridNum);
     
            
            

                

                std::cout << ir<< ", "<<iih; std::cout << std::endl;
            
                
                der_r.push_back(point1_value);
                std::cout << (point1_value)<< ", "; std::cout << std::endl;
                count1=count1+1;

                delete[] c0;
                delete[] w0;

            }


            
            
                
            
        }
        
        
        double pp;

        outname=("/data1/quijote/lzy_pos/res_cyls/4/qjthalocyls"+std::to_string(ih)+".txt");
        std::ofstream outFile3(outname, std::ios_base::binary);
        for(int i = 0; i < NUMTEST*NUMTEST1; ++i)
        {
            pp=(der_r[i]);
            outFile3.write(as_bytes(pp), sizeof(double));
        }
        outFile3.close();
        p.clear();

        der_r.clear();



    }


        
  

        


}


