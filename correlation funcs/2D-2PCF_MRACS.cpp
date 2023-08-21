// Second Order STAtistics (SOSTA)of MRACS code, which include:
// (1) two-point correlation (with or without volume-average) function
// (2) density variance with top-hat filter
#include"mracs.h"

#define R0 0.5           // Mpc/h
#define R1 50.           // Mpc/h
#define NUMTEST 20

int main()
{
    
    read_parameter();
    std::string ifname{"/data1/quijote/lzy_pos/dm/halo/real_space/qjt_halo_0.bin"};
    std::string outname{"/data1/quijote/lzy_pos/result/qjt_halo_0.txt"};

    std::vector<double> r_log;
    std::vector<double> theta_log;
    std::vector<double> rofh, xi_r, var_r;

    for(int i = 0; i < NUMTEST; ++i) 
    {
        r_log.push_back(25+5*i);
    }

    for(int i = 0; i < NUMTEST; ++i) 
    {
        theta_log.push_back(0.1+0.15*i);
    }


    for(int ih=0; ih<1000; ++ih)
    {
        ifname=("/data1/quijote/lzy_pos/dm/halo/red_space/qjthalorsd"+std::to_string(ih)+".bin");
        std::cout << ifname << ", "; std::cout << std::endl;
        std::vector<Particle> p   =  read_in_Halo_4vector(ifname);
        auto s = sfc(p);
        auto sc  = sfc_r2c(s);
        

        force_kernel_type(3);   //ring kernal
        int count1 = 0;
        for(int i=0; i<NUMTEST; ++i)
        {
            for(int j=0; j<NUMTEST; ++j)
            {
                double ir=r_log[i];
                double iih=theta_log[j];
                
                auto w = wfc(ir,iih,0);
                auto c = convol_c2r(sc, w);
                xi_r.push_back(inner_product(s, c, GridNum) * GridNum/pow(p.size(), 2) - 1);
                std::cout << ir<<iih; std::cout << std::endl;
                std::cout << inner_product(s, c, GridNum) * GridNum/pow(p.size(), 2) - 1 << ", "; std::cout << std::endl;
                std::cout << count1<< ", "; std::cout << std::endl;

                count1=count1+1;

                delete[] c;
                delete[] w;
            }



        }


        float pp;

        outname=("/data1/quijote/lzy_pos/result_rsd2pcf/qjthalorsd2pcf"+std::to_string(ih)+".txt");
        std::ofstream outFile3(outname, std::ios_base::binary);
        for(int i = 0; i < NUMTEST*NUMTEST; ++i)
        {
            pp=(xi_r[i]);
            outFile3.write(as_bytes(pp), sizeof(float));
        }
        outFile3.close();
  

        p.clear();
        xi_r.clear();




    }
}