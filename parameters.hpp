class Vector;
class Node;

#ifndef PARAMETERS_H
#define PARAMETERS_H

struct Parameters
{
        int integration_scheme;
        double eta;
        double dt_min;
        int random_seed;
        double f_rand_vel;
        
        double CONST_G;
        
        double SMBH_mass;
        double SMBH_pos_x;
        double SMBH_pos_y;
        double SMBH_pos_z;
        Vector SMBH_pos;
        
        double end_time;
        int N_snap;
        
        int Nx,Ny;
        
        int verbose_flag;
        
        std::string output_filename;
        
        Parameters();
        void sanity_checks();
};
#endif
