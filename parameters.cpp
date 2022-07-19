#include "tools.hpp"
#include "integrator.hpp"
#include "parameters.hpp"
#include "springs.hpp"

Parameters::Parameters()
{
    integration_scheme = 0;
    eta = 0.005;
    dt_min = 1.0e-10;
    random_seed = 0;
    f_rand_vel = 0.1;
    
    CONST_G = 1.0;
    
    SMBH_mass = 0*50.0;
    SMBH_pos_x = -10.0;
    SMBH_pos_y = -10.0;
    SMBH_pos_z = 0.0;
    SMBH_pos = Vector(SMBH_pos_x,SMBH_pos_y,SMBH_pos_z);
    
    Nx = Ny = 4;
    
    end_time = 10.0;
    N_snap = 100;
    
    output_filename = "node_data.txt";
    
    verbose_flag = 1;
};

void Parameters::sanity_checks()
{
    if (!(integration_scheme == 0 || integration_scheme == 1))
    {
        throw std::invalid_argument(std::string("Invalid integration scheme ") + std::to_string(integration_scheme) + "; should be either 0 (Euler) or 1 (RK4)");
    }
    if (eta <= 0)
    {
        throw std::invalid_argument(std::string("Eta = ") + std::to_string(eta) + " should be nonzero and positive");
    }
    if (dt_min <= 0)
    {
        throw std::invalid_argument(std::string("dt_min = ") + std::to_string(dt_min) + " should be nonzero and positive");
    }
    if (f_rand_vel < 0 || f_rand_vel > 1)
    {
        throw std::invalid_argument(std::string("f_rand_vel = ") + std::to_string(f_rand_vel) + " should lie between 0 and 1");
    }
    if (SMBH_mass < 0)
    {
        throw std::invalid_argument(std::string("SMBH_mass = ") + std::to_string(SMBH_mass) + " should be zero or positive");
    }
    if (Nx <= 0)
    {
        throw std::invalid_argument(std::string("Nx = ") + std::to_string(Nx) + " should be a nonzero and positive integer");
    }
    if (Ny <= 0)
    {
        throw std::invalid_argument(std::string("Ny = ") + std::to_string(Ny) + " should be a nonzero and positive integer");
    }
    if (N_snap <= 0)
    {
        throw std::invalid_argument(std::string("N_snap = ") + std::to_string(N_snap) + " should be a nonzero and positive integer");
    }
    if (end_time <= 0)
    {
        throw std::invalid_argument(std::string("end_time = ") + std::to_string(end_time) + " should be nonzero and positive");
    }
    if (verbose_flag < 0)
    {
        throw std::invalid_argument(std::string("verbose_flag = ") + std::to_string(verbose_flag) + " should be a positive integer");
    }
}

