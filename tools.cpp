#include "tools.hpp"
#include "integrator.hpp"
#include "parameters.hpp"
#include "springs.hpp"

double Vector::norm()
{
    return sqrt(x*x + y*y + z*z);
};
double Vector::norm_squared()
{
    return x*x + y*y + z*z;
};
double Vector::dot(Vector &v)
{
    return x*v.x + y*v.y + z*v.z;
};
Vector Vector::operator = (const Vector &v)
{
    if (this == &v)
    {
        return *this;
    }
    (*this).x = v.x;
    (*this).y = v.y;
    (*this).z = v.z;
    return *this;
}
Vector Vector::operator + (const Vector &v)
{
    Vector res;
    res.x = x + v.x;
    res.y = y + v.y;
    res.z = z + v.z;
    return res;
};
Vector Vector::operator += (const Vector &v)
{
    if (this == &v)
    {
        this->x *= 2;
        this->y *= 2;
        this->z *= 2;
        return *this;
    }
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return *this;
};
Vector Vector::operator - (const Vector &v)
{
    Vector res;
    res.x = this->x - v.x;
    res.y = this->y - v.y;
    res.z = this->z - v.z;
    return res;
};
Vector Vector::operator * (const double s)
{
    Vector res;
    res.x = this->x * s;
    res.y = this->y * s;
    res.z = this->z * s;
    return res;
};
Vector Vector::operator / (const double s)
{
    if (s == 0.0)
    {
        std::cout << "Error: dividing by zero!" << std::endl;
    }
    Vector res;
    res.x = x / s;
    res.y = y / s;
    res.z = z / s;
    return res;
};
void Vector::print()
{
    std::cout << "x = " << x << "; y = " << y << "; z = " << z << std::endl;
};

Node::Node()
{
    index = max_index;
    max_index++;
    mass = 0.0;
    pos = Vector();
    vel = Vector();
    acc = Vector();
};
Node::Node(Vector pos_, Vector vel_, double mass_)
{
    index = max_index;
    max_index++;
    mass = mass_;
    pos = pos_;
    vel = vel_;
};
Node::~Node()
{
    //cout << "Destructing Node" << endl;
}
void Node::print()
{
    std::cout << std::string(100, '=') << std::endl;
    std::cout << "Node -- index " << index << " -- mass " << mass << std::endl;
    std::cout << std::string(100, '-') << std::endl;
    std::cout << "Pos: " << "x = " << pos.x << "; y = " << pos.y << "; z = " << pos.z << std::endl;
    std::cout << "Vel: " << "x = " << vel.x << "; y = " << vel.y << "; z = " << vel.z << std::endl;
    std::cout << "Connected nodes (N = " << connected_nodes.size() << ") indices: ";
    int i=0;
    for (auto np = connected_nodes.begin(); np != connected_nodes.end(); np++,i++)
    {
        if (i != connected_nodes.size()-1)
        {
            std::cout << (**np).index << "; ";
        }
        else
        {
            std::cout << (**np).index << std::endl;
        }
    }
};
void Node::connect_node(Node *n, Spring *spring)
{
    if (n==this)
    {
        std::cout << " ERROR  Node::connect_node index" << index << std::endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        connected_nodes.push_back(n);
        connecting_springs.push_back(spring);
    }
};
int Node::get_number_of_neighbours()
{
    return connected_nodes.size();
};

Spring::Spring()
{
    k = 0.0;
    b = 0.0;
    u0 = 1.0;
}
Spring::Spring(double k_, double b_, double u0_)
{
    k = k_;
    b = b_;
    u0 = u0_;
}
Spring::~Spring()
{
    //std::cout << "Destructing Spring" << std::endl;
}

Node_Collection::Node_Collection(Parameters *parameters)
{
    this->set_parameters(parameters);
    this->set_integrator(this->parameters);
    
    t = 0.0;
    output_file_created = false;
    initial_energy_calculated = false;
}
Node_Collection::~Node_Collection()
{
    //std::cout << "Destructing Node_Collection" << std::endl;
    
    for (int i=0; i<this->nodes.size(); i++)
    {
        for (int j=0; j<this->nodes[i]->connecting_springs.size(); j++)
        {
            delete this->nodes[i]->connecting_springs[j];
        }
        delete this->nodes[i];
    }
    delete this->integrator;
    
};

void Node_Collection::set_parameters(Parameters *parameters)
{
    try
    {
        parameters->sanity_checks();        
    }
    catch(std::invalid_argument &e)
    {
        std::cerr << e.what() << std::endl;
        exit(-1);
    }
    this->parameters = parameters;
}

void Node_Collection::add_node(Node* n)
{
    this->nodes.push_back(n);
};

void Node_Collection::set_integrator(Parameters *parameters)
{
    this->integrator = new Integrator_Euler();
    if (this->parameters->integration_scheme == 1)
    {
        this->integrator = new Integrator_RK4();
    }
}

void Node_Collection::print()
{
    std::cout << std::string(100, '=') << std::endl;
    std::cout << std::string(100, '=') << std::endl;
    std::cout << "Printing Node Collection" << std::endl;
    for (int i=0; i<nodes.size(); i++)
    {
        nodes[i]->print();
    }
    std::cout << std::string(100, '=') << std::endl;
};

void Node_Collection::file_dump()
{
    const int N_nodes = this->nodes.size();

    if (this->output_file_created == false)
    {
        std::ofstream file(this->parameters->output_filename.c_str(), std::ios::out); // create new output file
        
        if (!file.fail())
        {
            std::ostringstream oss;
            
            oss << "N_nodes" << "\t" << std::to_string(N_nodes) << "\t" << "m" << "\t" << this->nodes[0]->mass << "\t" << "k" << "\t" << this->nodes[0]->connecting_springs[0]->k << "\t" << "b" << "\t" << this->nodes[0]->connecting_springs[0]->b << std::endl;
            
            file << oss.str();
        }
        this->output_file_created = true;
        file.close();
    }
    else
    {
        std::ofstream file(this->parameters->output_filename.c_str(), std::ios::app); // append to previously-created output file

        if (!file.fail())
        {
            std::ostringstream oss;
            oss.precision(OUTPUT_PRECISION);
            oss.setf(std::ios::scientific);
            
            oss << t << "\t";
            for (int i=0; i<N_nodes; i++)
            {
                oss << this->nodes[i]->pos.x << "\t" << this->nodes[i]->pos.y << "\t" << this->nodes[i]->pos.z << "\t";
            }
            
            oss << std::endl;
            file << oss.str();
        }
        file.close();
    }
        
}

void Node_Collection::set_up_2d_box_of_nodes(const double mass)
{
    
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    std::mt19937 rng(this->parameters->random_seed); //Standard mersenne_twister_engine (default seed 0)

    const int Nx = this->parameters->Nx;
    const int Ny = this->parameters->Ny;
    
    const double fx = 1.0/((double) Nx);
    const double fy = 1.0/((double) Ny);
    const double f_rand_vel = this->parameters->f_rand_vel;
    
    Vector vel_av(0,0,0);
    for (int ix=0; ix<Nx; ix++)
    {
        for (int iy=0; iy<Ny; iy++)
        {
            Vector pos((1+ix) * fx,(1+iy) * fy,0.0);
            Vector vel(f_rand_vel*dis(rng),f_rand_vel*dis(rng),0);
            vel_av += vel;

            Node *n = new Node(pos,vel,mass);
            this->add_node(n);        
        }
    }
    
    vel_av / ((double) (Nx*Ny));
    std::cout << "Average velocity of initiated nodes: " << std::endl;
    vel_av.print();
}

void Node_Collection::connect_all_nodes_in_node_collection(const double k, const double b)
{
    /* Double-loop over all nodes in the system and connect all pairs (except self-pairs),
     * taking k and b to be the same for all nodes. */
    for (auto &np1 : this->nodes)
    {
        for (auto &np2 : this->nodes)
        {
            if (np1 != np2) // avoid self-pairs
            {
                double u0 = ((*np2).pos - (*np1).pos).norm(); // rest distance for which there is no spring force
                Spring *spring = new Spring(k, b, u0);

                (*np1).connect_node(np2, spring);
            }
        }
    }
    this->initial_energy = this->calculate_energy();
}

double Node_Collection::get_relative_energy_error()
{
    if (this->initial_energy_calculated == false)
    {
        this->initial_energy = this->calculate_energy();
        this->initial_energy_calculated = true;
    }
    this->current_energy = this->calculate_energy();
    return fabs( (this->initial_energy - this->current_energy) / this->initial_energy);
}

double Node_Collection::calculate_energy()
{
    /* Compute total energy of system.
     * Currently does not take into account dissipative forces (spring dampening). */
     
    double energy = 0.0;
    int i1 = 0, i2 = 0;
    for (auto np1 = this->nodes.begin(); np1 != this->nodes.end(); np1++, i1++)
    {
        const double m1 = (**np1).mass;
        
        energy += 0.5 * m1 * (**np1).vel.norm_squared(); // kinetic energy
        
        /* Take into account gravity from a supermassive black hole (fixed at user-specified position, i.e., without back-reaction onto the SMBH. */
        Vector r_vec = (**np1).pos - parameters->SMBH_pos;
        const double r = r_vec.norm();
        energy += -parameters->CONST_G * parameters->SMBH_mass / r;

        /* Spring-related forces */
        i2 = 0;
        for (auto np2 = (**np1).connected_nodes.begin(); np2 != (**np1).connected_nodes.end(); np2++, i2++)
        {
            double m2 = (**np2).mass;
            
            Vector d_vec = (**np2).pos - (**np1).pos; // `distance vector': points from node 1 to node 2
            Vector v_vec = (**np2).vel - (**np1).vel; // relative velocity vector
            const double d = d_vec.norm(); // distance between nodes
            const double v = v_vec.norm();
            Vector d_vec_hat = d_vec / d; // normalised distance vector
            const double u0 = (**np1).connecting_springs[i2]->u0; // rest distance (for which there is no spring force)
            const double u = d - u0; // spring extension distance

            const double k = (**np1).connecting_springs[i2]->k; // spring constant
            const double b = (**np1).connecting_springs[i2]->b; // dampening factor

            if (i2 < i1) 
            {
                energy += 0.5 * k * u * u; // calculate spring potential energy only once for each pair
            }
       }
    }
    
    return energy;
}

