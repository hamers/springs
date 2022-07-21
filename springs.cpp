#include "tools.hpp"
#include "integrator.hpp"
#include "parameters.hpp"
#include "springs.hpp"

int Node::max_index = 0;

void parse_args(int argc, char* argv[], double &m, double &k, double &b, Parameters &parameters)
{
    std::cout << argc << std::endl;
    if (argc == 1)
    {
        std::cerr << "Usage: please provide command line arguments: node mass (-m), the spring constant (-k), spring dampening factor (-b), and integration scheme (-i; 0: Euler; 1: fourth-order Runge-Kutta)" << std::endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        int opt = 1;
        
        while ((opt = getopt(argc, argv, "m:k:b:i:o:")) != -1)
        {
            switch (opt)
            {
                case 'm': m = atof(optarg); break;
                case 'k': k = atof(optarg); break;
                case 'b': b = atof(optarg); break;
                case 'i': parameters.integration_scheme = atoi(optarg); break;
                case 'o': parameters.output_filename = std::string(optarg); break;
            }
        }
    }
}
    
void ClientCode(int argc, char *argv[])
{
    
    double m = 1, k = 0, b = 0; // Mass m, spring constant k, and spring dampening coefficient b are passed through the command line, whereas other parameters are assumed to be defined in parameters.cpp
    Parameters parameters;
    parse_args(argc, argv, m, k, b, parameters);
    std::cout << "m " << m << " k " << k << " b " << b << " integration scheme " << parameters.integration_scheme << " output filename " << parameters.output_filename << std::endl;

    Node_Collection node_collection = Node_Collection(&parameters);
    Node_Factory_Box_of_Nodes node_factory = Node_Factory_Box_of_Nodes(m,k,b);
    node_collection.generate_nodes(node_factory);

    node_collection.file_dump();
    
    std::cout << "Node collection initialization complete" << std::endl;
    node_collection.print();

    /* Outer loop / high level client usage of Node_Collection */ 
    const double t_end = parameters.end_time;
    const int N_snap = parameters.N_snap;
    const double dt_snap = t_end / ((double) N_snap);
    
    double t = 0.0;
    while (t < t_end)
    {
        node_collection.integrate(t);
        node_collection.file_dump();
        
        if (parameters.verbose_flag > 0)
        {
            std::cout << "Outer loop t " << t << " E_err (conservative forces only) " << node_collection.get_relative_energy_error() << std::endl;
        }

        t += dt_snap;
    }

    std::cout << "Post integration" << std::endl;
    node_collection.print();

    std::cout << "Invoking Python3 for data visualization" << std::endl;
    std::string command = "python3 plot_node_data.py --file " + parameters.output_filename;
    system(command.c_str());
}

int main(int argc, char *argv[])
{
    ClientCode(argc, argv);

    return 0;
}
