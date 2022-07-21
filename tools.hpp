#include <iostream>
#include <math.h>
#include <vector>
#include <map>
#include <random>
#include <unistd.h> // for option parsing
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>

template <class T>
using owner = T;

#define OUTPUT_PRECISION 15
#define CONST_1DIV_2 (double) 0.5
#define CONST_1DIV_6 (double) 1.0/6.0

class Integrator;
class Parameters;
class Node_Factory;

#ifndef VECTOR_H
#define VECTOR_H
class Vector
{
    public:
        double x,y,z;
        
        Vector() : x(0), y(0), z(0) {}
        Vector(double x, double y, double z) : x(x), y(y), z(z) {}
        
        double norm();
        double norm_squared();
        double dot(Vector &v);
                
        Vector operator = (const Vector &v);
        Vector operator + (const Vector &v);
        Vector operator += (const Vector &v);
        Vector operator - (const Vector &v);
        Vector operator * (const double s);
        Vector operator / (const double s);

        void print();
    private:
    
};
#endif

#ifndef SPRING_H
#define SPRING_H
struct Spring
{
    double k,b,u0;
    Spring();
    Spring(double k_, double b_, double u0_);
    ~Spring();
};
#endif

#ifndef NODE_H
#define NODE_H
class Node
{
    public:
        Vector pos;
        Vector vel;
        Vector acc;

        double mass;
        int index;
        std::vector<Node*> connected_nodes;
        std::vector<Spring*> connecting_springs;

        Node();
        Node(Vector pos_, Vector vel_, double mass_);
        ~Node();
        void print();
        void connect_node(Node *n, Spring *spring);
        int get_number_of_neighbours();

    private:
        static int max_index;

};
#endif

#ifndef NODE_COLLECTION_H
#define NODE_COLLECTION_H
class Node_Collection
{
    public:
        Node_Collection(Parameters *parameters);
        ~Node_Collection();

        void set_parameters(Parameters *parameters);
        void generate_nodes(Node_Factory &node_factory); // Use Factory Method to set up the nodes
        void integrate(double t_end);

        double get_relative_energy_error();
        
        void file_dump();

        void print();

    private:
        Parameters *parameters; // NOT owned by Node_Collection object
        owner<Integrator> *integrator; // Owned by Node_Collection object
        
        void check_for_node_collection_initialization();
        void set_integrator();

        std::vector<Node*> nodes;
        
        double t,dt;
        
        bool output_file_created;

        double calculate_energy();
        double current_energy;
        double initial_energy;
        bool initial_energy_calculated;

};
#endif

#ifndef NODE_FACTORY_H
#define NODE_FACTORY_H
class Node_Factory
{
    /* Class that contains methods to set up specific configurations (specified in subclasses) */
    
    public:
        Node_Factory() {} ;

        virtual std::vector<Node*> generate_nodes(const Parameters &parameters) = 0;
        virtual std::vector<Node*> set_up_2d_box_of_nodes(const Parameters &parameters) = 0;
        void connect_all_nodes_in_node_collection(std::vector<Node*> &nodes, const Parameters &parameters); // General function
            
    protected:
        double mass_, k_, b_;

};

class Node_Factory_Box_of_Nodes : public Node_Factory
{
    public:
        Node_Factory_Box_of_Nodes(double mass, double k, double b)
        {
            mass_ = mass;
            k_ = k;
            b_ = b;
        }
    private:
        std::vector<Node*> generate_nodes(const Parameters &parameters);
        std::vector<Node*> set_up_2d_box_of_nodes(const Parameters &parameters);
};

#endif


