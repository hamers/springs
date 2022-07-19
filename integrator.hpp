class Vector;
class Node;
class Parameters;

#include <vector>

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

class Integrator
{
    public:
        virtual void integrate_step(std::vector<Node *> &nodes, double &dt, const Parameters &parameters) = 0;
        void compute_all_node_accelerations(std::vector<Node *> &nodes, const Parameters &parameters);
        void compute_node_vel_acc_at_given_pos_vel(std::vector<Node *> &nodes, const Node node, Vector pos, Vector vel, Vector &vel_out, Vector &acc_out, const Parameters &parameters);
};

class Integrator_Euler: public Integrator
{
    public:
        void integrate_step(std::vector<Node *> &nodes, double &dt, const Parameters &parameters);
};

class Integrator_RK4: public Integrator
{
    public:
        void integrate_step(std::vector<Node *> &nodes, double &dt, const Parameters &parameters);
};

#endif
