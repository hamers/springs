#include "tools.hpp"
#include "integrator.hpp"
#include "parameters.hpp"
#include "springs.hpp"

void Node_Collection::integrate(double t_end)
{
    this->check_for_node_collection_initialization();
    
    if (t_end - this->t >= this->dt) // Make sure the internal time does not overshoot the user-desired end time
    {
        this->dt = this->t - t_end;
    }

    while (this->t < t_end) // Internal time loop
    {
        if (this->t + this->dt >= t_end) // Adjust internal time step when t_end would be overshot
        {
            this->dt = t_end - this->t;
        }

        this->t += dt;
        this->integrator->integrate_step(this->nodes, dt, *(this->parameters)); // This will use the old dt, and update dt for the next loop iteration
                
        if (this->parameters->verbose_flag > 1)
        {
            std::cout << "internal loop t " << t << " dt " << dt << std::endl;
        }
    }
}

void Integrator::compute_node_vel_acc_at_given_pos_vel(const std::vector<Node *> &nodes, const Node node, Vector pos, Vector vel, Vector &vel_out, Vector &acc_out, const Parameters &parameters)
{
    const double m = node.mass;
    vel_out = vel;
    acc_out = Vector(0.0,0.0,0.0);

    Vector r_vec = pos - parameters.SMBH_pos;
    const double r = r_vec.norm();
    acc_out += r_vec * (-parameters.CONST_G * parameters.SMBH_mass) / (r*r*r);

    int i2 = 0;
    for (auto np2 = node.connected_nodes.begin(); np2 != node.connected_nodes.end(); np2++, i2++)
    {
       Vector d_vec = (**np2).pos - pos; // `distance vector': points from input node to node 2
       Vector v_vec = (**np2).vel - vel; // relative velocity vector
       const double v = v_vec.norm();
       const double d = d_vec.norm(); // distance between nodes
       Vector d_vec_hat = d_vec / d; // normalised distance vector
       const double u0 = node.connecting_springs[i2]->u0; // rest distance (for which there is no spring force)
       const double u = d - u0; // spring extension distance

       const double k = node.connecting_springs[i2]->k; // spring constant
       const double b = node.connecting_springs[i2]->b; // dampening factor

       acc_out += d_vec_hat * k * u / m; // acceleration due to spring force
       acc_out += v_vec * b / m; // acceleration due to spring dampening
    }
}

void Integrator_Euler::integrate_step(std::vector<Node *> &nodes, double &dt, const Parameters &parameters)
{
    /* Simple Euler scheme with adaptive time steps */
    
    for (auto &np : nodes)
    {
        Integrator::compute_node_vel_acc_at_given_pos_vel(nodes, *np, (*np).pos, (*np).vel, (*np).vel, (*np).acc, parameters);
    }
  
    double dt_min = 1e100;
        
    for (auto np = nodes.begin(); np != nodes.end(); np++)
    {        
        (**np).pos += (**np).vel * dt;
        (**np).vel += (**np).acc * dt;
        dt_min = std::min(dt_min, parameters.eta * (**np).vel.norm() / ((**np).acc.norm())); // simple time-step criterion based on velocity and acceleration
    }
    
    dt = dt_min;
}

void Integrator_RK4::integrate_step(std::vector<Node *> &nodes, double &dt, const Parameters &parameters)
{
    /* "Classic" fourth-order Runge-Kutta method with adaptive time steps */
    
    std::vector<Vector> k_vels; // temporary container for velocities
    std::vector<Vector> k_accs; // temporary container for accelerations
       
    const double half_dt = dt * CONST_1DIV_2;
    const double sixth_dt = dt * CONST_1DIV_6;

    int i=0;    
    for (auto np = nodes.begin(); np != nodes.end(); np++)
    {
        Node node = **np;

        Vector k1_vel,k1_acc;
        Integrator::compute_node_vel_acc_at_given_pos_vel(nodes, node, node.pos,                        node.vel,                       k1_vel, k1_acc, parameters);

        Vector k2_vel,k2_acc;
        Integrator::compute_node_vel_acc_at_given_pos_vel(nodes, node, node.pos + k1_vel * half_dt,     node.vel + k1_acc * half_dt,    k2_vel, k2_acc, parameters);

        Vector k3_vel,k3_acc;
        Integrator::compute_node_vel_acc_at_given_pos_vel(nodes, node, node.pos + k2_vel * half_dt,     node.vel + k2_acc * half_dt,    k3_vel, k3_acc, parameters);
        
        Vector k4_vel,k4_acc;
        Integrator::compute_node_vel_acc_at_given_pos_vel(nodes, node, node.pos + k3_vel * dt,          node.vel + k3_acc * dt,         k4_vel, k4_acc, parameters);

        Vector k_vel = (k1_vel + k2_vel*2.0 + k3_vel*2.0 + k4_vel)*CONST_1DIV_6;
        Vector k_acc = (k1_acc + k2_acc*2.0 + k3_acc*2.0 + k4_acc)*CONST_1DIV_6;
        
        k_vels.push_back(k_vel);
        k_accs.push_back(k_acc);
    }

    double dt_min = 1e100;    
    i = 0;
    for (auto np = nodes.begin(); np != nodes.end(); np++, i++)
    {
        (**np).pos += k_vels[i] * dt;
        (**np).vel += k_accs[i] * dt;
        (**np).acc = k_accs[i];
        
        dt_min = std::min(dt_min, parameters.eta * (**np).vel.norm() / ((**np).acc.norm())); // simple time-step criterion based on velocity and acceleration
    }
    
    dt = dt_min;
}
