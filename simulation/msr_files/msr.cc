#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

class MSRSimulation {
private:
    const int nx, ny, n_steps;
    const double dt, dx, dy;
    const double x_min, x_max, y_min, y_max;
    
    // Physical parameters
    const double U0 = 1.0;          // Base flow velocity
    const double kappa = 0.1;       // Thermal diffusivity
    const double beta = 0.001;      // Thermal expansion coefficient
    const double T0 = 600.0;        // Reference temperature
    const double g = 9.81;          // Gravity
    
    // Fields
    std::vector<std::vector<double>> T;     // Temperature
    std::vector<std::vector<double>> u;     // x-velocity
    std::vector<std::vector<double>> v;     // y-velocity
    std::vector<std::complex<double>> Z;    // Complex coordinates

public:
    MSRSimulation(int nx_=100, int ny_=100, int steps=200) 
        : nx(nx_), ny(ny_), n_steps(steps),
          dt(0.01), dx(10.0/nx), dy(10.0/ny),
          x_min(-5), x_max(5), y_min(-5), y_max(5),
          T(ny_, std::vector<double>(nx_, T0)),
          u(ny_, std::vector<double>(nx_, 0.0)),
          v(ny_, std::vector<double>(nx_, 0.0)),
          Z(nx_ * ny_) {
        initialize_fields();
    }

    void initialize_fields() {
        // Setup complex coordinate grid
        for(int i = 0; i < ny; i++) {
            for(int j = 0; j < nx; j++) {
                double x = x_min + j * dx;
                double y = y_min + i * dy;
                Z[i*nx + j] = std::complex<double>(x, y);
            }
        }

        // Initialize temperature field with heat source
        int center_x = nx/2;
        int center_y = ny/2;
        int radius = nx/10;
        
        for(int i = 0; i < ny; i++) {
            for(int j = 0; j < nx; j++) {
                double r2 = std::pow(j - center_x, 2) + std::pow(i - center_y, 2);
                if(r2 < radius*radius) {
                    T[i][j] = T0 + 100.0;
                }
            }
        }
    }

    std::complex<double> complex_velocity(std::complex<double> z, double temp) {
        std::complex<double> w = U0 * (1.0 + 0.001*(temp - T0));
        
        // Add circulation points
        std::complex<double> z1(2, 2), z2(-2, -2);
        if(abs(z - z1) > 0.1) w += 0.5/(z - z1);
        if(abs(z - z2) > 0.1) w += 0.5/(z - z2);
        
        // Add buoyancy
        w += std::complex<double>(0, g * beta * (temp - T0));
        
        return w;
    }

    void update_velocity() {
        for(int i = 0; i < ny; i++) {
            for(int j = 0; j < nx; j++) {
                auto w = complex_velocity(Z[i*nx + j], T[i][j]);
                u[i][j] = std::real(w);
                v[i][j] = std::imag(w);
            }
        }
    }

    void update_temperature() {
        std::vector<std::vector<double>> T_new = T;
        
        for(int i = 1; i < ny-1; i++) {
            for(int j = 1; j < nx-1; j++) {
                // Diffusion term
                double laplacian = (T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] - 4*T[i][j])/(dx*dx);
                
                // Advection term
                double dTdx = (T[i][j+1] - T[i][j-1])/(2*dx);
                double dTdy = (T[i+1][j] - T[i-1][j])/(2*dy);
                
                // Heat source
                double q = 0.0;
                double x = x_min + j*dx;
                double y = y_min + i*dy;
                if(x*x + y*y < 1.0) q = 10.0;
                
                T_new[i][j] = T[i][j] + dt*(kappa * laplacian - u[i][j]*dTdx - v[i][j]*dTdy + q);
            }
        }
        T = T_new;
    }

    void run_simulation(const std::string& filename) {
        std::ofstream outfile(filename, std::ios::binary);
        
        // Write grid dimensions and number of steps
        outfile.write(reinterpret_cast<const char*>(&nx), sizeof(int));
        outfile.write(reinterpret_cast<const char*>(&ny), sizeof(int));
        outfile.write(reinterpret_cast<const char*>(&n_steps), sizeof(int));
        
        // Write spatial bounds
        outfile.write(reinterpret_cast<const char*>(&x_min), sizeof(double));
        outfile.write(reinterpret_cast<const char*>(&x_max), sizeof(double));
        outfile.write(reinterpret_cast<const char*>(&y_min), sizeof(double));
        outfile.write(reinterpret_cast<const char*>(&y_max), sizeof(double));
        
        for(int step = 0; step < n_steps; step++) {
            update_velocity();
            update_temperature();
            
            // Write fields
            for(int i = 0; i < ny; i++) {
                for(int j = 0; j < nx; j++) {
                    outfile.write(reinterpret_cast<const char*>(&u[i][j]), sizeof(double));
                    outfile.write(reinterpret_cast<const char*>(&v[i][j]), sizeof(double));
                    outfile.write(reinterpret_cast<const char*>(&T[i][j]), sizeof(double));
                }
            }
        }
        outfile.close();
    }
};

int main() {
    MSRSimulation sim;
    sim.run_simulation("msr_data.bin");
    return 0;
}