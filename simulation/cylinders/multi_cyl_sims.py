import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle

class MultiCylinderFlow:
    def __init__(self, cylinders, U=1.0, nx=100, ny=100, x_range=(-5, 5), y_range=(-5, 5)):
        """
        Initialize flow field for multiple cylinders
        cylinders: list of tuples (x, y, radius) for each cylinder
        U: free stream velocity
        """
        self.cylinders = cylinders
        self.U = U
        self.nx, self.ny = nx, ny
        self.x_range = x_range
        self.y_range = y_range
        
        # Create grid
        self.x = np.linspace(x_range[0], x_range[1], nx)
        self.y = np.linspace(y_range[0], y_range[1], ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        self.Z = self.X + 1j*self.Y

    def complex_potential(self, z):
        """Compute complex potential using superposition method"""
        w = self.U * z  # Free stream
        
        # Add contribution from each cylinder
        for x0, y0, R in self.cylinders:
            z0 = x0 + 1j*y0
            zeta = z - z0
            w += self.U * (R**2/zeta)  # Add dipole for each cylinder
            
        return w

    def get_velocity(self, z):
        """Calculate velocity field using derivative of complex potential"""
        dw = self.U  # Free stream contribution
        
        # Add contribution from each cylinder
        for x0, y0, R in self.cylinders:
            z0 = x0 + 1j*y0
            zeta = z - z0
            dw += -self.U * R**2/zeta**2
            
        return dw

    def compute_flow(self):
        """Compute velocity field"""
        # Calculate complex velocity
        w_prime = self.get_velocity(self.Z)
        
        # Extract velocity components
        u = np.real(w_prime)
        v = -np.imag(w_prime)
        
        # Create mask for points inside cylinders
        mask = np.ones_like(u, dtype=bool)
        for x0, y0, R in self.cylinders:
            mask = mask & (((self.X - x0)**2 + (self.Y - y0)**2) >= R**2)
            
        # Apply mask
        self.u = np.ma.masked_where(~mask, u)
        self.v = np.ma.masked_where(~mask, v)
        
        return self.u, self.v

def animate_flow(cylinders, n_frames=200, interval=50):
    """Create animation of flow around multiple cylinders"""
    # Initialize flow field
    flow = MultiCylinderFlow(cylinders)
    u, v = flow.compute_flow()
    
    # Setup figure
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')
    ax.set_xlim(flow.x_range)
    ax.set_ylim(flow.y_range)
    
    # Add cylinders
    for x0, y0, R in cylinders:
        circle = Circle((x0, y0), R, color='gray', fill=True)
        ax.add_patch(circle)
    
    # Initialize particles
    n_particles = 50
    particles = np.zeros((n_particles, 2, n_frames))
    # Start particles from left side
    particles[:, 0, 0] = flow.x_range[0] + 0.2
    particles[:, 1, 0] = np.linspace(flow.y_range[0]/2, flow.y_range[1]/2, n_particles)

    # Calculate particle trajectories
    dt = 0.1
    for i in range(n_particles):
        for t in range(n_frames-1):
            x, y = particles[i, :, t]
            z = x + 1j*y
            
            # Skip if inside any cylinder
            inside_cylinder = False
            for x0, y0, R in cylinders:
                if ((x - x0)**2 + (y - y0)**2) <= R**2:
                    inside_cylinder = True
                    break
            
            if inside_cylinder:
                particles[i, :, t+1] = particles[i, :, t]
                continue
                
            # Get velocity at current position
            w_prime = flow.get_velocity(z)
            vx = np.real(w_prime)
            vy = -np.imag(w_prime)
            
            # Update position
            particles[i, 0, t+1] = x + vx*dt
            particles[i, 1, t+1] = y + vy*dt

    # Plot streamlines
    ax.streamplot(flow.X, flow.Y, u, v, density=2.0, color='lightblue', linewidth=0.5)
    

    def update(frame):
        # Store circle patches
        circles = [p for p in ax.collections if isinstance(p, Circle)]
        
        # Clear all collections except circles
        for collection in ax.collections[:]:
            if not isinstance(collection, Circle):
                collection.remove()
        
        # Plot current particle positions
        scatter = ax.scatter(particles[:, 0, frame], particles[:, 1, frame], 
                            c='navy', s=20, alpha=0.6)
        
        # Remove ticks and labels
        ax.tick_params(
            axis='both',       # changes apply to both axes
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            left = False,      # ticks along the left edge are off
            right = False,     # ticks along the right edge are off
            labelleft=False,   # labels along the left edge are off
            labelbottom=False) # labels along the bottom edge are off

        # Set ylimit
        ax.set_ylim((-2.5,2.5))

        # Plot trails
        trails = []
        for i in range(n_particles):
            if frame > 0:
                line, = ax.plot(particles[i, 0, max(0, frame-20):frame],
                            particles[i, 1, max(0, frame-20):frame],
                            '-', alpha=0.2,color = 'navy')
                trails.append(line)
    
        # Return all artists that need to be updated
        return circles + [scatter] + trails


    anim = FuncAnimation(fig, update, frames=n_frames,
                        interval=interval, blit=True)
    return anim

if __name__ == "__main__":
    # Example configurations
    
    # Two cylinders in a line
    cylinders_2 = [(-1, 0, 0.5), (1, 0, 0.5)]
    
    # Three cylinders in a line
    cylinders_3 = [(-2, 0, 0.5), (0, 0, 0.5), (2, 0, 0.5)]
    
    # Three cylinders in a triangle
    cylinders_triangle = [(0, 0.6, 0.7), (-0.55, -0.5, 0.15), (0.45, -0.5, 0.3)]
    
    # # Create animation (choose one configuration)
    anim = animate_flow(cylinders_2)
    anim.save('cylinder_2.gif', writer='pillow')

    anim = animate_flow(cylinders_3)
    anim.save('cylinder_3.gif', writer='pillow')

    # anim = animate_flow(cylinders_triangle)
    # anim.save('cylinder_triangle.gif', writer='pillow')