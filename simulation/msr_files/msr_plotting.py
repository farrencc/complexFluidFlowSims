import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import struct

class MSRVisualizer:
    def __init__(self, filename):
        self.load_data(filename)
        self.setup_figure()
        
    def load_data(self, filename):
        """Load simulation data from binary file"""
        with open(filename, 'rb') as f:
            # Read header information
            self.nx = struct.unpack('i', f.read(4))[0]
            self.ny = struct.unpack('i', f.read(4))[0]
            self.n_steps = struct.unpack('i', f.read(4))[0]
            
            # Read spatial bounds
            self.x_min = struct.unpack('d', f.read(8))[0]
            self.x_max = struct.unpack('d', f.read(8))[0]
            self.y_min = struct.unpack('d', f.read(8))[0]
            self.y_max = struct.unpack('d', f.read(8))[0]
            
            # Create coordinate grids
            self.x = np.linspace(self.x_min, self.x_max, self.nx)
            self.y = np.linspace(self.y_min, self.y_max, self.ny)
            self.X, self.Y = np.meshgrid(self.x, self.y)
            
            # Read velocity and temperature fields
            self.u = np.zeros((self.n_steps, self.ny, self.nx))
            self.v = np.zeros((self.n_steps, self.ny, self.nx))
            self.T = np.zeros((self.n_steps, self.ny, self.nx))
            
            for t in range(self.n_steps):
                for i in range(self.ny):
                    for j in range(self.nx):
                        self.u[t,i,j] = struct.unpack('d', f.read(8))[0]
                        self.v[t,i,j] = struct.unpack('d', f.read(8))[0]
                        self.T[t,i,j] = struct.unpack('d', f.read(8))[0]
    
    def setup_figure(self):
        """Initialize the figure and axes"""

        # Set up colormaps and normalizations
        print('Setting up colormaps and normalizations')
        print(f'Temperature range: {np.min(self.T):.2f} to {np.max(self.T):.2f}')
        self.temp_norm = plt.Normalize(vmin=np.min(self.T), vmax=np.max(self.T))
        print(f'Velocity range: {np.min(self.u**2 + self.v**2)} to {np.max(self.u**2 + self.v**2)}')
        self.speed_norm = plt.Normalize(vmin=0, vmax=np.max(np.sqrt(self.u**2 + self.v**2)))
        
        # Temperature
        print(f'Creating figure with {self.n_steps} frames')
        self.fig, (self.ax1, self.ax2) = plt.subplots(1, 2, figsize=(15, 6))
        temp_plot = self.ax1.pcolormesh(self.X, self.Y, self.T[0], cmap='hot', norm=self.temp_norm)
        self.fig.colorbar(temp_plot, ax=self.ax1, label='Temperature')

        # Velocity
        speed0 = np.sqrt(self.u[0]**2 + self.v[0]**2)
        self.strm = self.ax2.streamplot(self.X, self.Y, self.u[0], self.v[0],
                                        density=2, color=speed0,
                                        cmap='viridis', norm=self.speed_norm)
        self.fig.colorbar(self.strm.lines, ax=self.ax2, label='Velocity magnitude')

        self.ax1.tick_params(
            axis='both',       # changes apply to both axes
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            left = False,      # ticks along the left edge are off
            right = False,     # ticks along the right edge are off
            labelleft=False,   # labels along the left edge are off
            labelbottom=False) # labels along the bottom edge are off

        self.ax2.tick_params(
            axis='both',       # changes apply to both axes
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            left = False,      # ticks along the left edge are off
            right = False,     # ticks along the right edge are off
            labelleft=False,   # labels along the left edge are off
            labelbottom=False) # labels along the bottom edge are off
    #
        self.fig.suptitle('Molten Salt Reactor Simulation')
        
    
    def update(self, frame):
        """Update function for animation"""
        self.ax1.clear()
        self.ax2.clear()

        # Update min/max for this frame
        current_temp_min = np.min(self.T[frame])
        current_temp_max = np.max(self.T[frame])
        self.temp_norm.vmin = current_temp_min
        self.temp_norm.vmax = current_temp_max
        
        temp = self.ax1.pcolormesh(self.X, self.Y, self.T[frame],
                                cmap='hot', norm=self.temp_norm)
        
        self.ax1.set_title(f'Temperature Field (t={frame*0.01:.2f}s)')
        
        # Plot velocity field
        speed = np.sqrt(self.u[frame]**2 + self.v[frame]**2)
        strm = self.ax2.streamplot(self.X, self.Y, self.u[frame], self.v[frame],
                                density=2, color=speed,
                                cmap='viridis', norm=self.speed_norm)
        self.ax2.set_title(f'Velocity Field (t={frame*0.01:.2f}s)')
        
        
        for ax in [self.ax1, self.ax2]:
            ax.set_xlim(self.x_min, self.x_max)
            ax.set_ylim(self.y_min, self.y_max)
            ax.set_aspect('equal')
    
    def create_animation(self, filename=None):
        """Create and save the animation"""
        print('Creating animation...')
        anim = FuncAnimation(self.fig, self.update, frames=self.n_steps,
                           interval=25, blit=False)

        print(f'Saving animation to {filename}')
        if filename:
            anim.save(filename, writer='pillow')
        # plt.show()
        return anim

if __name__ == "__main__":
    viz = MSRVisualizer("msr_data.bin")
    viz.create_animation("msr_simulation.gif")
    print('Done!')