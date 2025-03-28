import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches

# Parameters
U = 1.0  # Free stream velocity
a = 1.0  # Cylinder radius
nx, ny = 50, 50  # Grid points
dt = 0.1  # Time step
n_particles = 74
n_frames = 100

# Create grid
x = np.linspace(-3, 3, nx)
y = np.linspace(-3, 3, ny)
X, Y = np.meshgrid(x, y)
Z = X + 1j*Y

def get_velocity(z):
    z = np.where(abs(z) < 1e-10, 1e-10, z)
    w_prime = U * (1 - a**2/z**2)
    return w_prime.real, -w_prime.imag

# Initialize particles
particles = np.zeros((n_particles, 2, n_frames))
particles[:, 0, 0] = -2.8  # Start from left side
particles[:, 1, 0] = np.linspace(-2.9, 2.9, n_particles)  # Spread vertically

# Calculate particle trajectories
for i in range(n_particles):
    for t in range(n_frames-1):
        z = particles[i, 0, t] + 1j*particles[i, 1, t]
        if abs(z) <= a:  # Skip if inside cylinder
            particles[i, :, t+1] = particles[i, :, t]
            continue
        u, v = get_velocity(z)
        particles[i, 0, t+1] = particles[i, 0, t] + u*dt
        particles[i, 1, t+1] = particles[i, 1, t] + v*dt

# Setup plot
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_aspect('equal')
circle = patches.Circle((0, 0), a, color='gray', fill=True)
ax.add_patch(circle)

# Calculate background streamlines
u, v = get_velocity(X + 1j*Y)
mask = np.abs(X + 1j*Y) >= a
u = np.ma.masked_where(~mask, u)
v = np.ma.masked_where(~mask, v)
ax.streamplot(X, Y, u, v, density=1.5, color='lightblue', linewidth=0.5)

# Animation function
def animate(frame):
    for collection in ax.collections[:]:
        collection.remove()
    # Use add_patch instead of add_collection for circle
    ax.add_patch(circle)
    # Plot current position of particles
    ax.scatter(particles[:, 0, frame], particles[:, 1, frame], 
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
    # Plot trails
    for i in range(n_particles):
        if frame > 0:
            ax.plot(particles[i, 0, max(0, frame-20):frame],
                   particles[i, 1, max(0, frame-20):frame],
                   '-', alpha=0.2, color = 'navy')
    return ax.collections

# Create animation
anim = FuncAnimation(fig, animate, frames=n_frames,
                    interval=50, blit=True)
anim.save('cylinder_flow.gif', writer='pillow')
plt.close()