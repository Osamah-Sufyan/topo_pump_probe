import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
from numba import njit, prange
import numpy as np
import math
from math import e
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy.linalg import eigh
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from numba import njit
from scipy.ndimage import gaussian_filter1d
import random
import json
from numba import njit, prange
from scipy.signal import savgol_filter, butter, filtfilt

from matplotlib.patches import Arc, RegularPolygon
from numpy import radians as rad


def node_dist(x, y, cx, cy):
    """Calculate the Manhattan distance of a node from the lattice center."""
    return abs(cx - x) + abs(cy - y)

def remove_unwanted_nodes(G, m):
    """Remove nodes that don't belong to an m-layer hexagonal ring."""
    cx, cy = m - 0.5, 2 * m - (m % 2)  # Center of the lattice
    
    unwanted = []
    for n in G.nodes:
        x, y = n
        if node_dist(x, y, cx, cy) > 2 * m:
            unwanted.append(n)
    
    for n in unwanted:
        G.remove_node(n)
    
    return G

def get_boundary_atoms(G):
    """
    Extract boundary atoms by identifying nodes on the perimeter.
    A node is on the boundary if it has any "missing" neighbors compared
    to what a complete hexagonal lattice would have.
    """
    boundary_nodes = set()
    
    expected_neighbor_vectors = [
        (0, 1),     # up
        (1, 0),     # right
        (1, -1),    # down-right
        (0, -1),    # down
        (-1, 0),    # left
        (-1, 1),    # up-left
        (-1, -1),
        (-1, 0),
        (0, -1),
        (1, 1)
    ]
    
    for node in G.nodes():
        x, y = node
        neighbor_positions = {(n[0]-x, n[1]-y) for n in G.neighbors(node)}
        
        for dx, dy in expected_neighbor_vectors:
            potential_neighbor = (x + dx, y + dy)
            if potential_neighbor not in G.nodes():
                boundary_nodes.add(node)
                break
    
    return list(boundary_nodes)

def get_vertical_edge1(G, pos):
    """
    Extract nodes that form the rightmost vertical edge of the hexagonal flake.
    """
    x_values = sorted(set(x for x, _ in pos.values()), reverse=True)
    second_max_x = x_values[1] if len(x_values) > 1 else None


    max_x = max(x for x, _ in pos.values())
    
    vertical_edge_nodes = []
    for node, position in pos.items():
        if abs(position[0] - second_max_x) < 0.0001 or abs(position[0] - max_x) < 0.0001:
            vertical_edge_nodes.append(node)
    
    vertical_edge_nodes.sort(key=lambda n: pos[n][1], reverse=True)
    return vertical_edge_nodes


def get_vertical_edge2(G, pos):
    """
    Extract nodes that form the rightmost vertical edge of the hexagonal flake.
    """
    x_values = sorted(set(x for x, _ in pos.values()), reverse=True)
    second_max_x = x_values[-1] if len(x_values) > 1 else None


    max_x = max(x for x, _ in pos.values())
    
    x_value_1 = x_values[-1]
    x_value_2 = x_values[-2]
    vertical_edge_nodes = []
    for node, position in pos.items():
        if abs(position[0] - x_value_1) < 0.0001 or abs(position[0] - x_value_2) < 0.0001:
            vertical_edge_nodes.append(node)
    
    vertical_edge_nodes.sort(key=lambda n: pos[n][1], reverse=True)
    return vertical_edge_nodes

def plot_hexagonal_flake(G, pos, m, vertical_edge_nodes):
    """Plot the hexagonal flake with highlighted vertical edge."""
    edge_G = G.subgraph(vertical_edge_nodes)
    
    plt.figure(figsize=(8, 8))
    
    # Plot all atoms in grey
    nx.draw_networkx_nodes(G, pos, node_color='lightgrey', node_size=100)
    nx.draw_networkx_edges(G, pos, alpha=0.2)
    
    # Highlight the vertical edge in red
    nx.draw_networkx_nodes(edge_G, pos, node_color='red', node_size=100)
    ax = plt.gca()
    
    def drawCirc(ax, radius, centX, centY, angle_, theta2_, color_='black'):
        #========Line
        arc = Arc([centX, centY], radius, radius, angle=angle_,
                  theta1=0, theta2=theta2_, capstyle='round', linestyle='-', lw=2, color=color_)
        ax.add_patch(arc)

        #========Create the arrow head
        endX = centX + (radius / 2) * np.cos(rad(theta2_ + angle_))  # Do trig to determine end position
        endY = centY + (radius / 2) * np.sin(rad(theta2_ + angle_))

        ax.add_patch(  # Create triangle as arrow head
            RegularPolygon(
                (endX, endY),  # (x,y)
                numVertices=3,  # number of vertices
                radius=radius / 9,  # radius
                orientation=rad(angle_ + theta2_),  # orientation
                color=color_
            )
        )
        ax.set_xlim([centX - radius, centX + radius])
        ax.set_ylim([centY - radius, centY + radius])
        # Make sure you keep the axes scaled or else arrow will distort

    # Draw the circular arrow on top of the hexagon plot
    drawCirc(ax, 30, 20, 20, 0, 250, color_='blue')
        
    plt.axis('equal')
    plt.title(f'Hexagonal Flake (m={m}) with Vertical Edge Highlighted')
    #plt.savefig('plots/hexagonal_flake.png')
    plt.show()




def plot_edge():
    # Parameters
    m = 5  # Flake with side size m=3
    
    # Create the hexagonal lattice and process it
    G = nx.hexagonal_lattice_graph(2 * m - 1, 2 * m - 1, periodic=False, with_positions=True)
    G = remove_unwanted_nodes(G, m)
    pos = nx.get_node_attributes(G, 'pos')
    
    # Scale positions
    scale = 2.68
    for node, position in pos.items():
        pos[node] = tuple(round(i * scale, 4) for i in position)
    
    # Get boundary atoms
    boundary_atoms = get_boundary_atoms(G)
    
    # Get vertical edge
    vertical_edge_nodes = get_vertical_edge2(G, pos)
    
    # Plot the result
    plot_hexagonal_flake(G, pos, m, vertical_edge_nodes)
    
    # Print information

m = 5 # Flake with side size m

# Create the hexagonal lattice and process it
G = nx.hexagonal_lattice_graph(2 * m - 1, 2 * m - 1, periodic=False, with_positions=True)
G = remove_unwanted_nodes(G, m)
pos = nx.get_node_attributes(G, 'pos')

# Scale positions
scale = 2.68
for node, position in pos.items():
    pos[node] = tuple(round(i * scale, 4) for i in position)

# Get boundary atoms with improved detection
boundary_atoms = get_boundary_atoms(G)

# Create a new graph with only boundary atoms
boundary_G = G.subgraph(boundary_atoms)
  
dict0= pos

vertical_edge_nodes2 = get_vertical_edge2(G, pos)
vertical_edge_nodes1 = get_vertical_edge1(G, pos)
# plt.figure(figsize=(8, 8))
# for node in vertical_edge_nodes:
#     x, y = pos[node]
#     plt.scatter(x, y, color='red', s=100)
#     plt.gca().set_aspect('equal', adjustable='box')
# plt.show()



q= list(dict0.values())



r1=[]
for i in np.arange(0,len(dict0)):
    s=q[i][0]
    r1.append(s)
r2=[]
for i in np.arange(0,len(dict0)):
    s=q[i][1]
    r2.append(s)



A=np.zeros((len(dict0),3))
for i in np.arange(0, len(dict0)):
    A[i][1]= r1[i]
    A[i][2]= r2[i]
    A[i] =i
for i in np.arange(0, len(dict0)):
    A[i][1]= r1[i]
    A[i][2]= r2[i]



coord_to_index = {tuple(A[i, 1:]): int(A[i, 0]) for i in range(len(A))}

# Extract indices of vertical_edge_nodes


vertical_edge_indices1 = [coord_to_index[tuple(pos[node])] for node in vertical_edge_nodes1]
print(len(vertical_edge_indices1))
vertical_edge_indices2 = [coord_to_index[tuple(pos[node])] for node in vertical_edge_nodes2]
print(len(vertical_edge_indices2))

vertical_edge_indices_all = vertical_edge_indices1 + vertical_edge_indices2
print(len(vertical_edge_indices1))












Asx= np.arange(0.5,102.5,1.5)
Asx2 = np.arange(1.34, 1.34*1000, 4.0200000000000005)

Acells=[]
Bcells =[]
for i in np.arange(0,len(dict0)):
    if A[i][1] in np.round(Asx2,4):
    
        Acells.append(A[i])
    else:
        Bcells.append(A[i])
        
Acellsin=[]
for i in np.arange(0, len(Acells)):
    x= Acells[i][0]
    Acellsin.append(x)
    
Bcellsin=[]
for i in np.arange(0, len(Bcells)):
    x= Bcells[i][0]
    Bcellsin.append(x)

a_1= (3/2 * scale , np.sqrt(3)/2*scale)
a_2= (-3/2* scale, np.sqrt(3)/2* scale)
a_3= (0,-1*scale*np.sqrt(3))
pAcells=[]
for i in np.arange(0,len(Acells)):
    for j in np.arange(0,len(Acells)):
        dif =[Acells[i][1]-Acells[j][1], Acells[i][2]-Acells[j][2]]
        if round(dif[0]-a_1[0],2)==0 and round(dif[1]-a_1[1])==0:
           x = Acells[i][0]
           y = Acells[j][0]
           pAcells.append((x,y))
        
        
for i in np.arange(0,len(Acells)):
    for j in np.arange(0,len(Acells)):
        dif =[Acells[i][1]-Acells[j][1], Acells[i][2]-Acells[j][2]]
        if round(dif[0]-a_2[0],2)==0 and round(dif[1]-a_2[1])==0:
           x = Acells[i][0]
           y = Acells[j][0]
           pAcells.append((x,y))
for i in np.arange(0,len(Acells)):
    for j in np.arange(0,len(Acells)):
        dif =[Acells[i][1]-Acells[j][1], Acells[i][2]-Acells[j][2]]
        if round(dif[0]-a_3[0],2)==0 and round(dif[1]-a_3[1])==0:
           x = Acells[i][0]
           y = Acells[j][0]
           pAcells.append((x,y))
        
        

pBcells=[]
a_1b= (-3/2*scale, -np.sqrt(3)/2*scale)
a_2b= (3/2*scale, -1*scale*np.sqrt(3)/2)
a_3b= (0,scale*np.sqrt(3))

for i in np.arange(0,len(Bcells)):
    for j in np.arange(0,len(Bcells)):
        dif =[Bcells[i][1]-Bcells[j][1], Bcells[i][2]-Bcells[j][2]]
        if round(dif[0]-a_1b[0],2)==0 and round(dif[1]-a_1b[1])==0:
           x = Bcells[i][0]
           y = Bcells[j][0]
           pBcells.append((x,y))
        
        
for i in np.arange(0,len(Bcells)):
    for j in np.arange(0,len(Bcells)):
        dif =[Bcells[i][1]-Bcells[j][1], Bcells[i][2]-Bcells[j][2]]
        if round(dif[0]-a_2b[0],2)==0 and round(dif[1]-a_2b[1])==0:
           x = Bcells[i][0]
           y = Bcells[j][0]
           pBcells.append((x,y))
for i in np.arange(0,len(Bcells)):
    for j in np.arange(0,len(Bcells)):
        dif =[Bcells[i][1]-Bcells[j][1], Bcells[i][2]-Bcells[j][2]]
        if round(dif[0]-a_3b[0],2)==0 and round(dif[1]-a_3b[1])==0:
           x = Bcells[i][0]
           y = Bcells[j][0]
           pBcells.append((x,y))
        
        
def NN():
    near= []
    for i in np.arange(0,len(dict0)):
        for j in np.arange(0,len(dict0)):
            if round(np.sqrt((A[i][1]-A[j][1])**2+(A[i][2]-A[j][2])**2), 4)==2.68:
                near.append((i,j))
    return near


def NNN():
    nnear= []
    for i in np.arange(0,len(dict0)):
        for j in np.arange(0,len(dict0)):
            if round(np.sqrt((A[i][1]-A[j][1])**2+(A[i][2]-A[j][2])**2), 4)==4.636:
                nnear.append((i,j))
    return nnear

def hamiltonian(M,t_1,t_2):
    hamiltonian = np.zeros((len(dict0),len(dict0)),dtype=np.complex128)
    j=complex(0,1)
    for k in np.arange(0,len(dict0)):
        for l in np.arange(0,len(dict0)):
            if k in Acellsin:
                hamiltonian[k,k] = -M
            
                
            elif k in Bcellsin:
               
                hamiltonian[k,k] = M


    for k in np.arange(0,len(dict0)):
        for l in np.arange(0,len(dict0)):        
            if (k,l) in pAcells:
                hamiltonian[k,l] = t_2*j
                hamiltonian[l,k] = -t_2*j
                
            if (k,l) in pBcells:
                hamiltonian[k,l] = t_2*j
                hamiltonian[l,k] = -t_2*j
            

    for (k,l) in NN():
        hamiltonian[k,l] = t_1
        hamiltonian[l,k] = t_1
    return hamiltonian



#############section 2: solving the time-indepednent Schrödinger equation in real space and plotting the eigenvalues vs space index and spatial distribution of the eigenvectors, this part has been tested so you can skip it as well#################
t1 = -0.1
Delta = 0.01
t2 = 0.03
H=hamiltonian(Delta,t1,t2)




# Confirm the lengths of both t and A_x

# Plot the results





amplitudes = [1e0, 2e-1,  2e-2, 2e-3]

t2_values = np.arange(0,1, 1/20)  # Angular frequency in atomic units
current_x = []
current_y = []

 # Angular frequency in atomic units
plt.figure(figsize=(12, 6))
for factor in t2_values:
    t1 = -0.1
    Delta = 0.01
    t2 = factor
    H=hamiltonian(Delta,t1,t2)
    A_0 = 2e-1



    ew, ev = LA.eigh(H)
    dt = 12000
    #omega_0 = 0.08 # 570 nm
    omega_0 = 0.057 # 800 nm

    n_0 = 10

# Original time setup
    t = np.arange(0, 2 * np.pi * n_0 / omega_0, 2 * np.pi * n_0 / omega_0 / dt)

    # Generate the A_x signal
    sin2_term = np.sin(omega_0 * t / (2 * n_0)) ** 2
    A_x = A_0 * sin2_term * np.sin(omega_0 * t)
    A_y = -A_0 * sin2_term * np.cos(omega_0 * t)

    # Extend time by 25% and generate corresponding extended A_x values
    extra_time = len(t) * 0.025  # 25% extra time
    t_after = np.linspace(t[-1], t[-1] + extra_time, 5000)  # New time points

    # Extend the time vector
    t = np.concatenate((t, t_after))

    # Extend A_x by zero padding for the extra time
    zero_padding = np.zeros(len(t_after))  # Padding zeros for the extra time
    A_x = np.concatenate((A_x, zero_padding))  # Append zero padding to A_x 
    
   
    A_y = np.concatenate((A_x, zero_padding))

    time_steps = t

    # plot A_x and A_y vs time
    # plt.figure(figsize=(12, 6))
    # plt.plot(t, A_x, label="A_x", alpha=0.7, linewidth=2)
    # plt.plot(t, A_y, label="A_y", alpha=0.7, linewidth=2)
    # plt.title("Field Amplitudes A_x and A_y vs Time")
    # plt.xlabel(r"$t[a.u.]$", fontsize=14)
    # plt.ylabel(r"$A_{x,y}$", fontsize=14)
    # plt.xticks(fontsize=14)
    # plt.yticks(fontsize=14)
    # plt.legend()
    # plt.grid(True)
    # plt.savefig('field_amplitudes.pdf')
    # plt.show()

    # Define timing parameters
    # T_first_pulse = total_time / 8
    # T_decay = T_first_pulse / 2
    # T_buffer = total_time / 4
    # T_buffer_start = T_first_pulse + T_decay
    # T_second_pulse_start = T_buffer_start + T_buffer
    # T_second_pulse_end = T_second_pulse_start + 4*T_first_pulse  # Define end of second pulse

    # # Create envelope for first pulse with smooth decay
    # envelope_first = np.zeros_like(t)
    # mask_first = t <= 2*T_first_pulse
    # envelope_first[mask_first] = np.sin(np.pi * t[mask_first] / (2 * T_first_pulse)) ** 2

    # Add decay period
    # mask_decay = (t > T_first_pulse) & (t <= T_buffer_start)
    # t_decay = t[mask_decay] - T_first_pulse
    # envelope_first[mask_decay] = (np.cos(np.pi * t_decay / (2 * T_decay)) ** 2) * envelope_first[mask_first][-1]

    # Create envelope for second pulse
    # Create envelope for the second pulse
    # envelope_second = np.zeros_like(t)
    # mask_second = (t >= T_second_pulse_start) & (t <= T_second_pulse_end)
    # t_shifted = t[mask_second] - T_second_pulse_start  # Shifted time for second pulse
    # envelope_second[mask_second] = np.sin(np.pi * t_shifted / ((T_second_pulse_end - T_second_pulse_start))) ** 2


    # First pulse (circular polarization)
    # A_x = A_0 * envelope_first * np.sin(omega_0 * t)
    # A_y = -A_0 * envelope_first * np.cos(omega_0 * t)

    # # Second pulse (linear polarization in y direction)
    # A_y_second = A_0 * envelope_second * np.sin(0.25*omega_0 * t)
    # A_x_second = np.zeros_like(t)

    # # Combine pulses
    # A_x_total = A_x + A_x_second
    # A_y_total = A_y + A_y_second
    # A_total = np.sqrt(A_x_total**2 + A_y_total**2)

    # Plotting
    # plt.figure(figsize=(12, 6))
    # plt.plot(t, A_x_total, label="A_x (total)", alpha=0.7)
    # plt.plot(t, A_y_total, label="A_y (total)", alpha=0.7)
    # plt.plot(t, A_total, label="Total Field Amplitude", linestyle="--", color='black')
    # plt.title("Complete Pulse Sequence with Buffer")
    # plt.xlabel("Time (a.u.)")
    # plt.ylabel("Field Amplitude (a.u.)")
    # plt.legend()
    # plt.grid(True)
    # plt.savefig('pulse_sequence.pdf')
    #plt.show()

    Ht=np.zeros((len(t),len(dict0),len(dict0)),dtype=np.complex128)

    @njit(parallel=True)
    def build_Ht(H, A, A_x, A_y, Ht):
        for t_i in prange(len(t)):
            # if t_i > 
            for k in range(len(H)):
                for j in range(len(H)):
                    Ht[t_i, k, j] = H[k, j] * np.exp(-1j * (A[k][1] - A[j][1]) * A_x[t_i])* np.exp(-1j * (A[k][2] - A[j][2]) *(-1)* A_y[t_i])
        return Ht

    Ht = build_Ht(H, A ,A_x,A_y,  Ht)

    wavelength_um = 0.04564 / omega_0
    intensity_Wcm = 3.51e13 * (omega_0 * A_0) ** 2
    print(f"Angular frequency: {omega_0} a.u. -> Wavelength: {wavelength_um:.2f} µm")
    print(f"intensity: {A_0} a.u. -> intensity: {intensity_Wcm / 1e9:.3f} 10^9 W/cm^2")
    @njit
    def Ha(t, time_steps, Ht):
        # Find the index of the closest time step
        t_index = np.argmin(np.abs(time_steps - t))
        return Ht[t_index]

    @njit
    def schrodinger(t, psi, time_steps, Ht):
        return -1j * Ha(t, time_steps, Ht) @ psi

    @njit
    def RK4_step(f, t, y, dt, time_steps, Ht):
        k1 = dt * f(t, y, time_steps, Ht)
        k2 = dt * f(t + 0.5 * dt, y + 0.5 * k1, time_steps, Ht)
        k3 = dt * f(t + 0.5 * dt, y + 0.5 * k2, time_steps, Ht)
        k4 = dt * f(t + dt, y + k3, time_steps, Ht)
        return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    # Time steps and initial conditions
    t_start = time_steps[0]
    t_end = time_steps[-1]
    dt = time_steps[1] - time_steps[0]
    solutions = []

    dim = len(ev)
    for i in range(dim):
        psi0 = ev[:, i]
        psi_t = psi0
        psi_t_list = [psi0]
        for t in time_steps[:-1]:
            psi_t = RK4_step(schrodinger, t, psi_t, dt, time_steps, Ht)
            psi_t_list.append(psi_t)
        solutions.append(np.array(psi_t_list).T)

    @njit
    def calculate_current(solutions, t_index, time_steps, Ht, vertical_edge_indices):
        J = np.zeros(2, dtype=np.complex128)
        H = Ha(time_steps[t_index], time_steps, Ht)
        for l in range(len(solutions)//2 - 1):
            for i, i1 in enumerate(vertical_edge_indices):
                for j, j1 in enumerate(vertical_edge_indices):
                    if j1 and i1 in vertical_edge_indices1:

                        J[0] += (A[i1][1] - A[j1][1]) * -1j * np.conj(solutions[l][i1][t_index]) * H[i1][j1] * solutions[l][j1][t_index]
                        J[1] += (A[i1][2] - A[j1][2]) * -1j * np.conj(solutions[l][i1][t_index]) * H[i1][j1] * solutions[l][j1][t_index]
                    elif j1 and i1 in vertical_edge_indices2:
                        J[0] += (A[i1][1] - A[j1][1]) * -1j * np.conj(solutions[l][i1][t_index]) * H[i1][j1] * solutions[l][j1][t_index]
                        J[1] += -1 * (A[i1][2] - A[j1][2]) * -1j * np.conj(solutions[l][i1][t_index]) * H[i1][j1] * solutions[l][j1][t_index]

        return J
    #second_pulse_indices = np.where((time_steps >= T_second_pulse_start) & (time_steps <= T_second_pulse_end))[0]

    J_t = [calculate_current(solutions, t_index, time_steps, Ht, vertical_edge_indices_all) for t_index in range(len(time_steps))]
    #J_second_pulse = [J_t[index] for index in second_pulse_indices]

    J_x = [J[0] for J in J_t]
    J_y = [J[1] for J in J_t]

    # Multiply by 6
    #J_x_multiplied = [6 * J.real for J in J_x]
    J_y_multiplied = [6 * J.real for J in J_y]
    num = random.randint(0, 10000)
    # Plot J_x and J_y multiplied by 6 vs time
    average_J_x = np.mean(J_x[-2000:])
    average_J_y = np.mean(J_y[-2000:])
    current_x.append(average_J_x)
    current_y.append(average_J_y)
    
    

    #plt.plot(time_steps, J_x_multiplied, label=r'$J_{\perp}$', color='red')
plt.plot(t2_values, current_x, label = r"$\overline{J_x}$")
plt.plot(t2_values, current_y, label = r"$\overline{J_y}$")
plt.xlabel(r'$t_2$ (a.u.)', fontsize=14)
plt.ylabel(r'$I_p$ (a.u.)', fontsize=14)
plt.title(f'parallel current components along the edge, left helical, flakesize= {m}, frequency = {wavelength_um:.2f}, A0 = {A_0}')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=17)
plt.grid(True)
plt.savefig(f'plots/edge_photocurrents_all_neg_coord_t2_{num}.pdf')
plt.show()