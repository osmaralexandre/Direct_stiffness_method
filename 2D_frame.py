# This program solves plane frames with hinges and support displacements.
# The structure is solved by using direct stiffness method.
import os
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym

os.system('cls')

sym.init_printing()

fig, axs = plt.subplots(2, 2)

# Matrix with nodal coordinates.
"""
    The first and the second rows represents, respectively, the x and y coordinates of each node.
    Make shure that you made a coherent structure.
"""
coord = np.matrix([[0, 0, 2, 4, 4],
                    [0, 2, 5, 2, 1]])

# Matrix that defines the support conditions. 0 for free and 1 for support.
"""
    Each node of the frame has 3 degrees of freedon (dof). The first and second rows represents the displacements
    in x and y directions and the third row is the rotational dof.
"""
support = np.matrix([[1, 0, 0, 0, 1],
                    [1, 0, 0, 0, 1],
                    [0, 0, 0, 0, 1]])
               
# Matrix that indicates the nodal loads.
"""
    If the load position coincides with 1 in cont, it will correspond to a support displacement.
"""
load = np.matrix([[0, 1, 0, 0, 0.1],
                    [0, 0, 0, 0, 0],
                    [0, 0, 20, 0, 0]])

# Incidence matrix.
"""
    Indicates the inicial and the final node of each element.
"""
inc = np.matrix([[0, 1, 2, 1, 3],
                    [1, 2, 3, 3, 4]])

# Matrix with the position of the hinges.
"""
    1 if the node is articulated in that element. Incidence matrix shows the correlation between the nodes and the elements.
"""
hinge = np.matrix([[0, 0, 0, 1, 1],
                    [0, 0, 0, 1, 0]])

# Stiffness of each element.
EA = np.matrix([1,1,1,1,1])*50
EJ = np.matrix([1,1,1,1,1])*60

# Some important constants.
nn = np.size(coord,1); nb = np.size(inc,1); C = 10**10; cdef = 5; cap = 0.2; cn = 0.2; cq = 0.2; cf = 0.2; d_hinge = 0.05

# Stiffness matrix assembly
ngl = 3*nn
knl = []; kll = []
K = np.zeros((ngl,ngl))
for i in range(nb):
    dx = coord[0,inc[1,i]] - coord[0,inc[0,i]]
    dy = coord[1,inc[1,i]] - coord[1,inc[0,i]]
    L = np.sqrt(dx**2+dy**2)
    c = dx/L; s = dy/L

    # Stiffness matrix in a natural system
    kn = np.matrix([[EA[0,i]/L, 0, 0],
                    [0, 4*EJ[0,i]/L, 2*EJ[0,i]/L],
                    [0, 2*EJ[0,i]/L, 4*EJ[0,i]/L]])
    knl.append(kn) # A list with each stiffness matrix in a natural system

    # Consideration of hinges
    if hinge[0,i] == 1:
        kn = kn - kn[0:3, 1]*kn[1, 0:3] / kn[1, 1]
    if hinge[1, i] == 1:
        kn = kn - kn[0:3, 2]*kn[2, 0:3] / kn[2, 2]

    # Transformation to the local system
    T = np.matrix([[-c, -s, 0, c, s, 0],
                    [-s/L, c/L, 1, s/L, -c/L, 0],
                    [-s/L, c/L, 0, s/L, -c/L, 1]])
    kl = (T.transpose().dot(kn)).dot(T)
    kll.append(kl)

    # Global stiffness matrix assembly
    Lm = [3*(inc[0, i]+1)-2, 3*(inc[0, i]+1)-1, 3*(inc[0, i]+1), 3*(inc[1, i]+1)-2, 3*(inc[1, i]+1)-1, 3*(inc[1, i]+1)]
    for j in range(6):
        for k in range(6):
            K[Lm[j]-1, Lm[k]-1] = K[Lm[j]-1, Lm[k]-1] + kl[j, k]

# Nodal loading vector assembly and consideration of boundary conditions
p = np.zeros((ngl,1))
for i in range(nn):
    for j in range(3):
        p[i * 3 + j] = load[j, i] * (1 + C * support[j, i])
        K[i * 3 + j, i * 3 + j] = K[i * 3 + j, i * 3 + j] + C * support[j, i]

# Solving  the problem
d = np.linalg.solve(K, p)
print(f'Displacements in each node:\n {d}')

# Post processing of the results
epsilonl = np.zeros((6,1)); p = np.zeros((ngl,1))
normal = []; shear = []; bending = []
for i in range(nb):
    Lm = [3*(inc[0,i]+1)-2, 3*(inc[0,i]+1)-1, 3*(inc[0,i]+1), 3*(inc[1,i]+1)-2, 3*(inc[1,i]+1)-1, 3*(inc[1,i]+1)]
    for j in range(6):
        epsilonl[j,0] = d[Lm[j]-1]
    sigmal = kll[i].dot(epsilonl)
    for j in range(6):
        p[Lm[j]-1] = p[Lm[j]-1] + sigmal[j]
    dx = coord[0, inc[1, i]] - coord[0, inc[0, i]]
    dy = coord[1, inc[1, i]] - coord[1, inc[0, i]]
    L = np.sqrt(dx**2+dy**2)
    c = dx/L; s = dy/L
    Ta = np.matrix([[-c, -s, 0, 0, 0, 0],
                    [-s, c, 0, 0, 0, 0],
                    [0, 0, -1, 0, 0, 0],
                    [0, 0, 0, c, s, 0],
                    [0, 0, 0, s, -c, 0],
                    [0, 0, 0, 0, 0, 1]])
    sigmaa = Ta.dot(sigmal)
    print(f'Forces and moments in element {i+1}:\n {sigmaa}')

    # Construction of tables to plot the diagrams
    normal.append([[coord[0, inc[0, i]], coord[1, inc[0, i]]],
                    [-cn * s * sigmaa[0] + coord[0, inc[0, i]], cn * c * sigmaa[0] + coord[1, inc[0, i]]],
                    [-cn * s * sigmaa[3] + coord[0, inc[1, i]], cn * c * sigmaa[3] + coord[1, inc[1, i]]],
                    [coord[0, inc[1, i]], coord[1, inc[1, i]]]])

    shear.append([[coord[0, inc[0, i]], coord[1, inc[0, i]]],
                    [-cq*s*sigmaa[1]+coord[0, inc[0, i]], cq*c*sigmaa[1]+coord[1, inc[0, i]]],
                    [-cq*s*sigmaa[4]+coord[0, inc[1, i]], cq*c*sigmaa[4]+coord[1, inc[1, i]]],
                    [coord[0, inc[1, i]], coord[1, inc[1, i]]]])

    bending.append([[coord[0, inc[0, i]], coord[1, inc[0, i]]],
                    [cf*s*sigmaa[2]+coord[0, inc[0, i]], -cf*c*sigmaa[2]+coord[1, inc[0, i]]],
                    [cf*s*sigmaa[5]+coord[0, inc[1, i]], -cf*c*sigmaa[5]+coord[1, inc[1, i]]],
                    [coord[0, inc[1, i]], coord[1, inc[1, i]]]])

    xi = sym.Symbol('xi')
    disp = np.matrix([[xi-1, 0, 0, xi, 0, 0],[0, 1-3*xi**2+2*xi**3, (-xi**2+2*xi-1)*xi*L, 0, 2*xi**3-3*xi**2, (xi-1)*L*xi**2]])
    T = np.matrix([[-c,-s,0,c,s,0],[-s/L,c/L,1,s/L,-c/L,0],[-s/L,c/L,0,s/L,-c/L,1]])
    epsilonn = T.dot(epsilonl)
    knli = knl[i]
    if hinge[0,i] == 1:
        if hinge[1,i] == 1:
            T = np.linalg.inv(knli[1:3,1:3]).dot(knli[1:3,0])
            epsilonl[2] = epsilonl[2] - epsilonn[1] - T[0] * epsilonn[0]
            epsilonl[5] = epsilonl[5] - epsilonn[2] - T[1] * epsilonn[0]
        else:
            epsilonl[2] = epsilonl[2] - knli[1,0:3].dot(epsilonn) / knli[1, 1]
    elif hinge[1, i] == 1:
        epsilonl[5] = epsilonl[5] - knli[2,0:3].dot(epsilonn) / knli[2, 2]
    uv = (disp.dot(Ta)).dot(epsilonl)*cdef
    xd = coord[0, inc[0, i]] + xi * dx + c * uv[0] - s * uv[1]
    yd = coord[1, inc[0, i]] + xi * dy + s * uv[0] + c * uv[1]

    # Displacements
    datad_x = np.zeros((11, 1)); datad_y = np.zeros((11, 1))
    for j in range(0, 11, 1):
        datad_x[j] = xd[0, 0].subs(xi, j / 10)
        datad_y[j] = yd[0, 0].subs(xi, j / 10)
    axs[0, 0].plot(datad_x, datad_y, 'r')

# Plot of the structure, displacements, normal forces, shear forces and bending moments
# Plot of the structure
for i in range(nb):
    xdata = np.matrix([[coord[0,inc[0,i]]],[coord[0,inc[1,i]]]])
    ydata = np.matrix([[coord[1,inc[0,i]]],[coord[1,inc[1,i]]]])
    axs[0, 0].plot(xdata,ydata,'black')
    axs[1, 0].plot(xdata, ydata, 'black')
    axs[0, 1].plot(xdata, ydata, 'black')
    axs[1, 1].plot(xdata, ydata, 'black')

# Plot of the supports
ax = np.matrix([[0,0],[-1,.5],[-1,-.5],[0,0]])
ay = np.matrix([[0,0],[-.5,-1],[.5,-1],[0,0]])
ar = 0.5*np.matrix([[-1,-1],[1,-1],[1,1],[-1,1],[-1,-1]])
ipx = 0; ipy = 0; ipr = 0; irot = 0
supportx_x = []; supporty_x = []; supportr_x = []
supportx_y = []; supporty_y = []; supportr_y = []
for i in range(nn):
    if support[0,i] == 1:
        ipx = ipx + 1
        for j in range(4):
            supportx_x.append(ax[j,0]*cap + coord[0, i])
            supportx_y.append(ax[j,1]*cap + coord[1, i])
        axs[0, 0].plot(supportx_x[4*ipx-4:4*ipx],supportx_y[4*ipx-4:4*ipx],'b')
        axs[1, 0].plot(supportx_x[4*ipx-4:4*ipx],supportx_y[4*ipx-4:4*ipx],'b')
        axs[0, 1].plot(supportx_x[4*ipx-4:4*ipx],supportx_y[4*ipx-4:4*ipx],'b')
        axs[1, 1].plot(supportx_x[4*ipx-4:4*ipx],supportx_y[4*ipx-4:4*ipx],'b')
    if support[1,i] == 1:
        ipy = ipy + 1
        for j in range(4):
            supporty_x.append(ay[j,0]*cap + coord[0, i])
            supporty_y.append(ay[j,1]*cap + coord[1, i])
        axs[0, 0].plot(supporty_x[4*ipy-4:4*ipy],supporty_y[4*ipy-4:4*ipy],'b')
        axs[1, 0].plot(supporty_x[4*ipy-4:4*ipy],supporty_y[4*ipy-4:4*ipy],'b')
        axs[0, 1].plot(supporty_x[4*ipy-4:4*ipy],supporty_y[4*ipy-4:4*ipy],'b')
        axs[1, 1].plot(supporty_x[4*ipy-4:4*ipy],supporty_y[4*ipy-4:4*ipy],'b')
    if support[2,i] == 1:
        ipr = ipr + 1
        for j in range(5):
            supportr_x.append(ar[j,0]*cap + coord[0, i])
            supportr_y.append(ar[j,1]*cap + coord[1, i])
        axs[0, 0].plot(supportr_x[5*ipr-5:5*ipr],supportr_y[5*ipr-5:5*ipr],'b')
        axs[1, 0].plot(supportr_x[5*ipr-5:5*ipr],supportr_y[5*ipr-5:5*ipr],'b')
        axs[0, 1].plot(supportr_x[5*ipr-5:5*ipr],supportr_y[5*ipr-5:5*ipr],'b')
        axs[1, 1].plot(supportr_x[5*ipr-5:5*ipr],supportr_y[5*ipr-5:5*ipr],'b')

# Plot of the hinges
t_hinge_x = []; t_hinge_y = []
for i in range(nb):
    if hinge[0,i] == 1:
        irot = irot + 1
        t_hinge_x.append(coord[0,inc[0,i]]*(1-d_hinge)+coord[0,inc[1,i]]*d_hinge)
        t_hinge_y.append(coord[1, inc[0, i]] * (1 - d_hinge) + coord[1, inc[1, i]] * d_hinge)
    if hinge[1, i] == 1:
        irot = irot + 1
        t_hinge_x.append(coord[0, inc[1, i]] * (1 - d_hinge) + coord[0, inc[0, i]] * d_hinge)
        t_hinge_y.append(coord[1, inc[1, i]] * (1 - d_hinge) + coord[1, inc[0, i]] * d_hinge)
for j in range(irot):
    axs[0, 0].plot(t_hinge_x[j],t_hinge_y[j],'gray',marker='o')
    axs[1, 0].plot(t_hinge_x[j],t_hinge_y[j],'gray',marker='o')
    axs[0, 1].plot(t_hinge_x[j],t_hinge_y[j],'gray',marker='o')
    axs[1, 1].plot(t_hinge_x[j],t_hinge_y[j],'gray',marker='o')

# Bending moment
bend_x = np.zeros((nb*4,1)); bend_y = np.zeros((nb*4,1)); count = 0
for i in range(nb):
    for j in range(4):
        bend_x[count] = bending[i][j][0]
        bend_y[count] = bending[i][j][1]
        count = count+1

# Shear forces
shear_x = np.zeros((nb*4,1));shear_y = np.zeros((nb*4,1)); count = 0
for i in range(nb):
    for j in range(4):
        shear_x[count] = shear[i][j][0]
        shear_y[count] = shear[i][j][1]
        count = count+1

# Normal forces
norm_x = np.zeros((nb*4,1));norm_y = np.zeros((nb*4,1)); count = 0
for i in range(nb):
    for j in range(4):
        norm_x[count] = normal[i][j][0]
        norm_y[count] = normal[i][j][1]
        count = count+1

axs[0, 0].set_title("Displacements")
axs[0, 0].axis('scaled')
axs[1, 0].plot(bend_x,bend_y,'lime')
axs[1, 0].set_title("Bending moment")
axs[1, 0].sharex(axs[0, 0])
axs[1, 0].axis('scaled')
axs[0, 1].plot(shear_x,shear_y,'y')
axs[0, 1].set_title("Shear")
axs[0, 1].axis('scaled')
axs[1, 1].plot(norm_x,norm_y,'mediumturquoise')
axs[1, 1].set_title("Normal")
axs[1, 1].axis('scaled')
fig.tight_layout()
plt.show()
