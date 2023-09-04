#!/usr/bin/env python3
#
from scipy.sparse import coo_matrix, csr_matrix, csc_matrix, spdiags, bmat
import argparse



import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)

# 装饰子：指明被装饰函数输入的是笛卡尔坐标点
from fealpy.decorator import cartesian
import time
# 网格工厂：生成常用的简单区域上的网格
from fealpy.mesh import MeshFactory as MF
from fealpy.mesh import HalfEdgeMesh2d

# 均匀剖分的时间离散
from fealpy.timeintegratoralg import UniformTimeLine

# AllenCahn_model
from AllenCahnequation import LShapeData as PDE

# Lagrange 有限元空间
from fealpy.functionspace import LagrangeFiniteElementSpace

# Dirichlet 边界条件
from fealpy.boundarycondition import DirichletBC
from fealpy.tools.show import showmultirate

# solver
from scipy.sparse.linalg import spsolve
from fealpy.tools.show import showmultirate
from fealpy.tools.show import show_error_table
#拷贝对象
import copy

## 参数解析
parser = argparse.ArgumentParser(description=
        """
        单纯形网格（三角形、四面体）网格上任意次有限元方法求解热传导方程
        """)

parser.add_argument('--degree',
        default=1, type=int,
        help='Lagrange 有限元空间的次数, 默认为 1 次.')

parser.add_argument('--ns',
        default=10, type=int,
        help='空间各个方向剖分段数， 默认剖分 10 段.')

parser.add_argument('--nt',
        default=1000, type=int,
        help='时间剖分段数，默认剖分 100 段.')

parser.add_argument('--tol',
        default=0.05, type=float,
        help='自适应加密停止限度，默认设定为 0.05 段.')


parser.add_argument('--T',
        default=1, type=float,
        help='自适应加密停止限度，默认设定为 0.05 段.')

parser.add_argument('--TOLtime',
        default=0.1, type=float,
        help='时间误差限度，默认设定为 0.05.')

args = parser.parse_args()

degree = 1 #只能用p1元
ns = args.ns
nt = args.nt
tol = args.tol
T = args.T
TOLtime = args.TOLtime
theta = 0.7
ctheta = 0.2
pde = PDE()
domain = pde.domain()
tmesh = UniformTimeLine(0, T, nt) # 均匀时间剖分

smesh = pde.init_mesh()
smesh = HalfEdgeMesh2d.from_mesh(smesh, NV=3)


def velocity_matrix_test3( u, q=None):
    '''
    uh^2 \\cdot \phii \\cdot \\phij
    '''
    space = u.space
    mesh = space.mesh
    GD = mesh.geo_dimension()
    qf = space.integrator if q is None else mesh.integrator(q, etype='cell')
    # bcs.shape == (NQ, TD+1)
    # ws.shape == (NQ, )
    bcs, ws = qf.get_quadrature_points_and_weights()
    cellmeasure = space.cellmeasure
    gdof = space.number_of_global_dofs()
    c2d = space.cell_to_dof()
    uvalue = u.value(bcs)
    uvalue1 = uvalue**2
    phi = space.basis(bcs)
    shape = c2d.shape + c2d.shape[1:]
    I = np.broadcast_to(c2d[:, :, None], shape=shape)
    J = np.broadcast_to(c2d[:, None, :], shape=shape)
    val1 = np.einsum('q, qc, qci, qcj, c->cij', ws,
                    uvalue1, phi, phi, cellmeasure)
    M = csr_matrix(
        (val1.flat, (I.flat, J.flat)),
        shape=(gdof, gdof)
    )
    return M


t0 = 0
dt = tmesh.current_time_step_length()
@cartesian
def source0(p):
    return pde.source(p,t0)
e = pde.e
index = 50*0.025**2

space = LagrangeFiniteElementSpace(smesh, p=degree)
uh0 = space.interpolation(pde.init_value)
uh00 = space.function()
uh00[:] = uh0
@cartesian
def source(p):
    return pde.source(p, t1)



t1 = 0
compute_time = 0
total_dof = 0
max_dof = space.number_of_global_dofs()
max_refine_count = 0
time_step_count = 0
start =time.time()
#高斯点和权重，为了计算\int_{t^{n-1}}^{t^{n}} f(x,t) dt
Gw = np.array([0.3478548451,0.3478548451,0.6521451549,0.6521451549])
Gp = np.array([0.8611363116,-0.8611363116,0.3399810436,-0.3399810436])
def plot_color(space, uh, save):
    ipoints = space.interpolation_points()
    xx = ipoints[:, 0]
    yy = ipoints[:, 1]
    zz = np.zeros_like(xx)
    zz[:] = uh[:]
    plt.scatter(xx, yy, c=zz)
    plt.colorbar()
    plt.savefig(save)
    plt.close()
plot_color(space, uh0, "eeeee.png")
gdof = space.number_of_global_dofs()
solution = np.zeros((nt, gdof))
i = 0
while t1 <= T:
    t1 = t0 + dt
    uh1 = space.function()
    L = space.stiff_matrix()*index
    M = space.mass_matrix()
    print("t0:",t0)
    print("t1:",t1)
    Glp = (t1 - t0)/2*Gp + (t1 + t0)/2
    Glw = (t1 - t0)/2*Gw
    @cartesian
    def source_int(p):
        return (pde.source(p, Glp[0])*Glw[0] + pde.source(p, Glp[1])*Glw[1] + pde.source(p, Glp[2])*Glw[2] + pde.source(p, Glp[3])*Glw[3])/dt
    F = space.source_vector(source_int)
    #计算代数系统的右端项b
    b = 1/dt*M@uh0
    uh_old = space.function()
    uh_old[:] = uh0
    #牛顿迭代
    for j in range(20):
        J = velocity_matrix_test3(uh_old)
        #非线性方程求解矩阵
        A = M*1/dt + L  - e**(-2)*M + e**(-2)*J
        #Jacobi
        J = M*1/dt + L + 3*e**(-2)*J - e**(-2)*M
        #牛顿迭代
        res = b - A@uh_old
        uh1[:] = spsolve(J,res)
        uh1[:] += uh_old
        if np.linalg.norm((uh1[:] - uh_old[:])) < 10**(-5):
            print('-----------------')
            print('牛顿迭代次数:',j+1)
            print('norm:', np.linalg.norm((uh1[:] - uh_old[:])) )
            break
        uh_old[:] = uh1
    uh0[:] = uh1[:]
    u = np.zeros(gdof)
    u[:] = uh0[:]
    solution[i] = u
    i +=1
    t0 = t1
    if (abs(t1 - 0) < 1e-5) | (abs(t1 - 0.2) < 1e-5 ) | (abs(t1 -0.5) < 1e-5) | (abs(t1 - 1) < 1e-5):
        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1, projection='3d')
        uh0.add_plot(axes, cmap='rainbow')
        plt.savefig("uh/test-" + str(j) + "time:" + str(t1) + '-dof-'+ str(space.number_of_global_dofs()) +  '.png')
        plt.close()
        save = "uh/heat" + str(j) + "time:" + str(t1) + '-dof-'+ str(space.number_of_global_dofs()) +  '.png'
        plot_color(space, uh0, save)	
        smesh.add_plot(plt)
        plt.savefig('mesh/test-' + str(j) + "time:" + str(t1) + '-dof-'+ str(space.number_of_global_dofs()) +  '.png')
        plt.close()
        

ipoints = space.interpolation_points()
np.save("./xy.npy", ipoints)
np.save('./solution_tu.npy', solution)