"""
Dataset1d  TODO
Dataset1dT

Dataset2d
Dataset2dT  TODO

Dataset3d  TODO
Dataset3dT  TODO
"""

import numpy as np
import torch


class _Dataset:
    def __init__(self) -> None:
        pass

    def update_all_data(self):
        """更新所有数据，根据任务重写此方法"""
        return NotImplementedError

    def external_data(self):
        """若数据集从外部导入，则重写此方法"""
        return NotImplementedError

    def func_res(self):
        """请重写控制方程右端项"""
        return NotImplementedError

    def func_bcs(self):
        """请重写边界条件右端项"""
        return NotImplementedError

    def func_ics(self):
        """请重写初始条件右端项"""
        return NotImplementedError

    def array2tensor(self):
        """将data dict中的array变为tensor，并转移到device"""
        for k, v in self.datad.items():
            self.datad[k] = torch.from_numpy(v).float().to(self.device)

    def tensor2array(self):
        """将data dict中的tensor变为array"""
        for k, v in self.datad.items():
            self.datad[k] = v.detach().cpu().numpy()

    def to(self, device):
        """将data dict中所有的tensor转移到device"""
        for k, v in self.datad.items():
            self.datad[k] = v.to(device)
    
    def statistic(self):
        """获取内部点的统计量，用以标准化"""
        self.datad["mean"] = self.datad["X_res"].mean(axis=0)
        self.datad["std"] = self.datad["X_res"].std(axis=0)
        self.datad["var"] = self.datad["X_res"].var(axis=0)


class Dataset1dT(_Dataset):
    def __init__(self, domain, device):
        self.xmin, self.xmax, self.tmin, self.tmax = domain
        self.device = device
        self.datad = {}  # data dictionary 将数据存放于字典

    def interior_grid(self, nx, nt):
        """对内部网格采样"""
        x = np.linspace(self.xmin, self.xmax, nx)
        t = np.linspace(self.tmin, self.tmax, nt)
        xx, tt = np.meshgrid(x, t)
        X_res = np.stack([xx.flatten(), tt.flatten()], axis=1)
        self.datad["X_res"] = X_res

    def interior_random(self, nr):
        """对内部随机采样"""
        x = np.random.rand(nr) * (self.xmax - self.xmin) + self.xmin
        t = np.random.rand(nr) * (self.tmax - self.tmin) + self.tmin
        X_res = np.stack([x, t], axis=1)
        self.datad["X_res"] = X_res

    def boundary_grid(self, nb):
        """对边界网格采样"""
        x = np.repeat(self.xmin, nb)
        t = np.linspace(self.tmin, self.tmax, nb)
        b_min = np.stack([x, t], axis=1)

        x = np.repeat(self.xmax, nb)
        b_max = np.stack([x, t], axis=1)

        X_bcs = np.concatenate([b_min, b_max], axis=0)
        self.datad["X_bcs"] = X_bcs

    def boundary_random(self, nb):
        """对边界随机采样"""
        x = np.repeat(self.xmin, nb)
        t = np.random.rand(nb) * (self.tmax - self.tmin) + self.tmin
        b_min = np.stack([x, t], axis=1)

        x = np.repeat(self.xmax, nb)
        b_max = np.stack([x, t], axis=1)

        X_bcs = np.concatenate([b_min, b_max], axis=0)
        self.datad["X_bcs"] = X_bcs

    def initial_grid(self, n0):
        """对初始网格采样"""
        x = np.linspace(self.xmin, self.xmax, n0)
        t = np.repeat(self.tmin, n0)
        X_ics = np.stack([x, t], axis=1)
        self.datad["X_ics"] = X_ics

    def initial_random(self, n0):
        """对初始随机采样"""
        x = np.random.rand(n0) * (self.xmax - self.xmin) + self.xmin
        t = np.repeat(self.tmin, n0)
        X_ics = np.stack([x, t], axis=1)
        self.datad["X_ics"] = X_ics


class Dataset2d(_Dataset):
    def __init__(self, domain, device):
        self.xmin, self.xmax, self.ymin, self.ymax = domain
        self.device = device
        self.datad = {}  # data dictionary 将数据存放于字典

    def interior_grid(self, nx, ny):
        """对内部网格采样"""
        x = np.linspace(self.xmin, self.xmax, nx)
        y = np.linspace(self.ymin, self.ymax, ny)
        xx, yy = np.meshgrid(x, y)
        X_res = np.stack([xx.flatten(), yy.flatten()], axis=1)
        self.datad["X_res"] = X_res

    def interior_random(self, nr):
        """对内部随机采样"""
        x = np.random.rand(nr) * (self.xmax - self.xmin) + self.xmin
        y = np.random.rand(nr) * (self.ymax - self.ymin) + self.ymin
        X_res = np.stack([x, y], axis=1)
        self.datad["X_res"] = X_res

    def boundary_grid(self, nb):
        """对边界网格采样"""
        b1 = self.line_points(self.xmin, self.xmax, self.ymin, nb)
        b2 = self.line_points(self.xmin, self.xmax, self.ymax, nb)
        b3 = self.line_points(self.ymin, self.ymax, self.xmin, nb, fory=True)
        b4 = self.line_points(self.ymin, self.ymax, self.ymax, nb, fory=True)
        X_bcs = np.concatenate([b1, b2, b3, b4], axis=0)
        self.datad["X_bcs"] = X_bcs

    def boundary_random(self, nb):
        """对边界网格采样"""
        b1 = self.line_points(self.xmin, self.xmax, self.ymin, nb, method="random")
        b2 = self.line_points(self.xmin, self.xmax, self.ymax, nb, method="random")
        b3 = self.line_points(self.ymin, self.ymax, self.xmin, nb, fory=True, method="random")
        b4 = self.line_points(self.ymin, self.ymax, self.ymax, nb, fory=True, method="random")
        X_bcs = np.concatenate([b1, b2, b3, b4], axis=0)
        self.datad["X_bcs"] = X_bcs

    @staticmethod
    def line_points(x1, x2, y0, n, fory=False, method="grid"):
        """获取一条线上的点，用以边界采样"""
        if method == "grid":
            x = np.linspace(x1, x2, n)
        elif method == "random":
            x = np.random.rand(n) * (x2 - x1) + x1
        else:
            return NotImplementedError

        y = np.repeat(y0, n)

        if not fory:
            X = np.stack([x, y], axis=1)
        else:
            X = np.stack([y, x], axis=1)
        return X


class Dataset2dT(_Dataset):
    def __init__(self, domain, device):
        self.xmin, self.xmax, self.ymin, self.ymax, self.tmin, self.tmax = domain
        self.device = device
        self.datad = {}  # data dictionary 将数据存放于字典

    def interior_grid(self, nx, ny, nt):
        """对内部网格采样"""
        x = np.linspace(self.xmin, self.xmax, nx)
        y = np.linspace(self.ymin, self.ymax, ny)
        t = np.linspace(self.tmin, self.tmax, nt)

        xx, yy = np.meshgrid(x, y)

        xxx = np.repeat(xx[:, :, None], nt, axis=2)
        yyy = np.repeat(yy[:, :, None], nt, axis=2)

        tt = np.repeat(t[None, :], nx, axis=0)
        ttt = np.repeat(tt[None, :, :], ny, axis=0)

        X_res = np.stack([xxx.flatten(), yyy.flatten(), ttt.flatten()], axis=1)
        self.datad["X_res"] = X_res

    def interior_random(self, nr):
        """对内部随机采样"""
        x = np.random.rand(nr) * (self.xmax - self.xmin) + self.xmin
        y = np.random.rand(nr) * (self.ymax - self.ymin) + self.ymin
        t = np.random.rand(nr) * (self.tmax - self.tmin) + self.tmin
        X_res = np.stack([x, y, t], axis=1)
        self.datad["X_res"] = X_res

    def boundary_random(self, nb):
        """对边界随机采样"""
        x = np.repeat(self.xmin, nb)
        y = np.random.rand(nb) * (self.ymax - self.ymin) + self.ymin
        t = np.random.rand(nb) * (self.tmax - self.tmin) + self.tmin
        b_xmin = np.stack([x, y, t], axis=1)

        x = np.repeat(self.xmax, nb)
        # y = np.random.rand(nb) * (self.ymax - self.ymin) + self.ymin
        # t = np.random.rand(nb) * (self.tmax - self.tmin) + self.tmin
        b_xmax = np.stack([x, y, t], axis=1)

        x = np.random.rand(nb) * (self.xmax - self.xmin) + self.xmin
        y = np.repeat(self.ymin, nb)
        t = np.random.rand(nb) * (self.tmax - self.tmin) + self.tmin
        b_ymin = np.stack([x, y, t], axis=1)

        # x = np.random.rand(nb) * (self.xmax - self.xmin) + self.xmin
        y = np.repeat(self.ymax, nb)
        # t = np.random.rand(nb) * (self.tmax - self.tmin) + self.tmin
        b_ymax = np.stack([x, y, t], axis=1)

        X_bcs = np.concatenate([b_xmin, b_xmax, b_ymin, b_ymax], axis=0)
        self.datad["X_bcs"] = X_bcs

    def initial_random(self, n0):
        """对初始随机采样"""
        x = np.random.rand(n0) * (self.xmax - self.xmin) + self.xmin
        y = np.random.rand(n0) * (self.ymax - self.ymin) + self.ymin
        t = np.repeat(self.tmin, n0)
        X_ics = np.stack([x, y, t], axis=1)
        self.datad["X_ics"] = X_ics
