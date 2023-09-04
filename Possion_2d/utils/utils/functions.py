import time
import os
import random

import torch

import numpy as np


def used_time(start_time):
    seconds = time.time() - start_time
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    print(f'used time: {hours:.0f}h {minutes:.0f}m {seconds:.0f}s')


def seed_torch(seed=1024):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)  # 为了禁止hash随机化，使得实验可复现
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    # torch.use_deterministic_algorithms(True)  # 有检查操作，看下文区别


def filter2d(X, X1, X2):
    def _filter2d(X, x1, x2):
        # 索引 1 个元素
        a = np.where(X[:, 0] == x1)[0]
        b = np.where(X[:, 1] == x2)[0]
        idx = a[np.isin(a, b)]
        return X[idx], idx

    # 循环索引若干个元素
    idxes = []
    for i in range(len(X1)):
        for j in range(len(X2)):
            _, idx = _filter2d(X, X1[i], X2[j])
            idxes.append(idx[0])
    idxes = np.array(idxes)

    return X[idxes], idxes
