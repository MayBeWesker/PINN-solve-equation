"""
包括MLP ResNet MultiHeadMLP及其变种
MLP
ModifiedMLP
ResNet
ModifiedResNet
MultiHeadMLP
MultiHeadModifiedMLP
"""
import torch.nn as nn


class MLP(nn.Module):
    def __init__(self, mlp_layers):
        super(MLP, self).__init__()
        self.model = nn.Sequential()
        for i in range(len(mlp_layers) - 2):
            self.model.add_module(f'fc{i + 1}', nn.Linear(mlp_layers[i], mlp_layers[i + 1], bias=True))
            self.model.add_module(f'act{i + 1}', nn.Tanh())
        self.model.add_module(f'fc{len(mlp_layers) - 1})', nn.Linear(mlp_layers[-2], mlp_layers[-1], bias=False))

    def forward(self, X):
        return self.model(X)


class ModifiedMLP(nn.Module):
    def __init__(self, mlp_layers):
        super(ModifiedMLP, self).__init__()

        self.encoder_u = nn.Sequential()
        self.encoder_u.add_module('fc_u', nn.Linear(
            mlp_layers[0], mlp_layers[1], bias=True))
        self.encoder_u.add_module('act_u', nn.Tanh())

        self.encoder_v = nn.Sequential()
        self.encoder_v.add_module('fc_v', nn.Linear(
            mlp_layers[0], mlp_layers[1], bias=True))
        self.encoder_v.add_module('act_v', nn.Tanh())

        self.model = nn.Sequential()
        for i in range(len(mlp_layers)-2):
            layer = nn.Sequential()
            layer.add_module(f'fc{i}', nn.Linear(
                mlp_layers[i], mlp_layers[i+1], bias=True))
            layer.add_module(f'act{i}', nn.Tanh())
            self.model.add_module(f'layer{i}', layer)

        last_layer = nn.Sequential()
        last_layer.add_module(
            f'fc{len(mlp_layers)-2}', nn.Linear(mlp_layers[-2], mlp_layers[-1], bias=False))
        self.model.add_module(f'layer{len(mlp_layers)-2}', last_layer)

        # for param in self.parameters():
        #     if len(param.shape) > 1:
        #         nn.init.kaiming_normal_(param)

    def forward(self, X):
        u = self.encoder_u(X)
        v = self.encoder_v(X)

        for i in range(len(self.model) - 1):
            X = self.model[i](X)
            X = X / 2.
            X = (1 - X) * u + X * v
        return self.model[-1](X)


class ResNet(nn.Module):
    def __init__(self, nn_layers):
        super(ResNet, self).__init__()
        self.model = nn.Sequential()

        first_layer = nn.Sequential()
        first_layer.add_module(f'fc0', nn.Linear(
            nn_layers[0], nn_layers[1], bias=True))
        first_layer.add_module(f'act0', nn.Tanh())
        self.model.add_module(f'first', first_layer)

        for i in range(1, len(nn_layers)-2):
            block = nn.Sequential()
            block.add_module(f'fc{i}_0', nn.Linear(
                nn_layers[i], nn_layers[i+1], bias=True))
            block.add_module(f'act{i}_0', nn.Tanh())
            block.add_module(f'fc{i}_1', nn.Linear(
                nn_layers[i], nn_layers[i+1], bias=True))
            block.add_module(f'act{i}_1', nn.Tanh())
            self.model.add_module(f'block{i}', block)

        last_layer = nn.Sequential()
        last_layer.add_module(
            f'fc{len(nn_layers)-2}', nn.Linear(nn_layers[-2], nn_layers[-1], bias=False))
        self.model.add_module(f'last', last_layer)

        # for param in self.parameters():
        #     if len(param.shape) > 1:
        #         nn.init.kaiming_normal_(param)

    def forward(self, X):
        X = self.model[0](X)
        for i_block in range(1, len(self.model) - 1):
            # 残差连接
            X_ = self.model[i_block](X)
            X = X_ + X
        return self.model[-1](X)


class ModifiedResNet(nn.Module):

    def __init__(self, mlp_layers):
        super(ModifiedResNet, self).__init__()

        self.encoder_u = nn.Sequential()
        self.encoder_u.add_module('fc_u', nn.Linear(
            mlp_layers[0], mlp_layers[1], bias=True))
        self.encoder_u.add_module('act_u', nn.Tanh())

        self.encoder_v = nn.Sequential()
        self.encoder_v.add_module('fc_v', nn.Linear(
            mlp_layers[0], mlp_layers[1], bias=True))
        self.encoder_v.add_module('act_v', nn.Tanh())

        self.model = nn.Sequential()

        first_layer = nn.Sequential()
        first_layer.add_module(f'fc0', nn.Linear(
            mlp_layers[0], mlp_layers[1], bias=True))
        first_layer.add_module(f'act0', nn.Tanh())
        self.model.add_module(f'first', first_layer)

        for i in range(1, len(mlp_layers)-2):
            block = nn.Sequential()
            block.add_module(f'fc{i}_0', nn.Linear(
                mlp_layers[i], mlp_layers[i+1], bias=True))
            block.add_module(f'act{i}_0', nn.Tanh())
            block.add_module(f'fc{i}_1', nn.Linear(
                mlp_layers[i], mlp_layers[i+1], bias=True))
            block.add_module(f'act{i}_1', nn.Tanh())
            self.model.add_module(f'block{i}', block)

        last_layer = nn.Sequential()
        last_layer.add_module(
            f'fc{len(mlp_layers)-2}', nn.Linear(mlp_layers[-2], mlp_layers[-1], bias=False))
        self.model.add_module(f'last', last_layer)

#         for param in self.parameters():
#             if len(param.shape) > 1:
#                 nn.init.kaiming_normal_(param)

    def forward(self, X):
        u = self.encoder_u(X)
        v = self.encoder_v(X)

        X = self.model[0](X)
        for i_block in range(1, len(self.model) - 1):
            X_ = self.model[i_block](X)
            X = X + X_
            X = X / 2.
            X = (1 - X) * u + X * v
        return self.model[-1](X)


class MultiHeadMLP(nn.Module):
    def __init__(self, mlp_layers):
        super(MultiHeadMLP, self).__init__()

        self.model = nn.Sequential()
        for i in range(len(mlp_layers)-2):
            layer = nn.Sequential()
            layer.add_module(f'fc{i}', nn.Linear(
                mlp_layers[i], mlp_layers[i+1], bias=True))
            layer.add_module(f'act{i}', nn.Tanh())
            self.model.add_module(f'layer{i}', layer)

        # 利用ModuleList构造多头
        self.heads = nn.ModuleList()
        for i in range(mlp_layers[-1]):
            h = nn.Sequential()
            h.add_module(f'head{i}', nn.Linear(mlp_layers[-2], 1, bias=False))
            self.heads.append(h)

#         for param in self.parameters():
#             if len(param.shape) > 1:
#                 nn.init.kaiming_normal_(param)

    def output(self, X):
        """分别经过各个头 拼接输出"""
        out = []
        for h in self.heads:
            out.append(h(X))
        return torch.cat(out, dim=1)

    def forward(self, X):
        X = self.model(X)
        return self.output(X)


class MultiHeadModifiedMLP(nn.Module):
    def __init__(self, mlp_layers):
        super(MultiHeadModifiedMLP, self).__init__()

        self.encoder_u = nn.Sequential()
        self.encoder_u.add_module('fc_u', nn.Linear(
            mlp_layers[0], mlp_layers[1], bias=True))
        self.encoder_u.add_module('act_u', nn.Tanh())

        self.encoder_v = nn.Sequential()
        self.encoder_v.add_module('fc_v', nn.Linear(
            mlp_layers[0], mlp_layers[1], bias=True))
        self.encoder_v.add_module('act_v', nn.Tanh())

        self.model = nn.Sequential()
        for i in range(len(mlp_layers)-2):
            layer = nn.Sequential()
            layer.add_module(f'fc{i}', nn.Linear(
                mlp_layers[i], mlp_layers[i+1], bias=True))
            layer.add_module(f'act{i}', nn.Tanh())
            self.model.add_module(f'layer{i}', layer)

        # 利用ModuleList构造多头
        self.heads = nn.ModuleList()
        for i in range(mlp_layers[-1]):
            h = nn.Sequential()
            h.add_module(f'head{i}', nn.Linear(mlp_layers[-2], 1, bias=False))
            self.heads.append(h)

#         for param in self.parameters():
#             if len(param.shape) > 1:
#                 nn.init.kaiming_normal_(param)

    def output(self, X):
        """分别经过各个头 拼接输出"""
        out = []
        for h in self.heads:
            out.append(h(X))
        return torch.cat(out, dim=1)

    def forward(self, X):
        u = self.encoder_u(X)
        v = self.encoder_v(X)

        for i in range(len(self.model) - 1):
            X = self.model[i](X)
            X = X / 2.
            X = (1 - X) * u + X * v
        return self.output(X)
