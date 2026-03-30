import numpy as np

class Linear:
    def __init__(self, module) -> None:
        self.weight = module.weight.detach().numpy()
        self.bias = module.bias.detach().numpy() if module.bias is not None else None

    def __call__(self, x):
        x = x @ self.weight.T
        
        if self.bias is not None:
            x += self.bias

        return x

    def __repr__(self) -> str:
        return f'Linear({self.weight.shape[1]}, {self.weight.shape[0]})'
    
class BatchNorm1d:
    def __init__(self, module) -> None:
        self.weight = module.weight.detach().numpy()
        self.bias = module.bias.detach().numpy()
        self.eps = module.eps
        self.running_mean = module.running_mean.detach().numpy()
        self.running_var = module.running_var.detach().numpy()

    def __call__(self, x):
        return (x - self.running_mean) / np.sqrt(self.running_var + self.eps) * self.weight + self.bias
    
    def __repr__(self) -> str:
        return f'BatchNorm1d({self.weight.shape[0]})'
    
class Sigmoid:
    def __init__(self, module) -> None:
        pass

    def __call__(self, x):
        return 1 / (1 + np.exp(-x))
    
    def __repr__(self) -> str:
        return f'Sigmoid()'
    
class Dropout:
    def __init__(self, module) -> None:
        self.p = module.p

    def __call__(self, x):
        return x
    
    def __repr__(self) -> str:
        return f'Dropout({self.p:g})'
    
modules_names = {
    'Linear': Linear,
    'BatchNorm1d': BatchNorm1d,
    'Sigmoid': Sigmoid,
    'Dropout': Dropout
}

class Sequential:
    def __init__(self, module) -> None:
        self.layers = [modules_names[layer.__class__.__name__](layer) for layer in module]

    def __call__(self, x):
        for layer in self.layers:
            x = layer(x)
        return x