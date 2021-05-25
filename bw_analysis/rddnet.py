
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import matplotlib.pyplot as plt

class RDDWeights(nn.Module):
    def __init__(self, layers):
        super().__init__()
        d_in = [1] + layers
        d_out = layers + [1]
        self.layers = nn.ModuleList([nn.Linear(i, o) for i, o in zip(d_in, d_out)])

    def weightNet(self, x):
        for layer in self.layers[:-1]:
            x = F.relu(layer(x))
        return self.layers[-1](x)

    def forward(self, X):
        rlist = [self.weightNet(x) for x in torch.unbind(X, dim=1) ]
        return torch.stack(rlist, dim=1)

#weights are (n, dw)
#y is (n, dw)
#target is (n, 1)
def rddCriterion(weights, y, target, lamda):
    loss = target - torch.sum(weights*y, 1)/torch.sum(weights, 1)
    loss = torch.square(loss) + lamda*torch.norm(weights)
    return loss


def getRandomCutoff(y, x):
    if x[0] > 0:
        xc = np.random.choice(x[x<np.median(x)])
    else:
        xc = np.random.choice(x[x>np.median(x)])
    target = np.reshape(np.mean(y[x == xc]), (1, 1)).astype(np.float32)
    input = x - xc
    if x[0] > 0:
        input = input[x > xc]
        y = y[x>xc]
    else:
        input = input[x < xc]
        y = y[x<xc]
    nw = input.shape[0]
    input = np.reshape(input, (1, nw)).astype(np.float32)
    y = np.reshape(y, (1, nw)).astype(np.float32)
    return torch.from_numpy(input), torch.from_numpy(y), torch.from_numpy(target)

def trainWeights(y, x, layers):
    model = RDDWeights(layers)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    for epoch in range(0, 10000):
        input, values, target = getRandomCutoff(y, x)
        optimizer.zero_grad()
        weights = model(input)
        loss = rddCriterion(weights, values, target, 0.1)
        loss.backward()
        optimizer.step()
    xtorch = np.reshape(x, (1, x.shape[0]))
    xtorch = torch.from_numpy(xtorch.astype(np.float32))
    return model(xtorch).detach().numpy()

def rddNetEstimate(y, x, c=0, layers=[10, 10, 10]):
    y1 = y[x>0]
    x1 = x[x>0]
    y0 = y[x<0]
    x0 = x[x<0]

    weights1 = trainWeights(y1, x1, layers)
    weights0 = trainWeights(y0, x0, layers)
    print(x0)
    print(weights0)
    return np.sum(weights1*y[x>c])/np.sum(weights1) - np.sum(weights0*y[x<c])/np.sum(weights0)
