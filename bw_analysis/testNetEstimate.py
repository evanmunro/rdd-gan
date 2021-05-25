from rddnet import *
import pandas as pd

df = pd.read_csv("../data/cleaned/lee.csv")
y = df.y.values[:500]
x = df.x.values[:500]

y1 = y[x>0]
x1 = x[x>0]
y0 = y[x<0]
x0 = x[x<0]

#input, y, target = getRandomCutoff(y1, x1)
#model = RDDWeights([10, 10, 10])
#model(input)
#trainWeights(y1, x1, [10, 10, 10])

rddNetEstimate(y, x)
