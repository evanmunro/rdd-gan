from rdgan import GanRDD

import pandas as pd


lee_data = pd.read_csv("data/cleaned/lee.csv")
leeGAN = GanRDD("lee", lee_data, xbound = [-1, 1], ybound = [0, 1])
leeGAN.train()
leeGAN.save_models()
#leeGAN.load()

print("GROUND TRUTH: ")
print(leeGAN.groundTruth())
leeGAN.evaluate_results()
leeGAN.generate_data(save=True, sample_size = 1e7)


lee_data = pd.read_csv("data/cleaned/jl_math.csv")
leeGAN = GanRDD("lee", lee_data, xbound = [-1, 1], ybound = [0, 1])
leeGAN.train()
leeGAN.save_models()
#leeGAN.load()

print("GROUND TRUTH: ")
print(leeGAN.groundTruth())
leeGAN.evaluate_results()
leeGAN.generate_data(save=True, sample_size = 1e7)


#jl_data = pd.read_csv("")
#jlGAN =
