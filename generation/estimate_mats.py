from rdgan import GanRDD
import pandas as pd


mats_data = pd.read_csv("data/cleaned/mats_math.csv")
matsGAN = GanRDD("mats_math", mats_data, xbound = [mats_data.x.min(), mats_data.x.max()],
                                   ybound = [mats_data.y.min(), mats_data.y.max()],
                                   epochs=750)
#matsGAN.load()
matsGAN.train()
matsGAN.save_models()

print("GROUND TRUTH: ")
print(matsGAN.groundTruth())
matsGAN.evaluate_results()
matsGAN.generate_data(save=True, sample_size = 1e7)
