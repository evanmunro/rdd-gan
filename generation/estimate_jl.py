from rdgan import GanRDD
import pandas as pd

jl_data = pd.read_csv("data/cleaned/jl_math.csv")
jlGAN = GanRDD("jl_math", jl_data, xbound = [jl_data.x.min(), jl_data.x.max()],
                                   ybound = [jl_data.y.min(), jl_data.y.max()],
                                   epochs = 750)
#jlGAN.load()
jlGAN.train()
jlGAN.save_models()

print("GROUND TRUTH: ")
print(jlGAN.groundTruth())
jlGAN.evaluate_results()
jlGAN.generate_data(save=True, sample_size = 1e7)
