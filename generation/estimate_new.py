from rdgan import GanRDD
import pandas as pd

def run_training_loop(dname, xbound, ybound):
    data = pd.read_csv("data/cleaned/"+dname+".csv")
    print(dname)
    ganModel = GanRDD(dname, data, xbound = xbound, ybound = ybound, epochs=4000)
    #leeGAN.load()
    ganModel.train()
    ganModel.save_models()

    print("GROUND TRUTH: ")
    print(ganModel.groundTruth())
    ganModel.evaluate_results()
    ganModel.generate_data(save=True, sample_size = 1e7)

#run_training_loop("meyersson", [-100, 100], [0, 100])
#run_training_loop("senate", [-100, 100], [0, 100])
#run_training_loop("brazil", [-100, 100], [-100, 100])
run_training_loop("curved", [-1, 1], [-1000, 1000])
