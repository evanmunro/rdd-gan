
import wgan
import pandas as pd
import torch
import numpy as np

class GanRDD(object):

    def __init__(self, name, df, xbound, ybound):
        dfa = df[df['x']>0].copy()
        dfb = df[df['x']<=0].copy()
        self.name = name
        self.y1_GAN = GanWrapper(name + "_y1", dfa, outcome=['y'], context = ['x'],
                                 lbound = {'y': ybound[0]}, ubound = {'y':ybound[1]})
        self.y0_GAN = GanWrapper(name + "_y0", dfb, outcome=['y'], context = ['x'],
                                 lbound = {'y': ybound[0]}, ubound = {'y':ybound[1]})
        self.x_GAN = GanWrapper(name + "_x", df, outcome = ['x'], context = [],
                                 lbound = {'x': ybound[0]}, ubound = {'x':ybound[1]})
        self.data = df.copy()

    def train(self):
        self.y1_GAN.train()
        self.y0_GAN.train()
        self.x_GAN.train()
    def load(self, path="trained_models"):
        self.y1_GAN.load(path)
        self.y0_GAN.load(path)
        self.x_GAN.load(path)

    def save_models(self, path="trained_models"):
        self.y1_GAN.save(path)
        self.y0_GAN.save(path)
        self.x_GAN.save(path)

    def groundTruth(self):
        df = pd.DataFrame(np.zeros((1000000, 2)), columns=['x', 'y'])
        dfa = self.y1_GAN.generate(df)
        dfb = self.y0_GAN.generate(df)
        return dfa['y'].mean() - dfb['y'].mean()

    def evaluate_results(self):
        df_fake = self.generate_data(sample_size = self.data.shape[0])
        scatterplot = dict(x= ["x"],
                     y= ["y"],
                     samples = 5000, smooth = 0)
        histogram = dict(variables=['x','y','x','y'],
                       nrow=2, ncol=2)
        wgan.compare_dfs(self.data, df_fake, figsize=5, histogram=histogram, scatterplot=scatterplot,
                         save = True, path = "output")
        return 0

    def generate_data(self, save = False, sample_size=1e4):
        s = self.data.sample(int(sample_size), replace=True)
        s = self.x_GAN.generate(s)
        dfa = s[s['x']>0]
        dfb = s[s['x']<=0]
        dfa = self.y1_GAN.generate(dfa)
        dfb = self.y0_GAN.generate(dfb)
        df_fake = pd.concat([dfa,dfb], axis=0, ignore_index=True)
        if save:
             df_fake.drop(columns=['source']).to_feather("data/generated/{}_generated.feather".format(self.name))
        else:
            return df_fake

class GanWrapper(object):
    def __init__(self, name, df, outcome, context, lbound, ubound):
        self.name = name
        self.dwrapper = wgan.DataWrapper(df, outcome, continuous_lower_bounds=lbound,
                                              continuous_upper_bounds=ubound)
        outcome_scaled, context_scaled = self.dwrapper.preprocess(df)
        self.outcome = outcome_scaled
        self.context = context_scaled

        self.specs = wgan.Specifications(self.dwrapper,
                                       batch_size=64,
                                       critic_steps=2,
                                       optimizer=wgan.OAdam,
                                       generator_optimizer=wgan.OAdam,
                                       critic_lr = 1e-4,
                                       generator_lr = 1e-4,
                                       critic_d_hidden = [128, 128, 128],
                                       generator_d_hidden = [128, 128, 128],
                                       critic_gp_factor = 5,
                                       max_epochs = 2000,
                                       generator_d_noise = 2,
                                       print_every=100)
        self.critic = wgan.Critic(self.specs)
        self.generator = wgan.Generator(self.specs)

    def generate(self, df):
        data = self.dwrapper.apply_generator(self.generator, df)
        return data

    def load(self, path = "trained_models"):
        state_dict_critic = torch.load(path + "/" + self.name + "_C.pth")
        self.critic.load_state_dict(state_dict_critic)
        state_dict_gen = torch.load(path + "/" + self.name + "_G.pth")
        self.generator.load_state_dict(state_dict_gen)

    def train(self):
        wgan.train(self.generator, self.critic, self.outcome, self.context, self.specs)

    def save(self, path="trained_models"):
        torch.save(self.critic.state_dict(), path + "/" + self.name + "_C.pth")
        torch.save(self.generator.state_dict(), path + "/" + self.name + "_G.pth")
