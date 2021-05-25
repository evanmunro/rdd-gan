
import wgan

class GanRDD(object):

    def __init__(self, df):
        dfa = df[df['x']>0].copy()
        dfb = df[df['x']<=0].copy()
        self.y1_GAN = ModelWGAN(dfa, context = ['x'])
        self.y0_GAN = ModelWGAN(dfb, context = ['x'])
        self.x_GAN = ModelWGAN(x, context=[])
    def train_models(self):
        return 0
    def load_models(self, path):
        return 0
    def save_models(self):
        return 0
    def calculate_gt(self):
        return gt
    def generate_data(self, size=1e4):
        return data

class ModelWGAN(object):
    def __init__(self, name, data, context):
        self.name = name
        self.dwrapper = 0
        self.critic = 0
        self.generator = 0
        self.specs = wgan.Specifications(self.dwrapper,
                                       batch_size=bsize,
                                       critic_steps=2,
                                       optimizer=wgan.OAdam,
                                       generator_optimizer=wgan.OAdam,
                                       critic_lr = c_lr,
                                       generator_lr = g_lr,
                                       critic_d_hidden = cd,
                                       generator_d_hidden = gd,
                                       critic_gp_factor = factor,
                                       max_epochs = epochs,
                                       generator_d_noise = noise,
                                       print_every=100)
        return 0

    def generate_data(self, df, size=1e4):
        return 0

    def train_model(self):
        wgan.train(self.generator, self.critic, self.y, self.x, self.specs)

    def save_model(self):
        return 0
