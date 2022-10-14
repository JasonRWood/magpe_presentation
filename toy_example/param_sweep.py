import runner
import pandas as pd
import imageio as im
import matplotlib.pyplot as plt
import os

sol = runner.PySolver()

beta_max = 5.0
alpha_max = 5.0
bs = [2.0 + i/20 for i in range(21)]
q = 0.1
d = 0.5
gamma = 0.5
seed = 100
alpha_init = 10
fig_base_string = lambda x: f"../figures/temp_figure_{x}.png"
for i, b in enumerate(bs):
    seed = 100 + i
    sol.alpha_ad_dyn(beta_max, alpha_max, b, q, d, gamma, seed, alpha_init)

    df = pd.read_csv(f"../data/evo_sims/data_set{seed}.csv")
    alpha_vals = df["alpha_val"].values
    evo_steps = df["Evolutionary_step"].values

    evo_unique = []
    S_density = []
    I_density = []
    for step in set(evo_steps):
        dft = df[df["Evolutionary_step"]==step]
        S_density.append(dft["Density_of_Hosts"].iloc[0])
        I_density.append(sum(dft["Density_of_parasite"].values))
        evo_unique.append(step)

    fig = plt.figure(figsize = (10, 10))
    gs = fig.add_gridspec(2,2)
    ax = [fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0,1]), fig.add_subplot(gs[1,1])]

    ax[0].scatter(alpha_vals,evo_steps)
    ax[0].set_xlim([0, alpha_max])
    ax[0].set_ylim([0, evo_unique[-1]])
    ax[0].set_xlabel(r"$\alpha$ value", fontsize = 18)
    ax[0].set_ylabel("Evolutionary Time Step", fontsize = 18)

    ax[1].plot(evo_unique, S_density)
    
    ax[1].set_xlim([0, evo_unique[-1]])
    ax[1].set_ylim([0, max(max(S_density), 2)*1.5])
    ax[1].set_ylabel(r"Host Density", fontsize = 18)

    
    ax[2].plot(evo_unique, I_density)
    ax[2].set_xlim([0, max(evo_steps)])
    ax[2].set_ylim([0, 10])
    ax[2].set_ylabel(r"Parasite Density", fontsize = 18)
    fig.suptitle(f"Evolutionary and Ecological dynamics when b = {b}", fontsize = 22)
    plt.savefig(fig_base_string(i), bbox_inches = "tight")
    plt.close()
    # print(df)


with im.get_writer("../figures/Varying_b.gif", mode='I') as writer:
    for i in range(len(bs)):
        image = im.imread(fig_base_string(i))
        writer.append_data(image)
        os.remove(fig_base_string(i))

b = 2.0

gammas = [0.0 + i/20 for i in range(21)]

for i, gamma in enumerate(gammas):
    seed = 100 + i
    sol.alpha_ad_dyn(beta_max, alpha_max, b, q, d, gamma, seed, alpha_init)

    df = pd.read_csv(f"../data/evo_sims/data_set{seed}.csv")
    alpha_vals = df["alpha_val"].values
    evo_steps = df["Evolutionary_step"].values

    evo_unique = []
    S_density = []
    I_density = []
    for step in set(evo_steps):
        dft = df[df["Evolutionary_step"]==step]
        S_density.append(dft["Density_of_Hosts"].iloc[0])
        I_density.append(sum(dft["Density_of_parasite"].values))
        evo_unique.append(step)

    fig = plt.figure(figsize = (10, 10))
    gs = fig.add_gridspec(2,2)
    ax = [fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0,1]), fig.add_subplot(gs[1,1])]

    ax[0].scatter(alpha_vals,evo_steps)
    ax[0].set_xlim([0, alpha_max])
    ax[0].set_ylim([0, evo_unique[-1]])
    ax[0].set_xlabel(r"$\alpha$ value", fontsize = 18)
    ax[0].set_ylabel("Evolutionary Time", fontsize = 18)

    ax[1].plot(evo_unique, S_density)
    
    ax[1].set_xlim([0, evo_unique[-1]])
    ax[1].set_ylim([0, max(max(S_density), 2)*1.5])
    ax[1].set_ylabel(r"Host Density", fontsize = 18)
    ax[1].set_xlabel("Evolutionary Time", fontsize = 18)

    
    ax[2].plot(evo_unique, I_density)
    ax[2].set_xlim([0, max(evo_steps)])
    ax[2].set_ylim([0, 10])
    ax[2].set_ylabel(r"Parasite Density", fontsize = 18)
    ax[2].set_xlabel("Evolutionary Time", fontsize = 18)
    fig.suptitle(fr"Evolutionary and Ecological dynamics when $\gamma$ = {gamma}", fontsize = 22)
    plt.savefig(fig_base_string(i), bbox_inches = "tight")
    plt.close()
    # print(df)


with im.get_writer("../figures/Varying_gamma.gif", mode='I') as writer:
    for i in range(len(gammas)):
        image = im.imread(fig_base_string(i))
        writer.append_data(image)
        os.remove(fig_base_string(i))