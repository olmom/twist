from pylab import*
import pandas as pd
import seaborn as sns
from sklearn.linear_model import LinearRegression
plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = '13'
plt.rcParams['legend.fontsize'] = '12'
plt.rcParams['xtick.labelsize'] = '12'
plt.rcParams['ytick.labelsize'] = '12'
current_palette = sns.color_palette()
colors = ['#1B9E77', '#D95F02', '#7570B3']
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

##########################################

figure(figsize=(11, 5.5/2.*1.))

# We have chosen to show results for k4 of Gonze model but any
# the bifurcation diagrams of any parameter of Gonze/Goodwin
# could be plotted -- just adjust here
data = pd.read_csv("./results/bifurcations/gonze_k4_allinfo.dat", 
                   sep=" ", header=None)

print(data.keys())

# panel A
subplot(131)
subset=data.loc[data[0]==1]
x1, stFP = subset[3].values.reshape(-1,1), subset[6].values
model1 = LinearRegression()
model1.fit(x1, stFP)
scatter(subset[3], stFP/stFP, s=2, color="k")

subset2=data.loc[data[0]==2]
x2, ustFP = subset2[3].values[1:].reshape(-1,1), subset2[6].values[1:]
model2 = LinearRegression()
model2.fit(x2, ustFP)
scatter(subset2[3].values[1:], ustFP/ustFP, s=2, color="gray")

subset3=data.loc[data[0]==3]
unstFP_fit = model2.predict(subset3[3].values.reshape(-1,1)) 
max_norm, min_norm = subset3[6]/unstFP_fit, subset3[9]/unstFP_fit
plot(subset3[3], max_norm, lw=3.5, color=colors[1])
plot(subset3[3], min_norm, lw=3.5, color=colors[1])

idx = np.where(subset3[3].values[1:] < 0.35)[0][0]
plot(subset3[3].values[idx], max_norm.values[idx], marker='*',
        color="k", markersize=15)
plot(subset3[3].values[idx], min_norm.values[idx], marker='*', 
        color="k", markersize=15)

xlabel("$k_4$ (nM/h)")
ylabel("$x$min,max")
xlim(0.15,0.5)
ylim(0.40, 1.60)

# panel B
subplot(132)
plot(subset3[3],  subset3[5], lw=3.5, color=colors[1])
plot(subset3[3].values[idx], subset3[5].values[idx], marker='*',
        color="k", markersize=15)
xlabel("$k_4$ (nM/h)")
ylabel("period (h)")

# panel C
subplot(133)
plot(subset3[5],  max_norm-min_norm, lw=3.5, color=colors[1])
plot(subset3[5].values[idx], max_norm.values[idx]-min_norm.values[idx],
     marker='*', color="k", markersize=15)
ylabel("$x$max $-x$min")
xlabel("period (h)")
title("increasing $k_4$", x=0.65, y=0.1)

subplots_adjust(top=0.844,
bottom=0.246,
left=0.072,
right=0.982,
hspace=0.2,
wspace=0.414)
#savefig("./figures/fig4.pdf", format='pdf')
