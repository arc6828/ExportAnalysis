import matplotlib.pyplot as plt
# from matplotlib.pyplot import figure, show
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import difflib
from sklearn.neighbors import KDTree


class ZoomPan:
    def __init__(self):
        self.press = None
        self.cur_xlim = None
        self.cur_ylim = None
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.xpress = None
        self.ypress = None
        self.zoom_level = 0
        self.k = 5


    def zoom_factory(self, ax, base_scale = 2.):
        def zoom(event):
            cur_xlim = ax.get_xlim()
            cur_ylim = ax.get_ylim()

            xdata = event.xdata # get event x location
            ydata = event.ydata # get event y location

            if event.button == 'up':
                # deal with zoom in
                scale_factor = 1 / base_scale
                self.zoom_level = self.zoom_level + 1
            elif event.button == 'down':
                # deal with zoom out
                scale_factor = base_scale  
                self.zoom_level = self.zoom_level - 1
            else:
                # deal with something that should never happen
                scale_factor = 1
                print (event.button)

            new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
            new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor

            relx = (cur_xlim[1] - xdata)/(cur_xlim[1] - cur_xlim[0])
            rely = (cur_ylim[1] - ydata)/(cur_ylim[1] - cur_ylim[0])

            ax.set_xlim([xdata - new_width * (1-relx), xdata + new_width * (relx)])
            ax.set_ylim([ydata - new_height * (1-rely), ydata + new_height * (rely)])
            ax.figure.canvas.draw()
            ax.set_title(title + "- Zoom Level : [" + str(self.zoom_level) +"]")
            

        fig = ax.get_figure() # get the figure of interest
        fig.canvas.mpl_connect('scroll_event', zoom)

        return zoom

    def pan_factory(self, ax):
        def onPress(event):
            if event.inaxes != ax: return
            self.cur_xlim = ax.get_xlim()
            self.cur_ylim = ax.get_ylim()
            self.press = self.x0, self.y0, event.xdata, event.ydata
            self.x0, self.y0, self.xpress, self.ypress = self.press

        def onRelease(event):
            self.press = None
            ax.figure.canvas.draw()

        def onMotion(event):
            # print(event,event.inaxes)
            if event.inaxes != ax: return
            if self.press is None: 
                #ONHOVER CASE
                # vis = annot.get_visible()
                exist, obj = sc.contains(event)

                if exist:                    
                    index = obj["ind"][0]
                    position = sc.get_offsets()[index]
                    print(exist, obj,  position)            
                    if self.k > len(sc.get_offsets()) : 
                        self.k = len(sc.get_offsets())
                    distance_list , indexes_list = tree.query([position] , k=self.k)   
                    print(indexes_list)
                    for i in range(len(indexes_list[0])) : 
                        index = indexes_list[0][i]
                        position = sc.get_offsets()[index]
                        annotation_list[i].xy = position
                        annotation_list[i].set_text(C[index] + " (" + str(T[index])+")")
                        annotation_list[i].get_bbox_patch().set_facecolor("green")
                        annotation_list[i].get_bbox_patch().set_alpha(0.4)

                        # update_annot(ind)
                        annotation_list[i].set_visible(True)
                    fig.canvas.draw_idle()
                # else:
                #     annotation_list[0].set_visible(False)
                #     fig.canvas.draw_idle()

                # pos = sc.get_offsets()[ind["ind"][0]]
                # annot.xy = pos
                # text = "{}, {}".format(" ".join(list(map(str,ind["ind"]))), 
                #                     " ".join([names[n] for n in ind["ind"]]))
                # annot.set_text(text)
                # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
                # annot.get_bbox_patch().set_alpha(0.4)

                
                return
            dx = event.xdata - self.xpress
            dy = event.ydata - self.ypress
            self.cur_xlim -= dx
            self.cur_ylim -= dy
            ax.set_xlim(self.cur_xlim)
            ax.set_ylim(self.cur_ylim)

            ax.figure.canvas.draw()

        fig = ax.get_figure() # get the figure of interest

        # attach the call back
        fig.canvas.mpl_connect('button_press_event',onPress)
        fig.canvas.mpl_connect('button_release_event',onRelease)
        fig.canvas.mpl_connect('motion_notify_event',onMotion)

        #return the function
        return onMotion

fig = plt.figure(figsize =(8,8))
ax = fig.add_subplot(111, xlim=(-100,100), ylim=(-100,100), autoscale_on=False)
# ax.set_title('Click to zoom')
title = 'export every country in 2019'
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title(title)

# DATA START HERE
df_SU = pd.read_csv("df_SU.csv")
df_SU = df_SU.loc[df_SU['1'] == 2019]
# df_SU = df_SU.loc[df_SU['0'] == "Thailand"]
SU = df_SU.to_numpy()

X = SU[:,2].astype(float) 
Y = SU[:,3].astype(float) 
X = list(np.nan_to_num( X ))
Y = list(np.nan_to_num( Y ))
XY = np.transpose(np.asarray([X,Y]))
tree = KDTree(XY) 
C = list(SU[:,0])
T = list(SU[:,1])
sc = ax.scatter(X , Y , color='r', alpha=0.5)


# ax.scatter(x,y,s,c)
scale = 1.1
zp = ZoomPan()
figZoom = zp.zoom_factory(ax, base_scale = scale)
figPan = zp.pan_factory(ax)


annotation_list = []
for index in range(zp.k):
    # text = C[index] + " (" + str(T[index])+")"
    temp = ax.annotate("", xy=(100,100), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annotation_list.append(temp)

for i in range(len(C)) :
    if C[i] == "Thailand":
        plt.text(X[i], Y[i],   "*"  , size=10, rotation=0.,
                ha="right", va="top",
                bbox=None
                )

# manager = plt.get_current_fig_manager()
# fig.canvas.manager.full_screen_toggle()
plt.show()